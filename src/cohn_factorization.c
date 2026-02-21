# include <stdio.h>
# include <stdbool.h>
# include <stdlib.h>
# include <signal.h>
# include <omp.h>
# include "flint/flint.h"
# include "flint/ulong_extras.h"
# include "flint/fmpz_poly.h"
# include "flint/fmpz_poly_factor.h"


volatile sig_atomic_t keep_running = true;


typedef struct {
    int integer;
    int divisors;
    int poly_divisors;
} data_collection_t;


void interrupt(int sig) {
    keep_running = false;
}


int count_poly_divisors(fmpz_poly_t poly, fmpz_poly_factor_t poly_fac, ulong n) {
    
    int counter = 0;

    ulong temp = n;

    slong exp = 0;

    ulong digit = 0;
    int running_gcd = 0;
    int factored = 0;


    for (ulong b = 2; b <= n_sqrt(n); b++) {
        
        running_gcd = temp % b;
        temp /= b;
        
        fmpz_poly_set_coeff_ui(poly, exp, running_gcd);
        exp++;
        
        while (temp > 0) {
            digit = temp % b;
            temp /= b;
            running_gcd = n_gcd(running_gcd, digit);

            fmpz_poly_set_coeff_ui(poly, exp, digit);
        
            exp++;
        }

        if (running_gcd > 1) {
            counter++;
        } else {
            fmpz_poly_factor(poly_fac, poly);

            // Did it factor?
            factored = (poly_fac->num > 1) || (poly_fac->num == 1 && poly_fac->exp[0] >1);

            if (factored) {
                counter++;
            }     
        }
        
        // Reset for next loop
        fmpz_poly_zero(poly);
        temp = n;
        exp = 0;
    }
    
    return counter;
}


int main(int argc, char *argv[]) {

    signal(SIGINT, interrupt);

    if (argc != 2) {
        printf("Please provide (exactly) one integer argument.");
        return -1;
    }

    int run_length = atoi(argv[1]);

    data_collection_t *data = (data_collection_t *)malloc((run_length - 1) * sizeof(data_collection_t));
    if (data == NULL) {
        perror("Failed to allocate memory");
        return -1;
    }

    ulong items_processed = 0; 

    // Thread local code block
    #pragma omp parallel
    {

        n_factor_t factors;
        n_factor_init(&factors);

        fmpz_poly_t poly;
        fmpz_poly_factor_t poly_fac;
        fmpz_poly_init(poly);
        fmpz_poly_factor_init(poly_fac);
        
        #pragma omp for
        for (ulong i = 2; i <= (ulong)run_length; i++) {

            if (!keep_running) {
                continue;
            }

            int divisors_counter = 1;
            n_factor(&factors, i, 0);
            for (int j = 0; j < factors.num; j++) {
                divisors_counter *= (factors.exp[j] + 1);
            }
            
            data[i - 2].integer = i;
            data[i - 2].divisors = divisors_counter;
            
            if (divisors_counter == 2) {
                data[i - 2].poly_divisors = 0;
            } else {
                data[i - 2].poly_divisors = count_poly_divisors(poly, poly_fac, i);
            }

            #pragma omp atomic
            items_processed++;
        }
    
        // Cleanup
        fmpz_poly_clear(poly);
        fmpz_poly_factor_clear(poly_fac);

    }

    FILE *file = fopen("data.bin", "wb");
    if (file != NULL) {
        fwrite(data, sizeof(data_collection_t), items_processed, file);
        fclose(file);
        printf("Wrote %lu items to file \n", items_processed);
    } else {
        perror("Failed to open file");
    }

    free(data);
 
    return 0;
}
