import numpy as np
import matplotlib.pyplot as plt

#####################################
#
# Auxillary functions
#
#####################################

def sieve():
    # Idea: filter through semi-primes, constructing a sieve.
    # For each semiprime, check if it factors among the primes in
    # the sieve; if not, start checking the next multiples of 2...
    # Only need to check up to 6481 since data is only 42,000,000

    return



#####################################
#
# Data set-up
#
#####################################

data_collection = np.dtype([
    ('num', 'i4'),
    ('divisors', 'i4'),
    ('poly_divisors', 'i4')
])

data = np.fromfile('data.bin', dtype=data_collection)
data = data[data['num'] != 0]

is_semiprime = (data['divisors'] == 4)

semiprimes = data['num'][is_semiprime]
poly_div_for_semis = data['poly_divisors'][is_semiprime]

###### Random samples #######
full_sample_size = 1000000
semi_prime_sample_size = 100000
random_indices_full = np.random.choice(len(data['num']), size = full_sample_size, replace = False)
random_indices_semiprimes = np.random.choice(len(semiprimes), size = semi_prime_sample_size, replace = False)

max_full = np.max(data['num'][random_indices_full])
max_sp = np.max(semiprimes[random_indices_semiprimes])
max_tb = np.max(data['poly_divisors'][random_indices_full])
max_sp_tb = np.max(poly_div_for_semis[random_indices_semiprimes])

# Variables for plot limits
temp = max_tb
hundred_ceil_full = 0
while temp > 0:
    hundred_ceil_full += 100
    temp -= 100

temp = max_sp_tb
hundred_ceil_sp = 0
while temp > 0:
    hundred_ceil_sp += 100
    temp -= 100

#####################################
#
# graphs below
#
#####################################

num_tests = len(data['num'])
num_semiprimes = len(semiprimes)

####### Scatter plots #########

# Integers vs divisors
plt.xlim(0, max_full)
plt.ylim(0, hundred_ceil_full)
plt.title(fr"{full_sample_size:,} integers $n$ with $2\leq n \leq {len(data['num'])+1:,}$ vs $\tau(n)$")
plt.scatter(data['num'][random_indices_full], data['divisors'][random_indices_full],
            marker=".", color="xkcd:blood")
plt.show()

# Integers vs poly_divisors
plt.xlim(0, max_full)
plt.ylim(0, hundred_ceil_full)
plt.title(fr"{full_sample_size:,} integers $n$ with $2 \leq n \leq {len(data['num'])+1:,}$ vs $\tau_{{\beta,2}}(n)$")
plt.scatter(data['num'][random_indices_full], data['poly_divisors'][random_indices_full],
            marker=".", color="xkcd:blood")
plt.show()

# Semi-primes vs poly_divisors
plt.xlim(0, max_sp)
plt.ylim(0, hundred_ceil_sp)
plt.title(fr"{semi_prime_sample_size:,} semiprimes $n$ with $2\leq n \leq {len(data['num'])+1:,}$ vs $\tau_{{\beta,2}}(n)$")
plt.scatter(semiprimes[random_indices_semiprimes], poly_div_for_semis[random_indices_semiprimes],
            marker=".", color = "xkcd:blood")
plt.show()

# ?(p-q)/n vs poly_divisors
# Have to write a function to find log(|p-q|)
log_poly_divs_sp = np.log(poly_div_for_semis)
log_sp = np.log(semiprimes)


####### Cumulative Averages #########

# Cumulative avg. for 'Integers vs poly_divisors'
cumulative_avg = np.cumsum(data['poly_divisors']) / np.arange(1, len(data['poly_divisors'])+1)
plt.xlim(0, len(data['num']))
plt.ylim(0, max_tb//4)
plt.title(fr"Cumulative average of $\tau_{{\beta,2}}(n)$ for all $2\leq n \leq {len(data['num'])+1}$")
plt.plot(data['num'], cumulative_avg, color='xkcd:blood')
plt.show()

# Log-log-plot of avg.'s 'Integers vs poly_divisors'
log_num = np.log(data['num'][4:])
log_poly_divs_avg = np.log(cumulative_avg[4:])
plt.xlim(1, np.ceil(np.log(len(data['num'])-2)))
plt.ylim(-1.3, np.log(max_tb))
plt.title(fr"Cumulative average of $\tau_{{\beta,2}}$ with logarithmic axes")

m, c = np.polyfit(log_num, log_poly_divs_avg, 1)
best_fit_line = m * log_num + c

plt.plot(log_num, best_fit_line, color='xkcd:blood', linestyle='--', 
         label=f'Best Fit: $y = ({m:.4f})x + ({c:.4f})$')

plt.plot(log_num, log_poly_divs_avg)
plt.legend()
plt.show()
