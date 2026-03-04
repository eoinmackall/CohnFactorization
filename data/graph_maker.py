import numpy as np
import matplotlib.pyplot as plt

#####################################
#
# Auxillary functions
#
#####################################

def lowest_prime_factor_array(n, semiprimes):
    
    # Setting-up the sieve
    max_val = int(n**0.5) +1
    sieve = np.ones(max_val, dtype=bool)
    
    sieve[:2] = False
    for i in range(2, int(max_val**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    
    primes = np.nonzero(sieve)[0]

    # Finding the smallest factor
    unfactored_idx = np.arange(len(semiprimes))
    p_factors = np.zeros_like(semiprimes)    
    
    for p in primes:
        if len(unfactored_idx) == 0:
            break  

        divisible = (semiprimes[unfactored_idx] % p == 0)
        factored = unfactored_idx[divisible]
        p_factors[factored] = p

        unfactored_idx = unfactored_idx[~divisible]

    return p_factors


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

# Variables for plot y-axis limits
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

num_tests = len(data['num'])
num_semiprimes = len(semiprimes)


#####################################
#
# Graphs below here
#
#####################################

# Plot sizes
plt.rcParams["figure.figsize"] = (4, 4) 

####### Scatter plots #########

# Integers vs divisors
plt.xlim(0, max_full)
plt.ylim(0, hundred_ceil_full)
plt.xlabel("n")
plt.ylabel(r"$\tau(n)$")
plt.title(fr"{full_sample_size:,} integers $n$ with $2\leq n \leq {len(data['num'])+1:,}$ vs $\tau(n)$")
plt.scatter(data['num'][random_indices_full], data['divisors'][random_indices_full],
            marker=".", color="xkcd:blood")
plt.show()

# Integers vs poly_divisors
# plt.xlim(0, max_full)
# plt.ylim(0, hundred_ceil_full)
# plt.title(fr"{full_sample_size:,} integers $n$ with $2 \leq n \leq {len(data['num'])+1:,}$ vs $\tau_{{\beta,2}}(n)$")
# plt.scatter(data['num'][random_indices_full], data['poly_divisors'][random_indices_full],
#             marker=".", color="xkcd:blood")
# plt.show()

# Semi-primes vs poly_divisors
plt.xlim(0, max_sp)
plt.ylim(0, hundred_ceil_full)
plt.xlabel("n")
plt.ylabel(r"$\tau_{\geq 2}\beta(n)$")
plt.title(fr"{semi_prime_sample_size:,} semiprimes $n$ with $2\leq n \leq {len(data['num'])+1:,}$ vs $\tau_{{\geq 2}}\beta(pq)$")
plt.scatter(semiprimes[random_indices_semiprimes], poly_div_for_semis[random_indices_semiprimes],
            marker=".", color = "xkcd:blood")
plt.show()

# log(p-q)/log(pq) vs poly_divisors
p_factors = lowest_prime_factor_array(len(data['num']),semiprimes)
q_factors = semiprimes//p_factors

num = np.log(q_factors-p_factors)
denom = np.log(semiprimes)

ratio = num/denom

plt.xlim(0,1)
plt.xlabel(r"$\log(|p-q|)/\log(pq)$")
plt.ylabel(r"$\tau_{\geq 2}\beta(pq)$")
plt.title(fr"$\log(|p-q|)/\log(pq)$ vs $\tau_{{\geq 2}}\beta(pq)$ for {semi_prime_sample_size:,} semiprimes $n=pq$ with $6\leq n \leq {len(data['num'])+1:,}$")
plt.scatter(ratio[random_indices_semiprimes], poly_div_for_semis[random_indices_semiprimes],
            marker=".", color="xkcd:blood")
plt.show()

####### Cumulative Averages #########

# Cumulative avg. for 'Integers vs poly_divisors'
# cumulative_avg = np.cumsum(data['poly_divisors']) / np.arange(1, len(data['poly_divisors'])+1)
# plt.xlim(0, len(data['num']))
# plt.ylim(0, max_tb//4)
# plt.title(fr"Cumulative average of $\tau_{{\beta,2}}(n)$ for all $2\leq n \leq {len(data['num'])+1:,}$")
# plt.plot(data['num'], cumulative_avg, color='xkcd:blood')
# plt.show()

# Cumulative avg. for 'Integers vs poly_divisors' with 'Integers vs poly_divisors'
cumulative_avg = np.cumsum(data['poly_divisors']) / np.arange(2, len(data['poly_divisors'])+2)
plt.xlim(0, max_full)
plt.ylim(0, hundred_ceil_full)
plt.xlabel("n")
plt.ylabel(r"$\tau_{\geq 2}\beta(n)$")
plt.title(fr"${full_sample_size:,}$ integers $n$ with $2\leq n \leq {len(data['num'])+1:,}$ vs $\tau_{{\geq 2}}\beta(n)$")
plt.scatter(data['num'][random_indices_full], data['poly_divisors'][random_indices_full],
            marker=".", color="xkcd:blood")
plt.plot(data['num'], cumulative_avg, color='xkcd:charcoal', label=r"Cum. Avg. of $\tau_{\geq 2}\beta(n)$")
plt.legend()
plt.show()

# Log-log-plot of avg.'s 'Integers vs poly_divisors'
log_num = np.log(data['num'][4:])
log_poly_divs_avg = np.log(cumulative_avg[4:])
plt.xlim(1, np.ceil(np.log(len(data['num'])-2)))
plt.ylim(-1.3, np.log(max_tb))
plt.xlabel(r"$\log(n)$")
plt.ylabel(r"$\log(\tau_{\geq 2}\beta(n))$")
plt.title(r"Cumulative average of $\tau_{\geq 2}\beta(n)$ with logarithmic axes")

m, c = np.polyfit(log_num, log_poly_divs_avg, 1)
best_fit_line = m * log_num + c

plt.plot(log_num, best_fit_line, color='xkcd:blood', linestyle='--', 
         label=f'Least Squares Regression: $y = {m:.4f}x + ({c:.4f})$')

plt.plot(log_num, log_poly_divs_avg, color="xkcd:charcoal")
plt.grid(True, zorder=0)
plt.axhline(y=0, color="black", linewidth=1)
plt.legend()
plt.show()
