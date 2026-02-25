import numpy as np
import matplotlib.pyplot as plt

data_collection = np.dtype([
    ('num', 'i4'),
    ('divisors', 'i4'),
    ('poly_divisors', 'i4')
])

data = np.fromfile('data.bin', dtype=data_collection)

is_semiprime = (data['divisors'] == 4)

semiprimes = data['num'][is_semiprime]
poly_div_for_semis = data['poly_divisors'][is_semiprime]

###### Random samples #######
sample_size = 10000
random_indices_full = np.random.choice(len(data['num']), size = sample_size, replace = False)
random_indices_semiprimes = np.random.choice(len(semiprimes), size = sample_size, replace = False)

#####################################
#
# graphs below
#
#####################################

num_tests = len(data['num'])
num_semiprimes = len(semiprimes)

####### Scatter plots #########

# Integers vs divisors
plt.scatter(data['num'][random_indices_full], data['divisors'][random_indices_full])
plt.show()

# Integers vs poly_divisors
plt.scatter(data['num'][random_indices_full], data['poly_divisors'][random_indices_full])
plt.show()

# Semi-primes vs poly_divisors
plt.scatter(semiprimes[random_indices_semiprimes], poly_div_for_semis[random_indices_semiprimes])
plt.title(f"Semiprimes {len(semiprimes)}")
plt.show()

# ?(p-q)/n vs poly_divisors
# Have to write a function to find log(|p-q|)
log_poly_divs = np.log(poly_div_for_semis)

####### Cumulative Averages #########

# Cumulative avg. for 'Integers vs poly_divisors'

# Log-plot of avg.'s 'Integers vs poly_divisors'

