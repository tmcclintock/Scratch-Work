import numpy as np

# Generate some fake noisy data.
x = 10 * np.sort(np.random.rand(10))
yerr = 0.0 * 0.2 * np.ones_like(x)
y = np.sin(x) + yerr * np.random.randn(len(x))

import george
from george.kernels import ExpSquaredKernel

# Set up the Gaussian process.
kernel = ExpSquaredKernel(1.0)
gp = george.GP(kernel)

# Pre-compute the factorization of the matrix.
gp.compute(x, yerr)

# Compute the log likelihood.
print(gp.lnlikelihood(y))

t = np.linspace(0, 10, 500)
mu, cov = gp.predict(y, t)
std = np.sqrt(np.diag(cov))

import matplotlib.pyplot as plt
plt.plot(t, mu)
plt.fill_between(t, mu+std, mu-std, alpha=0.2)
plt.errorbar(x, y, yerr, c='k', ls='', marker='o')
plt.show()
