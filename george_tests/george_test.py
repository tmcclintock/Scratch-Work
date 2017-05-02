import numpy as np
import sys

# Generate some fake noisy data.
#x = 10 * np.sort(np.random.rand(10))
#yerr = 0.0 * 0.2 * np.ones_like(x)
#y = np.sin(x) + yerr * np.random.randn(len(x))

#Try using the e0 data 
yind = 5
y = np.genfromtxt("efg_means.txt").T[yind]
yerr = np.sqrt(np.genfromtxt("efg_vars.txt").T[yind])

x = np.loadtxt("cosmos.txt")
x = x[:-1] #Remove box 39, it's broken
x = x[:,1:]#Remove the boxnum

cosind = 0
s8 = x[:,cosind] #try
p = np.var(x, 0)
meanx = np.mean(x, 0) #The mean cosmo
varx = np.var(x, 0)
for i in range(len(meanx)):
    print meanx[i], np.sqrt(varx[i])

print y.shape, yerr.shape, x.shape

import george
from george.kernels import ExpSquaredKernel
from george.kernels import ExpKernel

# Set up the Gaussian process.
#kernel = ExpKernel(np.sqrt(varx), ndim=len(meanx))
kernel = ExpSquaredKernel(varx, ndim=len(meanx))
#kernel = ExpSquaredKernel(1.0, ndim=len(meanx))
gp = george.GP(kernel)

# Pre-compute the factorization of the matrix.
gp.compute(x, yerr)

# Compute the log likelihood.
print(gp.lnlikelihood(y))

s8t = np.linspace(min(s8), max(s8), 500)
t = np.array([meanx*np.ones_like(meanx) for i in xrange(500)])
t[:,cosind] = s8t
mu, cov = gp.predict(y, t)
std = np.sqrt(np.diag(cov))

import matplotlib.pyplot as plt
plt.plot(s8t, mu)
plt.fill_between(s8t, mu+std, mu-std, alpha=0.2)
plt.errorbar(s8, y, yerr, c='k', ls='', marker='o')
plt.show()
