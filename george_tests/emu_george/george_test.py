import numpy as np
import sys
import scipy
import scipy.spatial

#Try using the ROTATED e0 data 
yind = 0
y = np.genfromtxt("rotated_efg_means.txt").T[yind]
yerr = np.sqrt(np.genfromtxt("rotated_efg_vars.txt").T[yind])
#print y, yerr

x = np.loadtxt("cosmos.txt")
x = x[:-1] #Remove box 39, it's broken
x = x[:,1:]#Remove the boxnum
x = np.delete(x, 4, 1) #Remove ln10As

xstar = x[-1]
x = x[:-1] #TEST ON BOX 38
ystartrue = y[-1]
ystarvartrue = yerr[-1]**2
y = y[:-1]
yerr = yerr[:-1]
print x.shape, y.shape, yerr.shape, xstar.shape

ndim = len(x[0])
Ncos = len(x)
lguess = (np.max(x,0) - np.min(x,0))/Ncos #Av sep between points in each dimension
#k0guess = np.var(y)

import george
from george.kernels import ExpSquaredKernel

# Set up the Gaussian process.
kernel = ExpSquaredKernel(metric=lguess, ndim=ndim)
gp = george.GP(kernel)

# Pre-compute the factorization of the matrix.
gp.compute(x, yerr)
gp.optimize(x, y, yerr)
mu, cov = gp.predict(y, np.atleast_2d(xstar))
print mu, cov
print ystartrue, ystarvartrue
