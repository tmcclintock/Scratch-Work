import numpy as np
import sys

# Generate some fake noisy data.
#x = 10 * np.sort(np.random.rand(10))
#yerr = 0.0 * 0.2 * np.ones_like(x)
#y = np.sin(x) + yerr * np.random.randn(len(x))

#Try using the e0 data 
yind = 0
y = np.genfromtxt("efg_means.txt").T[yind]
yerr = np.sqrt(np.genfromtxt("efg_vars.txt").T[yind])
ind = np.argmax(y)

x = np.loadtxt("cosmos.txt")
x = x[:-1] #Remove box 39, it's broken
x = x[:,1:]#Remove the boxnum
x = np.delete(x, 4, 1)

p = np.var(x, 0)
meanx = np.mean(x, 0) #The mean cosmo
varx = np.var(x, 0)
for i in range(len(meanx)):
    continue
    print meanx[i], np.sqrt(varx[i]), varx[i]

print y.shape, yerr.shape, x.shape
#sys.exit()

import george
from george.kernels import ExpSquaredKernel
from george.kernels import ExpKernel

# Set up the Gaussian process.
guess = varx#np.ones_like(varx)
kernel = ExpSquaredKernel(metric=guess.copy(), ndim=len(meanx))
gp = george.GP(kernel, mean=np.mean(y))

# Pre-compute the factorization of the matrix.
gp.compute(x, yerr)

# Now find the best hyperparameters
import scipy.optimize as op
print guess

print "prepars = ",kernel.pars
def nll(p):
    gp.kernel[:] = p
    ll = gp.lnlikelihood(y, quiet=False)
    return ll if np.isfinite(ll) else np.inf

def grad_nll(p):
    gp.kernel[:] = p
    return -gp.grad_lnlikelihood(y, quiet=True)

p0 = gp.kernel.vector
results = op.minimize(nll, p0, jac=grad_nll, method='CG')
print results
print results['success']
if results['success']:
    gp.kernel[:] = results.x
else:
    gp.kernel.pars = guess
    gp.kernel[:] = np.log(gp.kernel.pars)
print "pars:",gp.kernel.pars
print "pars?",gp.kernel[:]
print np.log(gp.kernel.pars)
print "varx:",varx
gp.compute(x, yerr)


# Compute the log likelihood.
print(gp.lnlikelihood(y))

cosind = 0
s8 = x[:,cosind]
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
