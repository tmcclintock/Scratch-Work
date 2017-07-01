import numpy as np
import sys

# Generate some fake noisy data.
#x = 10 * np.sort(np.random.rand(10))
#yerr = 0.0 * 0.2 * np.ones_like(x)
#y = np.sin(x) + yerr * np.random.randn(len(x))

#Try using the ROTATED e0 data 
yind = 1
#y = np.genfromtxt("efg_means.txt").T[yind]
#yerr = np.sqrt(np.genfromtxt("efg_vars.txt").T[yind])
y = np.genfromtxt("rotated_efg_means.txt").T[yind]
yerr = np.sqrt(np.genfromtxt("rotated_efg_vars.txt").T[yind])
meany = np.mean(y)


x = np.loadtxt("cosmos.txt")
x = x[:-1] #Remove box 39, it's broken
x = x[:,1:]#Remove the boxnum
x = np.delete(x, [0,1,2,3,4,5,6], 1) #Remove ln10As

meanx = np.mean(x, 0) #The mean cosmo
varx = np.var(x, 0)

print y.shape, yerr.shape, x.shape

import george
from george.kernels import ExpSquaredKernel
from george.kernels import ExpKernel

# Set up the Gaussian process.
guess = varx
kernel = ExpSquaredKernel(metric=guess.copy(), ndim=len(meanx))
#gp = george.GP(kernel, mean=meany)

# Pre-compute the factorization of the matrix.
#gp.compute(x, yerr)

# Now find the best hyperparameters
import scipy.optimize as op

def nll(p):
    #amp = np.exp(p[0])
    #amp = p[0]
    #if 4.0 < amp or amp < 0.0: 
    #    return np.inf
    taus = p.copy()
    gp = george.GP(ExpSquaredKernel(metric=taus, ndim=len(meanx)), 
                   mean=meany)
    try:
        gp.compute(x, yerr)
    except (np.linalg.linalg.LinAlgError, ValueError) as e:
        return np.inf
    ll = gp.lnlikelihood(y, quiet=False)
    return ll if np.isfinite(ll) else 1e25

#def grad_nll(p):
#    gp.kernel[:] = p[1:]
#    return -gp.grad_lnlikelihood(y, quiet=True)

#p0 = np.insert(guess.copy(), 0, 1.0)
p0 = guess.copy()

print "guess",guess.shape
#results = op.minimize(nll, p0, jac=grad_nll, method='Nelder-Mead')
results = op.minimize(nll, p0, method='Nelder-Mead')

p = results['x']
gp = george.GP(ExpSquaredKernel(metric=p, ndim=len(meanx)), mean=meany)
#gp = george.GP(p[0] * ExpSquaredKernel(metric=p[1:].copy(), ndim=len(meanx)), mean=meany)

if results['success']:
    print results
else:
    print "didn't work"
    sys.exit()
gp.compute(x, yerr)
#print "pars:",gp.kernel.pars
#print "pars?",gp.kernel[:]
#print np.log(gp.kernel.pars)
#print "varx:",varx
#gp.compute(x, yerr)


# Compute the log likelihood.
print(gp.lnlikelihood(y))

cosind = -1
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
