"""
This is an example of how the boost factor data will be fit.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text',usetex=True,fontsize=24)

boost_filename = "boost_data.dat"
data = np.genfromtxt(boost_filename)
R,boost,err = data.T

lam_pivot = 30.0
z_pivot = 0.5
R_pivot = 0.5 #Mpc physical; NOT Mpc/h and NOT comoving
pivots = [lam_pivot,z_pivot,R_pivot]

def get_model(params,R,lam,z,pivots):
    B0,E_R = params
    lp,zp,Rp = pivots
    return 1.0-B0*(lam/lp)**1.0*((1+z)/(1+zp))**1.0*(R/Rp)**E_R

def lnprior(params):
    if any(params>10.0): return -np.inf
    return 0.0

def lnprob(params,boost,err,R,lam,z,pivots):
    prior = lnprior(params)
    if np.isinf(prior): return -np.inf
    model = get_model(params,R,lam,z,pivots)
    return -0.5*np.sum((boost-model)**2/err**2 + np.log(2*np.pi*err**2))

import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
lnprob_args = (boost,err,R,30.0,0.5,pivots)
result = op.minimize(nll,[-1.0,-1.0],args=lnprob_args,method='Nelder-Mead')
print result

run_mcmc = True
ndim,nwalkers = 2, 32
if run_mcmc:
    import emcee
    pos = [result['x']+1e-3*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=lnprob_args)

    print "running MCMC"
    sampler.run_mcmc(pos,500)
    chain = sampler.chain
    flatchain = chain.reshape((-1,ndim))
    np.savetxt("flatchain.txt",flatchain)


flatchain = np.loadtxt("flatchain.txt")
print flatchain[:,1]
print flatchain.shape
import corner
fig = corner.corner(flatchain,labels=[r'$B_0$',r'$E_R$'])
fig.savefig('corner.png')
plt.show()
plt.clf()

means = np.mean(flatchain,0)
B0,E_R = means
my_boost = get_model([B0,E_R],R,30.0,0.5,pivots)
plt.errorbar(R,boost,err,ls='',marker='o',c='k')
plt.xscale('log')
plt.plot(R,my_boost)
plt.xlabel(r"$R\ {\rm [Mpc]}$")
plt.ylabel(r"$1-f_{\rm cl}$")
plt.subplots_adjust(bottom=0.15)
plt.show()
