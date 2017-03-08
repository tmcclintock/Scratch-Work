"""
This does the boost factor fitting to all the data.

This is the file that will be turned into a cosmosis module.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text',usetex=True,fontsize=24)
labels = [r'$B_0$',r'$C_\lambda$',r'$D_z$',r'$E_R$']
zlabels = [r"$z\in(0.2,0.35)$",r"$z\in(0.35,0.5)$",r"$z\in(0.5,0.65)$"]

#Input filenames
dataname = "real_data/y1clust_l%d_z%d_corr_boost.dat"
covname  = "real_data/y1clust_l%d_z%d_corr_boost_cov.dat"

#Richnesses and redshifts
lams = np.loadtxt("real_data/mean_richnesses.txt")
zs = np.loadtxt("real_data/mean_redshifts.txt")

#Prepare the data vectors
boost_array = []
err_array = []
cov_array = []
icov_array = []
det_array = []
R_array = []
N_z, N_lam = 3,7
for i in xrange(0,N_z):
    for j in xrange(0,N_lam):
        lo = 3
        R,boost,err = np.loadtxt(dataname%(j,i)).T
        cov = np.loadtxt(covname%(j,i))
        R = R[lo:]
        boost = boost[lo:]
        err = err[lo:]
        cov = cov[lo:]
        cov = cov[:,lo:]
        boost_array.append(boost)
        err_array.append(err)
        R_array.append(R)
        cov_array.append(cov)
        icov_array.append(np.linalg.inv(cov))
        det_array.append(np.linalg.det(cov))

#Pivot values
lam_pivot = 30.0
z_pivot = 0.5
R_pivot = 0.5 #Mpc physical; NOT Mpc/h and NOT comoving
pivots = [lam_pivot,z_pivot,R_pivot]

lnprob_args = (boost_array,icov_array,R_array,lams,zs,pivots)

def get_model(params,R,lam,z,pivots):
    B0,C_L,D_Z,E_R = params
    lp,zp,Rp = pivots
    return 1.0-B0*(lam/lp)**C_L*((1+z)/(1+zp))**D_Z*(R/Rp)**E_R

def lnprior(params):
    if any(np.fabs(params)>15.0): return -np.inf
    if params[3] > 0.0: return -np.inf
    if params[0] > 0.0: return -np.inf
    return 0.0

def lnprob(params,boost_array,icov_array,R_array,lams,zs,pivots):
    prior = lnprior(params)
    if np.isinf(prior): return -np.inf
    llike = 0.0
    for i in xrange(0,N_z):
        for j in xrange(0,N_lam):
            lam = lams[i,j]
            z = zs[i,j]
            index = i*N_lam + j
            boost = boost_array[index]
            R = R_array[index]
            icov = icov_array[index]
            model = get_model(params,R,lam,z,pivots)
            X = boost-model
            llike += -0.5*np.dot(X,np.dot(icov,X))
    return llike + prior

#Flow control
do_single_test = True
do_maximization = True
do_mcmc = True
see_corner = True
see_fits = True

if do_single_test:
    print "Single test"
    params = [-0.14,-1.0,-1.0,-0.7]
    print lnprob(params,boost_array,icov_array,R_array,lams,zs,pivots)
            
if do_maximization:
    print "Finding most likely values"
    params = [-0.1,-0.5,-0.5,-0.5]
    import scipy.optimize as op
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll,params,args=lnprob_args,method='Powell')
    print result
    np.savetxt("results/maximization.txt",result['x'])

best_params = np.loadtxt("results/maximization.txt")
Ndim = 4
Nwalkers = 32
Nburn,Nsteps = 1500,3000
if do_mcmc:
    print "Performing MCMC on boost factors"
    import emcee
    pos = [best_params+1e-4*np.random.randn(Ndim) for i in range(Nwalkers)]
    sampler = emcee.EnsembleSampler(Nwalkers,Ndim,lnprob,args=lnprob_args)

    sampler.run_mcmc(pos,Nsteps)
    chain = sampler.chain
    flatchain = chain.reshape((-1,Ndim))
    #np.savetxt("results/chain.txt",chain)
    np.savetxt("results/flatchain.txt",flatchain)

if see_corner:
    print "Analyzing the chain"
    flatchain = np.loadtxt("results/flatchain.txt")
    print flatchain.shape
    usechain = flatchain[Nburn*Nwalkers:]
    import corner
    fig = corner.corner(usechain,labels=labels)
    plt.show()

if see_fits:
    print "Making fits"
    flatchain = np.loadtxt("results/flatchain.txt")
    usechain = flatchain[Nburn*Nwalkers:]
    #means = np.mean(usechain,0)
    means = best_params
    #means = [-0.04,1.0,0.0,-0.3]
    print means
    cmap = plt.get_cmap('jet')
    colors = [cmap(float(N_lam-i)/(N_lam-1)) for i in range(0,N_lam)]
    for i in xrange(0,N_z):
        for j in xrange(0,N_lam):
            index = i*N_lam + j
            R = R_array[index]
            plt.errorbar(R,boost_array[index],err_array[index],c=colors[j],marker='.',alpha=0.3)
            model = get_model(means,R,lams[i,j],zs[i,j],pivots)
            plt.plot(R,model,c=colors[j])
        plt.xscale('log')
        plt.xlabel(r"$R\ {\rm[Mpc]}$")
        plt.ylabel(r"$1-f_{\rm cl}$")
        plt.title(zlabels[i])
        plt.subplots_adjust(wspace=0.001,bottom=0.15)
        plt.ylim(0.9,1.6)
        plt.xlim(0.1,100)
        plt.show()
