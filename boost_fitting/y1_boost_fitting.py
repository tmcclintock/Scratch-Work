"""
Fit the Y1 boost factors given by tamas.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
from chainconsumer import ChainConsumer
plt.rc("text", usetex=True)
plt.rc("font", size=14)

nwalkers = 16
nsteps = 5000

#Boost model
def model(params, l, z, R):
    if len(params)==4:
        b0,c,d,e = params
    elif len(params)==5:
        b0,c,d,e,sigma = params
    return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(R/0.5)**e

def scatter_model(sigma, R):
    return (sigma/R)**2 #R is in Mpc, pivot is 1 Mpc

#Define the likelihoods
def lnprior(params):
    if len(params)==4:
        b0,c,d,e = params
    elif len(params)==5:
        b0,c,d,e,sigma = params
    #if any(np.fabs(params)>20.0): return -np.inf
    return 0.0

def lnlike(params, lams, zs, R, Bp1, Berr):
    if len(params)==4:
        b0,c,d,e = params
        sigma=0
    elif len(params)==5:
        b0,c,d,e,sigma = params
    LL = 0
    for i in range(len(Bp1)): #Loop over all boost files
        for j in xrange(0,len(Bp1[i])):
            Bmodel = model(params, lams[i,j], zs[i,j], R[i][j])
            scatter = scatter_model(sigma, R[i][j])
            inds = R[i][j] > 0.2 #Mpc; small-scale cutoff
            CHI2 = (Bp1[i][j]-Bmodel)**2/(Berr[i][j]**2+scatter)
            CHI2 = CHI2[inds]
            OTHER_TERM = np.log(Berr[i][j]**2+scatter)
            OTHER_TERM = OTHER_TERM[inds]
            LL += -0.5*CHI2.sum()
            LL += -0.5*OTHER_TERM.sum()
            #LL += np.sum(-0.5*(Bp1[i][j]-Bmodel)**2/(Berr[i][j]**2+scatter))
            #LL += np.sum(-0.5*np.log(Berr[i][j]**2+scatter))
    return LL

def lnprob(params, lams, zs, R, Bp1, Berr):
    lnp = lnprior(params)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, lams, zs, R, Bp1, Berr)

use_blue = True
def get_data(zs):
    datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/data_files/blinded_tamas_files/full-mcal-raw_y1clust_l%d_z%d_pz_boost.dat"
    if use_blue:
        datapath = "/home/tmcclintock/Desktop/boost_files/bluecurves/blue_z%d_l%d.txt"
        covpath  = "/home/tmcclintock/Desktop/boost_files/bluecurves/cov_z%d_l%d.txt"
    else:
        datapath = "/home/tmcclintock/Desktop/boost_files/redcurves/red_z%d_l%d.txt"
        covpath  = "/home/tmcclintock/Desktop/boost_files/redcurves/cov_z%d_l%d.txt"

    #Read in all data
    Bp1  = []
    Berr = []
    R    = []
    for i in range(len(zs)):
        Bp1i  = []
        Berri = []
        Ri    = []
        for j in xrange(0,len(zs[i])):
            Rij, Bp1ij, Berrij = np.loadtxt(datapath%(i, j), unpack=True)
            Bp1ij  = Bp1ij[Berrij > 1e-8]
            Rij    = Rij[Berrij > 1e-8]
            Berrij = Berrij[Berrij > 1e-8]
            Bp1i.append(Bp1ij)
            Berri.append(Berrij)
            Ri.append(Rij)
        Bp1.append(Bp1i)
        Berr.append(Berri)
        R.append(Ri)
    return Bp1, Berr, R

def bestfit(Bp1, Berr, lams, zs, R):
    #Parameters: B0, C, D, E, sigma(R=1 MPC)
    guess = [-1.0, 1.0, 2.0, -1.0, -1.0]
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, guess, args=(lams, zs, R, Bp1, Berr), 
                         method='Nelder-Mead')
    return result

def do_mcmc(bf, Bp1, Berr, lams, zs, R):
    ndim = len(bf)
    pos = [bf + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(lams ,zs, R, Bp1, Berr))
    print "Starting MCMC"
    sampler.run_mcmc(pos, nsteps)
    chain = sampler.flatchain
    print chain.shape
    np.savetxt("chain.txt", chain)
    print "Chain saved"
    return

def plot_BF(bf, zs, lams):
    Nz = len(zs)
    Nl = len(zs[0])
    fig, axarr = plt.subplots(Nz, Nl, sharex=True, sharey = True)
    for i in range(Nz):
        for j in range(Nl):
            Bmodel = model(bf, lams[i,j], zs[i,j], R[i][j])
            if use_blue: color='b'
            else: color='r'
            axarr[i,j].fill_between(R[i][j], Bp1[i][j]-Berr[i][j], Bp1[i][j]+Berr[i][j], color=color)
            axarr[i,j].plot(R[i][j], Bmodel, c='k')

            axarr[i,j].fill_between([0.03,0.2],[1.8,1.8], color="lightgray", alpha=0.7, zorder=-1)
            axarr[i,j].set_xscale('log')
            axarr[i,j].set_xticks([0.1,1.0,10])
            axarr[i,j].set_yticks([1.0, 1.2, 1.4, 1.6, 1.8])
            axarr[i,j].set_xlim(0.03, 40)
            axarr[i,j].set_ylim(0.8, 1.8)
            axarr[i,j].grid(ls=':')
    axarr[1,0].set_ylabel("Boost Factor")#r"$R\ [{\rm Mpc}]$")
    axarr[2,3].set_xlabel(r"$R\ [{\rm Mpc}]$")
    fig.set_size_inches(10, 5)
    plt.subplots_adjust(hspace=0.01, wspace=0.01, left=0.15, bottom=0.15)
    plt.show()
    return

def see_chain():
    fullchain = np.loadtxt("chain.txt")
    nburn = 1000
    chain = fullchain[nwalkers*nburn:]
    if len(fullchain[0])==4:
        labels = [r"$B_0$", r"$C_\lambda$", r"$D_z$", r"$E_R$"]
    elif len(fullchain[0])==5:
        labels = [r"$B_0$", r"$C_\lambda$", r"$D_z$", r"$E_R$",r"$\sigma_{1\ {\rm Mpc}}$"]
    cc = ChainConsumer().add_chain(chain, parameters=labels, name="boost")
    cc.plot()
    plt.subplots_adjust(bottom=0.2, top=0.9, left=0.25, right = 0.85)
    plt.show()

if __name__ == "__main__":
    zs = np.loadtxt("/home/tmcclintock/Desktop/des_wl_work/Y1_work/data_files/Y1_meanz.txt")
    lams = np.loadtxt("/home/tmcclintock/Desktop/des_wl_work/Y1_work/data_files/Y1_meanl.txt")

    Bp1, Berr, R = get_data(zs)
    res = bestfit(Bp1, Berr, lams, zs, R)
    plot_BF(res['x'], zs, lams)
    #do_mcmc(res['x'], Bp1, Berr, lams, zs, R)
    #see_chain()
