"""
Fit the Y1 boost factors given by tamas.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
from chainconsumer import ChainConsumer
from models import *
from likelihoods import *
plt.rc("text", usetex=True)
plt.rc("font", size=14)

nwalkers = 16
nsteps = 5000

use_blue = False
def get_data(zs, full_data=False):
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
    cov  = []
    for i in range(len(zs)):
        Bp1i  = []
        Berri = []
        Ri    = []
        covi  = []
        for j in xrange(0,len(zs[i])):
            Rij, Bp1ij, Berrij = np.loadtxt(datapath%(i, j), unpack=True)
            if full_data: cut = (Berrij > 1e-8)
            else: cut = (Berrij > 1e-8)*(Rij > 0.2)
            Bp1ij  = Bp1ij[cut]
            Rij    = Rij[cut]
            Berrij = Berrij[cut]
            covij  = np.loadtxt(covpath%(i,j))
            covij  = covij[cut]
            covij  = covij[:,cut]
            Bp1i.append(Bp1ij)
            Berri.append(Berrij)
            Ri.append(Rij)
            covi.append(covij)
        Bp1.append(Bp1i)
        Berr.append(Berri)
        R.append(Ri)
        cov.append(covi)
    return Bp1, Berr, R, cov

model_name = "bpl"
with_scatter = True

def bestfit(Bp1, Berr, lams, zs, R, cov):
    #Params:   B0,   C,   D,    E, sigma(R=1 MPC)
    if model_name is "simple": guess = [-0.2, 0.2, -2.0, -1.0]# B0, CL, Dz, Er
    elif model_name is "r1r2": guess = [-0.2, 0.2, -2.0, 1.0, 1.0]# B0, CL, Dz, E1r, E2r
    elif model_name is "bpl" : guess = [-0.2, 0.2, -2.0, -1.0, -1.0, 10.]# B0, CL, Dz, E1r, E2r, break
    if with_scatter: guess.append(1e-1)
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, guess, args=(lams, zs, R, Bp1, Berr, cov, model_name), 
                         method='Nelder-Mead')
    print result
    return result

def do_mcmc(bf, Bp1, Berr, lams, zs, R, cov):
    ndim = len(bf)
    pos = [bf + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(lams ,zs, R, Bp1, Berr, cov))
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
            Bmodel = model(bf, lams[i,j], zs[i,j], R[i][j], model_name)
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

    Bp1, Berr, R, cov = get_data(zs)
    res = bestfit(Bp1, Berr, lams, zs, R, cov)
    Bp1, Berr, R, cov = get_data(zs, True)
    plot_BF(res['x'], zs, lams)
    #do_mcmc(res['x'], Bp1, Berr, lams, zs, R, cov)
    #see_chain()
