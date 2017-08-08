"""
Fit one bin at a time.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
from chainconsumer import ChainConsumer
plt.rc("text", usetex=True)
plt.rc("font", size=14)

fit_blue = True
pname = "nfw_single"
with_scatter = True

#Define the likelihoods
def lnprior(params, pname):
    lnp = 0
    if pname == "nfw_single":
        if len(params)==2: b0,rs = params
        if len(params)==3: b0,rs,sigma = params
        if rs < 0: return -np.inf
        if b0 > 0: return -np.inf
    return lnp

def lnlike(params, lam, z, R, Bp1, Berr, cov, pname):
    X = Bp1-model(params, lam, z, R, pname)
    scatter = scatter_model(params, R, pname)
    cov = np.diag(np.diagonal(cov))
    cov2 = np.diag(scatter)
    icov = np.linalg.inv(cov + cov2)
    CHI2 = np.dot(X, np.dot(icov, X))
    return -0.5*(CHI2+np.log(np.linalg.det(cov+cov2)))

def lnprob(params, lam, z, R, Bp1, Berr, cov, pname):
    lnp = lnprior(params, pname)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, lam, z, R, Bp1, Berr, cov, pname)

#Boost model
def model(params, l, z, R, pname="nfw_single"):
    if pname == "nfw_single":
        if len(params)==2: b0,rs = params
        if len(params)==3: b0,rs,sigma = params
        x = R/rs
        i1 = np.where(x<1)[0]
        i2 = np.where(x>1)[0]
        Fx = np.ones_like(x)
        Fx[i2] *=  np.arctan(np.sqrt(x[i2]**2-1))/np.sqrt(x[i2]**2-1)
        Fx[i1] *= np.arctanh(np.sqrt(1-x[i1]**2))/np.sqrt(1-x[i1]**2)
        #return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d * (1-Fx)/(x**2-1)
        return 1.0 - b0 * (1-Fx)/(x**2-1)

def scatter_model(params, R, pname):
    if pname == "nfw_single":
        if len(params)==2: sigma = 0
        if len(params)==3: sigma = params[-1]
    return (sigma/R)**2 #R is in Mpc, pivot is 1 Mpc


def get_data(zs, full_data=False):
    datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/data_files/blinded_tamas_files/full-mcal-raw_y1clust_l%d_z%d_pz_boost.dat"
    if fit_blue:
        datapath = "/home/tom/Desktop/boost_files/bluecurves/blue_z%d_l%d.txt"
        covpath  = "/home/tom/Desktop/boost_files/bluecurves/cov_z%d_l%d.txt"
    else:
        datapath = "/home/tom/Desktop/boost_files/redcurves/red_z%d_l%d.txt"
        covpath  = "/home/tom/Desktop/boost_files/redcurves/cov_z%d_l%d.txt"
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

def bestfit(Bp1, Berr, lams, zs, R, cov):
    if pname is "nfw_single":guess= [-0.2, 1.0] #B0, Rs
    if with_scatter: guess.append(1e-1)
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, guess, args=(lams, zs, R, Bp1, Berr, cov, pname))
#                         method='Nelder-Mead')
    return result['x']

def plot_all(results, zs, lams):
    fig, axarr = plt.subplots(len(zs), len(zs[0]), sharex=True, sharey = True)
    for i in range(len(zs)):
        for j in range(len(zs[0])):
            Bmodel = model(results[i][j], lams[i,j], zs[i,j], R[i][j], pname)
            if fit_blue: color='b'
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

def plot_resid(results, zs, lams):
    fig, axarr = plt.subplots(len(zs), len(zs[0]), sharex=True, sharey = True)
    for i in range(len(zs)):
        for j in range(len(zs[0])):
            Bmodel = model(results[i][j], lams[i,j], zs[i,j], R[i][j], pname)
            if fit_blue: color='b'
            else: color='r'
            axarr[i,j].errorbar(R[i][j], Bp1[i][j]-Bmodel, Berr[i][j], color=color)
            axarr[i,j].fill_between([0.03,0.2],[1.8,1.8], color="lightgray", alpha=0.7, zorder=-1)
            axarr[i,j].set_xscale('log')
            axarr[i,j].set_xticks([0.1,1.0,10])
            axarr[i,j].set_xlim(0.2, 40)
            axarr[i,j].set_ylim(-0.12, 0.12)
            axarr[i,j].grid(ls=':')
            axarr[i,j].axhline(0.015, c='k', ls='--', zorder=-1)
            axarr[i,j].axhline(-0.02, c='k', ls='--', zorder=-1)
    axarr[1,0].set_ylabel("$\Delta(1+B)$")#r"$R\ [{\rm Mpc}]$")
    axarr[2,3].set_xlabel(r"$R\ [{\rm Mpc}]$")
    fig.set_size_inches(10, 5)
    plt.subplots_adjust(hspace=0.01, wspace=0.01, left=0.15, bottom=0.15)
    plt.show()
    return



if __name__ == "__main__":
    zs = np.loadtxt("/home/tom/Desktop/boost_files/Y1_meanz.txt")
    lams = np.loadtxt("/home/tom/Desktop/boost_files/Y1_meanl.txt")

    Bp1, Berr, R, cov = get_data(zs)
    results = []
    for i in range(len(zs)):
        resi = []
        for j in range(len(zs[0])):
            """
            if i is not 2:
                resi.append([0,1,0])
                continue
            if j is not 5:
                resi.append([0,1,0])
                continue
            """
            res = bestfit(Bp1[i][j], Berr[i][j], lams[i,j], zs[i,j], R[i][j], cov[i][j])
            #print res
            #print lnlike(res, lams[i][j], zs[i][j], R[i][j], Bp1[i][j], Berr[i][j], cov[i][j], pname)
            #res[0] -=0.02
            #print lnlike(res, lams[i][j], zs[i][j], R[i][j], Bp1[i][j], Berr[i][j], cov[i][j], pname)
            resi.append(res)
                

        print "Best fit z%d done"%(i)
        results.append(resi)
    Bp1, Berr, R, cov = get_data(zs)#, True)
    #plot_all(results, zs, lams)
    plot_resid(results, zs, lams)
