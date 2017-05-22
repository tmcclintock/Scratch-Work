"""
Fit the Y1 boost factors given by tamas.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
plt.rc("text", usetex=True, fontsize=20)

#Boost model
def model(params, l, z, R):
    b0,c,d,e = params
    return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(R/0.5)**e

#Define the likelihoods
def lnprior(params):
    b0,c,d,e = params
    #if any(np.fabs(params)>20.0): return -np.inf
    return 0.0

def lnlike(params, lams, zs, R, Bp1, Berr):
    b0,c,d,e = params
    LL = 0
    for i in range(len(Bp1)): #Loop over all boost files
        for j in xrange(0,len(Bp1[i])):
            Bmodel = model(params, lams[i,j], zs[i,j], R[i][j])
            LL += np.sum(-0.5*(Bp1[i][j]-Bmodel)**2/Berr[i][j]**2)
    return LL

def lnprob(params, lams, zs, R, Bp1, Berr):
    lnp = lnprior(params)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, lams, zs, R, Bp1, Berr)

def get_data(zs):
    datapath = "/home/tmcclintock/Desktop/des_wl_work/Y1_work/data_files/tamas_files/full-mcal-raw_y1clust_l%d_z%d_pz_boost.dat"
    #Read in all data
    Bp1  = []
    Berr = []
    R    = []
    for i in range(len(zs)):
        Bp1i  = []
        Berri = []
        Ri    = []
        for j in xrange(0,len(zs[i])):
            Rij, Bp1ij, Berrij = np.loadtxt(datapath%(j, i), unpack=True)
            Bp1ij  = Bp1ij[Berrij > 1e-3]
            Rij    = Rij[Berrij > 1e-3]
            Berrij = Berrij[Berrij > 1e-3]
            Bp1i.append(Bp1ij)
            Berri.append(Berrij)
            Ri.append(Rij)
        Bp1.append(Bp1i)
        Berr.append(Berri)
        R.append(Ri)
    return Bp1, Berr, R

def bestfit(Bp1, Berr, lams, zs, R):
    guess = [-1.0, 1.0, 1.0, -1.0]
    nll = lambda *args: -lnprob(*args)
    result = op.minimize(nll, guess, args=(lams, zs, R, Bp1, Berr), 
                         method='Nelder-Mead')
    return result

if __name__ == "__main__":
    zs = np.loadtxt("/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/meanz.txt")
    lams = np.loadtxt("/home/tmcclintock/Desktop/des_wl_work/Y1_work/blinded_data/meanl.txt")

    Bp1, Berr, R = get_data(zs)
    res = bestfit(Bp1, Berr, lams, zs, R)
    print res
    params = res['x']
    print params
    c = np.linspace(1.0, 0.3, len(zs[0]))
    cmaps = ["Blues", "Greens", "Reds"]
    f, axarr = plt.subplots(len(zs), sharex=True)
    for i in range(len(zs)):
        for j in range(len(zs[i])):
            Bmodel = model(params, lams[i,j], zs[i,j], R[i][j])
            col = plt.get_cmap(cmaps[i])(c[j])
            axarr[i].errorbar(R[i][j], Bp1[i][j], Berr[i][j], c=col, ls='', marker='o')
            axarr[i].plot(R[i][j], Bmodel, c=col, label="l%d"%j)
            axarr[i].set_ylabel(r"$B(z=%.2f)$"%zs[i,0], fontsize=14)
    plt.subplots_adjust(hspace=0.05, left=0.15, bottom=0.15)
    plt.xlabel(r"$R\ [{\rm Mpc}]$")
    plt.xscale('log')
    #plt.legend(loc=0)
    plt.show()
