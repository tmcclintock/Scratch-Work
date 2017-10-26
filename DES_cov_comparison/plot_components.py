import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.rc("text", usetex=True)
plt.rc("font", size=18)
plt.rc("font", family="serif")

Nbins = 15
binmin = 0.0323 #Mpc phys
binmax = 30.0 #Mpc phys
Redges = np.logspace(np.log(binmin), np.log(binmax), num=Nbins+1, base=np.e)
R = (Redges[:-1]+Redges[1:])/2.

def get_cLSS(zi, lj):
    cov = np.loadtxt("all_covs/cLSSell/cLSSell_z%d_l%d.txt"%(zi, lj))
    return np.sqrt(np.diag(cov))

def get_mcmlmis(zi, lj):
    cov = np.loadtxt("all_covs/mcmlmis/tom_covariance_z%d_l%d.txt"%(zi, lj))
    return np.sqrt(np.diag(cov))

def get_shape(zi, lj):
    cov = np.loadtxt("all_covs/shape_noise/full-unblind-mcal-zmix_y1subtr_l%d_z%d_n300_v1_shapecov.dat"%(lj, zi))
    return np.sqrt(np.diag(cov))

def get_uLSS(zi, lj):
    cov = np.loadtxt("all_covs/uLSS/full-unblind-mcal-zmix_y1subtr_l%d_z%d_lss_uncorr_y1mask_dst_cov.dat"%(lj, zi))
    return np.sqrt(np.diag(cov))

def get_jk(zi, lj):
    cov = np.loadtxt("all_covs/jks/full-unblind-mcal-zmix_y1subtr_l%d_z%d_dst_cov.dat"%(lj, zi))
    return np.sqrt(np.diag(cov))

if __name__ == "__main__":
    zi, lj = 0, 2
    clss = get_cLSS(zi, lj)
    ml = get_mcmlmis(zi, lj)
    shape = get_shape(zi, lj)
    ulss = get_uLSS(zi, lj)
    total = np.sqrt(clss**2 + ml**2 + shape**2 + ulss**2)

    jk = get_jk(zi, lj)
    pd = (jk - total)/total

    gs = gridspec.GridSpec(3, 3)
    axes = [plt.subplot(gs[0:2, :]), plt.subplot(gs[2, :])]
    #fig, axes = plt.subplots(2, sharex=True)
    
    axes[0].loglog(R, clss, label=r"cLSS+ell")
    axes[0].loglog(R, ml, label=r"M-c+M-$\lambda$+Mis")
    axes[0].loglog(R, shape, label=r"Shape")
    axes[0].loglog(R, ulss, label=r"uLSS")
    axes[0].loglog(R, total, label=r"Total")
    axes[0].loglog(R, jk, ls = '', marker='.', label=r"JK")
    axes[0].legend(loc=0, fontsize=10, frameon=False)

    axes[1].plot(R, pd, c='k', ls='-')
    axes[1].axhline(0, c='k', ls='--')

    axes[1].set_xlabel(r"$R\ [{\rm Mpc}]$")
    axes[1].set_ylabel(r"(JK-Total)/Total")
    axes[0].set_ylabel(r"$\delta\Delta\Sigma\ [{\rm M_\odot/pc^2}]$")
    axes[0].set_title("z%dl%d"%(zi, lj))
    plt.subplots_adjust(bottom=0.15, left=0.18, hspace=0)
    plt.show()
