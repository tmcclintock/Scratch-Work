import numpy as np
import matplotlib
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
    return cov

def get_mcmlmis(zi, lj):
    cov = np.loadtxt("all_covs/mcmlmis/tom_covariance_z%d_l%d.txt"%(zi, lj))
    return cov

def get_shape(zi, lj):
    cov = np.loadtxt("all_covs/shape_noise/full-unblind-mcal-zmix_y1subtr_l%d_z%d_n300_v1_shapecov.dat"%(lj, zi))
    return cov

def get_uLSS(zi, lj):
    cov = np.loadtxt("all_covs/uLSS/full-unblind-mcal-zmix_y1subtr_l%d_z%d_lss_uncorr_y1mask_dst_cov.dat"%(lj, zi))
    return cov

def get_jk(zi, lj):
    cov = np.loadtxt("all_covs/jks/full-unblind-mcal-zmix_y1subtr_l%d_z%d_dst_cov.dat"%(lj, zi))
    return cov

if __name__ == "__main__":
    zi, lj = 2, 6
    clss = get_cLSS(zi, lj)
    ml = get_mcmlmis(zi, lj)
    shape = get_shape(zi, lj)
    ulss = get_uLSS(zi, lj)
    total = np.sqrt(clss**2 + ml**2 + shape**2 + ulss**2)
    jk = get_jk(zi, lj)

    diff = jk - total

    fig, axes = plt.subplots(2, 4)#, sharex=True, sharey=True)


    extent = np.log10([0.0323, 30, 0.0323, 30])
    def plot_cov(cov, i, j, name="test"):
        ax = axes[i,j]
        ax.imshow(np.abs(cov[2:, 2:]), aspect='equal', interpolation='none', extent=extent, origin='lower', cmap="plasma", vmin=0.01, vmax=400, norm=matplotlib.colors.LogNorm())
        ax.set_xlabel(r"$R$ [Mpc]", fontsize=8)
        ax.set_ylabel(r"$R$ [Mpc]", fontsize=8)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(name)

    plot_cov(shape, 0, 0, "Shape")
    plot_cov(clss, 0, 1, "cLSS+ell")
    plot_cov(ml, 1, 0, r"M-c+M-$\lambda$+Mis")
    plot_cov(ulss, 1, 1, "uLSS")
    plot_cov(jk, 0, 2, "JK")
    plot_cov(total, 1, 2, "Total")
    plot_cov(jk-total, 0, 3, "JK-Total")
    fig.delaxes(axes[1,3])

    #axes[1].set_xlabel(r"$R\ [{\rm Mpc}]$")
    #axes[0].set_title("z%dl%d"%(zi, lj))
    plt.subplots_adjust(hspace=0)
    plt.show()
