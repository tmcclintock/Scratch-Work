"""
Let's do a single cluster. First we plot it.
"""
import fitsio
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text",usetex=True,fontsize=24)

fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
lam_trues = data['LAMBDA_CHISQ']
z_trues = data['Z_LAMBDA']
z_bests = np.loadtxt("z_best_all.txt")

names = data.dtype.names
lam_arrays = data['LAMBDA_CHISQS']
lam_arrays = data['LAMBDA_CHISQS']
print lam_arrays.shape
zname = 'redshift_list.txt'
zs = np.genfromtxt(zname).flatten()

def see_comparison(index):
    lam = lam_arrays[index]
    lam[lam<0.0] = 0.0
    lam_true = lam_trues[index]
    z_true = z_trues[index]
    z_best = z_bests[index]
    sigma_z = 0.0291
    plt.plot(zs,lam)
    plt.scatter(z_true,lam_true)
    plt.plot([z_true,z_true],[0,lam_true])
    def lam_model(z,sigmaz,z_fit,lam):
        return lam*np.exp(-0.5*(z_fit-z)**2/sigmaz**2)
    plt.plot(zs,lam_model(zs,sigma_z,z_best,lam_true))
    plt.xlim(0,0.45)
    plt.ylim(-10,250)
    plt.ylabel(r"Richness $\lambda$")
    plt.xlabel(r"Redshift $z$")
    plt.subplots_adjust(bottom=0.15)
    plt.show()

see_comparison(3901)
