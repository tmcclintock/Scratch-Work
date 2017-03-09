"""
Make a scatter plot of z vs sigma_z
with color coding by richness bins.
"""
import fitsio, sys, os
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rc("text",usetex=True,fontsize=24)

fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
lam_trues = data['LAMBDA_CHISQ']
z_trues = data['Z_LAMBDA']
lam_best = np.loadtxt("lam_best_all.txt")
z_best = np.loadtxt("z_best_all.txt")
sigma_z_all = np.loadtxt("sigma_z_all.txt")
sigma_z_all[sigma_z_all > 1.0] = -1.0

lo = (lam_trues >= 20) * (lam_trues < 30)
mid = (lam_trues >= 30) * (lam_trues < 60)
hi = lam_trues >= 60

use_best = False

if use_best:
    plt.scatter(z_best[lo],sigma_z_all[lo],c='g',label=r'$\lambda\in(20,30)$',alpha=0.9)
    plt.scatter(z_best[mid],sigma_z_all[mid],c='r',label=r'$\lambda\in(30,60)$',alpha=0.4)
    plt.scatter(z_best[hi],sigma_z_all[hi],c='b',label=r'$\lambda>60$',alpha=0.3)
    plt.xlabel(r"$z_{\rm best}$",fontsize=24)
else:
    plt.scatter(z_trues[lo],sigma_z_all[lo],c='g',label=r'$\lambda\in(20,30)$',alpha=0.9)
    plt.scatter(z_trues[mid],sigma_z_all[mid],c='r',label=r'$\lambda\in(30,60)$',alpha=0.4)
    plt.scatter(z_trues[hi],sigma_z_all[hi],c='b',label=r'$\lambda>60$',alpha=0.3)
    plt.xlabel(r"$z_{\rm true}$",fontsize=24)
plt.ylim(0.0,0.12)#max(sigma_z_all)*1.1)
plt.xlim(0.07,0.4)
plt.legend(loc='lower right',fontsize=16)
plt.ylabel(r"$\sigma_z$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
plt.clf()

pdz = np.fabs(z_best - z_trues)/z_trues
plt.scatter(z_trues[lo],pdz[lo],c='g',label=r'$\lambda\in(20,30)$',alpha=0.9)
plt.scatter(z_trues[mid],pdz[mid],c='r',label=r'$\lambda\in(30,60)$',alpha=0.4)
plt.scatter(z_trues[hi],pdz[hi],c='b',label=r'$\lambda>60$',alpha=0.3)
plt.xlabel(r"$z_{\rm true}$",fontsize=24)
#plt.legend(loc='upper left')
plt.ylabel(r"$|z_{\rm true}-z_{\rm best}|/z_{\rm true}$",fontsize=24)
plt.show()
plt.clf()

pdlam = np.fabs(lam_best - lam_trues)/lam_trues
plt.scatter(lam_trues[lo],pdlam[lo],c='g',label=r'$\lambda\in(20,30)$',alpha=0.9)
plt.scatter(lam_trues[mid],pdlam[mid],c='r',label=r'$\lambda\in(30,60)$',alpha=0.4)
plt.scatter(lam_trues[hi],pdlam[hi],c='b',label=r'$\lambda>60$',alpha=0.3)
plt.xlabel(r"$\lambda_{\rm true}$",fontsize=24)
#plt.legend(loc='upper left')
plt.ylabel(r"$|\lambda_{\rm true}-\lambda_{\rm best}|/\lambda_{\rm true}$",fontsize=24)
plt.show()
plt.clf()
