import numpy as np
import matplotlib.pyplot as plt
import fastcorr

Pklin = np.loadtxt("pklin.txt")
k = np.loadtxt("k.txt")

xidata = np.loadtxt("final_hmcf.txt")
print xidata.shape
R = xidata[:,0]
xi = xidata[:,1]
err = xidata[:,2]

xi_mm = fastcorr.calc_corr(R,k,Pklin)
bias = xi/xi_mm
berr = err/xi_mm

plt.loglog(R,xi)
plt.loglog(R,xi_mm)
#plt.show()
plt.clf()

plt.errorbar(R,bias,berr)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r"$\xi_{\rm hm}/\xi_{\rm mm}$",fontsize=24)
plt.xlabel(r"${\rm R}\ [{\rm Mpc/h}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
