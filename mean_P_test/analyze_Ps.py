"""
Now to analyze the power spectra I have.
"""
import numpy as np
import matplotlib.pyplot as plt

#Read in the clusters
cpath = "cluster_files/clusters_z0l4.txt"
zs, lams = np.genfromtxt(cpath, unpack=True)
zmean = np.mean(zs)
zerr = np.std(zs)
print "Working with %d clusters from the file %s"%(len(zs),cpath)
print "Mean redshift = %.3f +- %.3f"%(zmean,zerr)

meanzPnlpath  = "PK_files/Pmean_nl.txt"
meanzPlinpath = "PK_files/Pmean_lin.txt"
Pnlpath  = "PK_files/Pnl_%d.txt"
Plinpath = "PK_files/Plin_%d.txt"

k = np.loadtxt("PK_files/k.txt")
Pmeanz = np.loadtxt(meanzPnlpath)
Plinmeanz = np.loadtxt(meanzPlinpath)
Pnls = []
Plins = []
for i in range(len(zs)):
    Pnls.append(np.loadtxt(Pnlpath%i))
    Plins.append(np.loadtxt(Plinpath%i))
Pnls = np.array(Pnls)
Plins = np.array(Plins)

Pexp = np.mean(Pnls,0)
Plinexp = np.mean(Plins,0)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].loglog(k, Pmeanz)#, label=r"$P_{\rm nl}(k,\bar{z})$")
axarr[0].loglog(k, Pexp)#, label=r"$\bar{P}_{\rm nl}$")
axarr[0].loglog(k, Plinmeanz)#, label=r"$P_{\rm lin}(k,\bar{z})$")
axarr[0].loglog(k, Plinexp)#, label=r"$\bar{P}_{\rm lin}$")
axarr[1].plot(k, Pmeanz/Pexp, label=r"NL")
axarr[1].plot(k, Plinmeanz/Plinexp, label=r"LIN")
plt.legend(loc=0)
axarr[1].set_xlabel(r"$k\ [h/{\rm Mpc}]$", fontsize=24)
axarr[1].set_ylabel(r"$P(\bar{z})/\bar{P}$", fontsize=24)
axarr[0].set_ylabel(r"$P\ [h/{\rm Mpc}]^3$", fontsize=24)
plt.subplots_adjust(bottom=0.15, left=0.15)
plt.show()
