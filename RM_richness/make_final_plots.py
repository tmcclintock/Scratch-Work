"""
Make the final plots, which are 
histograms if sigma(z) inside z bins
as well as sigma(z) vs z.
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

N_bins = 10
edges = np.linspace(min(z_trues),max(z_trues)+0.00001,N_bins+1)

make_cluster_files = True
cpath = "cluster_files/clusters_%d.txt"
if make_cluster_files:
    cluster_files = []
    for i in range(N_bins):
        cluster_files.append(open(cpath%i,"w"))
        cluster_files[i].write("#zmin=%e\tzmax=%e\n"%(edges[i],edges[i+1]))
        cluster_files[i].write("#z lambda sigma_z\n")

    for i in xrange(0,len(z_best)):
        for j in range(N_bins):
            if edges[j] <= z_trues[i] < edges[j+1]: 
                cluster_files[j].write("%e %e %e\n"%(z_trues[i],lam_trues[i],sigma_z_all[i]))
                continue
        continue

    for i in range(N_bins):
        cluster_files[i].close()
#end if make_cluster_files

"""#sigma_z_all = sigma_z_all[sigma_z_all < 1.0]
n, bins, patches = plt.hist(sigma_z_all,50)
plt.ylabel("Number",fontsize=24)
plt.xlabel(r"$\sigma_z$",fontsize=24)
plt.title("All clusters")
#plt.show()
plt.clf()"""

lo = (lam_trues >= 20) * (lam_trues < 30)
mid = (lam_trues >= 30) * (lam_trues < 60)
hi = lam_trues >= 60
labels = [r"$\lambda\in(20,30)$",r"$\lambda\in(30,60)$",r"$\lambda>60$"]
cs = ['b','g','r']

zmeans = np.zeros((3,N_bins))
sigmazmeans = np.zeros((3,N_bins))
szerrs = np.zeros((3,N_bins))
for i,g in zip(xrange(0,3),[lo,mid,hi]):
    z = z_trues[g]
    sig = sigma_z_all[g]
    z = z[sig < 0.12]
    sig = sig[sig < 0.12]
    zmin,zmax = 0.08,0.35#min(z),max(z)
    dz = (zmax - zmin)/N_bins
    for j in range(N_bins):
        inds = (z > j*dz+zmin) * (z < (j+1)*dz+zmin)
        zmeans[i,j] = np.median(z[inds])
        sigmazmeans[i,j] = np.median(sig[inds])
        szerrs[i,j] = np.std(sig[inds])
        #print zmeans[i,j],sigmazmeans[i,j]
    plt.plot(zmeans[i],sigmazmeans[i],marker='o',label=labels[i],c=cs[i])
    plt.scatter(z,sig,alpha=0.05,marker='.',c=cs[i])
plt.ylabel(r"$\sigma_z$",fontsize=24)
plt.xlabel(r"$z_{\rm true}$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.legend(loc="lower right",fontsize=16)
plt.ylim(0,0.12)
plt.xlim(0.07,0.4)
plt.show()
plt.clf()

#Final curve with no richness binning
zmeans = np.zeros(N_bins)
sigmazmeans = np.zeros(N_bins)
szerrs = np.zeros(N_bins)
sig = sigma_z_all[sigma_z_all < 0.12]
z = z_trues[sigma_z_all < 0.12]
zmin,zmax = 0.08,0.35#min(z),max(z)
dz = (zmax - zmin)/N_bins
for j in range(N_bins):
    inds = (z > j*dz+zmin) * (z < (j+1)*dz+zmin)
    zmeans[j] = np.median(z[inds])
    sigmazmeans[j] = np.median(sig[inds])
    szerrs[j] = np.std(sig[inds])
output = np.array([zmeans,sigmazmeans]).T
np.savetxt("projection_effects.txt",output,header="z sigma(z)")
#plt.errorbar(zmeans,sigmazmeans,szerrs,marker='.')
plt.errorbar(zmeans,sigmazmeans,marker='o')
plt.ylabel(r"$\sigma_z$",fontsize=24)
plt.xlabel(r"$z_{\rm true}$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.legend(loc="lower right",fontsize=16)
plt.ylim(0,0.12)
plt.xlim(0.07,0.4)
plt.show()
plt.clf()
