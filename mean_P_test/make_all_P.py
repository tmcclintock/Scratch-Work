"""
Call class to calculate the power spectra at each redshift
for the clusters in the z0l4 file, which has 62 clusters.
"""
import numpy as np
from classy import Class

#Read in the clusters
cpath = "cluster_files/clusters_z0l4.txt"
zs, lams = np.genfromtxt(cpath, unpack=True)
zmean = np.mean(zs)
zerr = np.std(zs)
print "Working with %d clusters from the file %s"%(len(zs),cpath)
print "Mean redshift = %.3f +- %.3f"%(zmean,zerr)

#Fox cosmology
Ob = 0.049017
Om = 0.31834
Ocdm = Om - Ob
params = {
    'output': 'mPk transfer',
    "h":0.670435,
    "A_s":2.1e-9,
    "n_s":0.96191,
    "omega_b":Ob,
    "omega_cdm":Ocdm,
    'YHe':0.24755048455476272,#By hand, default value
    'P_k_max_1/Mpc':2000,
    'z_pk':0.5,
    'non linear':'halofit'}

cosmo = Class()
cosmo.set(params)
cosmo.compute()

k = np.logspace(-5, 3, base=10, num=50)
np.savetxt("PK_files/k.txt",k)
meanPnlpath  = "PK_files/Pmean_nl.txt"
meanPlinpath = "PK_files/Pmean_lin.txt"
Pmean    = np.array([cosmo.pk(ki, zmean) for ki in k])
Pmeanlin = np.array([cosmo.pk_lin(ki, zmean) for ki in k])
np.savetxt(meanPnlpath, Pmean)
np.savetxt(meanPlinpath, Pmeanlin)
print "Mean power spectra saved"

Pnlpath  = "PK_files/Pnl_%d.txt"
Plinpath = "PK_files/Plin_%d.txt"
for i in range(len(zs)):
    z = zs[i]
    print "Working on cluster %d at z=%.3f"%(i,z)
    P    = np.array([cosmo.pk(ki, z) for ki in k])
    Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])
    np.savetxt(Pnlpath%i, P)
    np.savetxt(Plinpath%i, Plin)
