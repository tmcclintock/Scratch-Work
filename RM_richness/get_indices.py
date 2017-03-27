"""
Get a set of indices to look at for the plots.
"""
import fitsio, sys, os
import numpy as np
fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
sigz = np.loadtxt("sigma_z_all.txt")
z = data['Z_LAMBDA']

inds = np.where((z>0.3)*(z<0.31))[0]
z = z[inds]
sigz = sigz[inds]
isort = np.argsort(sigz)
z = z[isort]
sigz = sigz[isort]
inds = inds[isort]
print sigz.shape
print inds[::30]
