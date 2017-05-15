"""
Compare DeltaSigma to Build_Delta_Sigma.
"""
import numpy as np
import os, sys
sys.path.insert(0, "../../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pyDS
sys.path.insert(0, "../../Build-Delta-Sigma/src/wrapper/")
import py_Build_Delta_Sigma
from colossus.halo import concentration as conc
from colossus.cosmology import cosmology as col_cosmology
import matplotlib.pyplot as plt
plt.rc("text", usetex=True, fontsize=24)

#Get the power spectra
z = 1.0
k = np.loadtxt("k.txt")
Plin = np.loadtxt("Plin_z%.2f.txt"%z)
Pnl  = np.loadtxt("Pnl_z%.2f.txt"%z)

#This is the fox sim cosmology
h = 0.670435
cosmo = {"h":h,"om":0.31834,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]
colcos = {"H0":cosmo['h']*100.,"Om0":cosmo['om'], 
          'Ob0': 0.049017, 'sigma8': 0.83495, 'ns': 0.96191, 'flat':True}
col_cosmology.addCosmology('fiducial_cosmology', colcos)
col_cosmology.setCosmology('fiducial_cosmology')
params = {"NR":300,"Rmin":0.01,
          "Rmax":200.0,"Nbins":15,"delta":200,
          "Rmis":0.2, "fmis":0.0,
          "miscentering":0,"averaging":0}

mass = 10**13
params['Mass'] = mass
params["concentration"] = conc.concentration(mass, 
                                             '200m', z, 
                                             model='diemer15')
params["R_bin_min"] = 0.0323*h*(1+z)
params["R_bin_max"] = 30.0*h*(1+z)

result = pyDS.calc_Delta_Sigma(k, Plin, k, Pnl, cosmo, params)
DSm = result['delta_sigma']
xihm = result['xi_hm']
R = result['R']
print "Have model"
print R.shape, xihm.shape, DSm.shape
#print result.keys()

bds_params = {"Mass": params['Mass'], "delta":200, 
              "timing":0, "miscentering":0,
              "concentration":params["concentration"]}
result2 = py_Build_Delta_Sigma.build_Delta_Sigma(R, xihm, cosmo, bds_params)
DSb = result2['delta_sigma']
print "Have built version"

plt.loglog(R, DSm, c='b')
plt.loglog(R, DSb, c='r', ls='--')
plt.show()
