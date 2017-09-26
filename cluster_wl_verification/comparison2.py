"""
Calculate stuff from my DeltaSigma repo and my Cluster_wl repo and compare the results.
"""
import sys
sys.path.insert(0, "../../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pDS
import numpy as np
import matplotlib.pyplot as plt

knl = np.loadtxt("./knl.txt")
Pnl = np.loadtxt("./pnl.txt")
klin = np.loadtxt("./klin.txt")
Plin = np.loadtxt("./plin.txt")
cosmo = {"h":0.7,"om":0.3}
om = cosmo['om']

NR = 1000
R = np.logspace(-2, 3, NR, base=10) #Xi_hm MUST be evaluated to higher than BAO
M = 1e14 
c = 5

input_params = {"Mass": M, "concentration":c, "NR":NR,"Rmin":0.01,
                "Rmax":400.,"Nbins":15,"R_bin_min":0.0323,"R_bin_max":30.0,
                "delta":200,"averaging":1}
outdict = pDS.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)

R = outdict['R']
xi_nfw = outdict['xi_1halo']
xi_hm = outdict['xi_hm']
Sigma = outdict['sigma']
DS = outdict['delta_sigma']
aDS = outdict['ave_delta_sigma']
Rb = outdict['Rbins']


from scipy.interpolate import InterpolatedUnivariateSpline as IUS
dsspl = IUS(R, DS)


pdiff = np.fabs(aDS-dsspl(Rb))/dsspl(Rb)
print max(pdiff)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].loglog(R, DS)
axarr[0].loglog(Rb, aDS, ls='--')

axarr[1].loglog(Rb, pdiff)
axarr[1].set_yscale('log')
#axarr[1].set_ylim(0.01, 2.0)
plt.show()
