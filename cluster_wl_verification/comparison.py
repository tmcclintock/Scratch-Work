"""
Calculate stuff from my DeltaSigma repo and my Cluster_wl repo and compare the results.
"""
import sys
sys.path.insert(0, "../../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma as pDS
import clusterwl
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

R_pds = outdict['R']
xi_nfw_pds = outdict['xi_1halo']
xi_hm_pds = outdict['xi_hm']
Sigma_pds = outdict['sigma']
DSpds = outdict['delta_sigma']
aDSpds = outdict['ave_delta_sigma']
Rbins_pds = outdict['Rbins']

xi_nfw   = clusterwl.xi.xi_nfw_at_R(R, M, c, om)
xi_mm    = clusterwl.xi.xi_mm_at_R(R, knl, Pnl)
bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)
Rp = np.logspace(-2, 2.7, NR, base=10)
Sigma  = np.zeros_like(Rp)
DeltaSigma = np.zeros_like(Rp)
Nbins = 15
Redges = np.logspace(np.log(0.0323), np.log(30.), num=Nbins+1, base=np.e)
aDS = np.zeros((Nbins))

clusterwl.deltasigma.calc_Sigma_at_R(Rp, R, xi_hm, M, c, om, Sigma)
clusterwl.deltasigma.calc_DeltaSigma_at_R(Rp, Rp, Sigma, M, c, om, DeltaSigma)
clusterwl.averaging.average_profile_in_bins(Redges, Rp, DeltaSigma, aDS)

x1 = DeltaSigma
x2 = DSpds

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
spl = IUS(Rp, x1)

print aDSpds - aDS

#x1 = DeltaSigma
#x2 = DSpds




f, axarr = plt.subplots(2, sharex=True)
axarr[0].loglog(Rp, x1)
axarr[0].loglog(R_pds, x2, ls='--')
axarr[0].loglog(Rbins_pds, aDSpds, ls=':')
axarr[0].loglog(Rbins_pds, aDS, ls=':')

pd2 = np.fabs(aDSpds - spl(Rbins_pds))/spl(Rbins_pds)
pd3 = np.fabs(aDS - spl(Rbins_pds))/spl(Rbins_pds)

axarr[1].loglog(R_pds, np.fabs(x2-spl(R_pds))/spl(R_pds))
axarr[1].loglog(Rbins_pds, pd2,ls="--")
axarr[1].loglog(Rbins_pds, pd3,ls=":")
axarr[1].set_yscale('log')
#axarr[1].set_ylim(0.01, 1.0)
plt.show()
