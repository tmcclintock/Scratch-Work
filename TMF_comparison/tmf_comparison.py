"""
A comparison between mine and Eduardo's TMF plots.
"""
import numpy as np
import matplotlib.pyplot as plt

ed = "mf_z0.0.dat"
ca = "mf_camb_z0.0.dat"
me, emf = np.genfromtxt(ed, unpack=True)
mc, cmf = np.genfromtxt(ca, unpack=True)
#Eduardo and camb give M*dn/dM = dn/dlnM

#plt.loglog(me, emf)
#plt.loglog(mc, cmf)
plt.plot(me, emf/cmf)
plt.xscale('log')
plt.ylabel("dn/dlnM")
plt.xlabel("Mass")
plt.title("Eduardo code vs CAMB")
#plt.show()
plt.clf()

cosmo_dict = {"om":0.3,"ob":0.05,"ol":0.68,"h":0.7,"s8":0.8,"ns":0.96,"w0":-0.8,"wa":0.2}
cosmo_dict['ok'] = 1.0 - cosmo_dict['om'] - cosmo_dict['ol']
import tinker_mass_function as TMFmod
TMF = TMFmod.TMF_model(cosmo_dict, redshift=0.0)
TMF.set_parameters(1.97,1.0,0.51,1.228)
params = np.array([1.97,1.0,0.51,1.228])
lM = np.log(me)
print dir(TMF)
dndlM = np.array([TMF.dndlM_at_lM(lMi,params) for lMi in lM])
plt.plot(me, dndlM/emf,label="Eduardo")
plt.plot(me, dndlM/cmf,label="camb")
plt.xscale('log')
plt.ylabel("dn/dlnM")
plt.xlabel("Mass")
plt.title("My code vs others")
plt.legend(loc=0)
plt.show()
