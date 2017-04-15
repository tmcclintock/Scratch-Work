"""
A comparison between mine and Eduardo's TMF plots.
"""
import cosmocalc as cc
import numpy as np
import matplotlib.pyplot as plt
import tinker_mass_function as TMF
cosmo_dict = {"om":0.3,"ob":0.05,"ol":0.68,"h":0.7,"s8":0.8,"ns":0.96,"w0":-0.8,"wa":0.2}
cosmo_dict['ok'] = 1.0 - cosmo_dict['om'] - cosmo_dict['ol']
cc.set_cosmology(cosmo_dict)

zs = np.array([0.0,1.0])

for z in zs:
    ed = "mf_z%.1f.dat"%z
    me, emf = np.genfromtxt(ed, unpack=True)
    #Eduardo has M*dn/dM = dn/dlnM
    sf = 1.0/(1+z)
    t08 = me*np.array([cc.tinker2008_mass_function(Mi,sf,200) for Mi in me])
    t10 = me*np.array([cc.tinker2010_mass_function(Mi,sf,200) for Mi in me])

    TMF = TMF.tinker_mass_function(cosmo_dict, redshift=z)
    TMF.set_parameters(1.97,1.0,0.51,1.228)
    params = np.array([1.97,1.0,0.51,1.228])
    lM = np.log(me)
    dndlM = np.array([TMF.dndlM_at_lM(lMi,params) for lMi in lM])
    plt.plot(me, dndlM/emf,label="Eduardo")
    plt.plot(me, dndlM/t08,label="T08")
    plt.plot(me, dndlM/t10,label="T10")
    plt.xscale('log')
    plt.ylabel("dn/dlnM")
    plt.xlabel("Mass")
    plt.title("My code vs others at z=%.1f"%z)
    #plt.ylim(0.95, 1.05)
    plt.legend(loc=0)
    plt.show()
