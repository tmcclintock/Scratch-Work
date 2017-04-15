"""
A test of g(sigma), where sigma(M) is calculated from Matt Becker's code.
"""
import numpy as np
import matplotlib.pyplot as plt
import cosmocalc as cc

zs = np.array([0.0,1.0])
sf = 1.0/(1+zs)

cosmo_dict = {"om":0.3,"ob":0.05,"ol":0.68,"h":0.7,"s8":0.8,"ns":0.96,"w0":-0.8,"wa":0.2}
cc.set_cosmology(cosmo_dict)

M,emf = np.genfromtxt("mf_z0.0.dat", unpack=True)

"""
for i in range(len(zs)):
    z = zs[i]
    a = sf[i]
    sigma = np.array([cc.sigmaMtophat_exact(Mi,a) for Mi in M])
    plt.plot(M,sigma,label="z=%.1f"%z)
plt.xscale('log')
plt.legend(loc=0)
#plt.show()
plt.clf()
"""

t08 = M*np.array([cc.tinker2008_mass_function(Mi,sf[0],200) for Mi in M])
t10 = M*np.array([cc.tinker2010_mass_function(Mi,sf[0],200) for Mi in M])

plt.plot(M, emf/t08)
plt.plot(M, emf/t10)
plt.xscale('log')
plt.show()
