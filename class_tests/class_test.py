"""
This is a test for using the CLASS module.
"""
from classy import Class
import numpy as np

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
k = np.logspace(-5,3,base=10,num=50)
z = np.array([0.5])
Pklin = np.array([cosmo.pk_lin(ki,z) for ki in k])
#Pklin = cosmo.get_pk_lin(k, z, len(k), 1, 1)
Pk = np.array([cosmo.pk(ki,z) for ki in k])
#Pk = cosmo.get_pk(k, z, len(k), 1, 1)
import matplotlib.pyplot as plt
plt.loglog(k,Pklin,label='linear')
plt.loglog(k,Pk,label='HF')
plt.legend()
plt.title("Fox power spectrum at z=%.1f"%z)
plt.show()

np.savetxt("pklin.txt",Pklin)
np.savetxt("pknl.txt",Pk)
np.savetxt("k.txt",k)
