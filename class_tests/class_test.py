"""
This is a test for using the CLASS module.
"""
from classy import Class
import numpy as np

Ob = 0.047
Om = 0.286
Ocdm = Om - Ob

"""params = {
    'output': 'mPk transfer',
    "h":0.670435,
    "A_s":2.1e-9,
    "n_s":0.96191,
    "omega_b":Ob,
    "omega_cdm":Ocdm,
    'YHe':0.24755048455476272,#By hand, default value
    'P_k_max_1/Mpc':2000,
    'z_pk':0.5,
    'non linear':'halofit'}"""
params = {
    'output': 'mPk transfer',
    "h":0.7,
    "ln10^{10}A_s":3.065525,
    "n_s":0.96,
    "Omega_b":Ob,
    "Omega_cdm":Ocdm,
    'P_k_max_1/Mpc':100,
    'z_pk':0.0}

cosmo = Class()

cosmo.set(params)
print cosmo.sigma8()
exit()


cosmo.compute()
k = np.logspace(-5,3,base=10,num=50)
zs = np.array([0.0,0.5, 1.0])
import matplotlib.pyplot as plt
for z in zs:
    Pklin = np.array([cosmo.pk_lin(ki,z) for ki in k])
    Pk = np.array([cosmo.pk(ki,z) for ki in k])
    plt.loglog(k,Pklin,label='linear z=%.1f'%z)
    plt.loglog(k,Pk,label='HF z=%.1f'%z)
plt.legend()
plt.show()
