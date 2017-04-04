"""
This is a test for using the CLASS module.
"""
from classy import Class
import numpy as np

params = {
    'output': 'mPk',
    'P_k_max_1/Mpc':2000,
    'non linear':'halofit'}

cosmo = Class()

cosmo.set(params)

cosmo.compute()
k = np.logspace(-5,3,base=10,num=50)
Pklin = np.array([cosmo.pk_lin(ki,0.0) for ki in k])
Pk = np.array([cosmo.pk(ki,0.0) for ki in k])
print dir(cosmo)
import matplotlib.pyplot as plt
plt.loglog(k,Pklin,label='linear')
plt.loglog(k,Pk,label='HF')
plt.legend()
plt.show()
