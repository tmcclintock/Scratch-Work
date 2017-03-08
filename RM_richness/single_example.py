"""
Let's do a single cluster. Then we will do the calculation for the width.
"""
import fitsio
import numpy as np
from scipy.optimize import minimize

fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
lam_true = data['LAMBDA_CHISQ'][0]
z_true = data['Z_LAMBDA'][0]
lam_arrays = data['LAMBDA_CHISQS']
zname = 'redshift_list.txt'
zs = np.genfromtxt(zname).flatten()

lam_data = lam_arrays[0]
#Fix the data
lam_data[lam_data<0.0] = 0.0
lam_n_data = lam_data/max(lam_data) #Normalized
print lam_n_data.shape,zs.shape

def get_lam_model(z,sigmaz,z_true,lam_true):
    return lam_true*np.exp(-0.5*(z_true-z)**2/sigmaz**2)

def total_diff(params,zs,lams):
    z_true,sigmaz,lam_true = params
    if any(params < 0.0): return np.inf
    if z_true > 2.0: return np.inf
    if lam_true > 300: return np.inf
    if sigmaz < 0.0: return np.inf
    lam_model = get_lam_model(zs,sigmaz,z_true,lam_true)
    return np.sum(np.fabs(lams-lam_model)/lam_model)
    #return np.sum(np.fabs(lam_n-lam_n_model)**2/lam_n_model**2)

x0 = np.array([z_true,0.2,lam_true])
print total_diff(x0,zs,lam_data)
print minimize(total_diff,x0=x0,args=(zs,lam_data),method='Nelder-Mead')
