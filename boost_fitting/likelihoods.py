"""
The likelihoods for the boost factor fitting.
"""
import numpy as np
from models import *

#Define the likelihoods
def lnprior(params, pname="simple"):
    lnp = 0
    if pname == "simple":
        if len(params)==4:
            b0,c,d,e = params
        elif len(params)==5:
            b0,c,d,e,sigma = params
        if e > 0: return -np.inf
    elif pname == "r1r2":
        if len(params)==5: b0,c,d,e1,e2 = params
        if len(params)==6: b0,c,d,e1,e2,sigma = params
        #if e1 < 0 or e2 < 0: return -np.inf
    elif pname == "bpl": #Broken power law
        if len(params)==6: b0,c,d,e1,e2,K = params
        if len(params)==7: b0,c,d,e1,e2,K,sigma = params
        if K < 0.2: return -np.inf
        if e1 > 0 or e2 > 0: return -np.inf
    elif pname == "nfw":
        if len(params)==6: b0,c,d,e1,e2,e3 = params
        if len(params)==7: b0,c,d,e1,e2,e3,sigma = params
        if e1 < 0: return -np.inf
    elif pname == "nfw_single":
        if len(params)==2: b0,rs = params
        if len(params)==3: b0,rs,sigma = params
        c,d = -1, 1 #garbage
        if rs < 0: return -np.inf
    if b0 > 0: return -np.inf
    if c < 0: return -np.inf
    if d > 0: return -np.inf
    return lnp

def lnlike(params, lams, zs, R, Bp1, Berr, covs, pname):
    LL = 0
    for i in range(len(Bp1)): #Loop over all boost files
        for j in xrange(3,len(Bp1[i])):
            LL += lnlike_single(params, lams[i,j], zs[i,j], R[i][j], Bp1[i][j],\
                                Berr[i][j], covs[i][j], pname)
    return LL

def lnlike_single(params, lam, z, R, Bp1, Berr, cov, pname):
    Bmodel = model(params, lam, z, R, pname)
    X = Bp1-Bmodel
    scatter = scatter_model(params, R, pname)
    cov2 = np.diag(scatter)
    icov = np.linalg.inv(cov + cov2)
    CHI2 = np.dot(X, np.dot(icov, X))
    return -0.5*(CHI2+np.log(np.linalg.det(cov+cov2)))

def lnprob(params, lams, zs, R, Bp1, Berr, cov, pname):
    lnp = lnprior(params, pname)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, lams, zs, R, Bp1, Berr, cov, pname)

