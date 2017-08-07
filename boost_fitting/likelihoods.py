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
        if e1 < 0: return - np.inf
    if b0 > 0: return -np.inf
    if c < 0: return -np.inf
    if d > 0: return -np.inf
    #if any(np.fabs(params)>20.0): return -np.inf
    return 0.0

def lnlike(params, lams, zs, R, Bp1, Berr, covs, pname):
    LL = 0
    for i in range(len(Bp1)): #Loop over all boost files
        for j in xrange(0,len(Bp1[i])):
            cov = np.copy(covs[i][j])
            Bmodel = model(params, lams[i,j], zs[i,j], R[i][j], pname)
            X = Bp1[i][j]-Bmodel
            scatter = scatter_model(params, R[i][j], pname)
            
            #cov = np.diag(np.diagonal(cov))
            #cov2 = np.outer(np.sqrt(scatter), np.sqrt(scatter))
            cov2 = np.diag(scatter)
            #print cov[0], cov[2]
            #print params
            #print cov+cov2
            icov = np.linalg.inv(cov + cov2)
            CHI2 = np.dot(X, np.dot(icov, X))
            OTHER = np.log(np.linalg.det(cov+cov2))
            LL += CHI2.sum()
            LL += OTHER
            #LL += np.sum(X**2/(Berr[i][j]**2+scatter)+\
            #                 np.log(Berr[i][j]**2+scatter))
    return -0.5*LL

def lnprob(params, lams, zs, R, Bp1, Berr, cov, pname):
    lnp = lnprior(params, pname)
    if not np.isfinite(lnp): return -np.inf
    return lnp + lnlike(params, lams, zs, R, Bp1, Berr, cov, pname)

