"""
The boost factor models.
"""
import numpy as np

#Boost model
def model(params, l, z, R, pname="simple"):
    if pname == "simple":
        if len(params)==4:
            b0,c,d,e = params
        elif len(params)==5:
            b0,c,d,e,sigma = params
        return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(R/0.5)**e
    elif pname == "r1r2":
        if len(params)==5:
            b0,c,d,e1,e2 = params
        if len(params)==6:
            b0,c,d,e1,e2,sigma = params
        return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(e1/R+e2/R**2)
    elif pname == "bpl": #Broken power law
        if len(params)==6:
            b0,c,d,e1,e2,K = params
        if len(params)==7:
            b0,c,d,e1,e2,K,sigma = params
        boost1 = 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(R/0.5)**e1
        boost2 = 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d*(R/0.5)**e2
        return boost1*(R<K) + boost2*(R>=K)
    elif pname == "nfw":
        if len(params)==5: b0,c,d,e1,e2 = params
        if len(params)==6: b0,c,d,e1,e2,sigma = params
        x = R/e1
        i1 = np.where(x<1)[0]
        i2 = np.where(x>1)[0]
        Fx = np.ones_like(x)
        Fx[i2] *=  np.arctan(np.sqrt(x[i2]**2-1))/np.sqrt(x[i2]**2-1)
        Fx[i1] *= np.arctanh(np.sqrt(1-x[i1]**2))/np.sqrt(1-x[i1]**2)
        return 1.0 - b0*(l/30.0)**c*((1.+z)/1.5)**d * e2*(1-Fx)/(x**2-1)

    
def scatter_model(params, R, pname):
    if pname == "simple":
        if len(params)==4: sigma = 0
        elif len(params)==5: sigma = params[-1]
    elif pname == "r1r2":
        if len(params)==5: sigma = 0
        if len(params)==6: sigma = params[-1]
    elif pname == "bpl": #Broken power law
        if len(params)==6: sigma = 0
        if len(params)==7: sigma = params[-1]
    elif pname == "nfw":
        if len(params)==5: sigma = 0
        if len(params)==6: sigma = params[-1]
    return (sigma/R)**2 #R is in Mpc, pivot is 1 Mpc
