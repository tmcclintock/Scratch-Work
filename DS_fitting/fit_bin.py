"""
Fitting some SV data but including halo concentration
as a parameter. Doing this outside of the 
CosmoSIS framework to see if it will work.
"""
import numpy as np
import emcee, sys
sys.path.insert(0,"../../Delta-Sigma/src/wrapper/")
import py_Delta_Sigma
import scipy.optimize as op
import corner

single_test = True
maxlike = True
domc = True
seecorner = True

ndim = 2
nwalkers, nsteps = 8, 1000
nburn = 200
corner_labels = [r"$\log_{10}M$",r"$c$"]

h = 0.7 #Hubble constant

#Get in the data
data = np.genfromtxt("profile_z0_l4.dat")[2:]
cov = np.genfromtxt("cov_t_z0_l4.dat")
cov = cov[2:] #limit is missing in data file
cov = cov[:,2:]#this is only a feature of the cov
boost = np.genfromtxt("boost_erinweight.txt")
boost = boost[boost[:,0]==0]
boost = boost[boost[:,1]==4]
boost = boost[:,3]
ds = data[:,1]
R = data[:,0] #Mpc physical
boost = boost[:-1] #
ds = ds[:-1]
print R
sys.exit()
R = R[:-1]
cov = cov[:-1]
cov = cov[:,:-1]
#Figure out the cut index. This will be passed through into the likelihood
cut_scale = 0.5 #Mpc physical
cut = min(np.where(R>cut_scale)[0])
R = R[cut:]
ds = ds[cut:]
boost = boost[cut:]
cov = cov[cut:]
cov = cov[:,cut:]
print R.shape,cov.shape
#Add 2 to the cut to signify the lower two bins
#that were missing from the profile.dat file but WERE
#present in the cov file. thanks, tamas
cut += 2

#Invert the cov matrix
icov = np.linalg.inv(cov)

redshift = 3.147470467649999826e-01
richness = 5.372359058220000350e+01
Rlam = 0.88314702 #Mpc/h
#Note to convert the DS model to be like the data
#I multiply by (1+redshift)**2 * h
binmin = 0.02*h*(1+redshift)
binmax = 30.0*h*(1+redshift)

#DS parameters
A = 1.026
fmis = 0.24
lnc = -1.06
Rmis = np.exp(lnc)*Rlam
#Boost parameters
B0 = 10**-1.399
Clam = 0.92
Dz = -4.0
ER = -0.98
zpivot = 0.5
lampivot = 30.0
Rpivot = 0.5 #Mpc phsyical
boost_model = 1+B0*((1+redshift)/(1+zpivot))**Dz * (richness/lampivot)**Clam * (R/Rpivot)**ER

#Create a dictionary with the cosmology
cosmo = {"h":h,"om":0.3,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]

#Load in the power spectra
klin = np.genfromtxt("klin.txt")
Plin = np.genfromtxt("Plin.txt")
knl  = np.genfromtxt("knl.txt")
Pnl  = np.genfromtxt("Pnl.txt")

def lnlike(params,data_args,param_args,pk_args):
    log10M,c = params
    ds_data,icov,boost,cut = data_args
    A,fmis,Rmis = param_args
    klin,Plin,knl,Pnl = pk_args

    input_params = {"Mass": 10**log10M,"NR":300,"Rmin":0.01,
                    "Rmax":100.0,"Nbins":15,
                    "R_bin_min":binmin,"R_bin_max":binmax,
                    "delta":200,"Rmis":Rmis,"fmis":fmis,
                    "timing":0,"miscentering":1,"averaging":1,
                    "concentration":c}
    
    return_dict = py_Delta_Sigma.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)
    ads = return_dict['ave_delta_sigma']
    mds = return_dict['ave_miscentered_delta_sigma']
    ds = (1-fmis)*ads + fmis*mds
    ds = ds[cut:-1] #chop off lower bins and upper bin (just for z0)
    ds *= (1+redshift)**2*h
    X = ds_data - A*ds/boost
    LL = -0.5 * np.dot(X,np.dot(icov,X))
    return LL

def lnprior(params):
    log10M,c = params
    if log10M < 10.0: return -np.inf
    if log10M > 17.0: return -np.inf
    if c < 0.05: return -np.inf
    return 0.0

def lnprob(params,data_args,param_args,pk_args):
    LP = lnprior(params)
    if np.isinf(LP): return LP
    return LP + lnlike(params,data_args,param_args,pk_args)

guess = [14.5,3.7]
data_args = [ds,icov,boost_model,cut]
param_args = [A,fmis,Rmis]
pk_args = [klin,Plin,knl,Pnl]

if single_test:
    print lnprob(guess,data_args,param_args,pk_args)

if maxlike:
    nll = lambda *args:-lnprob(*args)
    result = op.minimize(nll,guess,args=(data_args,param_args,pk_args))
    print "Best fit:",result
    lM = result['x'][0]
    print "True best mass (includes correction):",np.log10(10**lM*1.02/0.7)
    np.savetxt("outputs/500kpc_best_params.txt",result['x'])

if domc:
    start = np.loadtxt("outputs/500kpc_best_params.txt")
    pos = [start + 1e-2*np.random.randn(ndim) for k in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(data_args,param_args,pk_args))
    print "Starting emcee with %d steps"%nsteps
    sampler.run_mcmc(pos,nsteps)
    fullchain = sampler.flatchain
    likes = sampler.flatlnprobability
    np.savetxt("outputs/500kpc_fullchain.txt",fullchain)
    np.savetxt("outputs/500kpc_likes.txt",likes)
    print "Done with mcmc"

if seecorner:
    fullchain = np.loadtxt("outputs/500kpc_fullchain.txt")
    print fullchain.shape
    fullchain[:,0] = np.log10(10**fullchain[:,0]*1.02/0.7)
    chain = fullchain[nwalkers*nburn:]
    stds = np.sqrt(np.var(chain,0))
    print stds
    import matplotlib.pyplot as plt
    plt.rc('text',usetex=True, fontsize=20)

    fig = corner.corner(chain,labels=corner_labels,plot_datapoints=False,truths = [14.592,None])
    plt.show()
