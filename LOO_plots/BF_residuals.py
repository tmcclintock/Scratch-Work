"""
Show the residuals of all BFs, or just the bottom plots.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=24)
import tinker_mass_function as TMF
import pickle

xlabel  = r"$\log_{10}M\ [{\rm M_\odot}/h]$"
y0label = r"$N/[{\rm Gpc}^3\  \log_{10}{\rm M_\odot}/h]$"
y1label = r"$\%\ {\rm Diff}$"

base = "/home/tmcclintock/Desktop/Github_stuff/Mass-Function-Emulator/test_data/"
datapath = base+"N_data/Box%03d_full/Box%03d_full_Z%d.txt"
covpath  = base+"covariances/Box%03d_cov/Box%03d_cov_Z%d.txt"

N_cosmos = 39
N_z      = 10
c = np.linspace(1.0, 0.0, N_z) #Colors
cmap = plt.get_cmap('seismic') 
#cmap = plt.get_cmap('RdBu') 

scale_factors = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 
                          0.645161, 0.714286, 0.8, 0.909091, 1.0])
redshifts = 1./scale_factors - 1.0
volume = 1050.**3 #[Mpc/h]^3

cosmos = np.genfromtxt("cosmos.txt")

#This contains our parameterization
name = 'dfgB'
base_dir = "../../fit_mass_functions/output/%s/"%name
base_save = base_dir+"%s_"%name
best_fit_models = np.loadtxt(base_save+"bests.txt")

def get_params(model, sf):
    Tinker_defaults = {'d':1.97, 'e':1.0, "f": 0.51, 'g':1.228}
    B=None
    if name is 'dfgB':
        d0,d1,f0,f1,g0,g1,B = model
        e0 = Tinker_defaults['e']
        e1 = 0.0
    if name is 'defg':
        d0,d1,e0,e1,f0,f1,g0,g1 = model
    if name is 'dfg':
        d0,d1,f0,f1,g0,g1 = model
        e0 = Tinker_defaults['e']
        e1 = 0.0
    if name is 'efg':
        e0,e1,f0,f1,g0,g1 = model
        d0 = Tinker_defaults['d']
        d1 = 0.0
    if name is 'fg':
        f0,f1,g0,g1 = model
        d0 = Tinker_defaults['d']
        d1 = 0.0
        e0 = Tinker_defaults['e']
        e1 = 0.0
    k = sf - 0.5
    d = d0 + k*d1
    e = e0 + k*e1
    f = f0 + k*f1
    g = g0 + k*g1
    return d,e,f,g,B

make_BFs = True
BF_path = "saved_files/%s_bf_curves.p"%name
BFs = []
if make_BFs:
    for i in range(0, N_cosmos):
        #Get in the cosmology and create a cosmo_dict
        num, ombh2, omch2, w0, ns, ln10As, H0, Neff, sigma8 = cosmos[i]
        h = H0/100.
        Ob,Om = ombh2/(h**2), (ombh2+omch2)/(h**2)
        cosmo_dict = {"om":Om, "ob":Ob, "ol":1-Om, "ok":0.0, "h":h, 
                      "s8":sigma8, "ns":ns, "w0":w0, "wa":0.0}
        BFs_cos = []
        for j in range(N_z):
            data = np.loadtxt(datapath%(i, i, j))
            lM_bins = data[:,:2]

            #Now get the B.
            TMF_model = TMF.tinker_mass_function(cosmo_dict, redshifts[j])
            d,e,f,g,B = get_params(best_fit_models[i], scale_factors[j])
            TMF_model.set_parameters(d,e,f,g,B)
            N_bf = volume * TMF_model.n_in_bins(lM_bins)
            BFs_cos.append(N_bf)
        BFs.append(BFs_cos)
    pickle.dump(BFs, open(BF_path, "wb"))
BFs = pickle.load(open(BF_path, "rb"))

fig, axis = plt.subplots(1, sharex=True)

for i in range(0, N_cosmos):
    #Get in the cosmology and create a cosmo_dict
    num, ombh2, omch2, w0, ns, ln10As, H0, Neff, sigma8 = cosmos[i]
    h = H0/100.
    Ob,Om = ombh2/(h**2), ombh2/(h**2)+omch2/(h**2)
    cosmo_dict = {"om":Om, "ob":Ob, "ol":1-Om, "ok":0.0, "h":h, 
                  "s8":sigma8, "ns":ns, "w0":w0, "wa":0.0}

    for j in range(0,N_z):
        #First get the data.
        data = np.loadtxt(datapath%(i, i, j))
        lM_bins = data[:,:2]
        lM = np.mean(data[:, :2], 1)
        N = data[:,2]
        cov = np.loadtxt(covpath%(i, i, j))
        err = np.sqrt(np.diagonal(cov))

        #Now get the B.
        N_bf = BFs[i][j]

        #Plot the %difference
        pd  = 100.*(N-N_bf)/N_bf
        pde = 100.*err/N_bf
        for k in [0,-1]:
            axis.errorbar(lM[k], pd[k], pde[k], marker='.', ls='', c=cmap(c[j]), alpha=0.2)
        axis.plot(lM[1:-1], pd[1:-1], 
                  marker='.', ls='', c=cmap(c[j]), alpha=0.2)

    axis.axhline(0, c='k', ls='-', zorder=-1)

#Show
axis.set_xlabel(xlabel)
axis.set_ylabel(y1label)
axis.set_ylim(-20, 20)
plt.subplots_adjust(bottom=0.15, left=0.15, hspace=0.0)
fig.savefig("fig_BF_residuals.pdf")
plt.show()
