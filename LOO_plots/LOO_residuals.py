"""
Show the residuals of all regular emulators, or just the bottom panels.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=24)
import tinker_mass_function as TMF
import sys, os, pickle, emulator

def train(training_cosmos, training_data, training_errs):
    N_cosmos = len(training_cosmos)
    N_emulators = training_data.shape[1]
    emulator_list = []
    for i in range(N_emulators):
        y = training_data[:, i]
        yerr = training_errs[:, i]
        emu = emulator.Emulator(name="emu%d"%i, xdata=training_cosmos, 
                                ydata=y, yerr=yerr)
        emu.train()
        emulator_list.append(emu)
    return emulator_list

def predict_parameters(cosmology, emu_list):
    return np.array([emu.predict_one_point(cosmology)[0] for emu in emu_list])

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

scale_factors = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 
                          0.645161, 0.714286, 0.8, 0.909091, 1.0])
redshifts = 1./scale_factors - 1.0
volume = 1050.**3 #[Mpc/h]^3

cosmos = np.genfromtxt("cosmos.txt")

building_cosmos = np.delete(cosmos, 0, 1) #Delete boxnum
building_cosmos = np.delete(building_cosmos, 4, 1) #Delete ln10As
building_cosmos = np.delete(building_cosmos, -1, 0)#39 is broken

#This contains our parameterization
name = 'dfg'
base_dir = "../../fit_mass_functions/output/%s/"%name
base_save = base_dir+"%s_"%name
mean_models = np.loadtxt(base_save+"means.txt")
err_models = np.sqrt(np.loadtxt(base_save+"vars.txt"))

def get_params(model, sf):
    Tinker_defaults = {'d':1.97, 'e':1.0, "f": 0.51, 'g':1.228}
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
    return d,e,f,g

make_LOObad = False
LOObad_path = "saved_files/%s_LOObad_curves.p"%name
LOObad = []
if make_LOObad:
    for i in range(0, N_cosmos):
        #Get in the cosmology and create a cosmo_dict
        num, ombh2, omch2, w0, ns, ln10As, H0, Neff, sigma8 = cosmos[i]
        h = H0/100.
        Ob,Om = ombh2/(h**2), (ombh2+omch2)/(h**2)
        cosmo_dict = {"om":Om, "ob":Ob, "ol":1-Om, "ok":0.0, "h":h, 
                      "s8":sigma8, "ns":ns, "w0":w0, "wa":0.0}
        test_cosmo = building_cosmos[i]
        test_data  = mean_models[i]
        test_err   = err_models[i]
        training_cosmos = np.delete(building_cosmos, i, 0)
        training_data   = np.delete(mean_models, i, 0)
        training_errs   = np.delete(err_models, i, 0)
        #Train the emulators
        emu_list = train(training_cosmos, training_data, training_errs)
        emu_model = predict_parameters(test_cosmo, emu_list)

        LOObad_cos = []
        for j in range(N_z):
            data = np.loadtxt(datapath%(i, i, j))
            lM_bins = data[:,:2]

            #Now get the B.
            TMF_model = TMF.tinker_mass_function(cosmo_dict, redshifts[j])
            d,e,f,g = get_params(emu_model, scale_factors[j])
            TMF_model.set_parameters(d,e,f,g)
            N_emu = volume * TMF_model.n_in_bins(lM_bins)
            LOObad_cos.append(N_emu)
        LOObad.append(LOObad_cos)
    pickle.dump(LOObad, open(LOObad_path, "wb"))
LOObad = pickle.load(open(LOObad_path, "rb"))

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
        N_bf = LOObad[i][j]

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
fig.savefig("fig_LOObad_residuals.pdf")
plt.show()
