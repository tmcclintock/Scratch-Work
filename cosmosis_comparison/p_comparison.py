import numpy as np
from classy import Class
import os, sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

#Get the mean redshifts
zs = [1.0, 0.5, 0.25, 0.04]

#Fox cosmology
Ob = 0.049017
Om = 0.31834
Ocdm = Om - Ob
Omega_Lambda = 1. - Ocdm
h = 0.670435
print Ob, Ob, Ob, Ocdm, Ocdm, Ocdm
params = {
        'output': 'mPk',
        "h":h,
        "A_s":2.14013e-9,
        "n_s":0.96191,
        "tau_reio":0.08,
        "Omega_b":Ob,
        "Omega_cdm":Ocdm,
        'P_k_max_h/Mpc':2000.,
        'z_max_pk':1.0,
        'non linear':'halofit',
        'input_verbose':1,
        'background_verbose':1,
        'thermodynamics_verbose':1,
        'perturbations_verbose':1,
        'transfer_verbose':1,
        'primordial_verbose':1,
        'spectra_verbose':1,
        'nonlinear_verbose':1,
        'lensing_verbose':1,
        'output_verbose':1,
        'reio_parametrization':'reio_camb',
        #"k_pivot":50
        }

cosmo = Class()
cosmo.set(params)
#cosmo.set_default()
#cosmo.compute()
#print "pars:",cosmo.pars
"""
print "\n"
print "Neff:",cosmo.Neff()
print "H?:  ",cosmo.Hubble(1.0)
print "s8:  ",cosmo.sigma8()
print "Ob:  ",cosmo.Omega_b()
print "Om:  ",cosmo.Omega_m()
print "O0m: ",cosmo.Omega0_m()
print "Tcmb:",cosmo.T_cmb()
"""


k = np.logspace(-5, 3, base=10, num=4000)
for i in range(len(zs)):
    z = zs[i]

    k = np.loadtxt("../../fox_calibration/txt_files/P_files/k.txt")
    Pmm  = np.loadtxt("../../fox_calibration/txt_files/P_files/Pnl_z%.2f.txt"%(z))
    Plin = np.loadtxt("../../fox_calibration/txt_files/P_files/Plin_z%.2f.txt"%(z))
    #Pmm  = np.array([cosmo.pk(ki, z) for ki in k])
    #Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])

    kmmc = np.loadtxt("cosmosis_outputs/knl_z%.2f.txt"%z)
    Pmmc = np.loadtxt("cosmosis_outputs/pnl_z%.2f.txt"%z)
    klinc = np.loadtxt("cosmosis_outputs/klin_z%.2f.txt"%z)
    Plinc = np.loadtxt("cosmosis_outputs/plin_z%.2f.txt"%z)

    #s8class = cosmo.sigma8()
    #s8camb  = 0.83495
    #fix = s8camb**2/s8class**2
    
    cosspl = IUS(klinc, Plinc)
    claspl = IUS(k, Plin)
    kratio = 0.001
    ratio = cosspl(kratio)/ claspl(kratio)
#    ratio = 0.756
    #print "sigma8 = ", cosmo.sigma8()
    print "z=%.2f  ratio=%f  "%(z, ratio)

    #plt.loglog(k/h, Pmm*h**3, label="class nl")
    #plt.loglog(k/h, Plin*h**3, label="class lin")
    plt.loglog(k, Pmm, label="class nl")
    plt.loglog(k, Plin, label="class lin")
    plt.loglog(kmmc, Pmmc, label="camb nl")
    plt.loglog(klinc, Plinc, label="camb lin")
    plt.legend(loc=0)
    plt.title("z=%.2f"%z)
    plt.show()

    print "Done with z%.1f"%z
