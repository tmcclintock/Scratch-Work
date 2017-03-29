import sys,os
import numpy as np
sys.path.insert(0,"../../../Mass-Function-Emulator/")
import rot_mf_emulator as MFE
sys.path.insert(0,'../../../Mass-Function-Emulator/visualization/')
import visualize


#Scale factors and redshifts
scale_factors = np.array([0.25,0.333333,0.5,0.540541,0.588235,0.645161,0.714286,0.8,0.909091,1.0])
redshifts = 1./scale_factors - 1.0
volume = 1050.**3 #[Mpc/h]^3

#Read in the input cosmologies
cosmologies = np.genfromtxt("../test_data/building_cosmos.txt")
cosmologies = np.delete(cosmologies,5,1) #Delete ln10As
cosmologies = np.delete(cosmologies,0,1) #Delete boxnum
cosmologies = np.delete(cosmologies,-1,0)#39 is broken
N_cosmologies = len(cosmologies)
N_z = len(redshifts)

#Read in the input data
means = np.loadtxt("../full_training_data/rotated_mean_models.txt")
variances = np.loadtxt("../full_training_data/rotated_var_models.txt")
data = np.ones((N_cosmologies,len(means[0]),2)) #Last column is for mean/errs
data[:,:,0] = means
data[:,:,1] = np.sqrt(variances)

chi2s = np.zeros((N_cosmologies,N_z))
Nfp = np.zeros((N_cosmologies,N_z))

MF_path = "../../../../all_MF_data/"

percent = 10
per = percent/100.

#Loop over all boxes
for i in xrange(0,N_cosmologies):
    test_cosmo = cosmologies[i]
    test_data = data[i]
    training_cosmos = np.delete(cosmologies,i,0)
    training_data = np.delete(data,i,0)

    emu = MFE.mf_emulator("box%03d"%i)
    emu.train(training_cosmos,training_data)
    
    for j in xrange(0,N_z):
        MF_data = np.genfromtxt(MF_path+"/building_MF_data/full_mf_data/Box%03d_full/Box%03d_full_Z%d.txt"%(i,i,j))
        lM_bins = MF_data[:,:2]
        N_data = MF_data[:,2]
        cov_data = np.genfromtxt(MF_path+"/building_MF_data/covariances/Box%03d_cov/Box%03d_cov_Z%d.txt"%(i,i,j))
        N_err = np.sqrt(np.diagonal(cov_data))

        n = emu.predict_mass_function(test_cosmo,redshift=redshifts[j],lM_bins=lM_bins)
        N_emu = n*volume

        add_uncertainty = True
        if add_uncertainty:
            for ii in xrange(0,len(N_data)):
                for jj in xrange(0,len(N_data)):
                    cov_data[ii,jj] += per**2 * N_emu[ii] * N_emu[jj]
        
        chi2 = np.dot((N_data-N_emu),np.dot(np.linalg.inv(cov_data),(N_data-N_emu)))
        print "%d,%d chi2 = %f / %d"%(i,j,chi2,len(n))
        chi2s[i,j] = chi2
        Nfp[i,j] = len(n)

        lM = np.log10(np.mean(10**lM_bins,1))
        #visualize.NM_plot(lM,N_data,N_err,lM,N_emu)
np.savetxt("chi2s_p%dpc.txt"%percent,chi2s)
np.savetxt("Nfp.txt",Nfp)


import matplotlib.pyplot as plt
from scipy.stats import chi2
plt.rc('text',usetex=True, fontsize=20)

chi2s = np.loadtxt("chi2s_p%dpc.txt"%percent).flatten()
Nfp = np.loadtxt("Nfp.txt")

plt.hist(chi2s,20,normed=True) #Make the histogram
df = np.mean(Nfp)
mean,var,skew,kurt = chi2.stats(df,moments='mvsk')
x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99,df),100)
plt.plot(x,chi2.pdf(x,df))
plt.xlabel(r"$\chi^2$",fontsize=24)
plt.xlim(0,80)
plt.ylim(0,0.1)
plt.subplots_adjust(bottom=0.15)
plt.show()
