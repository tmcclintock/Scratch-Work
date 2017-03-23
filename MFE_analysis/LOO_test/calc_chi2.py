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
        
        chi2 = np.dot((N_data-N_emu),np.dot(np.linalg.inv(cov_data),(N_data-N_emu)))
        print "%d,%d chi2 = %f / %d"%(i,j,chi2,len(n))
        chi2s[i,j] = chi2
        Nfp[i,j] = len(n)

        lM = np.log10(np.mean(10**lM_bins,1))
        visualize.NM_plot(lM,N_data,N_err,lM,N_emu)
#np.savetxt("chi2s.txt",chi2s)
#np.savetxt("Nfp.txt",Nfp)
