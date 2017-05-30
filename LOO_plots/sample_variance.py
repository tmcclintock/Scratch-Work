import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=20)
import cosmocalc as cc

xlabel  = r"$\log_{10}M\ [{\rm M_\odot}/h]$"
ylabel  = r"$b(M, z)G(z)$"

base = "../../Mass-Function-Emulator/test_data/"
datapath = base+"N_data/Box%03d_full/Box%03d_full_Z%d.txt"
covpath  = base+"covariances/Box%03d_cov/Box%03d_cov_Z%d.txt"

cosmologies = np.genfromtxt("cosmos.txt")
scale_factors = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 
                          0.645161, 0.714286, 0.8, 0.909091, 1.0])

N_cosmos = 39
N_z      = len(scale_factors)
c = np.linspace(1.0, 0.0, N_z) #Colors
cmap = plt.get_cmap('seismic')

def get_bG(cosmo_dict, a, Masses):
    return cc.growth_function(a)*np.array([cc.tinker2010_bias(Mi, a, 200) for Mi in Masses])

def plot_tinker_bias(cosmo_dict, i):
    cc.set_cosmology(cosmo_dict)
    
    for j in range(N_z):
        a = scale_factors[j]
        #First plot the data.
        data = np.loadtxt(datapath%(i, i, j))
        lM_bins = data[:,:2]
        lM = np.mean(data[:, :2], 1)
        M = 10**lM
        N = data[:,2]
        cov = np.loadtxt(covpath%(i, i, j))
        err = np.sqrt(np.diagonal(cov))
        bG = get_bG(cosmo_dict, a, M)
        plt.plot(M, bG, c=cmap(c[j]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('log')
    plt.subplots_adjust(bottom=0.15)
    plt.show()


for i in range(0,1):
    #Get in the cosmology and create a cosmo_dict
    num,ombh2,omch2,w0,ns,ln10As,H0,Neff,sigma8 = cosmologies[i]
    h = H0/100.
    Ob,Om = ombh2/(h**2), ombh2/(h**2)+omch2/(h**2)
    cosmo_dict = {"om":Om, "ob":Ob, "ol":1-Om, "ok":0.0, "h":h, 
                  "s8":sigma8, "ns":ns, "w0":w0, "wa":0.0}

    plot_tinker_bias(cosmo_dict, i)
