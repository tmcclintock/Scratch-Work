"""
Create the mass function 'mosaic' plot. This is a plot of
N(z,M) for all boxes all overplotted on top of one another.

This is FIGURE 1.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=24)

xlabel  = r"$\log_{10}M\ [{\rm M_\odot}/h]$"
y0label = r"$N/[{\rm Gpc}^3\  \log_{10}{\rm M_\odot}/h]$"
y1label = r"$\%\ {\rm Diff}$"

base = "/home/tmcclintock/Desktop/Github_stuff/Mass-Function-Emulator/test_data/"
datapath = base+"N_data/Box%03d_full/Box%03d_full_Z%d.txt"
covpath  = base+"covariances/Box%03d_cov/Box%03d_cov_Z%d.txt"

N_cosmos = 39
N_z      = 10
#c = np.linspace(0.0, 1.0, N_z) #Colors
#cmap = plt.get_cmap('RdBu') 
c = np.linspace(1.0, 0.0, N_z) #Colors
cmap = plt.get_cmap('seismic') 


fig, axis = plt.subplots(1)

for i in range(N_cosmos):
    for j in range(N_z):
        data = np.loadtxt(datapath%(i, i, j))
        lM = np.mean(data[:, :2], 1)
        N = data[:,2]
        cov = np.loadtxt(covpath%(i, i, j))
        err = np.sqrt(np.diagonal(cov))
        #axis.errorbar(lM, N, err, marker='', c=cmap(c[j]), alpha=0.2)
        axis.plot(lM, N, marker='', c=cmap(c[j]), alpha=0.2)

plt.xlabel(xlabel)
plt.ylabel(y0label)
plt.yscale('log')
plt.subplots_adjust(bottom=0.15, left=0.15)
fig.savefig("fig1_allMFs.png")
plt.show()
