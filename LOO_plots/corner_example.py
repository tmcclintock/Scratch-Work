"""
Make an example corner plot.

This is Figure 3 and also Figure 5 if the rotated version is made.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=24)
import corner
from chainconsumer import ChainConsumer

name = 'dfg_rotated'
rotated = True
corner_labels = []
for i,l in zip(range(len(name)), name.split("_")[0]):
    if rotated:
        corner_labels.append(r"$%s_0'$"%l)
        corner_labels.append(r"$%s_1'$"%l)
    else:
        corner_labels.append(r"$%s_0$"%l)
        corner_labels.append(r"$%s_1$"%l)
print corner_labels
N_parameters = len(corner_labels)

base_dir = "../../fit_mass_functions/output/%s/"%name
base_save = base_dir+"%s_"%name

N_cosmos = 39

#MCMC configuration
nwalkers, nsteps = 16, 5000
nburn = 2000

for i in range(1): #N_cosmos):
    if not rotated:
        fullchain = np.loadtxt(base_dir+"chains/Box%03d_chain.txt"%(i))
    else:
        fullchain = np.loadtxt(base_dir+"rotated_chains/Rotated_Box%03d_chain.txt"%(i))
    chain = fullchain[nwalkers*nburn:]

    #Now with chainconsumer
    fig = ChainConsumer().add_chain(chain, parameters=corner_labels).plot()
    plt.subplots_adjust(bottom=0.15, left=0.15)
    if not rotated:
        fig.savefig("fig3_corner.png")
    else:
        fig.savefig("fig5_Rcorner.png")
    plt.show()
