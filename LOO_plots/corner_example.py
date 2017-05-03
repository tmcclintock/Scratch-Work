"""
Make an example corner plot.

This is Figure 3.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True, fontsize=24)
import corner
from chainconsumer import ChainConsumer

name = 'dfg'
corner_labels = [r"$d0$",r"$d1$",r"$e0$",r"$e1$",
                 r"$f0$",r"$f1$",r"$g0$",r"$g1$"]
if name is 'dfg':
    corner_labels = [r"$d0$",r"$d1$",r"$f0$",r"$f1$",r"$g0$",r"$g1$"]
if name is 'deg':
    corner_labels = [r"$d0$",r"$d1$",r"$e0$",r"$e1$",r"$g0$",r"$g1$"]
N_parameters = len(corner_labels)

base_dir = "../../fit_mass_functions/output/%s/"%name
base_save = base_dir+"%s_"%name

N_cosmos = 39

#MCMC configuration
nwalkers, nsteps = 16, 5000
nburn = 1000

for i in range(1): #N_cosmos):
    fullchain = np.loadtxt(base_dir+"chains/Box%03d_chain.txt"%(i))
    chain = fullchain[nwalkers*nburn:]

    #First with corner
    #fig = corner.corner(chain, labels=corner_labels, plot_datapoints=False)
    #plt.subplots_adjust(bottom=0.15, left=0.15)
    #plt.show()

    #Now with chainconsumer
    fig = ChainConsumer().add_chain(chain, parameters=corner_labels).plot()
    plt.subplots_adjust(bottom=0.15, left=0.15)
    fig.savefig("fig3_corner.png")
    plt.show()
