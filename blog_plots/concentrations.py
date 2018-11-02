import numpy as np
import cluster_toolkit as ct
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=14, family='serif')

r = np.logspace(-3, 0) #Mpc
M = 1e13 #Msun
cs = np.array([3,5,7,9])
Omega_m = 0.3

for c in cs:
    rho = ct.density.rho_nfw_at_R(r, M, c, Omega_m)
    plt.loglog(r, rho, label=r"$c=%d$"%c)
    continue
plt.xlabel(r"$r\ [{\rm Mpc}]$")
plt.ylabel(r"$\rho(r,M,c)\ [{\rm M_\odot/Mpc^3}]$")
plt.legend(frameon=False)
plt.savefig("conc_example.png", dpi=300, bbox_inches='tight')
plt.show()
