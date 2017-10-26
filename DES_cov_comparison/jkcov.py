import numpy as np
import matplotlib.pyplot as plt
import matplotlib

C = np.loadtxt("/Users/tmcclintock/Data/DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1subtr_l3_z0_dst_cov.dat")

x = np.linspace(0.0323, 30, 16)
extent = np.log10([0.0323, 30, 0.0323, 30])
fig = plt.figure(facecolor='white')
ax = fig.add_subplot(111)
ax.imshow(np.abs(C[2:, 2:]), aspect='equal', interpolation='none',
           extent=extent, origin='lower', cmap="plasma", vmin=0.01, vmax=400, norm=matplotlib.colors.LogNorm())
ax.set_xlabel("$log\mathrm{R}$ [Mpc]")
ax.set_ylabel("$log\mathrm{R}$ [Mpc]")
ax.set_xticklabels([])
ax.set_yticklabels([])
fig.savefig("testcov.png", bbox_inches="tight", dpi=300)
