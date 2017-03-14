import angcorr
import numpy as np
import matplotlib.pyplot as plt

#Load in a power spectrum
k = np.genfromtxt("./test_data/klin.txt")
P = np.genfromtxt("./test_data/plin.txt")

#Define angles
Nt = 1000
theta = np.linspace(0.01,10.0,Nt)

#Call angcorr
w = angcorr.calc_angcorr(theta,k,P)

plt.loglog(theta,w)
plt.show()
