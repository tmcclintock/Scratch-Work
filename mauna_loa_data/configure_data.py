"""
Rip open the mauna loa data from the maunaloa_c.dat
file and reconfigure so that it's useful.
"""
import numpy as np
import statsmodels.api as sm

#data = sm.datasets.get_rdataset("co2").data
data = np.genfromtxt("monthly_in_situ_co2_mlo.csv", delimiter=',', skip_header=57, usecols=(0,1,4))
print data.shape
print data[:3]
t,m, c = data.T
m = m/12.
t = t + m
good = c>0
t = t[good]
c = c[good]
import matplotlib.pyplot as plt
plt.scatter(t,c,marker='.')
out = np.array([t,c]).T
np.savetxt("co2.txt", out)
plt.show()
