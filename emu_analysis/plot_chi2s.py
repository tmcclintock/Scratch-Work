import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
plt.rc('text',usetex=True, fontsize=20)


chi2s = np.loadtxt("chi2s.txt").flatten()
Nfp = np.loadtxt("Nfp.txt")

plt.hist(chi2s,20,normed=True) #Make the histogram
df = np.mean(Nfp)
mean,var,skew,kurt = chi2.stats(df,moments='mvsk')
x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99,df),100)
plt.plot(x,chi2.pdf(x,df))
plt.xlabel(r"$\chi^2$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
