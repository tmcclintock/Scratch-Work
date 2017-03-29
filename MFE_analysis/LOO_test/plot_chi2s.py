import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
plt.rc('text',usetex=True, fontsize=20)

#The percentage of the model we add to the uncertainty
percent = 5
per = percent/100.

chi2s = np.loadtxt("chi2s_p%dpc.txt"%percent).flatten()
Nfp = np.loadtxt("Nfp.txt")

plt.hist(chi2s,20,normed=True) #Make the histogram
df = np.mean(Nfp)
mean,var,skew,kurt = chi2.stats(df,moments='mvsk')
x = np.linspace(chi2.ppf(0.01,df),chi2.ppf(0.99,df),100)
plt.plot(x,chi2.pdf(x,df))
plt.xlabel(r"$\chi^2$",fontsize=24)
plt.xlim(0,80)
plt.ylim(0,0.1)
plt.subplots_adjust(bottom=0.15)
plt.show()
