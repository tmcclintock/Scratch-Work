import numpy as np
import matplotlib.pyplot as plt

mu, sigma = 3., 1. # mean and standard deviation
s = np.random.lognormal(mu, sigma, 1000)
fig, axes = plt.subplots(2, sharey= True)
count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')

x = np.linspace(min(bins), max(bins), 10000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))/ (x * sigma * np.sqrt(2 * np.pi)))

for i in range(len(axes)):
    axes[i].plot(x, pdf, linewidth=2, color='r')
    count, bins, ignored = axes[i].hist(s, 100, normed=True, align='mid')
    plt.axis('tight')
axes[1].set_xscale('log')
axes[1].set_xlim(0.4, 700)
plt.show()
