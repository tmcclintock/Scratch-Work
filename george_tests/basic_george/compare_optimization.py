import george
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)

Nx = 40
x = 10*np.sort(np.random.rand(Nx))
yerr = 0.2 * np.ones_like(x)
y = np.cos(x) + 1 + yerr * np.random.randn(len(x))

#kernel = george.kernels.CosineKernel(1)
kernel = george.kernels.ExpSquaredKernel(2*np.pi)
gp = george.GP(kernel)
gp.compute(x, yerr)

t = np.linspace(min(x)-2, max(x)+2, 400)
mu1, cov1 = gp.predict(y, t)
err1 = np.sqrt(np.diag(cov1))

gp.optimize(x, y, yerr)
mu2, cov2 = gp.predict(y, t)
err2 = np.sqrt(np.diag(cov2))

plt.errorbar(x,y,yerr,c='k',ls='')
plt.plot(t, mu1, c='r')
plt.fill_between(t, mu1+err1, mu1-err1, color='r', alpha=0.2)
plt.plot(t, mu2, c='b')
plt.fill_between(t, mu2+err2, mu2-err2, color='b', alpha=0.2)
plt.show()
