import george
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)

Nx = 40
x = 10*np.sort(np.random.rand(Nx))
yerr = 0.2 * np.ones_like(x)
y = np.cos(x) + 1 + yerr * np.random.randn(len(x))

kernel = george.kernels.ExpSquaredKernel(2*np.pi)
gp = george.GP(kernel)
gp.optimize(x, y, yerr)
t = np.linspace(min(x)-2, max(x)+2, 400)
mu, cov = gp.predict(y, t)
err = np.sqrt(np.diag(cov))
conds = gp.sample_conditional(y, t, size= 50)

plt.errorbar(x,y,yerr,c='k',ls='')
plt.plot(t, mu, c='r')
plt.fill_between(t, mu+err, mu-err, color='r', alpha=0.2)
for i in range(len(conds)):
    plt.plot(t, conds[i], c='b', alpha=0.2)
plt.show()
