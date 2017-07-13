import numpy as np
from pyhmc import hmc

def logprob(x, ivar):
    logp = -0.5 * np.sum(ivar * x**2)
    grad = -ivar * x
    return logp, grad

#ivar = 1. / np.random.rand(5)
#samples = hmc(logprob, x0=np.random.randn(5), args=(ivar,), n_samples=1e4)

ivar = 1. / np.random.rand(5)
samples = hmc(logprob, x0=np.random.randn(5), args=(ivar,), n_samples=10000)

import corner
import matplotlib.pyplot as plt
fig = corner.corner(samples)
plt.show()
