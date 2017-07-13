import numpy as np
from pyhmc import hmc

def logprob(x, inverse_variance):
    logp = -0.5 * np.sum(inverse_variance * x**2)
    gradient = -ivar * x
    return logp, gradient

ivar = 1. / np.random.rand(5)

