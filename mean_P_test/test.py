"""
This is a small test concerning how to take the mean of power spectra
for clusters at different redshifts. For the DES WL analysis,
we simply took the mean redshift of all clusters, and calculated the
power spectrum from CAMB at that redshift to do our model.

There is an alternative, though, which is that we could instead
average the power spectra for each cluster at each of those redshifts.
Mathematically, these two cases are:
P_1 = P(k,<z>); <z> = \sum_{i=0}^{N} z_i

P_2 = \sum_{i=0}^{N} P(k,z_i)

What is the difference between these two models,
and how might these models affect the 2-halo term in
our DeltaSigma profile?
"""
import numpy as np
