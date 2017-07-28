"""
The splines in scipy has a very fast integrate() function.
I want to see if I can get a 1halo lensing signal from this
at a comparable speed to the C code.
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

def rhonfw(R, M, rs):
    print "not implemented yet"
    return 0
