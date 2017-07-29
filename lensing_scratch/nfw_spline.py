"""
The splines in scipy has a very fast integrate() function.
I want to see if I can get a 1halo lensing signal from this
at a comparable speed to the C code.
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

def xi_nfw(R, M, c):
    R200 = (M/(4./3.*np.pi*rhom*200))**.33333333333
    Rs = R200/c
    fc = np.log(1+c)-c/(1+c)
    xi = M/(4*np.pi*Rs**3*fc)/(R/Rs*(1+R/Rs)*(1+R/Rs))/rhom - 1.
    print "not implemented yet"
    return 0
