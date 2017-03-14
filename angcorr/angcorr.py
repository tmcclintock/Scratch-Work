import numpy as np
import os, inspect
from ctypes import c_double,c_int,POINTER,cdll
#Find the library no matter where we are.
sopath = os.path.join(os.path.dirname(__file__),"_angcorr.so")

def calc_angcorr(theta,k,P):
    theta = theta.copy()
    k = k.copy()
    P = P.copy()
    cclib = cdll.LoadLibrary(sopath)
    cac = cclib.ang_corr
    cac.restype = c_int

    """Args for the c code:
    k, P, Nk, theta, w, Ntheta
    """
    Nk = len(k)
    if Nk != len(P):
        raise Exception("len(k)!=len(P)")
    Nt = len(theta)
    cac.argtypes=[POINTER(c_double), POINTER(c_double), c_int, POINTER(c_double), POINTER(c_double), c_int]

    k_in = k.ctypes.data_as(POINTER(c_double))
    P_in = P.ctypes.data_as(POINTER(c_double))
    theta_in = theta.ctypes.data_as(POINTER(c_double))

    w = np.zeros(Nt)
    w_in = w.ctypes.data_as(POINTER(c_double))

    result = cac(k_in, P_in, Nk, theta_in, w_in, Nt)
    if result != 0:
        raise Exception("Error in angcorr.py")
    return w
