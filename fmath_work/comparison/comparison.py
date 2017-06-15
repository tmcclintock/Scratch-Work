"""
A very short code comparing the output of fmath and cmath exp() functions.
"""
import ctypes
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=14)

lib = ctypes.cdll.LoadLibrary("lol.so")
fexp = lib.fexp
sexp = lib.sexp
fexp.argtypes = (ctypes.c_float)
sexp.argtypes = (ctypes.c_float)
fexp.restype  = (ctypes.c_float)
sexp.restype  = (ctypes.c_float)

print fexp(3)
print sexp(3)
