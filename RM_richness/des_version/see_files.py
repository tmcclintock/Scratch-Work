"""
Let's start by looking at the stuff in the SDSS file.
"""
import fitsio
import numpy as np
import matplotlib.pyplot as plt

fname = "y1a1_gold_1.0.3_wide+d10-mof-001b_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
for name in data.dtype.names: print name
z = data['Z_LAMBDA_IN']
zt = data['Z_LAMBDA']
print z.shape, zt.shape
print min(z),max(z)
print min(zt),max(zt)
