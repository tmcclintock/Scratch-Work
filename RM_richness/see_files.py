"""
Let's start by looking at the stuff in the SDSS file.
"""
import fitsio
import numpy as np
import matplotlib.pyplot as plt

fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
for name in data.dtype.names: print name
z = data['Z_LAMBDA_IN']
print z.shape
