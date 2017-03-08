"""temp file"""
import numpy as np
import corner as corner
import matplotlib.pyplot as plt
import sys, os

data = np.loadtxt("data.txt")
fig = corner.corner(data)
plt.show()
