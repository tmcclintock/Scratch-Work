"""temp file"""
import numpy as np
import corner as corner
import matplotlib.pyplot as plt
import sys, os

old_labels = [r"$d0$",r"$d1$",r"$f0$",r"$f1$",r"$g0$",r"$g1$"]

data = np.loadtxt("Box034_chain.txt")
#fig = corner.corner(data)
#plt.show()

#Expectation values
eZ = np.mean(data,0)
outD = np.copy(data)
Rs = []

for i in range(0,2):
    D = data[:,i::2]
    print D.shape

    #Covariance matrix of D
    C= np.cov(D,rowvar=False)
    #Get the eigenvalues and eigenvectors
    w,R = np.linalg.eig(C)
    rD = np.dot(D[:],R)
    outD[:,i::2] = rD[:]

fig = corner.corner(outD,labels = old_labels)
plt.show()
