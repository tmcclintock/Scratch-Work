"""temp file"""
import numpy as np
import corner as corner
import matplotlib.pyplot as plt
import sys, os

data = np.loadtxt("data.txt")
#fig = corner.corner(data)
#plt.show()

#Expectation values
eZ = np.mean(data,0)
D = data[:] - eZ
d,g = D[:,0],D[:,1]

#Covariance matrix of D
C= np.cov(D,rowvar=False)
#Get the eigenvalues and eigenvectors
w,v = np.linalg.eig(C)
i = np.argmax(w)
v = v[:,i]
x = np.linspace(min(d),max(d),1000)
y =  v[1]*x# + v[0]
print v

plt.scatter(d,g,marker=".",alpha=0.1)
plt.plot(x,y)
plt.show()
