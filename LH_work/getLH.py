import lhsmdu

import numpy as np
np.random.seed(12345)
dim = 1
num = 20
k = lhsmdu.sample(dim, num) #LH sampling
k2 = lhsmdu.sample(dim, num) #LH sampling
l = lhsmdu.createRandomStandardUniformMatrix(dim, num) #MC sampled
print k.shape
print l.shape

def llike(x):
    x = np.array(x)
    return np.exp(-(x-0.5)**2/0.05)


import matplotlib.pyplot as plt
plt.plot(k[0], llike(k[0]), c='b', ls='', marker='.')
plt.plot(l[0], llike(l[0]), c='r', ls='', marker='.')
#plt.scatter([k[0]], [k[1]], c='b', label="LHS-MDU")
#plt.scatter([k2[0]], [k2[1]], c='g', label="LHS-MDU")
#plt.scatter([l[0]], [l[1]], c='r', label="MC")
plt.show()
