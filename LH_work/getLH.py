import lhsmdu

dim = 2
num = 20
k = lhsmdu.sample(dim, num) #LH sampling
l = lhsmdu.createRandomStandardUniformMatrix(dim, num) #MC sampled
print k.shape
print l.shape

import matplotlib.pyplot as plt
plt.scatter([k[0]], [k[1]], c='b', label="LHS-MDU")
plt.scatter([l[0]], [l[1]], c='r', label="MC")
plt.show()
