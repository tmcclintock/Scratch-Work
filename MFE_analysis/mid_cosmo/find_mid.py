import numpy as np

data = np.genfromtxt("../test_data/building_cosmos.txt")
print data.shape

#boxnum ombh2 omch2 w0 ns ln10As H0 Neff sigma8
data = np.delete(data,5,1) #Delete ln10As
data = np.delete(data,0,1) #Delete boxnum
data = np.delete(data,-1,0)#39 is broken
fulldata = np.copy(data)

ind_order = [0,6,1,4,2,3,5]
steps = [5,4,3,2,2,1,1]

for i in range(0,len(ind_order)):
    step = steps[i]
    print step, " is step"
    ind = ind_order[i]
    par = data[:,i]
    indsort = np.argsort(par)
    data = data[indsort[step:-step]]
    print data.shape
best = data[1]
print best
print np.argwhere(fulldata[:,0] == best[0]) # I get 34
print fulldata[34]

#Thus, cosmo 34 is the best one, since it is the "middle most" cosmology,
#or approximately is.
