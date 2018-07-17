import numpy as np
import sys, os
import pandas as pd

inbase = "/nfs/slac/g/ki/ki22/cosmo/tmcclint/NM_Z_data/reduced_mf_data/Box%03d_reduced/Box%03d_reduced_Z%d.list"

outbase = "/nfs/slac/g/ki/ki22/cosmo/tmcclint/NEWMF/Box%03d_MF/"

sfs = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0])

def run(box, snap):
    
    try:
        data = pd.read_csv_path(inbase%(box, box, snap), dtype='float64', header=0, delim_whitespace=True)
        data = data.as_matrix()
        M = data[:,2]
        x, y, z = data[:, 8:11]
        #x, y, z, M = np.loadtxt(inbase%sf, unpack=True)
    except ValueError:
        return
    os.system("mkdir -p %s"%(outbase%box))

    SIDE_LENGTH = 1050. #Mpc/h side of the whole simulation

    lM = np.log10(M)
    print "log bin edges:",min(lM), max(lM)
    Nbins = 15
    edges = np.linspace(min(lM), 16.5, Nbins+1)[1:] #only right edges
    lMi = np.digitize(lM, edges, True)
    N = np.bincount(lMi)
    Mave = np.ones_like(N)
    for i in range(len(Mave)):
        Mhere = M[lMi==i]/N[i]
        Mhere = Mhere.astype("float64")
        Mave[i] = np.sum(Mhere)

    print "N truth:", N
    print "total number ", N.sum()
    print "total input ", M.size

    Ndivs = 8
    Njk = Ndivs*Ndivs*Ndivs
    print "N per jk:",N/float(Njk)
    Vfactor = Njk/(Njk-1.) #gotta increase everything by this when we do the LOO
    L = SIDE_LENGTH/Ndivs #Mpc/h per side of jk subregion
    jk = np.floor(x/L) + Ndivs*np.floor(y/L) + Ndivs*Ndivs*np.floor(z/L)
    jk = jk.astype(int)
    Nsub = np.ones((Njk, Nbins))
    print "Jk facts: ",Njk, min(jk), max(jk)

    #Find the number in each subregion
    for i in range(Njk):
        lMjk = lM[jk==i] #the halos in the ith subregion
        lMj = np.digitize(lMjk, edges, True)
        Nsub[i] = np.bincount(lMj, minlength=len(N))#number in this subregion

    print "Resum from subs: ",np.sum(Nsub, axis=0).astype(int)
    Nloo = np.ones_like(Nsub)
    Nloo[:] *= N
    Nloo -= Nsub

    Nloo *= Vfactor #Rescale for the appropriate volume
    Nmean = np.mean(Nloo, axis=0)
    print "Nmean - N: ",Nmean - N
    X = Nloo[:] - Nmean

    C = np.zeros((Nbins, Nbins))
    for i in range(Nbins):
        for j in range(Nbins):
            C[i,j] = np.sum(X[:,i]*X[:,j])

    err = np.sqrt(np.diag(C))

    edges = np.linspace(min(lM), max(lM), Nbins+1)
    lMbins = (edges[1:] + edges[:-1])/2.
    out = np.zeros((Nbins, 4 ))

    for i in range(Nbins):
        out[i, 0] = edges[i]
        out[i, 1] = edges[i+1]
        out[i, 2] = N[i]
        out[i, 3] = Mave[i]
    
    np.savetxt("Box%03d_Z%d.txt"%(box, snap), out, header = "logM_lo logM_hi N Mave")
    np.savetxt("Box%03d_Z%d_cov.txt"%(box, snap), C)

if __name__ == "__main__":
    box, snap = 0,0
    run(box, sf)
