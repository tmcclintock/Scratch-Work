"""
Split the DES SV file into individual files with specific cluster samples.
"""
import fitsio
#fname = "sva1_gold_1.0.2_run_redmapper_v6.3.3_lgt20_catalog.fit"
fname = "redmapper_sva1_public_v6.3_catalog.fits"
data,header = fitsio.read(fname,header=True)

print data.shape
lams = data['LAMBDA']
zs = data['Z_LAMBDA']
print min(zs),max(zs)
print min(lams),max(lams)

zlims = [0.2, 0.4, 0.6, 0.8]
llims = [20, 35, 180]

outpath = "cluster_files/clusters_z%dl%d.txt"
outfiles = []

for i in range(len(zlims)-1):
    outl = []
    for j in range(len(llims)-1):
        outl.append(open(outpath%(i,j+3),"w"))
        outl[j].write("# z lambda\n")
    outfiles.append(outl)

for i in range(len(zlims)-1):
    for j in range(len(llims)-1):
        inds = (zs > zlims[i])*(zs < zlims[i+1])*(lams>=llims[j])*(lams<llims[j+1])
        zij = zs[inds]
        lamsij = lams[inds]
        print "z%dl%d have %d clusters"%(i,j+3,len(zij))
        for k in range(len(zij)):
            outfiles[i][j].write("%f\t%f\n"%(zij[k],lamsij[k]))

for i in range(len(zlims)-1):
    for j in range(len(llims)-1):
        outfiles[i][j].close()
