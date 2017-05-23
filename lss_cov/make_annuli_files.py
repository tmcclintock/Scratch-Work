import numpy as np
import sys
from astropy.cosmology import FlatLambdaCDM
h = 0.70
cosmo = FlatLambdaCDM(H0=h*100, Om0=0.3)

arcmin_per_rad = 3437.75

"""
The annuli file is formatted like:
Nbins
bin0_left bin0_right
bin1_left bin1_right
...
"""
Rmin = 0.0323 #Mpc physical
Rmax = 30.0 #Mpc physical
Nbins = 15
R = np.logspace(np.log10(Rmin), np.log10(Rmax), Nbins+1)

meanz = np.loadtxt("meanz.txt")
print meanz.shape
print meanz

for i in range(len(meanz)):
    for j in range(len(meanz[i])):
        z = meanz[i,j]
        annuli = R/(cosmo.angular_diameter_distance(z).value/arcmin_per_rad)
        fout = open("thetas_ds_z%d_l%d.txt"%(i,j), "w")
        fout.write("%d\n"%(len(annuli)-1))
        for k in range(len(annuli)-1):
            fout.write("%.4f %.4f\n"%(annuli[k], annuli[k+1]))
        fout.close()
        
