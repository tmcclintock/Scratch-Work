import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.integrate import quad
import aemulus_extras as ae
import aemulus_data as AD

delta = 200. #overdensity
rhocrit = 2.77533742639e+11 # h^2 Msun/Mpc^3

Volume = 1.05**3 #[Mpc/h]^3 #Aemulus volume
box = 0
snap = 9 #z=0
extra = ae.Extras(0)
cos = AD.building_box_cosmologies()[box]
h = cos[5]/100.
Omega_m = (cos[0] + cos[1])/h**2
rhomdelta = rhocrit * Omega_m * delta

M = extra.M
lM = np.log(M)
dndlM = extra.dndlM[snap]

spl = IUS(lM, dndlM)

def integrand(lnM):
    M = np.exp(lnM)
    return spl(lnM) * M

M_above = np.zeros_like(M)
vol_in = np.zeros_like(M_above)
frac = np.zeros_like(M_above)
lMmax = lM[-1]
for i, Mmin in enumerate(M):
    lMmin = lM[i]
    M_above[i] = Volume * quad(integrand, lMmin, lMmax)[0]
    vol_in[i] = M_above[i] / rhomdelta
    frac[i] = vol_in[i] / Volume

import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family='serif', size=14)
plt.plot(M, frac)
plt.xscale('log')
plt.ylabel(r"$f_{\rm interior}$")
plt.xlabel(r"$M_{\rm min}\ [h^{-1}{\rm M}_\odot]$")
plt.title("Volume fraction inside halos")
plt.show()
