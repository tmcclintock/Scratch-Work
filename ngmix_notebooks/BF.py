# necessary imports
import ngmix
from ngmix.observation import Observation, ObsList, MultiBandObsList
from ngmix.fitting import LMSimple
import numpy
from numpy import array
from numpy.random import uniform as urand

# eps is a constant that we need and represents 1% accuracy
eps = 0.01
# Seed the RNG
numpy.random.seed(8381)

gal_jacob = ngmix.UnitJacobian(row = 16.0, col = 16.0)
psf_jacob = ngmix.UnitJacobian(row = 12.0, col = 12.0)

# Object is an exponential disk approximated by gaussians.
pars  = [0.0, 0.0, 0.2, -0.1, 16.0, 100.0]
gmix0 = ngmix.GMixModel(pars, "exp")
# PSF is a single gaussian
psf_pars = [0.0, 0.0, -0.03, 0.02, 4.0, 1.0]
psf_gmix = ngmix.GMixModel(psf_pars, "gauss")
# Convolve the two
gmix = gmix0.convolve(psf_gmix)

dimensions = [32, 32]
image0 = gmix.make_image(dimensions, npoints=10, jacobian=gal_jacob)
psf_dimensions = [24, 24]
psf_image = psf_gmix.make_image(psf_dimensions, npoints=10, jacobian=psf_jacob)
# Add noise to the galaxy image
sigma = 0.01
noise = numpy.random.normal(scale = sigma, size = image0.shape)
image = image0 + noise

# Make an observation of the psf image
psf_obs = Observation(psf_image, jacobian=psf_jacob)
# We use a 'simple' model fit with 6 parameters.
# For simplicity we will guess these parameters before pixelization
pfitter = LMSimple(psf_obs, "gauss")
guess = array(psf_pars)
guess[0] += urand(low=-eps, high=eps)
guess[1] += urand(low=-eps, high=eps)
guess[2] += urand(low=-eps, high=eps)
guess[3] += urand(low=-eps, high=eps)
guess[4] *= (1.0 + urand(low=-eps, high=eps))
guess[5] *= (1.0 + urand(low=-eps, high=eps))
# Kick off the fitter and get out the mixture of the fit
print "pfitter.go start"
pfitter.go(guess)
print "pfitter.go end"
psf_gmix_fit = pfitter.get_gmix()
# Set the mixture to the observation. This is needed for galaxy fitting later.
psf_obs.set_gmix(psf_gmix_fit)

weight = numpy.zeros(image.shape) + 1.0/sigma**2
obs = Observation(image, weight=weight, jacobian = gal_jacob, psf = psf_obs)
fitter = LMSimple(obs, "exp")
guess = array(pars)
guess[0] += urand(low=-eps, high=eps)
guess[1] += urand(low=-eps, high=eps)
guess[2] += urand(low=-eps, high=eps)
guess[3] += urand(low=-eps, high=eps)
guess[4] *= (1.0 + urand(low=-eps, high=eps))
guess[5] *= (1.0 + urand(low=-eps, high=eps))
# kick off the fitter and get out the result of the fit
print "fitter.go start"
fitter.go(guess)
print "fitter.go end"
result = fitter.get_result()

obs_list = ObsList()
obs_list.append(obs) # Object 1
obs_list.append(obs) # Object 2
obs_list.append(obs) # Object 3
fitter = LMSimple(obs_list, "exp")
# Now with multiple bands
mb_obs_list = MultiBandObsList()
mb_obs_list.append(obs_list) #Object list for g
mb_obs_list.append(obs_list) #Object list for r
mb_obs_list.append(obs_list) #Object list for i
mb_obs_list.append(obs_list) #Object list for z
fitter = LMSimple(mb_obs_list, "exp")
