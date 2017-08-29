import numpy as np
import math

# magnitude parameters (from example.py)

A_pixel = 525 # photonelectrons/s mm
sigma_pixel = 525 # photonelectrons/s mm
sigma_psf = 0.5 # pixel
t_exp = 0.2 # s
aperture = 15 # mm
base_photons = 19100 # photoelectrons per mm^2 and second of a magnitude 0 G2 star

gaussian_noise_sigma = 20e-6 # rad

# image parameters
IMG_X=1920
IMG_Y=1440
DEG_DIAG=10
PIXSCALE=3600*DEG_DIAG/np.sqrt(IMG_X**2+IMG_Y**2)
ARC_ERR=gaussian_noise_sigma*180*3600/(np.pi*PIXSCALE) #px
MAX_FALSE_STARS=10
MIN_MAG=-2.5*math.log10((A_pixel+20*sigma_pixel)/(math.pi*base_photons*t_exp*(aperture*math.erf(2**.5*sigma_pixel**.5))**2))

print "IMG_X="+str(IMG_X)
print "IMG_Y="+str(IMG_Y)
print "DEG_DIAG="+str(DEG_DIAG)
print "PIXSCALE="+str(PIXSCALE)
print "ARC_ERR="+str(ARC_ERR)
print "MAX_FALSE_STARS="+str(MAX_FALSE_STARS)
print "MIN_MAG="+str(MIN_MAG)
print "PSF_RADIUS="+str(sigma_psf)
