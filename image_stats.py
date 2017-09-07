import numpy as np
import math
#from example.py
res_x = 1920 # pixels
res_y = 1440 # pixels

# normalized focal length
f = 0.5 / np.tan(np.deg2rad(10) / 2)

# pixel aspect ratio
pixel_ar = 1

# normalized principal point
ppx = 0.5
ppy = 0.5

gaussian_noise_sigma = 20e-6 # rad

cam = 0
# magnitude parameters

A_pixel = 525 # photonelectrons/s mm
sigma_pixel = 525 # photonelectrons/s mm

sigma_psf = 0.5 # pixel
t_exp = 0.2 # s
aperture = 15 # mm

base_photons = 19100 # photoelectrons per mm^2 and second of a magnitude 0 G2 star

magnitude_gaussian = 0.01 # mag
# star count

min_true = 3
max_true = 100
min_false = 0
max_false = 80


IMG_X=res_x
IMG_Y=res_y
PIXSCALE=3600*np.rad2deg(2*np.arctan(1/(2*f)))/IMG_X
DEG_X=IMG_X*PIXSCALE/3600
DEG_Y=IMG_Y*PIXSCALE/3600

ARC_ERR=2.2*gaussian_noise_sigma*180*3600/np.pi #arc seconds
POS_VARIANCE=(ARC_ERR/PIXSCALE)**2
#divide by sqrt(2) to correct for a mistake in the code
ARC_ERR=ARC_ERR/np.sqrt(2)
MAX_FALSE_STARS=max_false
EXPECTED_FALSE_STARS=(min_false+max_false)/2.0 # assume uniform distribution
MIN_MAG=-2.5*math.log10((A_pixel+20*sigma_pixel)/(math.pi*base_photons*t_exp*(aperture*math.erf(2**.5*sigma_pixel**.5))**2))
PSF_RADIUS=sigma_psf

#TODO: image variance set to zero for contest
#Don't factor image variance into estimate of star position error
IMAGE_VARIANCE=0


print "IMG_X="+str(IMG_X)
print "IMG_Y="+str(IMG_Y)
print "DEG_X="+str(DEG_X)
print "DEG_Y="+str(DEG_Y)
print "PIXSCALE="+str(PIXSCALE)
print "ARC_ERR="+str(ARC_ERR)
print "MAX_FALSE_STARS="+str(MAX_FALSE_STARS)
print "EXPECTED_FALSE_STARS="+str(EXPECTED_FALSE_STARS)
print "MIN_MAG="+str(MIN_MAG)
print "PSF_RADIUS="+str(PSF_RADIUS)
print "POS_VARIANCE="+str(POS_VARIANCE)
print "IMAGE_VARIANCE="+str(IMAGE_VARIANCE)
print "PHOTONS="+str(base_photons * t_exp * aperture ** 2 * np.pi)

#_ = scene.render(False)
