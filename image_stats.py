import numpy as np
import math
#from example.py

# resolution
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

min_true = 0
max_true = 100
min_false = 0
max_false = 80


IMG_X=res_x
IMG_Y=res_y
PIXSCALE=3600*np.rad2deg(2*np.arctan(1/(2*f)))/IMG_X
DEG_X=IMG_X*PIXSCALE/3600
DEG_Y=IMG_Y*PIXSCALE/3600

MAX_FALSE_STARS=max_false
EXPECTED_FALSE_STARS=(min_false+max_false)/2.0 # assume uniform distribution

PSF_RADIUS=sigma_psf
PHOTONS=base_photons * t_exp * np.pi * aperture ** 2  #photons recieved from a mag 0 star
BRIGHT_THRESH=A_pixel+5*sigma_pixel/(math.erf((2*PSF_RADIUS)**-0.5)**2/4)
MIN_MAG=-2.5*math.log10(BRIGHT_THRESH/PHOTONS)+.6 #fudge it
BRIGHT_THRESH=PHOTONS*10**(MIN_MAG/-2.5)


#IMAGE_VARIANCE=sigma_pixel
IMAGE_VARIANCE=0
POS_ERR_SIGMA=2
POS_VARIANCE=(gaussian_noise_sigma*180*3600/(np.pi*PIXSCALE))**2
#TODO: image variance set to zero for contest
#Don't factor image variance into estimate of star position error


print "IMG_X="+str(IMG_X)
print "IMG_Y="+str(IMG_Y)
print "PIXSCALE="+str(PIXSCALE)
print "DEG_X="+str(DEG_X)
print "DEG_Y="+str(DEG_Y)
print "REQUIRED_STARS=4"
print "MAX_FALSE_STARS="+str(MAX_FALSE_STARS)
print "EXPECTED_FALSE_STARS="+str(EXPECTED_FALSE_STARS)
print "PSF_RADIUS="+str(PSF_RADIUS)
print "PHOTONS="+str(PHOTONS)
print "BRIGHT_THRESH="+str(BRIGHT_THRESH)
print "MIN_MAG="+str(MIN_MAG)
print "IMAGE_VARIANCE="+str(IMAGE_VARIANCE)
print "POS_ERR_SIGMA="+str(POS_ERR_SIGMA)
print "POS_VARIANCE="+str(POS_VARIANCE)

#_ = scene.render(False)
