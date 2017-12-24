import numpy as np
import math

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

IMAGE_VARIANCE=BRIGHT_THRESH
POS_VARIANCE=(gaussian_noise_sigma*180*3600/(np.pi*PIXSCALE))**2
#TODO: image variance set to zero for contest
#Don't factor image variance into estimate of star position error


f_calib=open(CALIBRATION_FILE, 'w')
f_calib.write("IMG_X="+str(IMG_X)+"\n")
f_calib.write("IMG_Y="+str(IMG_Y)+"\n")
f_calib.write("PIXSCALE="+str(PIXSCALE)+"\n")
f_calib.write("DEG_X="+str(DEG_X)+"\n")
f_calib.write("DEG_Y="+str(DEG_Y)+"\n")
f_calib.write("DB_REDUNDANCY=1\n")
f_calib.write("DOUBLE_STAR_PX=3.5\n")
f_calib.write("REQUIRED_STARS=4\n")
f_calib.write("MAX_FALSE_STARS="+str(MAX_FALSE_STARS)+"\n")
f_calib.write("EXPECTED_FALSE_STARS="+str(EXPECTED_FALSE_STARS)+"\n")
f_calib.write("PSF_RADIUS="+str(PSF_RADIUS)+"\n")
f_calib.write("PHOTONS="+str(PHOTONS)+"\n")
f_calib.write("BRIGHT_THRESH="+str(BRIGHT_THRESH)+"\n")
f_calib.write("IMAGE_VARIANCE="+str(IMAGE_VARIANCE)+"\n")
f_calib.write("POS_ERR_SIGMA=2\n")
f_calib.write("POS_VARIANCE="+str(POS_VARIANCE)+"\n")
f_calib.close()
