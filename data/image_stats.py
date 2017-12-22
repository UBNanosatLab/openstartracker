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

CALIBRATION_FILE="calibration.txt"

execfile("write_image_stats.py")

#_ = scene.render(False)
