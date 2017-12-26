from os import listdir
from os.path import isfile, join
import cv2
import numpy as np
import math
from scipy.stats import poisson
import sys

mypath=sys.argv[1]+"/samples/"
image_names = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
num_images=len(image_names)

images = np.asarray([cv2.imread( join(mypath,image_names[n]) ).astype(float) for n in range(0, num_images)])
median_image=np.median(images,axis=0)
cv2.imwrite(sys.argv[1]+"/median_image.png",median_image)

#filter the background image for astrometry - more important for starfield generator
for n in range(0, num_images):
	cv2.imwrite(sys.argv[1]+"/processed/"+image_names[n],np.clip(images[n]-median_image,a_min=0,a_max=255).astype(np.uint8))

#stars arent that big a deal, but leaving them out can't hurt (2% more accurate)
IMAGE_VARIANCE=np.ma.average((images-median_image)**2,weights=images<median_image)

#set a threshold so that we have a PROB_FALSE_STAR chance of a star apearing above the threshold per image
print "IMAGE_VARIANCE="+str(IMAGE_VARIANCE)
