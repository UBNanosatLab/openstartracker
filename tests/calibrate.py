from __future__ import print_function
from os import listdir,system,environ
from os.path import isfile, join
import cv2
import numpy as np
import math
from scipy.stats import poisson
from scipy.optimize import least_squares
from scipy.special import erf
import sys
from astropy.io import fits
from astropy import wcs
from scipy import spatial

## Environment variables:
try: EXPOSURE_TIME = float(environ['EXPOSURE_TIME'])
except KeyError: EXPOSURE_TIME = 0.05 # s
try: APERTURE = float(environ['APERTURE'])
except KeyError: APERTURE = 60.7 # mm
try: DOUBLE_STAR_PX = float(environ['DOUBLE_STAR_PX'])
except KeyError: DOUBLE_STAR_PX = 3.5 #pixels of seperation needed to distinguish stars from each other
try: POS_ERR_SIGMA = float(environ['POS_ERR_SIGMA'])
except KeyError: POS_ERR_SIGMA = 2 #Check all constellations which fall inside these bounds
### Note: Increasing this can actualy reduce the probability of finding a match
### as the true match has to stand out against a larger crowd 


### NOTE: all of the following options multiply runtime by (N+2)^2
try: MAX_FALSE_STARS = int(environ['MAX_FALSE_STARS'])
except KeyError: MAX_FALSE_STARS = 2 #maximum number of objects that can be brighter than the two brightest stars 
try: DB_REDUNDANCY = int(environ['DB_REDUNDANCY'])
except KeyError: DB_REDUNDANCY = 1 #of the brightest DB_REDUNDANCY+2 stars, we need at least 2
try: REQUIRED_STARS = int(environ['REQUIRED_STARS'])
except KeyError: REQUIRED_STARS = 5 #How many stars should we try to match?
### For ultrawide fov this may be set to 3 for faster matching
### For ultranarrow fov, it may be necessary to set this to 5 (Also send me an email and we'll talk)
### TODO: figure out how big "ultrawide" and "ultranarrow" are

def angles2xyz(ra,dec):
	x=np.cos(np.radians(ra))*np.cos(np.radians(dec))
	y=np.sin(np.radians(ra))*np.cos(np.radians(dec))
	z=np.sin(np.radians(dec))
	return list((x,y,z))

#load our star catalog, converting from id,ra,dec to x,y,z,id
def getstardb(year=1991.25,filename="hip_main.dat"):
	yeardiff=year-1991.25
	stardb={}
	starfile = open(filename)
	for line in starfile.readlines():
		fields=line.split('|')
		try:
			HIP_ID=int(fields[1]);
			MAG=float(fields[5]);
			DEC=yeardiff*float(fields[13])/3600000.0 + float(fields[9]);
			cosdec=np.cos(np.pi*DEC/180.0);
			RA=yeardiff*float(fields[12])/(cosdec*3600000.0) + float(fields[8]);
			X=np.cos(np.pi*RA/180.0)*cosdec;
			Y=np.sin(np.pi*RA/180.0)*cosdec;
			Z=np.sin(np.pi*DEC/180.0);
		except ValueError:
			continue
		try:
			f6=int(fields[6])
		except ValueError:
			f6=0
		if (int(fields[29])==0 or int(fields[29])==1) and f6!=3:
			UNRELIABLE=0
		else:
			UNRELIABLE=1
		stardb[HIP_ID]=[HIP_ID,MAG,DEC,RA,X,Y,Z,UNRELIABLE]
	return stardb


def basename(filename):
	if "." in filename:
		filename=".".join(filename.split(".")[0:-1])
	return filename


def ost_gray_from_bgr(img):
	# cv2 stores PNGs as BGR.  The native OST path uses R + 2*G + B.
	return img[:,:,2] + 2.0*img[:,:,1] + img[:,:,0]


def pixel_integrated_gaussian(xs, ys, x0, y0, flux, sigma):
	s=math.sqrt(2.0)*sigma
	ex=erf((xs-x0+0.5)/s)-erf((xs-x0-0.5)/s)
	ey=erf((ys-y0+0.5)/s)-erf((ys-y0-0.5)/s)
	return 0.25*flux*ex*ey


def fit_star_psf_sigma(gray, x0, y0, image_variance, radius, max_sigma):
	h,w=gray.shape
	xi=int(math.floor(x0+0.5))
	yi=int(math.floor(y0+0.5))
	if xi < radius or xi >= w-radius or yi < radius or yi >= h-radius:
		return None
	xgrid,ygrid=np.meshgrid(np.arange(xi-radius, xi+radius+1, dtype=float),
	                     np.arange(yi-radius, yi+radius+1, dtype=float))
	vals=gray[yi-radius:yi+radius+1, xi-radius:xi+radius+1].astype(float)
	border=(np.abs(xgrid-xi)==radius) | (np.abs(ygrid-yi)==radius)
	bg=float(np.median(vals[border]))
	obs=vals-bg
	positive=np.clip(obs, 0.0, None)
	flux0=float(np.sum(positive))
	if flux0 <= image_variance:
		return None
	min_sigma=math.sqrt(1.0/12.0)
	max_sigma=max(max_sigma, min_sigma*1.25)
	sigma0=min(max(0.5, min_sigma*1.05), max_sigma*0.8)
	weight=np.sqrt(np.maximum(np.abs(vals), image_variance))
	p0=np.array([x0, y0, flux0, sigma0], dtype=float)
	lo=np.array([x0-2.0, y0-2.0, 0.0, min_sigma], dtype=float)
	hi=np.array([x0+2.0, y0+2.0, max(flux0*10.0, image_variance*100.0), max_sigma], dtype=float)
	p0=np.minimum(np.maximum(p0, lo+1e-6), hi-1e-6)
	def residual(p):
		return ((pixel_integrated_gaussian(xgrid, ygrid, p[0], p[1], p[2], p[3])-obs)/weight).ravel()
	try:
		res=least_squares(residual, p0, bounds=(lo, hi), max_nfev=100)
	except Exception:
		return None
	if not res.success:
		return None
	sigma=float(res.x[3])
	if sigma <= min_sigma*1.001 or sigma >= max_sigma*0.999:
		return None
	if abs(res.x[0]-x0) > 1.75 or abs(res.x[1]-y0) > 1.75:
		return None
	return sigma


def estimate_psf_sigma(images_by_name, astrometry_results, image_variance):
	# Estimate a calibration-time PSF width from all matched astrometry stars.
	# Each star gets an independent circular, pixel-integrated Gaussian fit in a
	# small median-subtracted aperture.  The production tracker then uses the
	# robust aggregate PSF sigma instead of re-estimating a shared sigma per frame.
	radius=int(environ.get('PSF_SAMPLE_RADIUS', '3'))
	max_sigma=float(environ.get('PSF_MAX_SIGMA', str(max(2.0, DOUBLE_STAR_PX))))
	sigmas=[]
	for name, rows in astrometry_results.items():
		if name not in images_by_name or len(rows)==0:
			continue
		gray=ost_gray_from_bgr(images_by_name[name])
		pos=np.asarray(rows[:,9:11], dtype=float)
		if len(pos)>1:
			tree=spatial.cKDTree(pos)
			dist=tree.query(pos, k=2)[0][:,1]
			isolated=dist > max(2.0*radius+1.0, DOUBLE_STAR_PX)
		else:
			isolated=np.ones(len(pos), dtype=bool)
		for p, ok in zip(pos, isolated):
			if not ok:
				continue
			sigma=fit_star_psf_sigma(gray, p[0], p[1], image_variance, radius, max_sigma)
			if sigma is not None and np.isfinite(sigma):
				sigmas.append(sigma)
	if not sigmas:
		raise RuntimeError("no usable PSF sigma fits")
	sigmas=np.asarray(sigmas, dtype=float)
	med=float(np.median(sigmas))
	mad=float(np.median(np.abs(sigmas-med)))
	if mad > 0:
		keep=np.abs(sigmas-med) <= 3.0*1.4826*mad
		if np.any(keep):
			sigmas=sigmas[keep]
	psf_sigma=float(np.median(sigmas))
	print("PSF_SIGMA: ", psf_sigma, "from", len(sigmas), "matched stars")
	return psf_sigma


#only do this part if we were run as a python script
if __name__ == '__main__':
	samplepath=sys.argv[1]+"/samples"
	image_names = [ f for f in listdir(samplepath) if isfile(join(samplepath,f)) ]
	num_images=len(image_names)
	#NOTE: if you get NoneType error, delete non-png files
	images = np.asarray([cv2.imread( join(samplepath,image_names[n]) ).astype(float) for n in range(0, num_images)])
	median_image=np.median(images,axis=0)
	cv2.imwrite(sys.argv[1]+"/median_image.png",median_image)
	system("md5sum "+samplepath+"/* >"+sys.argv[1]+"/checksum.txt")
	if system("diff -q "+sys.argv[1]+"/checksum.txt "+sys.argv[1]+"/calibration_data/checksum.txt")!=0:
		print ("Clearing old calibration data:")
		system("rm -rfv "+sys.argv[1]+"/calibration_data/* ")
	
	system("mv "+sys.argv[1]+"/checksum.txt "+sys.argv[1]+"/calibration_data/checksum.txt")
		
	stardb=getstardb()
	
	astrometry_results={}
	#filter the background image for astrometry - more important for starfield generator
	for n in range(0, num_images):
		images[n]-=median_image
		image_name=sys.argv[1]+"/calibration_data/"+basename(image_names[n])+".png"
		img=np.clip(images[n],a_min=0,a_max=255).astype(np.uint8)
		cv2.imwrite(image_name,img)
		solve_cmd="solve-field --skip-solved --cpulimit 60 "+image_name
		#solve_cmd="solve-field --skip-solved --cpulimit 60 -v "+image_name #verbose
		print (solve_cmd)
		system(solve_cmd)
		if isfile(basename(image_name)+'.wcs'):
			print ('wcsinfo '+basename(image_name)+'.wcs  | tr [:lower:] [:upper:] | tr " " "=" | grep "=[0-9.-]*$" > '+basename(image_name)+'.solved')
			system('wcsinfo '+basename(image_name)+'.wcs  | tr [:lower:] [:upper:] | tr " " "=" | grep "=[0-9.-]*$" > '+basename(image_name)+'.solved')
			hdulist=fits.open(basename(image_name)+".corr")
			astrometry_results[image_names[n]]=np.array([[i['flux'],i['field_x'],i['field_y'],i['index_x'],i['index_y']]+angles2xyz(i['index_ra'],i['index_dec']) for i in hdulist[1].data])
		
	
	#Use only values below the median for variance calculation.
	#This is equivalent to calculating variance after having filtered out
	#stars and background light
	THRESH_FACTOR=5
	IMAGE_VARIANCE=np.ma.average(images**2,weights=images<0)
	
	bestimage=""
	maxstars=0
	#for stars over 5*IMAGE_VARIANCE, find the corresponding star in the db
	sd = np.array(list(stardb.values()), dtype = object)	#<SB> had to explicitly convert to list for python3
	
	star_kd = spatial.cKDTree(sd[:,4:7])
	for i in astrometry_results:
		astrometry_results[i]=astrometry_results[i][astrometry_results[i][:,0]>IMAGE_VARIANCE*THRESH_FACTOR]
		astrometry_results[i]=np.hstack((sd[star_kd.query(astrometry_results[i][:,5:8])[1]],astrometry_results[i]))
		if len(astrometry_results[i])>maxstars:
			bestimage=i
			maxstars=len(astrometry_results[i])
	images_by_name={image_names[n]: images[n] for n in range(num_images)}
	PSF_SIGMA=estimate_psf_sigma(images_by_name, astrometry_results, IMAGE_VARIANCE)
	astrometry_results_all=np.vstack(list(astrometry_results.values()))
	# Expicitly convert to a float array to prevent numpy error
	astrometry_results_all = astrometry_results_all.astype('float')
	
	#find the dimmest star
	dimmest_match = astrometry_results_all[np.argmax(astrometry_results_all[:,1]),:]

	BASE_FLUX=dimmest_match[8]/pow(10.0,-dimmest_match[1]/2.5)
	print ("BASE_FLUX: ",BASE_FLUX) 
	
	db_img_dist=np.linalg.norm(astrometry_results_all[:,9:11]-astrometry_results_all[:,11:13],axis=1)
	db_img_dist=db_img_dist-IMAGE_VARIANCE/(astrometry_results_all[:,8])
	
	POS_VARIANCE=np.mean(db_img_dist)
	
	#<SB> execfile went away in python3
	#https://stackoverflow.com/questions/6357361/alternative-to-execfile-in-python-3
	filename = sys.argv[1]+"/calibration_data/"+basename(bestimage)+".solved"
	exec(compile(open(filename, "rb").read(), filename, 'exec'))
	
	f_calib=open(sys.argv[1]+"/calibration.txt", 'w')
	f_calib.write("IMG_X="+str(IMAGEW)+"\n")
	f_calib.write("IMG_Y="+str(IMAGEH)+"\n")
	f_calib.write("PIXSCALE="+str(PIXSCALE)+"\n")
	f_calib.write("DB_REDUNDANCY="+str(DB_REDUNDANCY)+"\n")
	f_calib.write("DOUBLE_STAR_PX="+str(DOUBLE_STAR_PX)+"\n")
	f_calib.write("REQUIRED_STARS="+str(REQUIRED_STARS)+"\n")
	f_calib.write("MAX_FALSE_STARS="+str(MAX_FALSE_STARS)+"\n")
	f_calib.write("BASE_FLUX="+str(BASE_FLUX)+"\n")
	f_calib.write("THRESH_FACTOR="+str(THRESH_FACTOR)+"\n")
	f_calib.write("IMAGE_VARIANCE="+str(IMAGE_VARIANCE)+"\n")
	f_calib.write("POS_ERR_SIGMA="+str(POS_ERR_SIGMA)+"\n")
	f_calib.write("POS_VARIANCE="+str(POS_VARIANCE)+"\n")
	f_calib.write("PSF_SIGMA="+str(PSF_SIGMA)+"\n")
	f_calib.write("APERTURE="+str(APERTURE)+"\n")
	f_calib.write("EXPOSURE_TIME="+str(EXPOSURE_TIME)+"\n")
	f_calib.close()
	
	print ("Calibration finished")
	print ("calibration.txt and median_image.png are in "+sys.argv[1]+"\n")
	system("cat "+sys.argv[1]+"/calibration.txt")
