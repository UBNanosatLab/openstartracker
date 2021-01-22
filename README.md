# openstartracker
A fast, robust, open source startracker based on a new class of baysian startracker algorithms

Features:

* Fast lost in space identification
* Image to image matching
* Collect and store size, shape and color information of unknown objects
* Tracks unknown objects between images
* Programable python frontend / reusable C++ backend (BEAST-2) with no external dependencies 
* Uses astrometry.net for calibration (check if your camera is good enough by uploading your star images to nova.astrometry.net)
* Supports python 2 and 3 (see bottom)

### Basic setup:


##### From a fresh ubuntu linux install
```
sudo apt-get install python-scipy libopencv-dev python-opencv swig python-systemd
```

Additional packages needed for calibration and unit testing:
~~~~
sudo apt-get install git astrometry.net python-astropy

cd /usr/share/astrometry

Download fits files corresponding to your camera fov size (see astrometry.net for details
sudo wget http://data.astrometry.net/4100/index-4112.fits
sudo wget http://data.astrometry.net/4100/index-4113.fits
sudo wget http://data.astrometry.net/4100/index-4114.fits
sudo wget http://data.astrometry.net/4100/index-4115.fits
sudo wget http://data.astrometry.net/4100/index-4116.fits
sudo wget http://data.astrometry.net/4100/index-4117.fits
sudo wget http://data.astrometry.net/4100/index-4118.fits
sudo wget http://data.astrometry.net/4100/index-4119.fits

git clone https://github.com/UBNanosatLab/openstartracker.git

cd openstartracker/tests
./unit_test.sh -crei science_cam_may8_0.05sec_gain40

~~~~
##### To calibrate a new camera:
~~~~
cd openstartracker/
mkdir yourcamera
mkdir yourcamera/samples
mkdir yourcamera/calibration_data
~~~~
add 3-10 star images of different parts of the sky taken with your camera to yourcamera/samples

edit APERTURE and EXPOSURE_TIME in calibrate.py (you want to take images with the lowest exposure time that consistently solves)


run ./unit_test.sh -crei yourcamera to recalibrate and test

The ESA test should have a score of >70. If its worse than this, play around with exposure time (50ms is a good starting point)

##### Python 3 support:

To enable python 3, you will need to edit 2 lines in two files:

beast/Makefile: PYTHONHEADERS=... 

tests/unit_test.sh: PYTHON=...

for python 3 you may need to install the python 3 versions of the dependencies - ie

~~~~
sudo apt-get python3-scipy python3-systemd python3-pip
sudo -H pip3 install opencv-python
sudo -H pip3 install astropy
~~~~

##### Reference frames used:

For RA,DEC,Ori, openstartracker uses the same convention as astrometry.net, where RA and DEC are in the same frame as the star positions specified in the hipparcos catalogue (updated to the current year). Orientation is degrees east of north (ie orientation 0 means that up and down aligns with north-south)

