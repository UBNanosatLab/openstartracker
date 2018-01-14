# openstartracker
A fast, robust, open source startracker based on a new class of baysian startracker algorithms

Features:

* Fast lost in space identification
* Image to image matching
* Collect and store size, shape and color information of unknown objects
* Tracks unknown objects between images
* Programable python frontend / reusable C++ backend (BEAST-2) with no external dependencies 

### Basic setup:

##### From a fresh xubuntu 16.04 linux install
Development computer setup:
~~~~
sudo apt-get install git libvte-dev libtool ctags gdb meld nemiver libwebkit-dev vim geany geany-plugins astrometry.net python-astropy  g++-multilib

Note: The following also needed for actually building (ie. on flight computer)
sudo apt-get install python-scipy libopencv-dev python-opencv swig
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
./unit_test.sh -crei xmas
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
