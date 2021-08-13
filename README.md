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
sudo apt-get install python3-scipy libopencv-dev python3-opencv swig python3-systemd
```

Additional packages needed for calibration and unit testing:
~~~~
sudo apt-get install git astrometry.net python3-astropy

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

#### Using Docker
The included Dockerfile allows the build and test environment to be configured in a reproducible way with all the dependencies.
According to their website,
> Docker is an open platform for developing, shipping, and running applications. Docker enables you to separate your applications from your infrastructure so you can deliver software quickly.

In practice, it's a nice way to ensure you always have the same dependencies installed, independent of the rest of your system.

**Note:** the Docker environment here is best used for development purposed. When running openstartracker on a flight computer, it is most
likely better to run the program directly, without the Docker container. Docker is used here just to create a standard, reproducible dev and testing environment.

To use the Docker environment to run openstartracker,

0. Install [Docker](https://docs.docker.com/get-docker/).
1. Run `source setup.sh`. This will execute a `docker build ...` command to create the Docker environment from the Dockerfile, copy
in all the code config files from this repo, 
and create a command alias called `dstart` that can be run to start the Docker environment and access the pre-built Docker environment.
  - This may take awhile to run the first time as the Docker image is built, but will run much more quickly in subsequent builds
2. Run `dstart` to start and enter the Docker environment.
3. Proceed with the following instructions.

**Important note:** Any code and file changes should be made **outside** the Docker environment -- consider it a 'run only' space.
Any time you make a change to the code in this repo, you'll want to re-run the `source setup.sh` command so that your changes are
included in the Docker build.

### Run on example imagery:
You can use the included `science_cam_may8_0.05sec_gain40` images to test out the calibration process:

```
# if using docker, run `source setup.sh` and then `dstart` to enter the Docker shell
cd tests/
./unit_test.sh -crei science_cam_may8_0.05sec_gain40/
```

This command will **c**alibrate your image sensor, **r** regenerate the test data, run an **E**SA test, and finally run the **i**mage test where images are fed to the calibrated star tracker program to produce an attitude fix.

The usage message for `unit_test.py` is here:
```
./unit_test.sh -h
Usage: ./unit_test.sh [options] testdir [cmd]

	-c	Calibrate based on images in testdir/samples/
	-r	Regenerate ESA test
	-e	Run ESA test
	-i	Run image test
```

### To calibrate a new camera:
1. Create directories for you camera's imagery:
~~~~
cd openstartracker/tests/
mkdir yourcamera
mkdir yourcamera/samples
mkdir yourcamera/calibration_data
~~~~
2. Add 3-10 star images of different parts of the sky taken with your camera to `tests/yourcamera/samples`
3. Edit `APERTURE` and `EXPOSURE_TIME` in `calibrate.py` (you want to take images with the lowest exposure time that consistently solves)
4. (if using docker, run `dstart` to enter the Docker environment)
5. Run 
  ```
  cd tests/
  ./unit_test.sh -crei yourcamera
  ```
  to calibrate and test

The ESA test should have a score of >70. If its worse than this, play around with exposure time (50ms is a good starting point)

##### Reference frames used:

For RA,DEC,Ori, openstartracker uses the same convention as astrometry.net, where RA and DEC are in the same frame as the star positions specified in the hipparcos catalogue (updated to the current year). Orientation is degrees east of north (ie orientation 0 means that up and down aligns with north-south)

