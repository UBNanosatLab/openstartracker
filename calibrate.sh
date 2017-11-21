#!/bin/bash
awk -F '|' '
BEGIN {
	deg2rad = atan2(0, -1)/180;
	#yeardiff=strftime("%Y")+(strftime("%m")-1)/12-1991.25
	yeardiff=0
}
{
	HIP_ID=$2;
	MAG=$6;
	DEC=yeardiff*$14/3600000 + $10;
	cosdec=cos(deg2rad*DEC);
	RA=yeardiff*$13/(cosdec*3600000) + $9;
	X=cos(deg2rad*RA)*cosdec;
	Y=sin(deg2rad*RA)*cosdec;
	Z=sin(deg2rad*DEC);
	MAX_90P=$50;
	MIN_90P=$51;
	UNRELIABLE=($30==0||$30==1)&&$7!=3?0:1;
	printf("%d %.12g %.12g %.12g %.12g %d\n",HIP_ID,MAG,X,Y,Z,UNRELIABLE)
}' hip_main.dat >catalog.dat

#this also copies the image for use by astrometry
python image_stats.py > calibration/calibration.txt
#python gendb.py calibration/stars.txt calibration/constellations.txt calibration/dbsize.txt
#sort -n calibration/constellations.txt > calibration/constellations2.txt
#mv calibration/constellations2.txt calibration/constellations.txt
#
#source calibration/calibration.txt
#source calibration/dbsize.txt
#cat calibration/dbsize.txt
#
#PARAMSIZE=$[$NUMCONST*(4*4)]
##4 integers to hold the star ids, and 1 pointer ofset from the base to the next location
#echo -n "calulated size:  "
#echo "$PARAMSIZE beastdb.bin"
#
#./beastgen
#echo -n "actual size:     "
#wc -c beastdb.bin
#
#gzip -f beastdb.bin
#echo -n "conpressed size: "
#wc -c beastdb.bin.gz
