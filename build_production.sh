#!/bin/sh
[ "$1" == "" ] && exit
cd beast
./go&&
cd ..&&
mkdir startracker-production&&
cp -r beast/ startracker-production/&&
cp tests/$1/calibration.txt startracker-production&&
cp tests/$1/median_image.png startracker-production&&
cp hip_main.dat startracker-production&&
cp startracker.py startracker-production/&&
touch startracker-production.tgz&&
rm startracker-production.tgz &&
rsync -av startracker-production/ root@192.168.100.213:~/startracker-production/ &&
rm -rf startracker-production

