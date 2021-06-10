#!/bin/bash

CALIBRATE=0
REGENERATE=0
ESA_TEST=0
IMG_TEST=0

#PYTHON="/usr/bin/python2.7"
PYTHON="/usr/bin/python3.8"

while getopts ":crei" opt; do
  case $opt in
    c)
	  CALIBRATE=1
      ;;
    r)
	  REGENERATE=1
      ;;
    e)
	  ESA_TEST=1
      ;;
    i)
	  IMG_TEST=1
      ;;
   \?)
      echo "Usage: ./unit_test.sh [options] testdir [cmd]"
      echo -e ""
      echo -e "\t-c\tCalibrate based on images in testdir/samples/"
      echo -e "\t-r\tRegenerate ESA test"
      echo -e "\t-e\tRun ESA test"
      echo -e "\t-i\tRun image test"
      echo -e ""
      echo -e "Example cmd:"
	  echo -e "\tmassif-visualizer: valgrind --tool=massif"
	  echo -e "\tkcachegrind: valgrind --tool=cachegrind"
      echo -e ""
      exit
      ;;
  esac
done
shift "$[$OPTIND-1]"

pushd "`dirname $0`">/dev/null

TESTDIR="$1"
if [ ! -d "$TESTDIR" ]; then
	echo "'$TESTDIR' is not a valid directory "
	exit
fi

shift

KILLPID=""
if [[ $ESA_TEST == 1 ]]; then
	make || exit
fi
if [[ $IMG_TEST == 1 ]]; then
	pushd ../beast >/dev/null
	./go || exit
	popd>/dev/null
fi
if [[ $CALIBRATE == 1 ]]; then
	echo "Calibrating..."
	time $PYTHON calibrate.py $TESTDIR || exit
fi
if [[ $REGENERATE == 1 ]]; then
	echo "Regenerating..."
	time $PYTHON simulator.py $TESTDIR/calibration.txt $TESTDIR/input.csv $TESTDIR/result.csv || exit
fi

if [[ $ESA_TEST == 1 ]]; then
	echo "ESA test..."
	make &&
	time $@ ./test $TESTDIR/input.csv $TESTDIR/calibration.txt 1991.25 | tee $TESTDIR/result_real.csv &&
	gprof test | gprof2dot -s | dot -Tpdf -o test.pdf &&
	echo "camera coverage simulation percent:" &&
	echo "100-`diff --suppress-common-lines --speed-large-files -y $TESTDIR/result.csv $TESTDIR/result_real.csv | wc -l`/1" | bc -l &&
	$PYTHON score.py $TESTDIR/result.csv $TESTDIR/result_real.csv 
fi

if [[ $IMG_TEST == 1 ]]; then
	$@ $PYTHON startracker.py $TESTDIR/calibration.txt 1991.25 $TESTDIR/median_image.png &
	KILLPID="$!"
	sleep 10
	#make sure we dont crash when given an image w/ no stars
	echo "rgb.solve_image('$TESTDIR/median_image.png')" | nc -w1 127.0.0.1 8010
	sleep 0.5
	for i in $TESTDIR/samples/*; do
		echo "rgb.solve_image('$i')" | nc -w1 127.0.0.1 8010
		sleep 0.5
		echo "rgb.solve_image('$i')" | nc -w1 127.0.0.1 8010
		sleep 0.5
	done
  #sleep 0.5
  #echo 'exception test' | nc 127.0.0.1 8010
  sleep 0.5
	echo 'quit()' | nc -w1 127.0.0.1 8010
fi
if [ "$KILLPID" != "" ] ; then 
	kill $KILLPID
fi
popd>/dev/null

