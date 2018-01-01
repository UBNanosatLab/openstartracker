#!/bin/bash

TESTDIR="$1"
TESTCMD="./test"
if [[ $# -ge 2 ]]; then
	shift
	#example commands:
	#massif-visualizer: valgrind --tool=massif
	#kcachegrind: valgrind --tool=cachegrind
	TESTCMD="$@ $TESTCMD"
fi

pushd "`dirname $0`"
if [ "$TESTDIR" == "xmas" ] ; then
	make &&
	#time python calibrate.py $TESTDIR &&
	#time python3 simulator.py $TESTDIR/calibration.txt $TESTDIR/input.csv $TESTDIR/result.csv&&
	time $TESTCMD $TESTDIR/input.csv $TESTDIR/calibration.txt 1991.25 | tee $TESTDIR/result_real.csv &&
	echo "camera coverage simulation percent:" &&
	echo "100-`diff --suppress-common-lines --speed-large-files -y $TESTDIR/result.csv $TESTDIR/result_real.csv | wc -l`/1" | bc -l &&
	python score.py $TESTDIR/result.csv $TESTDIR/result_real.csv
fi
popd
echo $TESTCMD $TESTDIR/input.csv $TESTDIR/calibration.txt
