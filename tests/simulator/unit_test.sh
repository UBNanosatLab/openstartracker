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
fi
popd
echo $TESTCMD $TESTDIR/input.csv $TESTDIR/calibration.txt
