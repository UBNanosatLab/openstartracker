#!/bin/bash

CALIBRATE=0
REGENERATE=0
ESA_TEST=0
IMG_TEST=0

#PYTHON="/usr/bin/python2.7"
PYTHON="${PYTHON:-$(command -v python3)}"
if [ -z "$PYTHON" ]; then
	echo "python3 not found"
	exit 1
fi

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
	if [ -f ost_imgtest.c ]; then
		make imgtest pipeline test || exit
	else
		pushd ../beast >/dev/null
		./go || exit
		popd>/dev/null
	fi
fi
if [[ $CALIBRATE == 1 ]]; then
	echo "Calibrating..."
	if command -v solve-field >/dev/null && $PYTHON - <<'PY' >/dev/null 2>&1
import astropy
import cv2
PY
	then
		time $PYTHON calibrate.py $TESTDIR || exit
	elif [ -f "$TESTDIR/calibration.txt" ] && [ -f "$TESTDIR/median_image.png" ]; then
		echo "Calibration dependencies unavailable; using existing calibration fixture."
	else
		echo "Calibration dependencies unavailable and no calibration fixture exists."
		exit 1
	fi
fi
if [[ $REGENERATE == 1 ]]; then
	echo "Regenerating..."
	time $PYTHON simulator.py $TESTDIR/calibration.txt $TESTDIR/input.csv $TESTDIR/result.csv || exit
fi

if [[ $ESA_TEST == 1 ]]; then
	echo "ESA test..."
	make &&
	time $@ ./test $TESTDIR/input.csv $TESTDIR/calibration.txt 1991.25 | tee $TESTDIR/result_real.csv || exit
	if [ -x ./test_cpp ]; then
		$@ ./test_cpp $TESTDIR/input.csv $TESTDIR/calibration.txt 1991.25 > $TESTDIR/result_cpp.csv 2>/dev/null || exit
		if [ "${OST_STRICT_CPP:-0}" = "1" ]; then
			diff -q $TESTDIR/result_cpp.csv $TESTDIR/result_real.csv >/dev/null || {
				echo "C BEAST regression: output differs from C++ test_cpp"
				exit 1
			}
		elif [ -f "$TESTDIR/result.csv" ]; then
			$PYTHON - "$TESTDIR/result.csv" "$TESTDIR/result_cpp.csv" "$TESTDIR/result_real.csv" <<'PY' || exit
import sys

def score(expected_path, actual_path):
    expected = open(expected_path).read().splitlines()
    actual = open(actual_path).read().splitlines()
    total = 0.0
    for e, a in zip(expected, actual):
        ee = e.split(",")
        aa = a.split(",")
        t = c = w = 0
        for x, y in zip(ee, aa):
            if int(x) != -1:
                t += 1
            if int(y) != -1:
                if x == y:
                    c += 1
                else:
                    w += 1
        total += max((c - 2 * w) / float(max(t, 1)), -1) if t else 0.0
    return total

expected, cpp, c = sys.argv[1:]
cpp_score = score(expected, cpp)
c_score = score(expected, c)
if c_score + 1e-12 < cpp_score:
    print("C BEAST regression: C score %.12g < C++ score %.12g" % (c_score, cpp_score))
    sys.exit(1)
PY
		fi
	fi
	./test --relative-self $TESTDIR/input.csv $TESTDIR/calibration.txt 1991.25 >/dev/null || {
		echo "C BEAST relative star self-test failed"
		exit 1
	}
	if command -v gprof2dot >/dev/null && command -v dot >/dev/null; then
		gprof test | gprof2dot -s | dot -Tpdf -o test.pdf
	fi
	echo "camera coverage simulation percent:" &&
	echo "100-`diff --suppress-common-lines --speed-large-files -y $TESTDIR/result.csv $TESTDIR/result_real.csv | wc -l`/1" | bc -l &&
	$PYTHON score.py $TESTDIR/result.csv $TESTDIR/result_real.csv 
fi

if [[ $IMG_TEST == 1 ]]; then
	if [ -x ./ost_imgtest ]; then
		PIPE_OUT=""
		PIPE_STARS=""
		PIPE_UNIT=""
		# Make sure we do not crash when given an image with no stars,
		# then run each sample twice like the original socket test.
		$@ ./ost_imgtest $TESTDIR/calibration.txt $TESTDIR/median_image.png $TESTDIR/median_image.png >/dev/null || exit
		for i in $TESTDIR/samples/*; do
			$@ ./ost_imgtest $TESTDIR/calibration.txt $TESTDIR/median_image.png "$i" "$i" >/dev/null || exit
		done
		if [ -x ./ost_pipeline ] && [ -x ./test ]; then
			PIPE_OUT="$(mktemp)"
			PIPE_STARS="$(mktemp)"
			PIPE_UNIT="$(mktemp)"
			$@ ./ost_pipeline --stars-out "$PIPE_STARS" \
				$TESTDIR/calibration.txt 1991.25 $TESTDIR/samples/* \
				> "$PIPE_OUT" || {
				rm -f "$PIPE_OUT" "$PIPE_STARS" "$PIPE_UNIT"
				exit 1
			}
			$@ ./test "$PIPE_STARS" $TESTDIR/calibration.txt 1991.25 \
				> "$PIPE_UNIT" 2>/dev/null || {
				rm -f "$PIPE_OUT" "$PIPE_STARS" "$PIPE_UNIT"
				exit 1
			}
			diff -q "$PIPE_OUT" "$PIPE_UNIT" >/dev/null || {
				echo "C image pipeline regression: tracker API output differs from unit-test executable"
				rm -f "$PIPE_OUT" "$PIPE_STARS" "$PIPE_UNIT"
				exit 1
			}
			rm -f "$PIPE_OUT" "$PIPE_STARS" "$PIPE_UNIT"
		fi
	else
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
fi
if [ "$KILLPID" != "" ] ; then 
	kill $KILLPID
fi
popd>/dev/null
