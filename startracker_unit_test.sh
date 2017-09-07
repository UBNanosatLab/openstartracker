#!/bin/bash
./go>/dev/stderr&&
./calibrate.sh >/dev/stderr&&
#python3 example.py
python3 example.orig.py
time ./main >result_real.csv
echo "camera coverage simulation percent:"
echo "100-`diff --suppress-common-lines --speed-large-files -y result.csv result_real.csv | wc -l`/1" | bc -l
