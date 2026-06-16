#!/bin/sh
set -eu

[ $# -ge 1 ] || { echo "usage: $0 testdir [dest]" >&2; exit 1; }
TESTDIR=$1
DEST=${2:-root@192.168.100.213:~/startracker-production/}
OUT=$(mktemp -d startracker-production.XXXXXX)
trap 'rm -rf "$OUT"' EXIT

make -C ost
cp -r beast ost startracker.py startracker_ost.py hip_main.dat "$OUT"/
cp "tests/$TESTDIR/calibration.txt" "tests/$TESTDIR/median_image.png" "$OUT"/
rsync -av "$OUT"/ "$DEST"
