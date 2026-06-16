#!/usr/bin/env python3
"""Startracker command-line frontend using the native OST image pipeline.

Unlike ``startracker.py``, this does not use OpenCV and does not require a
median image.  It reads PNGs with the native OST/libpng wrapper, extracts stars
with ``ost_bg_*``, and matches with the native ``ost`` Python API.
"""

import argparse
import sys

def read_png_rgba(filename):
    """Read a PNG via the native OST/libpng wrapper."""
    import ost

    return ost.read_png_rgba(filename)


def write_stars(f, stars):
    vals = []
    for s in stars:
        vals.extend(s)
    f.write(",".join("%.17g" % v for v in vals) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", default="hip_main.dat")
    parser.add_argument("--stars-out")
    parser.add_argument("calibration")
    parser.add_argument("year", type=float)
    parser.add_argument("images", nargs="+")
    args = parser.parse_args(argv)

    import ost

    cfg = ost.load_config(args.calibration)
    tracker = ost.Tracker(cfg).prepare_catalog(args.catalog, args.year)
    pipeline = ost.ImagePipeline(cfg)
    stars_file = open(args.stars_out, "w") if args.stars_out else None
    try:
        for image in args.images:
            width, height, rgba = read_png_rgba(image)
            if width != cfg.IMG_X or height != cfg.IMG_Y:
                raise ValueError("%s: got %dx%d, expected %dx%d" %
                                 (image, width, height, cfg.IMG_X, cfg.IMG_Y))
            stars, info = pipeline.measure_rgba(rgba)
            if stars_file:
                write_stars(stars_file, stars)
            ids = tracker.match_catalog_stars(stars)
            print(",".join(str(i) for i in ids))
            print("%s: components=%d fitted=%d used=%d dropped=%d sigma=%.9g" %
                  (image, info["components"], info["fitted"], info["used"],
                   info["dropped"], info["sigma"]), file=sys.stderr)
    finally:
        if stars_file:
            stars_file.close()


if __name__ == "__main__":
    main()
