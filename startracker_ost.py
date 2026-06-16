#!/usr/bin/env python3
"""Startracker command-line frontend using the native OST image pipeline.

Unlike ``startracker.py``, this does not use OpenCV and does not require a
median image.  It reads PNGs with the Python standard library, extracts stars
with ``ost_bg_*``, and matches with the native ``ost`` Python API.
"""

import argparse
import struct
import sys
import zlib

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _paeth(a, b, c):
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def read_png_rgba(filename):
    """Read an 8-bit, non-interlaced PNG and return ``width, height, rgba``."""
    data = open(filename, "rb").read()
    if not data.startswith(PNG_SIGNATURE):
        raise ValueError("%s: not a PNG file" % filename)
    pos = len(PNG_SIGNATURE)
    width = height = bit_depth = color_type = interlace = None
    palette = None
    trans = b""
    compressed = []
    while pos + 8 <= len(data):
        n = struct.unpack(">I", data[pos:pos + 4])[0]
        typ = data[pos + 4:pos + 8]
        chunk = data[pos + 8:pos + 8 + n]
        pos += 12 + n
        if typ == b"IHDR":
            width, height, bit_depth, color_type, _, _, interlace = struct.unpack(">IIBBBBB", chunk)
        elif typ == b"PLTE":
            palette = chunk
        elif typ == b"tRNS":
            trans = chunk
        elif typ == b"IDAT":
            compressed.append(chunk)
        elif typ == b"IEND":
            break
    if bit_depth != 8 or interlace != 0:
        raise ValueError("%s: only 8-bit non-interlaced PNGs are supported" % filename)
    channels = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}.get(color_type)
    if channels is None:
        raise ValueError("%s: unsupported PNG color type %s" % (filename, color_type))
    raw = zlib.decompress(b"".join(compressed))
    stride = width * channels
    out = bytearray(width * height * 4)
    prev = bytearray(stride)
    src = 0
    dst = 0
    for _ in range(height):
        filt = raw[src]
        src += 1
        row = bytearray(raw[src:src + stride])
        src += stride
        for i, val in enumerate(row):
            left = row[i - channels] if i >= channels else 0
            up = prev[i]
            up_left = prev[i - channels] if i >= channels else 0
            if filt == 1:
                val += left
            elif filt == 2:
                val += up
            elif filt == 3:
                val += (left + up) >> 1
            elif filt == 4:
                val += _paeth(left, up, up_left)
            elif filt != 0:
                raise ValueError("%s: unsupported PNG filter %d" % (filename, filt))
            row[i] = val & 255
        if color_type == 0:
            for x in row:
                out[dst:dst + 4] = bytes((x, x, x, 255))
                dst += 4
        elif color_type == 2:
            for i in range(0, stride, 3):
                out[dst:dst + 4] = bytes((row[i], row[i + 1], row[i + 2], 255))
                dst += 4
        elif color_type == 3:
            if palette is None:
                raise ValueError("%s: indexed PNG missing palette" % filename)
            for x in row:
                pi = 3 * x
                alpha = trans[x] if x < len(trans) else 255
                out[dst:dst + 4] = bytes((palette[pi], palette[pi + 1], palette[pi + 2], alpha))
                dst += 4
        elif color_type == 4:
            for i in range(0, stride, 2):
                out[dst:dst + 4] = bytes((row[i], row[i], row[i], row[i + 1]))
                dst += 4
        elif color_type == 6:
            out[dst:dst + 4 * width] = row
            dst += 4 * width
        prev = row
    return width, height, out


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
