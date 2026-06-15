FROM ubuntu:20.04

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /home

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      astrometry.net bc build-essential git graphviz libopencv-dev libpng-dev \
      pkg-config python3-astropy python3-opencv python3-pandas \
      python3-pkgconfig python3-scipy python3-systemd tzdata \
 && rm -rf /var/lib/apt/lists/*

ADD http://data.astrometry.net/4100/index-4112.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4113.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4114.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4115.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4116.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4117.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4118.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4119.fits /usr/share/astrometry/
