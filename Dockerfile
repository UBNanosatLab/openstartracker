
FROM ubuntu:20.04

# deal with timezones (change to what suits your use case)
ENV TZ=America/New_York

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR home

# need all of this to work (was missing some packages)
RUN apt-get update \
&&  apt-get -y install software-properties-common \
&&  add-apt-repository ppa:deadsnakes/ppa \
&&  apt-get install -y python3-scipy libopencv-dev python3-opencv \
                     swig python3-systemd git astrometry.net \
                     python3-astropy python3-pkgconfig \
                     python3-dev libpython3.7-dev libpython3.8-dev \
                     build-essential python3.6 libpython3.6-dev \
                     python3-pandas \
                     graphviz \
                     bc \
                     netcat \
                     tzdata \
&& rm -rf /var/lib/apt/lists/* 

ADD http://data.astrometry.net/4100/index-4112.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4113.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4114.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4115.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4116.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4117.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4118.fits /usr/share/astrometry/
ADD http://data.astrometry.net/4100/index-4119.fits /usr/share/astrometry/
