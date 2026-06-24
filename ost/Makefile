.PHONY: all clean

CC = gcc
CFLAGS ?= -O3 -g
LIBPNG_CFLAGS ?= $(shell pkg-config --cflags libpng)
LIBPNG_LIBS ?= $(shell pkg-config --libs libpng)

all: _ost.so

ost.o: ost.c ost.h
	$(CC) -std=c99 -Wall -Wextra $(CFLAGS) -fPIC -c ost.c -o ost.o

png.o: png.c ost.h
	$(CC) -std=c99 -Wall -Wextra $(CFLAGS) $(LIBPNG_CFLAGS) -fPIC -c png.c -o png.o

_ost.so: ost.o png.o
	$(CC) -shared ost.o png.o -o _ost.so $(LIBPNG_LIBS) -lm

clean:
	rm -f _ost.so ost.o png.o
