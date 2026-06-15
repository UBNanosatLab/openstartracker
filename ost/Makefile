.PHONY: all clean

CC = gcc
CFLAGS ?= -Ofast -g

all: _ost.so

ost.o: ost.c ost.h
	$(CC) -std=c99 -Wall -Wextra $(CFLAGS) -fPIC -c ost.c -o ost.o

_ost.so: ost.o
	$(CC) -shared ost.o -o _ost.so -lm

clean:
	rm -f _ost.so ost.o
