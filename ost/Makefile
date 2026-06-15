CC = gcc

all: _ost.so

_ost.so: ost.c ost.h
	$(CC) -std=c99 -Wall -Wextra -Ofast -g -fPIC -shared ost.c -o _ost.so -lm

clean:
	rm -f _ost.so
