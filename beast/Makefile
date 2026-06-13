CC = gcc
PYTHON ?= python3

all: _beast_py.so

_beast_py.so: beast_py.c beast.h
	$(CC) -std=c99 -Wall -Wextra -Ofast -g -fPIC -shared beast_py.c -o _beast_py.so -lm

clean:
	rm -f _beast_py.so beast.so beast.o beast.pyc beast_wrap.cxx beast_wrap.o
