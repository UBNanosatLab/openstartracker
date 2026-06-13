CC = gcc
CXX = g++
LIBPNG = $$(pkg-config --cflags --libs libpng)

all: test test_cpp imgtest pipeline

test: test_beast_c.c ../c_beast/c_beast.h
	$(CC) -std=c99 -Wall -Wextra -Ofast -g -pg test_beast_c.c -o test -lm

test_cpp: test.c
	$(CXX) -std=c++11 -Wall -Ofast -g -pg test.c -o test_cpp -lm

imgtest: ost_imgtest.c ../c_beast/c_beast.h
	$(CC) -std=c99 -Wall -Wextra -O3 ost_imgtest.c -o ost_imgtest $(LIBPNG) -lm

pipeline: ost_pipeline.c ../c_beast/c_beast.h
	$(CC) -std=c99 -Wall -Wextra -O3 ost_pipeline.c -o ost_pipeline $(LIBPNG) -lm
