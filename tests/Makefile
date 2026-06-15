CC = gcc
CXX = g++
LIBPNG = $$(pkg-config --cflags --libs libpng)

all: test test_cpp test-beast imgtest pipeline

test: test.c ../ost/ost.h ost_test_work.h
	$(CC) -std=c99 -Wall -Wextra -Ofast -g -pg test.c -o test -lm

test_cpp: test.c ../ost/ost.h ost_test_work.h
	$(CXX) -std=c++11 -Wall -Wextra -Ofast -g -pg test.c -o test_cpp -lm

test-beast: test-beast.cpp ../beast/beast.h ../ost/ost.h
	$(CXX) -std=c++11 -Wall -Wextra -Ofast -g -pg test-beast.cpp -o test-beast -lm

imgtest: ost_imgtest.c ../ost/ost.h
	$(CC) -std=c99 -Wall -Wextra -O3 ost_imgtest.c -o ost_imgtest $(LIBPNG) -lm

pipeline: ost_pipeline.c ../ost/ost.h ost_test_work.h
	$(CC) -std=c99 -Wall -Wextra -O3 ost_pipeline.c -o ost_pipeline $(LIBPNG) -lm
