.PHONY: all clean

CC = gcc
CXX = g++
CFLAGS ?= -Ofast -g
CXXFLAGS ?= -Ofast -g
LIBPNG = $$(pkg-config --cflags --libs libpng)
OST_OBJ = ../ost/ost.o
TRACKER_DEPS = ../ost/ost.h ost_test_work.h ost_test_tracker.h $(OST_OBJ)

all: test test_cpp test-beast imgtest pipeline

$(OST_OBJ): ../ost/ost.c ../ost/ost.h
	$(MAKE) -C ../ost ost.o

test: test.c $(TRACKER_DEPS)
	$(CC) -std=c99 -Wall -Wextra $(CFLAGS) test.c $(OST_OBJ) -o test -lm

test_cpp: test.c $(TRACKER_DEPS)
	$(CXX) -std=c++11 -Wall -Wextra $(CXXFLAGS) test.c $(OST_OBJ) -o test_cpp -lm

test-beast: test-beast.cpp ../beast/beast.h ../ost/ost.h $(OST_OBJ)
	$(CXX) -std=c++11 -Wall -Wextra $(CXXFLAGS) test-beast.cpp $(OST_OBJ) -o test-beast -lm

imgtest: ost_imgtest.c ../ost/ost.h $(OST_OBJ)
	$(CC) -std=c99 -Wall -Wextra -O3 ost_imgtest.c $(OST_OBJ) -o ost_imgtest $(LIBPNG) -lm

pipeline: ost_pipeline.c $(TRACKER_DEPS)
	$(CC) -std=c99 -Wall -Wextra -O3 ost_pipeline.c $(OST_OBJ) -o ost_pipeline $(LIBPNG) -lm

clean:
	rm -f test test_cpp test-beast ost_imgtest ost_pipeline gmon.out test.pdf */result_cpp.csv
	$(MAKE) -C ../ost clean
