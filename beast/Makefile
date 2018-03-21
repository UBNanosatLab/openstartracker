CC = g++

all: $(OBJ)
	swig -python -c++ beast.i
	$(CC) -g -std=c++11 -Ofast -fPIC -c beast_wrap.cxx -o beast_wrap.o -lstdc++ -I/usr/include/python2.7
	$(CC) -g -shared -fPIC beast_wrap.o -o _beast.so
	
clean:
	rm beast.so beast.o	beast.py beast.pyc beast_wrap.cxx beast_wrap.o
