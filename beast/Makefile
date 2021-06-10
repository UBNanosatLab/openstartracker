CC = g++
PYTHONHEADERS=/usr/include/python3.8
#PYTHONHEADERS=/usr/include/python2.7

all: $(OBJ)
	swig -python -c++ beast.i
	$(CC) -g -std=c++11 -Ofast -fPIC -c beast_wrap.cxx -o beast_wrap.o -lstdc++ -I$(PYTHONHEADERS)
	$(CC) -g -shared -fPIC beast_wrap.o -o _beast.so
	
clean:
	rm beast.so beast.o	beast.py beast.pyc beast_wrap.cxx beast_wrap.o
