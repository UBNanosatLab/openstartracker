CC = g++

all: $(OBJ)
	$(CC) beastgen.cpp -o beastgen $(OBJ)
	$(CC) -Ofast -fPIC -c beast.cpp -o beast.o
	$(CC) -Ofast -fPIC -c beast_wrap.cxx -o beast_wrap.o -lstdc++ -I/usr/include/python2.7
	$(CC) -shared -fPIC beast_wrap.o beast.o -o _beast.so

debug:
	$(CC) -Wall -ffast-math -Og -g beastgen.cpp -o beastgen $(OBJ)
	$(CC) -Wall -ffast-math -Og -g beast.cpp -o beast
