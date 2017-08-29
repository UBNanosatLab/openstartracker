CC = g++

all: $(OBJ)
	$(CC) beastgen.cpp -o beastgen $(OBJ)
	$(CC) -Ofast -fPIC -c beast.cpp -o beast.o
	
debug:
	$(CC) -Wall -ffast-math -Og -g beastgen.cpp -o beastgen $(OBJ)
	$(CC) -Wall -ffast-math -Og -g beast.cpp -o beast
