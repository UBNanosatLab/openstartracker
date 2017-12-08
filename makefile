CC = g++

all: $(OBJ)
	#$(CC) beastgen.cpp -o beastgen -lm
	$(CC) -Ofast -fPIC test.c -o test -lm
	
debug:
	#$(CC) -Wall -Og -pg -g beastgen.cpp -o beastgen -lm
	$(CC) -Wall -O0 -pg -g test.c -o test -lm
