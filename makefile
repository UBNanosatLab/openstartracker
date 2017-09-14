CC = g++

all: $(OBJ)
	$(CC) beastgen.cpp -o beastgen -lm
	$(CC) -Ofast -fPIC main.c -o main -lm
	
debug:
	$(CC) -Wall -Og -pg -g beastgen.cpp -o beastgen -lm
	$(CC) -Wall -Og -pg -g main.c -o main -lm
