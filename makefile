CC = gcc

all: $(OBJ)
	$(CC) beastgen.c -o beastgen -lm
	$(CC) -Ofast -fPIC main.c -o main -lm
	
debug:
	$(CC) -Wall -O0 -pg -g beastgen.c -o beastgen -lm
	$(CC) -Wall -O0 -pg -g main.c -o main -lm
