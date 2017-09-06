CC = gcc

all: $(OBJ)
	$(CC) beastgen.c -o beastgen -lm
	$(CC) -Ofast -fPIC main.c -o main -lm
	
debug:
	$(CC) -Wall -ffast-math -Og -g beastgen.c -o beastgen -lm
	$(CC) -Wall -ffast-math -Og -g main.c -o main -lm
