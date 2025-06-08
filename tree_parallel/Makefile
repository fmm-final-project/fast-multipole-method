CC      := gcc
CFLAGS  := -O3 -march=native -g

BIN     := tree_parallel
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): tree_parallel.c
	$(CC) -o $(BIN) $(CFLAGS) tree_parallel.c $(LIB)

clean:
	rm -f $(BIN)