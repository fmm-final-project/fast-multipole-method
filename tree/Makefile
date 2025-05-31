CC      := gcc
CFLAGS  := -O3 -march=native -g

BIN     := tree
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): tree.c
	$(CC) -o $(BIN) $(CFLAGS) tree.c $(LIB)

clean:
	rm -f $(BIN)