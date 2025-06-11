CC      := gcc
CFLAGS  := -O3 -march=native -g

BIN     := naive_gravity_parallel
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): naive_gravity_parallel.c
	$(CC) -o $(BIN) $(CFLAGS) naive_gravity_parallel.c $(LIB)

clean:
	rm -f $(BIN)