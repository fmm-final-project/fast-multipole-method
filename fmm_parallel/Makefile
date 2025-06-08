CC      := gcc
CFLAGS  := -O3 -march=native -g

BIN     := tree2fmm_parallel
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): tree2fmm_parallel.c
	$(CC) -o $(BIN) $(CFLAGS) tree2fmm_parallel.c $(LIB)

clean:
	rm -f $(BIN)