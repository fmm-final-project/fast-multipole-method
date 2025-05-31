CC      := gcc
CFLAGS  := -O3 -march=native -g

BIN     := tree2fmm
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): tree2fmm.c
	$(CC) -o $(BIN) $(CFLAGS) tree2fmm.c $(LIB)

clean:
	rm -f $(BIN)