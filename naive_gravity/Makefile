CC      := gcc
CFLAGS  := -O2

BIN     := naive_gravity_omp
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): naive_gravity_omp.c
	$(CC) -o $(BIN) $(CFLAGS) naive_gravity_omp.c $(LIB)

clean:
	rm -f $(BIN)