CC      := gcc
CFLAGS  := -O2

BIN     := naive_gravity
LIB     := -lm -fopenmp

all: $(BIN)

$(BIN): naive_gravity.c
	$(CC) -o $(BIN) $(CFLAGS) naive_gravity.c $(LIB)

clean:
	rm -f $(BIN)