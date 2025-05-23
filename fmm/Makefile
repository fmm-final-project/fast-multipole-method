CC      := gcc
CFLAGS  := -O2 -g

BIN     := fmm
LIB     := -lm

all: $(BIN)

$(BIN): fmm.c
	$(CC) -o $(BIN) $(CFLAGS) fmm.c $(LIB)

clean:
	rm -f $(BIN)