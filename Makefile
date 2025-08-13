CC = gcc
CFLAGS = -Iinclude -O2
LIBS = -lpbc -lgmp

all: build_circuit interpolate pot

build_circuit: src/build_circuit.c src/circuit.c
	$(CC) $(CFLAGS) -o $@ src/build_circuit.c src/circuit.c $(LIBS)

interpolate: src/interpolate.c src/circuit.c
	$(CC) $(CFLAGS) -o $@ src/interpolate.c src/circuit.c $(LIBS)

pot: src/pot.c
	$(CC) $(CFLAGS) -o $@ src/pot.c $(LIBS)

clean:
	rm -f build_circuit interpolate pot
