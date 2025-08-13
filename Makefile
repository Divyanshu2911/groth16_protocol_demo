CC = gcc
CFLAGS = -Iinclude -O2
LIBS = -lpbc -lgmp

FMT = src/fmt.c

all: build_circuit interpolate pot keygen prover verifier

build_circuit: src/build_circuit.c src/circuit.c $(FMT)
	$(CC) $(CFLAGS) -o $@ src/build_circuit.c src/circuit.c $(FMT) $(LIBS)

interpolate: src/interpolate.c src/circuit.c  $(FMT)
	$(CC) $(CFLAGS) -o $@ src/interpolate.c src/circuit.c $(FMT) $(LIBS)

pot: src/pot.c $(FMT)
	$(CC) $(CFLAGS) -o $@ src/pot.c $(FMT) $(LIBS)

keygen: src/keygen.c src/circuit.c src/poly.c $(FMT)
	$(CC) $(CFLAGS) -o $@ src/keygen.c src/circuit.c src/poly.c $(FMT) $(LIBS)

prover: src/prover.c src/circuit.c src/poly.c $(FMT)
	$(CC) $(CFLAGS) -o $@ src/prover.c src/circuit.c src/poly.c $(FMT) $(LIBS)

verifier: src/verifier.c $(FMT)
	$(CC) $(CFLAGS) -o $@ src/verifier.c $(FMT) $(LIBS)

clean:
	rm -f build_circuit interpolate pot keygen prover verifier
