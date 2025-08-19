# Groth16 Polynomial Demo

## Modules & Usage

### Compilation

```bash
make clean
make
```

### Execution steps

1. **build_circuit**: Constructs R1CS and prints wires
   ```bash
   ./build_circuit path/to/a.param [degree of the polynomial y = f(x)] x y a0…ad
   ```
2. **interpolate**: Builds QAP-polynomials A_j,B_j,C_j
   ```bash
   ./interpolate path/to/a.param [degree of the polynomial y = f(x)] x y a0…ad
   ```
3. **pot**: Powers-of-Tau ceremony
   ```bash
   ./pot path/to/a.param [deg]
   ```

4. **keygen**: Generates prover's and verifier's keys using pairing
   ```bash
   ./keygen path/to/a.param [deg] x y a0....ad
   ```

5. **prove**: Generates proof as proof_demo.bin
   ```bash
   ./prover path/to/a.param [deg] x y a0....ad
   ```

6. **verify**: Verifies the generated proof
   ```bash
   ./verifier path/to/a.param proof_demo.bin
   ```