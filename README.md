# Groth16 Polynomial Demo

## Modules & Usage

```bash
make clean
make
```

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
*/
