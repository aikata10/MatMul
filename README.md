# MatMul

## Details
1. Refer to \url{https://fherma.io/challenges/652bf669485c878710fd020b/overview} for a setup guide to OpenFHE.
2. The `Mat_Mult` directory contains, by default, the proposed arbitrary packing for square matrices.
3. The `square` sub-directory contains the same files as the `Mat_mult` folder.
4. The `rec` folder contains the files for the proposed rectangular multiplication with case 1. These can be compied into Mat_mult folder for testing.

## HOW TO RUN

A simple make command in the Mat_mul folder starts the program run. It is currently running for 100 iterations, hecne takes a long time. Global variables for matrix dimensions, CKKS Ring Dimensionans, and number of runs can be modified using the file `matrix_mult_ckks.cpp`.

