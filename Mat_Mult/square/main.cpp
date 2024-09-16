#include <iostream>
#include <fstream>
#include <string>
#include "matrix_mult_ckks.h"

using namespace lbcrypto;

int main(int argc, char *argv[])
{
    std::string pubKeyLocation;
    std::string multKeyLocation;
    std::string rotKeyLocation;
    std::string ccLocation;
    std::string matrixALocation;
    std::string matrixBLocation;
    std::string outputLocation;

    MatrixMultCKKS matrixMultCKKS(ccLocation, pubKeyLocation, multKeyLocation, rotKeyLocation, matrixALocation, matrixBLocation,
                                  outputLocation);
    matrixMultCKKS.eval();
    matrixMultCKKS.deserializeOutput();
    return 0;
}
