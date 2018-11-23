#include <vector>
#include <iostream>
#include "OpenNL_psm.h"


int main(int argc, char** argv) {
    std::vector<float> x(32);
    // declare new solver
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 32);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    // "lock" 3 samples (via a quadratic penalty term)
    float scale = 100.;
    nlBegin(NL_ROW);
    nlCoefficient(0, 1.*scale);
    nlRightHandSide(1.*scale);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW);
    nlCoefficient(18, 1.*scale);
    nlRightHandSide(2.*scale);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW);
    nlCoefficient(31, 1.*scale);
    nlRightHandSide(1.*scale);
    nlEnd(NL_ROW);

    // laplacian smoothing
    for (int i=1; i<32; i++) {
        nlBegin(NL_ROW);
        nlCoefficient(i,  1);
        nlCoefficient(i-1, -1);
        nlEnd(NL_ROW);
    }

    // finish building the matrix and solve the system
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    // retrieve the solution
    for (int i=0; i<32; i++) {
        x[i] = nlGetVariable(i);
        std::cout << x[i] << std::endl;
    }
    return 0;
}

