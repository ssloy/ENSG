#include <vector>
#include <iostream>
#include "geometry.h"
#include "model.h"
#include "OpenNL_psm.h"


int main(int argc, char** argv) {
    if (argc<2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    // load the model
    Model m(argv[1]);

    // for each vertex determine if it is located on the boundary
    std::vector<bool> boundary_verts(m.nverts(), false);
    for (int i=0; i<m.nhalfedges(); i++) {
        if (m.opp(i)<0) {
            boundary_verts[m.from(i)] = true;
            boundary_verts[m.to  (i)] = true;
        }
    }

    // we solve three independent systems for x,y and z
    for (int d=0; d<3; d++) {
        // declare new solver
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, m.nverts());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        // "lock" the boundary vertices (via a quadratic penalty term)
        for (int v=0; v<m.nverts(); v++) {
            if (boundary_verts[v]) {
                float scale = 10.;
                nlBegin(NL_ROW);
                nlCoefficient(v, 1.*scale);
                nlRightHandSide(m.point(v)[d]*scale);
                nlEnd(NL_ROW);
            }
        }

        // laplacian smoothing
        for (int hi=0; hi<m.nhalfedges(); hi++) {
            if (m.opp(hi)<0 || m.opp(hi)<hi) continue; // we skip boundary halfedges
            int v1 = m.from(hi);
            int v2 = m.to  (hi);
            nlBegin(NL_ROW);
            nlCoefficient(v1,  1);
            nlCoefficient(v2, -1);
            nlEnd(NL_ROW);
        }

        // finish building the matrix and solve the system
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        // retrieve the solution
        for (int i=0; i<m.nverts(); i++) {
            m.point(i)[d] = nlGetVariable(i);
        }
    }
    std::cout << m << std::endl;
    return 0;
}

