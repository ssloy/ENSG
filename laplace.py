from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr

m = Mesh("shell/slice.obj")
nboundary = sum(m.boundary)

for dim in range(2): # solve for x first, then for y
    A = scipy.sparse.lil_matrix((nboundary + m.ncorners, m.nverts))
    b = [0] * A.shape[0]
    for row in range(m.ncorners):
        i = m.org(row)
        j = m.dst(row)
        A[row, j] =  1
        A[row, i] = -1

    for (i,v) in enumerate(m.V):
        if m.on_border(i):
            A[row, i] = 1 * 100 # quadratic penalty to lock boundary vertices
            b[row] = v[dim] * 100
            row += 1
    A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
    x = lsmr(A, b)[0] # call the least squares solver
    for i in range(m.nverts): # apply the computed flattening
        m.V[i][dim] = x[i]

m.write_vtk("output.vtk")
m.write_obj("output.obj")

