from mesh import Mesh
import scipy.sparse
from scipy.sparse.linalg import lsmr
import numpy as np

m = Mesh("rocker_profile.obj")

# construire A et b

f = lsmr(A, b)[0] # call the least squares solver
m.write_vtk("rocker_profile.vtk", f)

isos = [4.54, 11.723, 16, 19]

for iso in isos:
    ntri = m.ntriangles
    for i in range(ntri):
        case = int(f[m.T[i][0]]>iso) + 2*int(f[m.T[i][1]]>iso) + 4*int(f[m.T[i][2]]>iso)
        if 0==case or 7==case: continue # do not cut the triangle, no iso-line here
        if 4==case or 3==case:           # rotate triangle indices
            m.T[i] = np.roll(m.T[i], +1) # so the tip we cut
        if 2==case or 5==case:           # is centered on the first vertex
            m.T[i] = np.roll(m.T[i], -1)
        a = (f[m.T[i][1]]-iso)/(f[m.T[i][1]]-f[m.T[i][0]]) # barycentric coordinates
        b = (f[m.T[i][2]]-iso)/(f[m.T[i][2]]-f[m.T[i][0]]) # of new vertices on edges 0-1 and 0-2
        f = np.append(f, [iso,iso])

        m.add_vertex(m.V[m.T[i][1]]*(1-a) + m.V[m.T[i][0]]*a)
        m.add_vertex(m.V[m.T[i][2]]*(1-b) + m.V[m.T[i][0]]*b)
        m.add_triangle([m.T[i][1],  m.T[i][2], m.nverts-2]) # these two new triangles
        m.add_triangle([m.nverts-2, m.T[i][2], m.nverts-1]) # form the base quad
        m.T[i][1] = m.nverts-2 # reuse the original triangle
        m.T[i][2] = m.nverts-1 # to form the tip

m.write_vtk("rocker_profile2.vtk", f)



