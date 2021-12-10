from mesh import Mesh
import numpy as np

def vec3_equal_vec3(a,b):
    threshold = 1e-6
    return np.dot(a-b, a-b)<threshold

def is_edge_present(a, b, h):
    for c in range(h.ncorners):
        vi = h.V[h.org(c)]
        vj = h.V[h.dst(c)]
        if (vec3_equal_vec3(vi, a) and vec3_equal_vec3(vj, b)) or (vec3_equal_vec3(vi, b) and vec3_equal_vec3(vj, a)):
            return True
    return False

models = [("chevron", 4), ("ifp1", 2), ("ifp2", 2), ("shell", 6)]
for model,nhor in models:
    m = Mesh(model+"/slice.obj")
    horizon_meshes = [ Mesh(model+"/horizon"+str(i)+".obj") for i in range(nhor) ]
    faults_mesh = Mesh(model+"/faults.obj")

    horizon_id = [-1]*m.ncorners
    is_fault   = [ False ]*m.ncorners
    fault_opposite = [-1]*m.ncorners

    for c in range(m.ncorners): # for all half-edges
        vi = m.V[m.org(c)]
        vj = m.V[m.dst(c)]
        for (hid, mh) in enumerate(horizon_meshes):
            if is_edge_present(vi, vj, mh):
                horizon_id[c] = hid
        if is_edge_present(vi, vj, faults_mesh):
            is_fault[c] = True

    fault_corners = [ c for c in range(m.ncorners) if is_fault[c] ]
    for c1 in fault_corners:
        for c2 in fault_corners:
            if c1==c2: continue
            if vec3_equal_vec3(m.V[m.org(c1)], m.V[m.dst(c2)]) and vec3_equal_vec3(m.V[m.org(c2)], m.V[m.dst(c1)]):
                fault_opposite[c1] = c2
                fault_opposite[c2] = c1

    out = open(model+"/attributes.py", 'w')
    print("horizon_id = " + str(horizon_id), file = out)
    print("is_fault = " + str(is_fault), file = out)
    print("fault_opposite = " + str(fault_opposite), file = out)
    out.close()
