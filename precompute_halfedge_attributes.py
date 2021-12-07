from mesh import Mesh
import numpy as np

def is_edge_present(a, b, h):
    threshold = 1e-6
    for c in range(h.ncorners):
        vi = h.V[h.org(c)]
        vj = h.V[h.dst(c)]
        if (np.dot(vi-a, vi-a)<threshold and np.dot(vj-b, vj-b)<threshold) or (np.dot(vj-a, vj-a)<threshold and np.dot(vi-b, vi-b)<threshold):
            return True;
    return False

models = [("chevron", 4), ("ifp1", 2), ("ifp2", 2), ("shell", 6)]
for model,nhor in models:
	m = Mesh(model+"/slice.obj")
	horizon_meshes = [ Mesh(model+"/horizon"+str(i)+".obj") for i in range(1,nhor+1) ]
	faults_mesh = Mesh(model+"/faults.obj")
 
	horizon_id = [-1]*m.ncorners
	is_fault   = [ False ]*m.ncorners
 
	for c in range(m.ncorners): # for all half-edges
	    vi = m.V[m.org(c)]
	    vj = m.V[m.dst(c)]
	    for (hid, mh) in enumerate(horizon_meshes):
        	if is_edge_present(vi, vj, mh):
	            horizon_id[c] = hid
	    if is_edge_present(vi, vj, faults_mesh):
	        is_fault[c] = True

	out = open(model+"/attributes.py", 'w')
	print("horizon_id = " + str(horizon_id), file = out)
	print("is_fault = " + str(is_fault), file = out)
	out.close()
