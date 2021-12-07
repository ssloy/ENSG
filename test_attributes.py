from mesh import Mesh
import importlib

model = "ifp2"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")

for c in range(m.ncorners):
	if (attr.horizon_id[c]>=0):
		m.V[m.org(c)][2] = m.V[m.dst(c)][2] = (1+attr.horizon_id[c])/37.76 # divide by something nice

for c in range(m.ncorners):
	if (attr.is_fault[c]):
		m.V[m.org(c)][2] = m.V[m.org(c)][2] - 1/68.54 # arbitrary scaling coefficent
		m.V[m.dst(c)][2] = m.V[m.dst(c)][2] - 1/68.54 # to make the result look nice
print(m)
