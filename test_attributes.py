from mesh import Mesh
import importlib

model = "chevron"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")

for c in range(m.ncorners): # lift all vertices of all horizons
    if attr.horizon_id[c]>=0:
        height = (1+attr.horizon_id[c])/37.76 # arbitrarily chosen coeff to get a visually nice result
        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height

for c in range(m.ncorners): # lower vertices in faults
    if attr.is_fault[c]:
        m.V[m.org(c)][2] -= 0.00431 # arbitrary scaling coefficent
        m.V[m.dst(c)][2] -= 0.00431 # to make the result look nice
print(m)
