from mesh import Mesh
import importlib

model = "shell"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")

print(f"n_corners = {m.ncorners}")
corner = 5722
print(f"corner = {corner}, vert = {m.org(corner)}, horizon = {attr.horizon_id[corner]}")
print(f"previous_corner = {m.prev(corner)}, vert = {m.org(m.prev(corner))}, horizon = {attr.horizon_id[m.prev(corner)]}")
print(f"next_corner = {m.next(corner)}, vert = {m.org(m.next(corner))}, horizon = {attr.horizon_id[m.next(corner)]}")

