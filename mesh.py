import numpy as np
import sys

class Mesh():
    def __init__(self, filename):
        # parse the .obj file
        V, T = [], []
        with open(filename) as f:
           for line in f.readlines():
               if line.startswith('#'): continue
               values = line.split()
               if not values: continue
               if values[0] == 'v':
                   V.append([float(x) for x in values[1:4]])
               elif values[0] == 'f':
                   T.append([int(x) for x in values[1:4]])
        self.V, self.T = np.array(V), np.array(T)-1

        # compute the adjacency
        self.v2c = np.array([-1]*self.nverts)
        self.c2c = np.array([-1]*self.ncorners)
        for c in range(self.ncorners):
            v = self.T[c//3][c%3]
            self.v2c[v] = c
        for c in range(self.ncorners):
            v = self.T[c//3][c%3]
            self.c2c[c] = self.v2c[v]
            self.v2c[v] = c

        # speed up the computations
        self.opp = np.array([-1]*self.ncorners)
        for c in range(self.ncorners):
            c_org = self.T[c//3][c%3]
            c_dst = self.T[c//3][(c+1)%3]
            cir = c
            opp = -1
            while True:
                cand = (cir//3)*3 + (cir+2)%3
                cand_org = self.T[cand//3][cand%3]
                cand_dst = self.T[cand//3][(cand+1)%3]
                if (cand_org == c_dst and cand_dst == c_org):
                    opp = cand # we suppose manifold input
                cir = self.c2c[cir]
                if (cir==c): break
            self.opp[c] = opp

        self.boundary = np.array([False]*self.nverts)
        for v in range(self.nverts):
            cir = self.v2c[v]
            if cir<0: continue
            while (True):
                if self.opp[cir]<0:
                    self.boundary[v] = True
                    break
                cir = self.c2c[cir]
                if (cir==self.v2c[v]): break

    @property
    def nverts(self):
        return len(self.V)

    @property
    def ntriangles(self):
        return len(self.T)

    @property
    def ncorners(self):
        return len(self.T)*3

    def org(self, c):
        return self.T[c//3][c%3]

    def dst(self, c):
        return self.T[c//3][(c+1)%3]

    def prev(self, c):
        return (c//3)*3 + (c+2)%3

    def next(self, c):
        return (c//3)*3 + 1

    def opposite(self, c):
        return self.opp[c]

    def on_border(self, v):
        return self.boundary[v]

    def neighbors(self,v):
        out = []
        cir = self.v2c[v]
        if cir<0: return []
        cnt = 0
        while True:
            neigh = self.T[cir//3][(cir+1)%3]
            out.append(neigh)
            cir = self.c2c[cir]
            if (cir==self.v2c[v]): break
            cnt = cnt + 1
        return out

    def __str__(self):
        ret = ""
        for v in self.V:
            ret = ret + ("v %f %f %f\n" % (v[0], v[1], v[2]))
        for t in self.T:
            ret = ret + ("f %d %d %d\n" % (t[0]+1, t[1]+1, t[2]+1))
        return ret

    def write_vtk(self, filename, scalar_field=None):
        original_stdout = sys.stdout
        with open(filename, 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.


            print("""# vtk DataFile Version 3.0
ENSG geometry processing course
ASCII
DATASET UNSTRUCTURED_GRID
""")
            print("POINTS %d double" % self.nverts)
            for v in self.V:
                print("%f %f %f" % (v[0], v[1], v[2]))
            print("\nCELLS %d %d" % (self.ntriangles, 4*self.ntriangles))
            for t in self.T:
                print("3 %d %d %d" % (t[0], t[1], t[2]))
            print("\nCELL_TYPES %d" % self.ntriangles)
            print( "5 " * self.ntriangles)

            if scalar_field is not None:
                assert len(scalar_field)==self.nverts
                print("""\nPOINT_DATA %d
SCALARS f double 1
LOOKUP_TABLE default""" % self.nverts)
                for v in scalar_field:
                    print("%f " % v, end="")
            sys.stdout = original_stdout # Reset the standard output to its original value



