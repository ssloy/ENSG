"""
Author(s): Mohamed Aidahi, Abdrahmane Berete
Date: 2023-12-26
"""
import os

from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError
import importlib
import sys


def deform_mesh(model_path: Path, output_path: Path = None, horizon_const: float = 10, fault_const: float = 10):
    """
    Deform a geological mesh using a least-squares approach.
    Make faults vertical and horizons horizontal.

    :param model_path: Path, folder which contains the mesh slice.obj and attributes.py files
    :param output_path: Path, output folder
    :param horizon_const: float, value for horizons constraint in the least-squares system
    :param fault_const: float, value for faults constraint in the least-squares system
    :return: None
    """
    # Read the model data
    m = Mesh(str(model_path.joinpath("slice.obj")))
    attr = importlib.import_module(str(model_path)  +".attributes")

    # Initialise the linear system
    A = scipy.sparse.lil_matrix((m.ncorners * 3, m.nverts * 2))
    b = [m.V[m.dst(c)][dim] - m.V[m.org(c)][dim] for c in range(m.ncorners) for dim in [0, 1]] + [0] * m.ncorners

    # First constraint: to avoid the perverted case of all vertices squeezed into the same position.
    for c in range(m.ncorners):
        for dim in [0, 1]:
            A[c * 2 + dim, m.dst(c) * 2 + dim] = 1
            A[c * 2 + dim, m.org(c) * 2 + dim] = -1


    row = m.ncorners * 2
    for c in range(m.ncorners):
        # the horizons constraint
        if attr.horizon_id[c] >= 0:
            A[row, m.dst(c) * 2 + 1] = horizon_const
            A[row, m.org(c) * 2 + 1] = -horizon_const

        # the faults constraints
        if attr.is_fault[c]:
            A[row + 1, m.dst(c) * 2] = fault_const
            A[row + 1, m.org(c) * 2] = -fault_const

            # joinning connex components
            i = m.dst(c) * 2
            j = m.org(attr.fault_opposite[c]) * 2
            A[row, i + 1] = 1
            A[row, j + 1] = -1

            A[row // 2, i] = 1
            A[row // 2, j] = -1

        row += 1

    # Solving the least-squares system for x and
    A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector multiplications
    x = lsmr(A, b)[0]

    # writing the solution
    for v in range(m.nverts):
        for dim in [0, 1]:
            m.V[v][dim] = x[v * 2 + dim]

    for c in range(m.ncorners):  # lift all vertices of all horizons
        if attr.horizon_id[c] >= 0:
            height = (1 + attr.horizon_id[c]) / 37.76  # arbitrarily chosen coeff to get a visually nice result
            m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height/10

    for c in range(m.ncorners):  # lower vertices in faults
        if attr.is_fault[c]:
            m.V[m.org(c)][2] -= 0.01  # arbitrary scaling coefficent
            m.V[m.dst(c)][2] -= 0.01  # to make the result look nice

    m.write_vtk(str(output_path.joinpath("deformed.vtk")))


def positive_float(x):
    try:
        x = float(x)
    except ValueError:
        raise ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0:
        raise ArgumentTypeError("%r is not >= 0" % (x,))

    return x


if __name__ == "__main__":

    parser = ArgumentParser(description="deform a geological mesh using a least-squares method.")
    # model path
    parser.add_argument("model_path",
                        type=str,
                        help="path to folder containing slice.obj & attributes.py")
    # output path
    parser.add_argument("output_path",
                        type=str,
                        help="path to output folder")
    # horizons constraint
    parser.add_argument("horizon_const",
                        type=positive_float,
                        help="value of the horizons constraint in the least-squares system")
    # faults constraint
    parser.add_argument("fault_const",
                        type=positive_float,
                        help="value of the faults constraint in the least-squares system")

    args = parser.parse_args(sys.argv[1:])

    model_path = Path(args.model_path)
    output_path = Path(args.output_path)

    try:
        # check if the mesh file exists
        assert model_path.joinpath('slice.obj').is_file()

        # check if the attributes.py file exists
        assert model_path.joinpath('attributes.py').is_file()

        # check if the output path exists
        assert output_path.is_dir()
        # check if write permission to output path is granted
        assert os.access(output_path, os.W_OK)

        # call the function
        deform_mesh(model_path, output_path, args.horizon_const, args.fault_const)

    except AssertionError:
        message = ""
        if not model_path.joinpath("slice.obj").is_file():
            message += "Could not find slice.obj file in " + str(model_path) + ".\n"

        if not model_path.joinpath("attributes.py").is_file():
            message += "Could not find attributes.py file in" + str(model_path) + ".\n"

        if not output_path.is_dir():
            message += "Could not find output path " + str(output_path) + ".\n"

        if not os.access(output_path, os.W_OK):
            message += "Write permission to output path is not granted.\n"

        print(message)
