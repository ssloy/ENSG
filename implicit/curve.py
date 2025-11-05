import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr

n = 32
if False:
    A = np.matrix(np.zeros((n+1,n)))
    np.fill_diagonal(A[:n-2,  :], -1.) # second order finite differences
    np.fill_diagonal(A[:n-2, 1:],  2.)
    np.fill_diagonal(A[:n-2, 2:], -1.)
    A[n-2,   0] = 1000*1. # "lock" vertices
    A[n-1, n-1] = 1000*1.
    A[n,  n//3] = 1000*1.
    
    b = np.matrix([[1.]]*(n-2) + [[1000*0.],[1000*0.],[1000*0.]])
    x = np.linalg.inv(A.T*A)*A.T*b
else:
    A = scipy.sparse.lil_matrix((n+1, n))
    for i in range(n-2):
        A[i,i] = -1
        A[i,i+1] = 2
        A[i,i+2] = -1
    A[n-2,   0] = 1000*1. # "lock" vertices
    A[n-1, n-1] = 1000*1.
    A[n,  n//3] = 1000*1.
    A = A.tocsc() # sparse row matrix for fast matrix-vector multiplication

    b = [[1.]]*(n-2) + [[1000*0.],[1000*0.],[1000*0.]]
    x = lsmr(A, b, maxiter=10000)[0] # call the least squares solver

import matplotlib.pyplot as plt
plt.plot(x, 'gs-', label='x(i)')
plt.legend(frameon=False)
plt.show()
