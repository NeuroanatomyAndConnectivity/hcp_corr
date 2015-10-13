"""
Given the upper diagonal of a correlation matrix, this script generates the
full correlation matrix. This is much less memory consuming than np.corrcoef,
however, it is tricky since the vector must be 'resized' outside of functions.

    vector: 1D numpy array
        Upper diagonal of a correlation matrix reduced in to a 1D array

    # Step 1
    find out original dimension of correlation matrix, NxN
    # Step 2
    resize 1D array as a 2D array, fill out it with zeros for empty elements
    # Step 3
    fill in new 2D array 

    # here we go...    
    
    N = N_original(vector)
    vector.resize([N,N])
    upper_to_down(vector)

"""

from scipy.optimize import fsolve

def N_original(a):
    x = a.shape[0]
    def func(N):
        return N*(N-1.0) / 2.0 -x
    n = int(round(fsolve(func, [x])))
    return n
    
def upper_to_down(M):

    if not M.flags['C_CONTIGUOUS']:
        raise Exception("C_CONTIGUOUS required")
    n = M.shape[0]
    size = (n - 1) * n / 2
    U = M.reshape([n*n,])
    k = size
    for i in range(n-1, -1, -1):
        len = n - 1 - i
        M[i+1:n,i] = U[k-len:k]
        M[i,i+1:n] = U[k-len:k]
        M[i,i] = 1.0
        k -= len
    return M