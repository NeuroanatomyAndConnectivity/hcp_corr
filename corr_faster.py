import warnings
import numpy as np
from scipy.linalg import get_blas_funcs


def corrcoef_upper(X):
    """
    .....numpy.corrcoef 1.9.1 adapted for mem-usage optimization....
    X : numpy 2D array, MxN
        
        Correlation coefficients of a given matris is row-wise calculated.  
        The upper diagonal of the resultant correlation matrix (of shape MxM) 
        is returned as a  vector
    
    output, c: numpy 1D array, length  M.(M-1)/2 

    """
    c, d = my_cov(X)

    d = np.sqrt(d)
    # calculate "c / multiply.outer(d, d)" row-wise for mem & speed
    k = 0
    for i in range(0, d.size - 1):
        len = d.size - i - 1
        c[k:k + len] /= (d[-len:] * d[i])
        k += len

    return c

def my_cov(m):
    # numpy.cov 1.9.1 adapted for real arrays
    m = np.asarray(m)
    dtype = np.result_type(m, np.float64)
    X = np.array(m, ndmin=2, dtype=dtype)
    N = X.shape[1]
    
    fact = float(N - 1)
    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice", RuntimeWarning)
        fact = 0.0

    X -= X.mean(axis=1, keepdims=True)
    # This returns np.dot(X, X.T) / fact 
    return cool_syrk(1.0/fact, X)

def cool_syrk(fact, X): 
    syrk = get_blas_funcs("syrk", [X])
    R = syrk(fact, X)
    d = np.diag(R).copy()
    size = mat_to_upper_F(R)
    R.resize([size,])
    return R,d

def mat_to_upper_F(A):
    # check if input array Fortran style contiguous
    if not A.flags['F_CONTIGUOUS']:
        raise Exception("F_CONTIGUOUS required")
    # row-number of matrix A        
    n = A.shape[0]
    # size of upper-triangular of matrix A
    size = (n - 1) * n / 2
    # reshape matrix A column-wise into U
    U = A.reshape([1,n*n], order='F') 
    k = 0
    # fill out U by upper diagonal elements of A
    for i in range(0, n-1):
        len = n - 1 - i
        U[0,k:k+len] = A[i,i+1:n]
        k += len
    return size

# save upper diagonal correlation matrix as 1D array
def write_upper(file, A, fmt="%g"):
    count = A.size
    print "count", count
    A = A.reshape([count,])
    step = 10000
    k = 0
    f = open(file,'wb') #?
    while k < count:
        i = min(step,count-k)
        np.savetxt(f,A[k:k+i].reshape([1,i]),fmt=fmt, delimiter='\n', newline='\n')
        k += i
    f.close()