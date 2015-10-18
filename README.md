# hcp_corr
load hcp (human connectome project) data, get correlation matrices with less memory usage and faster

Python modules required:

NumPy 1.9.X -- preferred with high performance blas/lapack libraries (e.g. openblas or intel-mkl)

SciPy 0.14.1

NiBabel 2.1.0dev -> $ git clone --branch enh/cifti2 https://github.com/satra/nibabel.git

numexpr 2.4.4

load_hcp.py : loading *GIFTI formatted hcp data for a given subject and generating its time-series
matrix. The function is affective to load data for many subjects e.g. in a loop and especially good
with "N_user" option, which decides the length of time-series for fast-forward runs. 

corr_faster.py : given a matrix X, this function returns a 1D array, whose elements correspond to
the upper-triangular part (excluding diagonals) of correlation coefficient matrix. "write_upper"
function exports the 1D array more efficient than numpy.savetxt.

corr_full.py : given a 1D array of upper triangular correlation matrix, this function returns the
full correlation matrix

Step 1, load time-series as a numpy array

K = load_hcp.t_series(data_path, subject, template, cnt_files, N_user, subject_path=None, dtype=None)

Step 2, get upper triangular of correlation matrix of time-series

K = corr_faster.corrcoef_upper(K)

Step 3, convert 1D array to its full-symmetric 2D correlation matrix. 1D array must be re-sized 
outside of function.

N_orig = corr_full.N_original(K) 

K.resize([N_orig,N_orig])

corr_full.upper_to_down(K)

References:
        http://www.humanconnectome.org/
        https://github.com/satra/nibabel/tree/enh/cifti2 
