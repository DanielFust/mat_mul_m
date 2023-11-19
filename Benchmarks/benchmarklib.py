'''
Python bindings for mat_mul_m benchmarks

benchmark1: B=matmul(A,A)
benchmark2: A=matmul(A,A)
benchmark3: B=c1*matmul(At,A) + c2*matmul(Bt,B)
benchmark4: 


cache_blocking:    naive vs blocked matrix multiplication
matmul_benchmark1: naive vs. BLAS vs. MATMUL vs. cache-blocked for B=matmul(A,A)
matmul_benchmark2: naive vs. BLAS vs. MATMUL vs. cache-blocked for B=matmul(At,A)

'''
import numpy as np
import ctypes

# Select the correct 'shared object' file extension for OS
extension = 'dylib' #MacOS
#extension = 'so'    #Linux
#extension = 'dll'   #Windows

'''
  Data structure to hold mat_mul_m operation execution speeds and
  control test execution speeds to compare performance 
''' 
class Benchmark:
    execuction_times = np.array([], dtype=ctypes.c_double)
    control_times    = np.array([], dtype=ctypes.c_double)

    def __init__(self,                                               \
                 execution_times = np.array([], dtype=ctypes.c_double), \
                 control_times   = np.array([], dtype=ctypes.c_double)):
        self.execuction_times = execution_times
        self.control_times    = control_times



'''
  Data structure to compare execution times for raw matrix multiplications
  using cache-blocking, MATMUL Fortran intrinsic, BLAS, and naive matrix
  multiplication 
''' 
class MatmulBenchmark:
    blocked_times   = np.array([], dtype=ctypes.c_double)
    intrinsic_times = np.array([], dtype=ctypes.c_double)
    blas_times      = np.array([], dtype=ctypes.c_double)
    naive_times     = np.array([], dtype=ctypes.c_double)

    def __init__(self,                                                    \
                 blocked_times    = np.array([], dtype=ctypes.c_double),  \
                 intrinsic_times  = np.array([], dtype=ctypes.c_double),  \
                 blas_times       = np.array([], dtype=ctypes.c_double),  \
                 naive_times      = np.array([], dtype=ctypes.c_double)):
        self.blocked_times   = blocked_times
        self.intrinsic_times = intrinsic_times
        self.blas_times      = blas_times
        self.naive_times     = naive_times







'''
benchmark1

Tests the execution time of mat_mul_m routines vs. Fortran's MATMUL intrinsic
for statement:
    B=matmul(A,A) for A,B of shape m-by-m 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
'''
def benchmark1(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm1 = benchmarks.benchmark1
    bm1.argtypes=[ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_int),    \
                  ctypes.c_int,                    \
                  ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    times = np.zeros(n_tests, dtype=ctypes.c_double)
    ctrl  = np.zeros(n_tests, dtype=ctypes.c_double)

    # Pointers to data 
    ms_ptr   = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    times_ptr = times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ctrl_ptr  = ctrl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm1(times_ptr,ctrl_ptr,ms_ptr,n_tests,n_iter)

    # Package results
    return Benchmark(times,ctrl)





'''
benchmark2

Tests the execution time of mat_mul_m routines vs. Fortran's MATMUL intrinsic
for statement:
    A=matmul(A,A) for A of shape m-by-m 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
'''
def benchmark2(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm2 = benchmarks.benchmark2
    bm2.argtypes=[ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_int),    \
                  ctypes.c_int,                    \
                  ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    times = np.zeros(n_tests, dtype=ctypes.c_double)
    ctrl  = np.zeros(n_tests, dtype=ctypes.c_double)

    # Pointers to data 
    ms_ptr   = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    times_ptr = times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ctrl_ptr  = ctrl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm2(times_ptr,ctrl_ptr,ms_ptr,n_tests,n_iter)

    # Package results
    return Benchmark(times,ctrl)





'''
benchmark3

Tests the execution time of mat_mul_m routines vs. Fortran's MATMUL intrinsic
for statement:
    B=c1*matmul(At,A) + c2*matmul(Bt,B) for A,B of shape m-by-1 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
'''
def benchmark3(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm3 = benchmarks.benchmark3
    bm3.argtypes=[ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_int),    \
                  ctypes.c_int,                    \
                  ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    times   = np.zeros(n_tests, dtype=ctypes.c_double)
    ctrl    = np.zeros(n_tests, dtype=ctypes.c_double)

    # Pointers to data 
    ms_ptr    = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    times_ptr = times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ctrl_ptr  = ctrl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm3(times_ptr,ctrl_ptr,ms_ptr,n_tests,n_iter)

    # Package results
    return Benchmark(times,ctrl)





'''
benchmark4

Tests the execution time of mat_mul_m routines vs. Fortran's MATMUL intrinsic
for statement:
    B=matmul(A,A) for A,B of shape m-by-m 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
'''
def benchmark4(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm4 = benchmarks.benchmark4
    bm4.argtypes=[ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_double), \
                  ctypes.POINTER(ctypes.c_int),    \
                  ctypes.c_int,                    \
                  ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    times   = np.zeros(n_tests, dtype=ctypes.c_double)
    ctrl    = np.zeros(n_tests, dtype=ctypes.c_double)

    # Pointers to data 
    ms_ptr    = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    times_ptr = times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ctrl_ptr  = ctrl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm4(times_ptr,ctrl_ptr,ms_ptr,n_tests,n_iter)

    # Package results
    return Benchmark(times,ctrl)    





'''
cache_blocking

Tests the execution time of blocked matrix multiplication vs. naive matrix multiplication
    B=matmul(A,A) for A,B of shape m-by-m 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
Returns:
  Benchmark......Data structure with execution times and control times  
'''
def cache_blocking(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    cb = benchmarks.cache_blocking
    cb.argtypes=[ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_int),    \
                 ctypes.c_int,                    \
                 ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    times = np.zeros(n_tests, dtype=ctypes.c_double)
    ctrl  = np.zeros(n_tests, dtype=ctypes.c_double)

    # Pointers to data 
    ms_ptr    = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    times_ptr = times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ctrl_ptr  = ctrl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    cb(times_ptr,ctrl_ptr,ms_ptr,n_tests,n_iter)

    # Package results
    return Benchmark(times,ctrl)





'''
matmul_benchmark1

Tests the execution time of blocked matrix multiplication vs. naive matrix multiplication
vs. BLAS DGEMM, vs. MATMUL Fortran intrinsic
    B=matmul(A,A) for A,B of shape m-by-m 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
Returns:
  MatmulBenchmark......Data structure with execution times for each algorithm 
'''
def matmul_benchmark1(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm = benchmarks.matmul_benchmark1
    bm.argtypes=[ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_int),    \
                 ctypes.c_int,                    \
                 ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    blocked_times   = np.zeros(n_tests, dtype=ctypes.c_double)
    intrinsic_times = np.zeros(n_tests, dtype=ctypes.c_double)
    blas_times      = np.zeros(n_tests, dtype=ctypes.c_double)
    naive_times     = np.zeros(n_tests, dtype=ctypes.c_double)
   

    # Pointers to data 
    ms_ptr        = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    blocked_ptr   = blocked_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    intrinsic_ptr = intrinsic_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    blas_ptr      = blas_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    naive_ptr     = naive_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm(blocked_ptr, intrinsic_ptr, blas_ptr, naive_ptr, ms_ptr, n_tests, n_iter)

    # Package results
    return MatmulBenchmark(blocked_times, intrinsic_times, blas_times, naive_times)





'''
matmul_benchmark2

Tests the execution time of blocked matrix multiplication vs. naive matrix multiplication
vs. BLAS DGEMM, vs. MATMUL Fortran intrinsic
    B=matmul(At,A) for A of shape m-by-1 

Parameters:
  mat_sizes......NumPy array of matrix sizes 'm'
  n_iter.........Number of times to perform each test
Returns:
  MatmulBenchmark......Data structure with execution times for each algorithm 
'''
def matmul_benchmark2(mat_sizes,n_iter=1000):
    # Import binary
    benchmarks = ctypes.CDLL('benchmarks.'+extension) 

    # Alias benchmark routines 
    bm = benchmarks.matmul_benchmark2
    bm.argtypes=[ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_double), \
                 ctypes.POINTER(ctypes.c_int),    \
                 ctypes.c_int,                    \
                 ctypes.c_int]

    # Execution times and control times
    n_tests = len(mat_sizes)
    blocked_times   = np.zeros(n_tests, dtype=ctypes.c_double)
    intrinsic_times = np.zeros(n_tests, dtype=ctypes.c_double)
    blas_times      = np.zeros(n_tests, dtype=ctypes.c_double)
    naive_times     = np.zeros(n_tests, dtype=ctypes.c_double)
   

    # Pointers to data 
    ms_ptr        = mat_sizes.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    blocked_ptr   = blocked_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    intrinsic_ptr = intrinsic_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    blas_ptr      = blas_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    naive_ptr     = naive_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call benchmark subroutine
    bm(blocked_ptr, intrinsic_ptr, blas_ptr, naive_ptr, ms_ptr, n_tests, n_iter)

    # Package results
    return MatmulBenchmark(blocked_times, intrinsic_times, blas_times, naive_times)    










