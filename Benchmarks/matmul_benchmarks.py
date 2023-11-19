'''
Benchmarks for matrix multiplication execution times comparing:
  - Fortran intrinsic MATMUL
  - BLAS DGEMM
  - Basic cache-blocking
  - Naive matrix multiplication

Run this script in the 'mat_mul_m' root directory (preferentially) or
from the 'benchmarks' subdirectory.  
'''
import ctypes as ct
import numpy as np
from pathlib import Path
import os
import benchmarklib as bmlib
import matplotlib
matplotlib.use('Qt5Agg') #May be necessary for macOS
from matplotlib import pyplot as plt

# If running from 'benchmarks' directory, change to mat_mul_m root directory 
delim = '/'  #Linux and MacOS
#delim = '\\' #Windows
pwd  = os.getcwd()
cdir = pwd.split(delim)[-1] 
if (cdir == 'benchmarks'):
    os.chdir('..')    


# Create a double array, pass it to Fotran as a pointer
#mat_sizes = np.array([4,8,16,32,64,128,256,512], dtype=ct.c_int)
mat_sizes = np.array([4,8,16,32,64,128], dtype=ct.c_int)
n         = 1000

# Benchmark 1
bm1 = bmlib.matmul_benchmark1(mat_sizes,n)
blocked_speedup1   = 100*(bm1.naive_times-bm1.blocked_times)/bm1.naive_times
intrinsic_speedup1 = 100*(bm1.naive_times-bm1.intrinsic_times)/bm1.naive_times
blas_speedup1      = 100*(bm1.naive_times-bm1.blas_times)/bm1.naive_times

# Benchmark 2
bm2 = bmlib.matmul_benchmark2(mat_sizes,n)
blocked_speedup2   = 100*(bm2.naive_times-bm2.blocked_times)/bm2.naive_times
intrinsic_speedup2 = 100*(bm2.naive_times-bm2.intrinsic_times)/bm2.naive_times
blas_speedup2      = 100*(bm2.naive_times-bm2.blas_times)/bm2.naive_times

print("Matrix Multiplication Benchmark1:")
print(bm1.naive_times)
print(bm1.blocked_times) 
print(bm1.intrinsic_times)
print(bm1.blas_times)


print("\n\nMatrix Multiplication Benchmark2:")
print(bm2.naive_times)
print(bm2.blocked_times) 
print(bm2.intrinsic_times)
print(bm2.blas_times)



# Results for Benchmark1
fig1, axs = plt.subplots(2,1, figsize=(10, 8))#, layout='constrained') 
fig1.suptitle(r'Matrix Multiplication $AB: \ A,B \in \BbbR^{N \times N}$',fontsize=24)

# Execution times
ax = axs[0]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Execution Time [s]',fontsize=14)     
ax.plot(mat_sizes,bm1.naive_times,'k-o', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='Naive')
ax.plot(mat_sizes,bm1.blocked_times,'r--s', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k',
        label='Blocked')
ax.plot(mat_sizes,bm1.intrinsic_times,'g-.^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='MATMUL')         
ax.plot(mat_sizes,bm1.blas_times,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='BLAS')    
ax.set_yscale('log')
ax.set_xscale('log')   
ax.grid(True)                    
ax.legend(frameon=True)


# Speedup
ax = axs[1]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Speedup [%]',fontsize=14)     
ax.plot(mat_sizes,blocked_speedup1,'r--s', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k',
        label='Blocked')
ax.plot(mat_sizes,intrinsic_speedup1,'g-.^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='MATMUL')         
ax.plot(mat_sizes,blas_speedup1,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='BLAS')    
ax.set_yscale('linear')
ax.set_xscale('log')   
ax.grid(True)                    
ax.legend(frameon=True)

plt.show(block=False)




# Results for Benchmark2
fig2, axs = plt.subplots(2,1, figsize=(10, 8))#, layout='constrained') 
fig2.suptitle(r'Matrix Multiplication $AB: \ A^T,B \in \BbbR^{1 \times N}$',fontsize=24)

# Execution times
ax = axs[0]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Execution Time [s]',fontsize=14)     
ax.plot(mat_sizes,bm2.naive_times,'k-o', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='Naive')
ax.plot(mat_sizes,bm2.blocked_times,'r--s', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k',
        label='Blocked')
ax.plot(mat_sizes,bm2.intrinsic_times,'g-.^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='MATMUL')         
ax.plot(mat_sizes,bm2.blas_times,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='BLAS')    
ax.set_yscale('log')
ax.set_xscale('log')   
ax.grid(True)                    
ax.legend(frameon=True)


# Speedup
ax = axs[1]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Speedup [%]',fontsize=14)     
ax.plot(mat_sizes,blocked_speedup2,'r--s', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k',
        label='Blocked')
ax.plot(mat_sizes,intrinsic_speedup2,'g-.^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='MATMUL')         
ax.plot(mat_sizes,blas_speedup2,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='BLAS')    
ax.set_yscale('linear')
ax.set_xscale('log')   
ax.grid(True)                    
ax.legend(frameon=True)

plt.show(block=True)