import ctypes as ct
import numpy  as np
import os
import benchmarklib as bmlib
from pathlib  import Path
import matplotlib
from   matplotlib.lines import Line2D
matplotlib.use('Qt5Agg') #May be necessary for macOS
from matplotlib         import pyplot as plt
from matplotlib.legend  import Legend

# If running from 'benchmarks' directory, change to mat_mul_m root directory 
delim = '/'  #Linux and MacOS
#delim = '\\' #Windows
pwd  = os.getcwd()
cdir = pwd.split(delim)[-1] 
if (cdir == 'benchmarks'):
    os.chdir('..')    

# Create a double array, pass it to Fortran as a pointer
mat_sizes = np.array([4,8,16,32,64], dtype=ct.c_int)
#mat_sizes = np.array([4,8,16,32,64,128,256], dtype=ct.c_int)
n         = 1000

# Run Benchmarks
bm1 = bmlib.benchmark1(mat_sizes,n)
bm2 = bmlib.benchmark2(mat_sizes,n)
bm3 = bmlib.benchmark3(mat_sizes,n)
bm4 = bmlib.benchmark4(mat_sizes,n)

bm1_speedup = 100*(bm1.control_times-bm1.execuction_times)/bm1.control_times
bm2_speedup = 100*(bm2.control_times-bm2.execuction_times)/bm2.control_times
bm3_speedup = 100*(bm3.control_times-bm3.execuction_times)/bm3.control_times
bm4_speedup = 100*(bm4.control_times-bm4.execuction_times)/bm4.control_times

print("Benchmark1 execution time ratio:")
print(bm1.execuction_times/bm1.control_times)

print("Benchmark2 execution time ratio:")
print(bm2.execuction_times/bm2.control_times)

print("Benchmark3 execution time ratio:")
print(bm3.execuction_times/bm3.control_times)

print("Benchmark4 execution time ratio:")
print(bm4.execuction_times/bm4.control_times)


# Execution times
fig1, axs = plt.subplots(2,1, figsize=(10, 8))#, layout='constrained') 
fig1.suptitle(r'mat_mul_m versus MATMUL Instrinsic',fontsize=24)

ax = axs[0]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Execution Time [s]',fontsize=14)    
ax.set_title('$A,B \\in \\mathbb{R}^{N \\times N}, \\quad C,D \\in \\mathbb{R}^{1 \\times N}$')

ax.plot(mat_sizes,bm1.execuction_times,'k-o', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k')#,
        #label='mat_mul_m $B=AA$')
ax.plot(mat_sizes,bm1.control_times,'k--o', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k')#,
        #label='matmul $B=AA$')

ax.plot(mat_sizes,bm2.execuction_times,'r-s', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k')#,
        #label='mat_mul_m $A=AA$')
ax.plot(mat_sizes,bm2.control_times,'r--s', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k')#,
        #label='matmul $A=AA$') 

ax.plot(mat_sizes,bm3.execuction_times,'g-^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k')#,
        #label='mat_mul_m $B=c_1 A^T A + c_2 B^T B$')
ax.plot(mat_sizes,bm3.control_times,'g--^', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k')#,
        #label='matmul $B=c_1 A^T A + c_2 B^T B$')   

ax.plot(mat_sizes,bm4.execuction_times,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k')#,
        #label='mat_mul_m $B=c_1 A^T A + c_2 B^T B$')
ax.plot(mat_sizes,bm4.control_times,'b-->', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k')   
    
ax.set_yscale('log')
ax.set_xscale('log')   
ax.grid(True)      
leg1_elements = [Line2D([0], [0], color='k', marker='o', lw=1.5, markersize=6, markeredgecolor='k', label='$B=AA$'),
                 Line2D([0], [0], color='r', marker='s', lw=1.5, markersize=6, markeredgecolor='k', label='$A=AA$'),
                 Line2D([0], [0], color='g', marker='^', lw=1.5, markersize=6, markeredgecolor='k', label='$B=C^T C + D^T D + C^T D + D^TC$'),
                 Line2D([0], [0], color='b', marker='^', lw=1.5, markersize=6, markeredgecolor='k', label='$A=((AA)A)(A(AA))^T - (A(AA))^T((AA)A) + A$')]
                 
leg1 = Legend(ax, handles=leg1_elements, labels=['$B=AA$','$A=AA$','$B=C^T C + D^T D + C^T D + D^TC$','$A=((AA)A)(A(AA))^T - (A(AA))^T((AA)A) + A$'],
             loc='upper left', frameon=True)
leg2_elements = [Line2D([0], [0], color='gray', linestyle='-',  lw=1.5, label='mat_mul'),
                 Line2D([0], [0], color='gray', linestyle='--', lw=1.5, label='MATMUL')]
leg2 = Legend(ax, handles=leg2_elements, labels=['mat_mul','MATMUL'],
             loc='center left', frameon=True)                              
ax.add_artist(leg1)   
ax.add_artist(leg2)          
#ax.legend(handles=legend_elements)
#ax.legend(frameon=True)


# Speedup
ax = axs[1]
ax.set_xlabel(xlabel='$N$',fontsize=14)
ax.set_ylabel(ylabel='Speedup [%]',fontsize=14)     
ax.plot(mat_sizes,bm1_speedup,'k-o', 
        linewidth=1.5, 
        markersize=6,
        markeredgecolor='k',
        label='$B=AA$')
ax.plot(mat_sizes,bm2_speedup,'r-s', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='$A=AA$')  
              
ax.plot(mat_sizes,bm3_speedup,'g-^', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='$B=C^T C + D^T D + C^T D + D^TC$')    
        
ax.plot(mat_sizes,bm4_speedup,'b->', 
        linewidth=1.5, 
        markersize=6, 
        markeredgecolor='k',
        label='$A=((AA)A)(A(AA))^T - (A(AA))^T((AA)A) + A$')            
ax.set_yscale('linear')
ax.set_xscale('log')  
ax.set_ylim([-100,100])
ax.grid(True)                    
ax.legend(frameon=True)

#     ((AA)A)(A(AA))^T - (A(AA))^T((AA)A) + A

'''
 mat_mul(mat_mul(mat_mul(A,A),A),transpose(mat_mul(A,mat_mul(A,A)))) &
          - mat_mul(transpose(mat_mul(A,mat_mul(A,A))), mat_mul(mat_mul(A,A),A)) &
          + A
'''


plt.show(block=True)









