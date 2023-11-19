
# mat_mul_m
---
 `mat_mul_m` is a Fortran module that minimizes dynamic allocations during matrix multiplications while providing effectively identical syntax to Fortran intrinsic `MATMUL`

 ## Building mat_mul_m
 Call `make mat_mul_m.o` from the shell to compile the module as a static object or `make mat_mul_m` to generate a linkable shared object. Call `make benchmarks` to compile shared object containing benchmarks with C bindings that are callable from Python. Call `make all` to build `mat_mul_m` and `benchmarks` static and shared objects.

 ### OS Compatibility
 Before building the binaries, verify that the *shared object* file extension for your operating system in the **Makefile** and in **benchmarklib.py** (typically *so* for Linux, *dll* for Windows, and *dylib* for MacOS) is correct. In **Makefile**, ensure that `FC` is set to your desired Fortran compiler and that the appropriate optimization and debug flags are assigned to `FFLAGS` and that the correct *position independent code* flag is assigned to `PIC`. 
 *Note: this project was built on MacOS Venture 13.4.1 with the GNU compiler. Alternate OS and compilers have not yet been tested.*

 ### Linking with BLAS
 In order to link **BLAS** for optimal matrix multiplication, locate its file path on your machine and assign `BLAS_DIR` to this directory in **Makefile**. **BLAS** is included in **LAPACK** and is typically installed in Linux and MacOS by default.

 ### Single vs. Double Precision
 `mat_mul_m` uses working precision, `WP`, set from the `precision_m` module in **precision_m.f90**, which  is set to double precision by default. To use single precision, change `integer, parameter :: WP = DP` to `integer, parameter :: WP = SP` in **precision_m.f90**. In subroutine `matrix_multiply_blas`, located in **mat_mul_m.f90**, use **BLAS** subroutine `DGEMM` for double precision and `SGEMM` for single precision.


### Using mat_mul_m without BLAS
`mat_mul_m` includes a basic cache-blocked algorithm for matrix multiplication that does not require calls to external libraries. To use `mat_mul_m` without **BLAS**, modify subroutine `matrix_multiply`, located in **mat_mul_m.f90**, to use either subroutine `matrix_multiply_blocked` or intrinsic function `MATMUL` instead of `matrix_multiply_blas`. If no **BLAS** library is available, all references to **BLAS** in **Makefile** and **mat_mul_m.f90** must be removed. Additionally, some benchmarks will no longer work.

 