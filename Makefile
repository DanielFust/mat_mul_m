
# Fortran compiler
FC        := gfortran

#Optimization and debug flags
FFLAGS    := -O3 

#Position indepedent code flag
PIC       := -fPIC

#BLAS Directory
BLAS_DIR := /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework

#File Extension for shared objects
EXT := dylib #MacOS
#EXT := so   #Linux
#EXT := dll  #Windows

#Project directories
SRC       := ./src
EX        := ./Examples
DEBUG_DIR := ./Debug
BM_DIR    := ./Benchmarks





OBJS := precision_m.o linked_list_m.o mat_mul_m.o

mat_mul_m.o: precision_m.o linked_list_m.o $(SRC)/mat_mul_m.f90
	$(FC) -c  $(FFLAGS) $(PIC) -L$(BLAS_DIR) $(SRC)/mat_mul_m.f90 -o mat_mul_m.o

mat_mul_m: mat_mul_m.o precision_m.o linked_list_m.o
	$(FC) -shared  $(FFLAGS) $(PIC) -L$(BLAS_DIR) -o mat_mul_m.$(EXT) -lblas

linked_list_m.o: $(SRC)/linked_list_m.f90
	$(FC) -c $(FFLAGS) $(PIC) $(SRC)/linked_list_m.f90 -o linked_list_m.o

precision_m.o: $(SRC)/precision_m.f90
	$(FC) -c $(FFLAGS) $(PIC) $(SRC)/precision_m.f90 -o precision_m.o

mat_mul_debug: mat_mul_m.o precision_m.o
	$(FC) $(FFLAGS) $(OBJS) $(DEBUG_DIR)/mat_mul_debug.f90 -o mat_mul_debug -lblas

benchmarks.o: mat_mul_m.o precision_m.o
	$(FC) -c $(PIC) $(FFLAGS) $(BM_DIR)/benchmarks.f90 -o benchmarks.o			

benchmarks: benchmarks.o mat_mul_m.o precision_m.o
	$(FC) -shared $(FFLAGS) $(OBJS) benchmarks.o -o benchmarks.$(EXT) -lblas	

test_cache_blocking: mat_mul_m.o precision_m.o
	$(FC) $(FFLAGS) -Wall $(OBJS) $(DEBUG_DIR)/test_cache_blocking.f90 -o test_cache_blocking -lblas	

all: mat_mul_m benchmarks #test_cache_blocking mat_mul_debug		

clean:
	rm -f *.o *.$(EXT) debug.* benchmarks.* test_cache_blocking.*



