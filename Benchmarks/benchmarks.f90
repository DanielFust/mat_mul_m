
!------------------------------------------------------------------------------
! subroutine benchmark1
!
! Times mat_mul_m versus intrinsice for matrix multiplication statement:
!    B = matmul(A,A)
!
! Parameters:
!   exec_times.......execution times array for mat_mul_m routines
!   control_times....execution times array for matmul intrinsic
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine benchmark1(exec_times, control_times, mat_sizes, n_tests, N) bind(C,name="benchmark1")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: exec_times(n_tests)
  real(c_double), intent(out)         :: control_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer               :: i,j,m
  real(WP), allocatable :: A(:,:), B(:,:)
  real(WP)              :: start_time, stop_time
  !----------------------------------------------------------------

  do i=1,n_tests
    m = mat_sizes(i)

    allocate(A(m,m),B(m,m))
    !write(*,*) "setting scratch array size..."
    call set_scratch_array_size(m*m)
    call update_scratch_memory()
    
    !-- Using mat_mul_m operations --
    !write(*,*) "about to measure cpu time..."
    call cpu_time(start_time)
    !write(*,*) "just did it homie..."
    do j=1,N
        !B = mat_mul(A,A) + mat_mul(A,B) + mat_mul(B,B)
        !A = mat_mul(A,mat_mul(A,A)) - mat_mul(A,mat_mul(A,A))
      B = mat_mul(A,A)
    enddo
    call cpu_time(stop_time)
    exec_times(i) = stop_time - start_time

    !-- Using intrinsic matmul (Control test) --
    call cpu_time(start_time)
    do j=1,N
      !B = matmul(A,A) + matmul(A,B) + matmul(B,B)
      !A = matmul(A,matmul(A,A)) - matmul(A,matmul(A,A))
      B = matmul(A,A)
    enddo
    call cpu_time(stop_time)
    control_times(i) = stop_time - start_time

    deallocate(A,B)
  enddo

end subroutine benchmark1








!------------------------------------------------------------------------------
! subroutine benchmark2
!
! Times mat_mul_m versus intrinsice for matrix multiplication statement:
!    A = matmul(A,A)
!
! Parameters:
!   exec_times.......execution times array for mat_mul_m routines
!   control_times....execution times array for matmul intrinsic
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine benchmark2(exec_times, control_times, mat_sizes, n_tests, N) bind(C,name="benchmark2")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: exec_times(n_tests)
  real(c_double), intent(out)         :: control_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer               :: i,j,m
  real(WP), allocatable :: A(:,:), B(:,:)
  real(WP)              :: start_time, stop_time
  !----------------------------------------------------------------

  do i=1,n_tests
    m = mat_sizes(i)

    allocate(A(m,m),B(m,m))
    !write(*,*) "setting scratch array size..."
    call set_scratch_array_size(m*m)
    call update_scratch_memory()
    
    !-- Using mat_mul_m operations --
    !write(*,*) "about to measure cpu time..."
    call cpu_time(start_time)
    !write(*,*) "just did it homie..."
    do j=1,N
      A(1:m,1:m) = 0.0_WP 
      A = mat_mul(A,A)
    enddo
    call cpu_time(stop_time)
    exec_times(i) = stop_time - start_time

    !-- Using intrinsic matmul (Control test) --
    call cpu_time(start_time)
    do j=1,N
      A(1:m,1:m) = 0.0_WP 
      A = matmul(A,A)
    enddo
    call cpu_time(stop_time)
    control_times(i) = stop_time - start_time

    deallocate(A,B)
  enddo

end subroutine benchmark2










!------------------------------------------------------------------------------
! subroutine benchmark3
!
! Times mat_mul_m versus intrinsice for matrix multiplication statement:
!    B = c1*A'A + c2*B'B for A,B with dimension(m,1)
!
! Parameters:
!   exec_times.......execution times array for mat_mul_m routines
!   control_times....execution times array for matmul intrinsic
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine benchmark3(exec_times, control_times, mat_sizes, n_tests, N) bind(C,name="benchmark3")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: exec_times(n_tests)
  real(c_double), intent(out)         :: control_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer                               :: i,j,m
  real(WP), allocatable, dimension(:,:) :: K, A, At, B, Bt
  real(WP)                              :: start_time, stop_time
  real(WP)                              :: c1, c2
  !----------------------------------------------------------------

  !write(*,*) 'in benchmark3...'
  do i=1,n_tests
    m = mat_sizes(i)
    allocate(A(m,1),At(1,m),B(m,1),Bt(1,m),K(m,m))

    call set_scratch_array_size(m*m)
    call update_scratch_memory()
    
    A(1:m,1)  = 1.0_WP
    At(1,1:m) = 1.0_WP
    B(1:m,1)  = 1.0_WP
    Bt(1,1:m) = 1.0_WP
    c1 =  2.0_WP
    c2 = -1.0_WP

    !-- Using mat_mul_m operations --
    call cpu_time(start_time)
    do j=1,N
      K = mat_mul(At,A) + mat_mul(Bt,B) + mat_mul(At,B) + mat_mul(Bt,A)
      !K = mat_mul(mat_mul(At,A),mat_mul(Bt,B)) + mat_mul(mat_mul(At,B), mat_mul(Bt,A)) 
      !K = c1*mat_mul(At,A) + c2*mat_mul(Bt,B) 
    enddo
    call cpu_time(stop_time)
    exec_times(i) = stop_time - start_time

    !-- Using intrinsic matmul (Control test) --
    call cpu_time(start_time)
    do j=1,N
      K = matmul(At,A) + matmul(Bt,B) + matmul(At,B) + matmul(Bt,A) 
      !K = matmul(matmul(At,A),matmul(Bt,B)) + matmul(matmul(At,B), matmul(Bt,A))
      !K = c1*matmul(At,A) + c2*matmul(Bt,B)
    enddo
    call cpu_time(stop_time)
    control_times(i) = stop_time - start_time

    deallocate(A,At,B,Bt,K)
  enddo

end subroutine benchmark3









!------------------------------------------------------------------------------
! subroutine benchmark4
!
! Stress test. 
!
! Parameters:
!   exec_times.......execution times array for mat_mul_m routines
!   control_times....execution times array for matmul intrinsic
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine benchmark4(exec_times, control_times, mat_sizes, n_tests, N) bind(C,name="benchmark4")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: exec_times(n_tests)
  real(c_double), intent(out)         :: control_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer                               :: i,j,m
  real(WP), allocatable, dimension(:,:) :: A,B
  real(WP)                              :: start_time, stop_time
  real(WP)                              :: c1, c2
  !----------------------------------------------------------------

  !write(*,*) 'in benchmark4...'
  do i=1,n_tests
    m = mat_sizes(i)
    allocate(A(m,m),B(m,m))

    call set_scratch_array_size(m*m)
    call update_scratch_memory()
    
    A(m,m) = 1.0_WP

    !-- Using mat_mul_m operations --
    call cpu_time(start_time)
    do j=1,N
      A =   mat_mul(mat_mul(mat_mul(A,A),A),transpose(mat_mul(A,mat_mul(A,A)))) &
          - mat_mul(transpose(mat_mul(A,mat_mul(A,A))), mat_mul(mat_mul(A,A),A)) &
          + A
    enddo
    call cpu_time(stop_time)
    exec_times(i) = stop_time - start_time

    !-- Using intrinsic matmul (Control test) --
    call cpu_time(start_time)
    do j=1,N
      A =   matmul(matmul(matmul(A,A),A),transpose(matmul(A,matmul(A,A)))) &
          - matmul( transpose( matmul(A,matmul(A,A))), matmul(matmul(A,A),A)) &
          + A
    enddo
    call cpu_time(stop_time)
    control_times(i) = stop_time - start_time

    deallocate(A,B)
  enddo

end subroutine benchmark4











!------------------------------------------------------------------------------
! subroutine cache_blocking
!
! Test for speed of cache-blocked matrix multiplication versus naive matrix
! multiplication
!
! Parameters:
!   blocked_times....execution times array for cache-blocked algorithm
!   control_times....execution times array for naive algorithm
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine cache_blocking(blocked_times, &
                          control_times, &
                          mat_sizes,     &
                          n_tests,       &
                          N) bind(C,name="cache_blocking")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: blocked_times(n_tests)
  real(c_double), intent(out)         :: control_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer               :: i,j,m
  real(WP), allocatable :: A(:,:), B(:,:)
  real(WP)              :: start_time, stop_time
  !----------------------------------------------------------------

  do i=1,n_tests
    m = mat_sizes(i)

    allocate(A(m,m),B(m,m))
    A(1:m,1:m) = 1.0_WP
    B(1:m,1:m) = 1.0_WP

    call set_scratch_array_size(m*m)
    call update_scratch_memory()
    
    !-- Using mat_mul_m operations --
    call cpu_time(start_time)
 
    do j=1,N
      call matrix_multiply(B,A,A)
    enddo
    call cpu_time(stop_time)
    blocked_times(i) = stop_time - start_time

    !-- Using intrinsic matmul (Control test) --
    call cpu_time(start_time)
    do j=1,N
      call matrix_multiply_naive(B,A,A)
      !B=matmul(A,A)
    enddo
    call cpu_time(stop_time)
    control_times(i) = stop_time - start_time

    deallocate(A,B)
  enddo

end subroutine cache_blocking








!------------------------------------------------------------------------------
! subroutine matmul_benchmark1
!
! Test for speed of: 
!   -cache-blocked matrix multiplication versus naive matrix
!   -intrinsic matmul
!   -BLAS DGEMM matrix multiplication
!   -navie algorithm
!
! On statement: B=matmul(A,A) for A,B with dimension(m,m)
!
! Parameters:
!   blocked_times....execution times array for cache-blocked algorithm
!   intrinsic_times..execution times array for MATMUL intrinsic
!   blas_times.......execution times array for BLAS subroutine
!   naive_times......execution times array for naive algorithm
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine matmul_benchmark1(blocked_times,   &
                             intrinsic_times, &
                             blas_times,      &
                             naive_times,     &
                             mat_sizes,       &
                             n_tests,         &
                             N) &
                            bind(C,name="matmul_benchmark1")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: blocked_times(n_tests)
  real(c_double), intent(out)         :: intrinsic_times(n_tests)
  real(c_double), intent(out)         :: blas_times(n_tests)
  real(c_double), intent(out)         :: naive_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer               :: i,j,m
  real(WP), allocatable :: A(:,:), B(:,:)
  real(WP)              :: start_time, stop_time
  !----------------------------------------------------------------

  do i=1,n_tests
    m = mat_sizes(i)

    allocate(A(m,m),B(m,m))
    A(1:m,1:m) = 1.0_WP
    B(1:m,1:m) = 1.0_WP

    !-- Using mat_mul_m operations --
    call cpu_time(start_time)

    do j=1,N
      call matrix_multiply_blocked(B,A,A)
    enddo
    call cpu_time(stop_time)
    blocked_times(i) = stop_time - start_time

    !-- Using intrinsic matmul --
    call cpu_time(start_time)
    do j=1,N
      B=matmul(A,A)
    enddo
    call cpu_time(stop_time)
    intrinsic_times(i) = stop_time - start_time

    !-- Using BLAS matrix multiplication --
    call cpu_time(start_time)
    do j=1,N
      call matrix_multiply_blas(B,A,A)
    enddo
    call cpu_time(stop_time)
    blas_times(i) = stop_time - start_time

    !-- Using naive matrix multiplication (Control test) --
    call cpu_time(start_time)
    do j=1,N
      call matrix_multiply_naive(B,A,A)
      !B=matmul(A,A)
    enddo
    call cpu_time(stop_time)
    naive_times(i) = stop_time - start_time

    deallocate(A,B)
  enddo

end subroutine matmul_benchmark1








!------------------------------------------------------------------------------
! subroutine matmul_benchmark2
!
! Test for speed of: 
!   -cache-blocked matrix multiplication versus naive matrix
!   -intrinsic matmul
!   -BLAS DGEMM matrix multiplication
!   -navie algorithm
!
! On statement: B=matmul(A',A) for A with dimension(m,1)
!
! Parameters:
!   blocked_times....execution times array for cache-blocked algorithm
!   intrinsic_times..execution times array for MATMUL intrinsic
!   blas_times.......execution times array for BLAS subroutine
!   naive_times......execution times array for naive algorithm
!   mat_sizes........matrix size array
!   n_tests..........number of tests (length of arrays)
!   N................number of times to perform each benchmark operation
!------------------------------------------------------------------------------
subroutine matmul_benchmark2(blocked_times,   &
                             intrinsic_times, &
                             blas_times,      &
                             naive_times,     &
                             mat_sizes,       &
                             n_tests,         &
                             N) &
                             bind(C,name="matmul_benchmark2")
  use mat_mul_m 
  use iso_c_binding
  implicit none
  !------------------------ Parameters ----------------------------
  real(c_double), intent(out)         :: blocked_times(n_tests)
  real(c_double), intent(out)         :: intrinsic_times(n_tests)
  real(c_double), intent(out)         :: blas_times(n_tests)
  real(c_double), intent(out)         :: naive_times(n_tests)
  integer(c_int), intent(in)          :: mat_sizes(n_tests)
  integer(c_int), intent(in),   value :: n_tests
  integer(c_int), intent(in),   value :: N
  !---------------------- Local Variables -------------------------
  integer               :: i,j,m
  real(WP), allocatable :: A(:,:), At(:,:), B(:,:)
  real(WP)              :: start_time, stop_time
  !----------------------------------------------------------------

  do i=1,n_tests
    m = mat_sizes(i)

    allocate(A(1,m),At(m,1),B(m,m))
    A(1,1:m)  = 1.0_WP
    At(1:m,1) = 1.0_WP

    !-- Using mat_mul_m operations --
    call cpu_time(start_time)

    do j=1,N
      call matrix_multiply_blocked(B,At,A)
    enddo
    call cpu_time(stop_time)
    blocked_times(i) = stop_time - start_time

    !-- Using intrinsic matmul --
    call cpu_time(start_time)
    do j=1,N
      B=matmul(At,A)
    enddo
    call cpu_time(stop_time)
    intrinsic_times(i) = stop_time - start_time

    !-- Using BLAS matrix multiplication --
    call cpu_time(start_time)
    do j=1,N
      call matrix_multiply_blas(B,At,A)
    enddo
    call cpu_time(stop_time)
    blas_times(i) = stop_time - start_time

    !-- Using naive matrix multiplication (Control test) --
    call cpu_time(start_time)
    do j=1,N
      call matrix_multiply_naive(B,At,A) 
    enddo
    call cpu_time(stop_time)
    naive_times(i) = stop_time - start_time

    deallocate(A,At,B)
  enddo

end subroutine matmul_benchmark2





