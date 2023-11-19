program test_cache_blocking
    use mat_mul_m
    use precision_m
    implicit none
    integer, parameter :: n=20
    integer :: i,j
    real(WP) :: A(n,n), B(n,n), C(n,n)
    real(WP) :: err
  
    do j=1,n
      do i=1,n
        A(i,j) = 1.0_WP!i*j + i*i - j*j
      enddo
    enddo
  
    call matrix_multiply_blocked(B,A,A)
    C = matmul(A,A)
  
    err = 0.0_WP
    do j=1,n
      do i=1,n
        err = err + (B(i,j) - C(i,j))**2 
      enddo
    enddo
  
    write(*,*) 'squared error=',err

    if (.FALSE.) then
      write(*,*) 'cache-blocked:'
      do i=1,n
        write(*,*) B(i,1:n)  
      enddo

      write(*,*) 'intrinsic:'
      do i=1,n
        write(*,*) C(i,1:n)  
      enddo
    end if
  end program
  