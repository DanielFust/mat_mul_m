
  
  subroutine do_something()
    use mat_mul_m
    implicit none
    real(WP), allocatable, dimension(:,:) :: A,B,C
    integer :: rows,cols,i
    !type(mat_mul_t) :: mm1, mm2
    
    rows=3
    cols=3

    allocate(B(rows,cols),C(rows,cols),A(rows,cols))
    B(1:rows,1:cols) = 1.0_WP
    A(1:rows,1:cols) = 1.0_WP

    !mm1 = mat_mul(B,B)
    !mm2 = mat_mul(B,B)



    write(*,*) 'B='
    do i=1,rows
      write(*,*) B(i,1:cols)
    enddo

    !B = mat_mul(B,B)
   !B = mat_mul(mm1,mm2)
    !B = mat_mul(B, (B+B)) !
    !C = matmul(B, B+B) * matmul(B,B+1.0_WP) - matmul(B,B)
    !B = mat_mul(B, B+B) * mat_mul(B,B+1.0_WP) - mat_mul(B,B)
    !A = mat_mul(A,mat_mul(A,A)) - mat_mul(A,mat_mul(A,A))
    A = mat_mul(mat_mul(A,A),mat_mul(B,B)) + mat_mul(mat_mul(A,B), mat_mul(B,A)) 
    !B=B+B

    write(*,*) 'mulitplication completed'

    write(*,*) 'A='
    do i=1,rows
      write(*,*) A(i,1:cols)
    enddo

    write(*,*) 'B='
    do i=1,rows
      write(*,*) B(i,1:cols)
    enddo

    write(*,*) 'C='
    do i=1,rows
      write(*,*) C(i,1:cols)
    enddo

  end subroutine do_something  
  
  
  
  program mat_mul_debug
    use precision_m
    use mat_mul_m
    implicit none
    real(WP), pointer :: A(:,:)
    real(WP), allocatable :: B(:,:)
    integer :: rows,cols,i
    
    rows=3
    cols=3


    call do_something()
    
  end program mat_mul_debug
  



