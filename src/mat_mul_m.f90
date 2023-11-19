!------------------------------------------------------------------------------
! module mat_mul_m
!
! This module serves as a working proof-of-concept for syntactically elegant 
! matrix multiplication of arrays that reduces (or eliminates entirely) dynamic
! allocations.
!
! mat_mul_m uses internal memory management to store all temporary results on
! one of (no less than) three preallocated scratch arrays. The size of these
! arrays is set by the module variable SCRATCH_MEMORY_SIZE and may be increased
! dynamically if DYNAMIC_ARRAY_SIZE is set to 'true'.
!
! The data type of arrays are to be of precision WP (working precision), which
! is set by an accompanying module 'precision_m'
!
! The interface 'transpose' extends the Fortran intrinsic to accomodate 
! internal data types of the mat_mul_m module. The interface 'mat_mul' is named
! thusly to distinguish from the Fortran intrinsic MATMUL.
!
! To change which matrix multiplication algorithm is used, modify the
! 'matrix_multiplication' subroutine to select for one of the following 
! subroutines:
!   matrix_multiply_blas
!   MATMUL (intrinsic calls a function rather than subroutine)
!   matrix_multiply_blocked
!   matrix_multiply_blocked_alt
!   matrix_multiply_naive
!
! For single precision matrix multiplication with BLAS, change subroutine call
! DGEMM to SGEMM in 'matrix_multiply_blas'
!------------------------------------------------------------------------------
module mat_mul_m
  use precision_m  
  use linked_list_m
  implicit none   

  !------------------- module constants -------------------------
  integer, parameter          :: DEFAULT_SCRATCH_MEMORY_SIZE = 256
  integer, parameter          :: DEFAULT_NUM_SCRATCH_ARRAYS  = 4
  logical, parameter, private :: DEBUG = .FALSE.
  integer, parameter          :: BLOCK_SIZE = 16
  !------------------- module variables -------------------------
  integer, save             :: SCRATCH_MEMORY_SIZE = DEFAULT_SCRATCH_MEMORY_SIZE
  integer, save             :: NUM_SCRATCH_ARRAYS  = DEFAULT_NUM_SCRATCH_ARRAYS
  logical, save             :: DYNAMIC_ARRAY_SIZE  = .TRUE.
  type(linked_list_t), save :: scratch_memory
  !--------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! Internal container for scratch memory. scratch_array_t objects own a scatch
  ! array of heap memory and tracks whether the array is currently in use by a
  ! temporary object (it does not track the object or objects that point to the
  ! scratch array).
  !----------------------------------------------------------------------------
  type, private :: scratch_array_t
    real(WP), allocatable :: data(:)
    logical               :: in_use = .FALSE.
    contains
    final                 :: scratch_array_destroy
  end type scratch_array_t

  !----------------------------------------------------------------------------
  ! Contains a pointer to scratch memory. Wrapping the pointer in this derived
  ! type signifies to mat_mul_t constructors that the array is temporary and 
  ! located in scratch memory, and therefore may need to be freed during the
  ! construction procedure.
  !----------------------------------------------------------------------------
  type, private :: temp_array_t
    real(WP), pointer :: data(:,:)
    contains
    procedure :: rows => temp_array_rows
    procedure :: cols => temp_array_cols
    final :: temp_array_finalize 
  end type
    
  
  !----------------------------------------------------------------------------
  ! User-facing interface for matrix multiplication
  !----------------------------------------------------------------------------
  interface mat_mul
    module procedure :: mat_mul_array_array
    module procedure :: mat_mul_temp_temp
    module procedure :: mat_mul_temp_array
    module procedure :: mat_mul_array_temp      
  end interface mat_mul

  !----------------------------------------------------------------------------
  ! Extends the addition operator to internal mat_mul_m data types
  !----------------------------------------------------------------------------
  interface operator(+)
    module procedure :: temp_plus_temp    
    module procedure :: temp_plus_array
    module procedure :: array_plus_temp    
    module procedure :: temp_plus_scalar
    module procedure :: scalar_plus_temp    
  end interface


  !----------------------------------------------------------------------------
  ! Extends the subtraction operator to internal mat_mul_m data types
  !----------------------------------------------------------------------------
  interface operator(-)
    module procedure :: temp_minus_temp    
    module procedure :: temp_minus_array
    module procedure :: array_minus_temp      
    module procedure :: temp_minus_scalar
    module procedure :: scalar_minus_temp    
  end interface


  !----------------------------------------------------------------------------
  ! Extends the multiplication operator (element-by-element) to internal 
  ! mat_mul_m data types
  !----------------------------------------------------------------------------
  interface operator(*)
    module procedure :: temp_times_temp
    module procedure :: temp_times_array
    module procedure :: array_times_temp
    module procedure :: temp_times_scalar
    module procedure :: scalar_times_temp
  end interface

  !----------------------------------------------------------------------------
  ! Extends the division operator (element-by-element) to internal 
  ! mat_mul_m data types
  !----------------------------------------------------------------------------
  interface operator(/)
    module procedure :: temp_divided_by_temp
    module procedure :: temp_divided_by_array
    module procedure :: array_divided_by_temp
    module procedure :: temp_divided_by_scalar
    module procedure :: scalar_divided_by_temp
  end interface

  !----------------------------------------------------------------------------
  ! Extends the tranpose intrinsic to internal mat_mul_m data types
  !----------------------------------------------------------------------------
  interface transpose
    module procedure :: transpose_temp_array
  end interface

  !----------------------------------------------------------------------------
  ! Assignment for internal mat_mul_m data types
  !----------------------------------------------------------------------------
  interface assignment (=)
    module procedure :: temp_array_assign
  end interface

  !----------------------------------------------------------------------------
  ! Constructor for temp_array_t data types
  !----------------------------------------------------------------------------
  interface temp_array
    module procedure :: temp_array_constructor
  end interface



  !-------------- Scratch Memory Control Procedures ---------------
  public :: set_scratch_array_size
  public :: set_dynamic_array_size
  public :: restore_defaults
  public :: deallocate_scratch_memory
  public :: update_scratch_memory

  !-------- Accessible Matrix Multiplication Procedures -----------
  public :: matrix_multiply
  public :: matrix_multiply_blas
  public :: matrix_multiply_naive
  public :: matrix_multiply_blocked
  public :: matrix_multiply_blocked_alt
  public :: additive_matmul

  !---------------------- Abort Execution -------------------------
  public :: graceful_exit

  !------------------ Internal Procedures -------------------------
  !------- Memory Management --------
  private :: reserve_temporary_memory
  private :: scratch_array_destroy
  private :: free_scratch_array
  private :: remove_scratch_memory_lock
  private :: add_scratch_array
  !------ mat_mul procedures ------
  private :: mat_mul_array_array
  private :: mat_mul_array_temp
  private :: mat_mul_temp_array
  private :: mat_mul_temp_temp
  !-------- Matrix Addition --------
  private :: temp_plus_temp    
  private :: temp_plus_array
  private :: array_plus_temp    
  private :: temp_plus_scalar
  private :: scalar_plus_temp
  !------ Matrix Subtraction -------
  private :: temp_minus_temp    
  private :: temp_minus_array
  private :: array_minus_temp      
  private :: temp_minus_scalar
  private :: scalar_minus_temp
  !---- Matrix Multiplication ------ 
  private :: temp_times_temp
  private :: temp_times_array
  private :: array_times_temp
  private :: temp_times_scalar
  private :: scalar_times_temp   
  !-------- Matrix Division -------- 
  private :: temp_divided_by_temp
  private :: temp_divided_by_array
  private :: array_divided_by_temp
  private :: temp_divided_by_scalar
  private :: scalar_divided_by_temp 
  !------ Transpose ------
  private :: transpose_temp_array
  !------ temp_array_t bound procedures ------   
  private :: temp_array_constructor
  private :: temp_array_assign
  private :: temp_array_rows
  private :: temp_array_cols

  
  contains


    !--------------------------------------------------------------------------
    ! subroutine set_scratch_array_size
    !
    ! Sets the size of the scratch arrays used in mat_mul_m. If the size is 
    ! insufficient, dynamic allocations may be required during execution.
    !
    ! Parameters:
    !   n........number of real(WP) elements in each scratch array
    !--------------------------------------------------------------------------
    subroutine set_scratch_array_size(n)
      implicit none
      !---------- Parameters ----------
      integer, intent(in), value :: n
      !--------------------------------
      SCRATCH_MEMORY_SIZE = n 
    end subroutine set_scratch_array_size
    
    
    !--------------------------------------------------------------------------
    ! subroutine set_dynamic_array_size
    !
    ! If .TRUE. the scratch array size may increase dynamically. If .FALSE. the
    ! scratch array size will reset to the previous size after execution.
    !
    ! Parameters:
    !   bool........T/F
    !--------------------------------------------------------------------------
    subroutine set_dynamic_array_size(bool)
      implicit none
      !------------ Parameters ------------
      logical, intent(in), value :: bool
      !------------------------------------
      DYNAMIC_ARRAY_SIZE = bool 
    end subroutine set_dynamic_array_size



    !--------------------------------------------------------------------------
    ! subroutine reset_defaults
    !
    ! Reset scratch array parameters to their default values
    !--------------------------------------------------------------------------
    subroutine restore_defaults()
      implicit none
      !------------------------------------
      NUM_SCRATCH_ARRAYS  = DEFAULT_NUM_SCRATCH_ARRAYS 
      SCRATCH_MEMORY_SIZE = DEFAULT_SCRATCH_MEMORY_SIZE
    end subroutine restore_defaults






    !--------------------------------------------------------------------------
    ! subroutine scratch_array_destroy
    !
    ! Deallocates scratch data
    !--------------------------------------------------------------------------
    subroutine scratch_array_destroy(this)
      implicit none
      !----------------------------------------------
      type(scratch_array_t), intent(inout) :: this
      !----------------------------------------------
      if (allocated(this%data)) deallocate(this%data)
    end subroutine scratch_array_destroy  


    !--------------------------------------------------------------------------
    ! subroutine deallocate_scratch_memory
    !
    ! Deallocates scratch memory arrays to free memory for other tasks. This
    ! may be called if there is a large gap in execution time before the next 
    ! usage of mat_mul
    !--------------------------------------------------------------------------
    subroutine deallocate_scratch_memory()
      implicit none
      !----------------------------------------------

      call scratch_memory%traverse(deallocate_scratch_array)
      
      contains
      subroutine deallocate_scratch_array(node)
        implicit none
        !----------------------- Parameters ---------------------
        type(linked_list_node_t), pointer, intent(inout) :: node
        !--------------------------------------------------------

        !-- type guard --
        select type (scratch => node%value)
          class is (scratch_array_t)
            if (allocated(scratch%data)) then
              deallocate(scratch%data)
              scratch%in_use = .FALSE.
            endif     
        end select
      end subroutine
    end subroutine deallocate_scratch_memory


    

    !--------------------------------------------------------------------------
    ! function reserve_temporary_memory
    !
    ! Returns a pointer with shape (/ rows, cols /) associated with module 
    ! scratch memory
    !
    ! Parameters:
    !   rows.....rows of requested temporary array
    !   cols.....columns of requested temporary array
    ! Returns:
    !   arr_ptr....pointer with shape (/rows,cols/) to matrix_m module scratch
    !              memory 
    !--------------------------------------------------------------------------
    function reserve_temporary_memory(rows,cols) result(arr_ptr)
      use precision_m  
      use linked_list_m
      implicit none
      !------------------------ Parameters ----------------------
      integer,  intent(in) :: rows, cols
      real(WP), pointer    :: arr_ptr(:,:)
      !---------------------- Local Variables -------------------
      integer                           :: arr_size, alloc_size
      type(linked_list_node_t), pointer :: node_ptr
      integer                           :: count
      !----------------------------------------------------------
      !write(*,*) 'in reserve_temporary_memory...'

      call update_scratch_memory()

      arr_ptr => null()

      arr_size = rows*cols
      alloc_size = MAX(arr_size,DEFAULT_SCRATCH_MEMORY_SIZE)

      !-- Loop through scratch arrays to find free memory --
      count = 0
      node_ptr => scratch_memory%head
      do while(associated(node_ptr))
        count=count+1
        !-- dynamic type guard --
        select type (scratch => node_ptr%value)
        class is (scratch_array_t)
          if (.NOT.scratch%in_use) then
            scratch%in_use = .TRUE.

            !-- if somehow the array is not already allocated --
            if (.NOT.allocated(scratch%data)) then
              allocate(scratch%data(alloc_size))

            !-- if array is larger than reserved memory --
            else if (size(scratch%data).LT.arr_size) then
              deallocate(scratch%data)
              allocate(scratch%data(alloc_size))
            endif
            
            !-- allow larger scratch arrays to be reserved --
            if (DYNAMIC_ARRAY_SIZE) SCRATCH_MEMORY_SIZE = alloc_size

            if (DEBUG) write(*,*) 'scratch array [',count,'] reserved'
            arr_ptr(1:rows,1:cols) => scratch%data(1:arr_size)
            exit
          endif
        end select  
        node_ptr => node_ptr%next
      enddo  
        
      if (.NOT.associated(arr_ptr)) then
        write(*,*) 'Error in reserve_temporary_memory: no scratch arrays are available for association'
        call graceful_exit(0)
      endif  
    end function reserve_temporary_memory



    !--------------------------------------------------------------------------
    ! subroutine free_scratch_array(arr_ptr)
    !
    ! Frees scratch array associated with arr_ptr from use
    !
    ! Parameters:
    !   arr_ptr.....pointer to array
    !--------------------------------------------------------------------------
    subroutine free_scratch_array(arr_ptr)
      use precision_m
      use linked_list_m
      implicit none
      !------------------- Parameters --------------------
      real(WP), pointer, intent(inout) :: arr_ptr(:,:) 
      !----------------- Local Variables -----------------
      real(WP), pointer                 :: ptr(:,:)
      type(linked_list_node_t), pointer :: node_ptr
      integer                           :: rows,cols
      integer                           :: count
      !---------------------------------------------------
      !write(*,*) 'in free_scratch_array...'

      ! DEV NOTE: 'associated' intrinsic can only compare pointer and
      !           target (or other pointer) with the same rank. the
      !           dummy pointer is used as a local alias for the scratch
      !           array with a compatible rank.

      rows = size(arr_ptr,1)
      cols = size(arr_ptr,2)

      count=0

      node_ptr => scratch_memory%head
      do while (associated(node_ptr))
        count=count+1
        !-- dynamic type guard --
        select type (scratch => node_ptr%value)
        class is (scratch_array_t)
          ptr(1:rows,1:cols) => scratch%data(1:rows*cols)
          if (associated(arr_ptr,ptr)) then
            scratch%in_use = .FALSE.
            nullify(arr_ptr)
            if (DEBUG) write(*,*) 'scratch array [',count,'] freed'
            exit
          endif
        end select
        node_ptr => node_ptr%next      
      enddo  

    end subroutine free_scratch_array



    !--------------------------------------------------------------------------
    ! subroutine remove_scratch_memory_lock(arr_ptr)
    !
    ! Removes the lock associated with arr_ptr so the scratch array may be used
    ! by other objects. DOES NOT NULLIFY arr_ptr.
    !
    ! WARNING: ONLY CALL THIS PROCEDURE ON TEMPORARY (CREATED AND DESTROYED IN 
    !          THE SAME LINE) ARRAY POINTERS AND ONLY IF THEY WILL NO LONGER
    !          BE REFERENCED. USE PROCEDURE 'free_scratch_array' INSTEAD 
    !          WHENEVER POSSIBLE.
    !
    ! Parameters:
    !   arr_ptr.....pointer to array
    !--------------------------------------------------------------------------
    subroutine remove_scratch_memory_lock(arr_ptr)
      use precision_m
      use linked_list_m
      implicit none
      !------------------- Parameters --------------------
      real(WP), pointer, intent(in) :: arr_ptr(:,:) 
      !----------------- Local Variables -----------------
      real(WP), pointer                 :: ptr(:,:)
      type(linked_list_node_t), pointer :: node_ptr
      integer                           :: rows,cols
      integer                           :: count
      !---------------------------------------------------

      !write(*,*) 'in remove_scratch_memory_lock...'

      rows = size(arr_ptr,1)
      cols = size(arr_ptr,2)

      count=0

      node_ptr => scratch_memory%head
      do while (associated(node_ptr))
        count=count+1
        !-- dynamic type guard --
        select type (scratch => node_ptr%value)
        class is (scratch_array_t)
          ptr(1:rows,1:cols) => scratch%data(1:rows*cols)
          if (associated(arr_ptr,ptr)) then
            scratch%in_use = .FALSE.
            if (DEBUG) write(*,*) 'scratch array [',count,'] freed'
            exit
          endif
        end select
        node_ptr => node_ptr%next      
      enddo  

    end subroutine remove_scratch_memory_lock

  
  
  
    !--------------------------------------------------------------------------
    ! function mat_mul_array_temp
    !
    ! matrix multiplication for two temp_array_t objects
    !
    ! Parameters:
    !   temp1.....temporary array type in scratch memory
    !   temp2.....temporary array type in scratch memory
    ! Result:
    !   matrix multiplcation result in scratch memeory (temp_array_t)
    !--------------------------------------------------------------------------
    function mat_mul_array_array(mat1,mat2) result(this)
      implicit none
      !-------------------- Parameters --------------------
      real(WP),          intent(in), target :: mat1(:,:)
      real(WP),          intent(in), target :: mat2(:,:)
      type(temp_array_t)                    :: this
      !------------------ Local Variables -----------------
      integer :: m1,n1,m2,n2,i,j
      !----------------------------------------------------

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      m2 = size(mat2,1)
      n2 = size(mat2,2)

      if (DEBUG) then
        if (n1.NE.m2) then
          write(*,*) 'Error in mat_mul_array_array: matrix dimensions incompatible for matrix multiplication'
          call graceful_exit(0)
        endif  
      endif
  
      !write(*,*) 'mat_mul_constructor_array_array called...'

      this = temp_array(m1,n2)
      call matrix_multiply(this%data,mat1,mat2)

    end function mat_mul_array_array




    !--------------------------------------------------------------------------
    ! function mat_mul_temp_temp
    !
    ! matrix multiplication for two temp_array_t objects
    !
    ! Parameters:
    !   temp1.....temporary array type in scratch memory
    !   temp2.....temporary array type in scratch memory
    ! Result:
    !   matrix multiplcation result in scratch memeory (temp_array_t)
    !--------------------------------------------------------------------------
    function mat_mul_temp_temp(temp1,temp2) result(this)
      implicit none
      !-------------------- Parameters --------------------
      type(temp_array_t), intent(in), target :: temp1
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: this
      !------------------ Local Variables -----------------
      integer :: m1,n1,m2,n2
      !----------------------------------------------------
      
      m1 = temp1%rows()
      n1 = temp1%cols()
      m2 = temp2%rows()
      n2 = temp2%cols()

      if (DEBUG) then
        if (n1.NE.m2) then
          write(*,*) 'Error in mat_mul_temp_temp: matrix dimensions incompatible for matrix multiplication'
          call graceful_exit(0)
        endif  
      endif
  
      !write(*,*) 'mat_mul_temp_temp called...'
      
      this = temp_array(m1,n2)
      call matrix_multiply(this%data,temp1%data,temp2%data)
      call remove_scratch_memory_lock(temp1%data)
      call remove_scratch_memory_lock(temp2%data)

    end function mat_mul_temp_temp




    !--------------------------------------------------------------------------
    ! function mat_mul_temp_array
    !
    ! Matrix multiplication for temp_array_t object and raw array.
    !
    ! Parameters:
    !   temp1.....temporary array type in scratch memory
    !   mat2......array/matrix
    ! Result:
    !   matrix multiplcation result in scratch memory (temp_array_t)
    !--------------------------------------------------------------------------
    function mat_mul_temp_array(temp1,mat2) result(this)
      implicit none
      !-------------------- Parameters --------------------
      type(temp_array_t), intent(in), target :: temp1
      real(WP),           intent(in), target :: mat2(:,:)
      type(temp_array_t)                     :: this
      !------------------ Local Variables -----------------
      integer            :: m1,n1,m2,n2
      type(temp_array_t) :: mat_copy
      !----------------------------------------------------
      
      m1 = temp1%rows()
      n1 = temp1%cols()
      m2 = size(mat2,1)
      n2 = size(mat2,2)
      
      if (DEBUG) then
        if (n1.NE.m2) then
          write(*,*) 'Error in mat_mul_temp_array: matrix dimensions incompatible for matrix multiplication'
          call graceful_exit(0)
        endif  
      endif
  
      !write(*,*) 'mat_mul_temp_array called...'
      
      this = temp_array(m1,n2)
      call matrix_multiply(this%data,temp1%data,mat2)
      call remove_scratch_memory_lock(temp1%data)

    end function mat_mul_temp_array



    !--------------------------------------------------------------------------
    ! function mat_mul_array_temp
    !
    ! Constructor for mat_mul_t object for raw array and temp_array_t object.
    !
    ! Parameters:
    !   mat1......array/matrix
    !   temp2.....temporary array type in scratch memory
    ! Result:
    !   matrix multiplcation result in scratch memeory (temp_array_t)
    !--------------------------------------------------------------------------
    function mat_mul_array_temp(mat1,temp2) result(this)
      implicit none
      !-------------------- Parameters --------------------
      real(WP),           intent(in), target :: mat1(:,:)
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: this
      !------------------ Local Variables -----------------
      integer            :: m1,n1,m2,n2
      type(temp_array_t) :: mat_copy
      !----------------------------------------------------
      
      m2 = temp2%rows()
      n2 = temp2%cols()
      m1 = size(mat1,1)
      n1 = size(mat1,2)

      if (DEBUG) then
        if (n1.NE.m2) then
          write(*,*) 'Error in mat_mul_constructor_array_temp: matrix dimensions incompatible for matrix multiplication'
          call graceful_exit(0)
        endif  
      endif
  
      !write(*,*) 'mat_mul_constructor_array_temp called...'

      this = temp_array(m1,n2)
      call matrix_multiply(this%data,mat1,temp2%data)
      call remove_scratch_memory_lock(temp2%data)
    
    end function mat_mul_array_temp



    !************************* BEGIN ADDITION BLOCK ***************************

    

    !--------------------------------------------------------------------------
    ! function temp_plus_temp
    !
    ! Addition of two temporary array temp_array_t objects.
    !
    ! Parameters:
    !   temp1......scratch matrix
    !   temp2......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_plus_temp(temp1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      type(temp_array_t), intent(in)         :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      !if (DEBUG) write(*,*) 'temp_plus_temp called...'
      
      m1=temp1%rows()
      n1=temp1%cols()
      m2=temp2%rows()
      n2=temp2%cols()

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp1%data(1,1) + temp2
        call remove_scratch_memory_lock(temp1%data)
      elseif ((m2.EQ.1).AND.(n2.EQ.1)) then 
        mat = temp1 + temp2%data(1,1) 
        call remove_scratch_memory_lock(temp2%data)
      else
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_plus_temp: incompatible array sizes for matrix addition.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of first array --
        mat%data => temp1%data
        mat%data(1:m1,1:n1) = temp1%data(1:m1,1:n1) + temp2%data(1:m1,1:n1)
        !-- Reclaim memory of second array --
        call remove_scratch_memory_lock(temp2%data)
      endif
    end function temp_plus_temp  
    


    !--------------------------------------------------------------------------
    ! function temp_plus_array
    !
    ! Addition of temporary array temp_array_t and array.
    !
    ! Parameters:
    !   temp......scratch matrix
    !   arr.......matrix/array
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_plus_array(temp,arr) result(mat)
      implicit none
      !------------------- Parameters --------------------  
      type(temp_array_t), intent(in), target :: temp
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t)                     :: mat
      !----------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !---------------------------------------------------

      if (DEBUG) write(*,*) 'temp_plus_array called...'
      
      m1=temp%rows()
      n1=temp%cols()
      m2=size(arr,1)
      n2=size(arr,2)

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp_array(m2,n2)
        mat%data(1:m2,1:n2) = temp%data(1,1) + arr(1:m2,1:n2)
        !-- Remove locks from temporary arrays --
        call remove_scratch_memory_lock(temp%data)
      else 
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_plus_array: incompatible array sizes for matrix addition.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = temp%data(1:m1,1:n1) + arr(1:m1,1:n1)
      endif  
    end function temp_plus_array
      
      
  

    !--------------------------------------------------------------------------
    ! function array_plus_temp
    !
    ! Addition of array and temporary array temp_array_t.
    !
    ! Parameters:
    !   arr.......matrix/array
    !   temp......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function array_plus_temp(arr,temp) result(mat)
      implicit none
      !----------------- Parameters --------------------  
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t), intent(in), target :: temp
      type(temp_array_t)                     :: mat
      !--------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'array_plus_temp called...'
      
      m1=size(arr,1)
      n1=size(arr,2)
      m2=temp%rows()
      n2=temp%cols()

      if ((m2.EQ.1).AND.(n2.EQ.1)) then
        mat = temp_array(m1,n2)
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) + temp%data(1,1)
        !-- Remove locks from temporary array --
        call remove_scratch_memory_lock(temp%data)
      else  
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in array_plus_temp: incompatible array sizes for matrix addition.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) + temp%data(1:m1,1:n1)
      endif  
    end function array_plus_temp
  



    !--------------------------------------------------------------------------
    ! function temp_plus_scalar
    !
    ! Addition of temp_array_t object and a scalar.
    !
    ! Parameters:
    !   temp1......temporary array object
    !   val2.......scalar
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_plus_scalar(temp1,val2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      real(WP),           intent(in), value  :: val2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_plus_scalar called...'
      
      m=temp1%rows()
      n=temp1%cols()

      mat%data => temp1%data
      mat%data(1:m,1:n) = temp1%data(1:m,1:n) + val2
    end function temp_plus_scalar



    !--------------------------------------------------------------------------
    ! function scalar_plus_temp
    !
    ! Addition of a scalar and temp_array_t object.
    !
    ! Parameters:
    !   val1.......scalar
    !   temp2......temporary array object
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function scalar_plus_temp(val1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      real(WP),           intent(in), value  :: val1
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'scalar_plus_temp called...'
      
      m=temp2%rows()
      n=temp2%cols()

      mat%data => temp2%data
      mat%data(1:m,1:n) = val1 + temp2%data(1:m,1:n)
    end function scalar_plus_temp

    !************************** END ADDITION BLOCK ****************************
  


    !*********************** BEGIN SUBTRACTION BLOCK **************************
    

    !--------------------------------------------------------------------------
    ! function temp_minus_temp
    !
    ! Subtraction of two temporary array temp_array_t objects.
    !
    ! Parameters:
    !   temp1......scratch matrix
    !   temp2......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_minus_temp(temp1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      type(temp_array_t), intent(in)         :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_minus_temp called...'
      
      m1=temp1%rows()
      n1=temp1%cols()
      m2=temp2%rows()
      n2=temp2%cols()

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp1%data(1,1) - temp2
        call remove_scratch_memory_lock(temp1%data)
      elseif ((m2.EQ.1).AND.(n2.EQ.1)) then 
        mat = temp1 - temp2%data(1,1) 
        call remove_scratch_memory_lock(temp2%data)
      else
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_minus_temp: incompatible array sizes for matrix subtraction.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of first array --
        mat%data => temp1%data
        mat%data(1:m1,1:n1) = temp1%data(1:m1,1:n1) - temp2%data(1:m1,1:n1)
        !-- Reclaim memory of second array --
        call remove_scratch_memory_lock(temp2%data)
      endif
    end function temp_minus_temp  
    


    !--------------------------------------------------------------------------
    ! function temp_minus_array
    !
    ! Subtraction of temporary array temp_array_t and array.
    !
    ! Parameters:
    !   temp......scratch matrix
    !   arr.......matrix/array
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_minus_array(temp,arr) result(mat)
      implicit none
      !------------------- Parameters --------------------  
      type(temp_array_t), intent(in), target :: temp
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t)                     :: mat
      !----------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !---------------------------------------------------

      if (DEBUG) write(*,*) 'temp_minus_array called...'
      
      m1=temp%rows()
      n1=temp%cols()
      m2=size(arr,1)
      n2=size(arr,2)

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp_array(m2,n2)
        mat%data(1:m2,1:n2) = temp%data(1,1) - arr(1:m2,1:n2)
        !-- Remove locks from temporary arrays --
        call remove_scratch_memory_lock(temp%data)
      else 
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_minus_array: incompatible array sizes for matrix subtraction.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = temp%data(1:m1,1:n1) - arr(1:m1,1:n1)
      endif  
    end function temp_minus_array
      
      
  

    !--------------------------------------------------------------------------
    ! function array_minus_temp
    !
    ! Subtraction of array and temporary array temp_array_t.
    !
    ! Parameters:
    !   arr.......matrix/array
    !   temp......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function array_minus_temp(arr,temp) result(mat)
      implicit none
      !----------------- Parameters --------------------  
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t), intent(in), target :: temp
      type(temp_array_t)                     :: mat
      !--------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'array_minus_temp called...'
      
      m1=size(arr,1)
      n1=size(arr,2)
      m2=temp%rows()
      n2=temp%cols()

      if ((m2.EQ.1).AND.(n2.EQ.1)) then
        mat = temp_array(m1,n2)
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) - temp%data(1,1)
        !-- Remove locks from temporary array --
        call remove_scratch_memory_lock(temp%data)
      else  
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in array_minus_temp: incompatible array sizes for matrix subtraction.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) - temp%data(1:m1,1:n1)
      endif  
    end function array_minus_temp
  




    !--------------------------------------------------------------------------
    ! function temp_minus_scalar
    !
    ! Subtraction of temp_array_t object and a scalar.
    !
    ! Parameters:
    !   temp1......temporary array object
    !   val2.......scalar
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_minus_scalar(temp1,val2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      real(WP),           intent(in), value  :: val2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_minus_scalar called...'
      
      m=temp1%rows()
      n=temp1%cols()

      mat%data => temp1%data
      mat%data(1:m,1:n) = temp1%data(1:m,1:n) - val2
    end function temp_minus_scalar



    !--------------------------------------------------------------------------
    ! function scalar_minus_temp
    !
    ! Subtraction of a scalar and temp_array_t object.
    !
    ! Parameters:
    !   val1.......scalar
    !   temp2......temporary array object
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function scalar_minus_temp(val1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      real(WP),           intent(in), value  :: val1
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'scalar_minus_temp called...'
      
      m=temp2%rows()
      n=temp2%cols()

      mat%data => temp2%data
      mat%data(1:m,1:n) = val1 - temp2%data(1:m,1:n)
    end function scalar_minus_temp
    !************************ END SUBTRACTION BLOCK ***************************

    !********************* BEGIN MULTIPLICATION BLOCK *************************
    

    !--------------------------------------------------------------------------
    ! function temp_times_temp
    !
    ! Multiplication of two temporary array temp_array_t objects.
    !
    ! Parameters:
    !   temp1......scratch matrix
    !   temp2......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_times_temp(temp1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      type(temp_array_t), intent(in)         :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_times_temp called...'
      
      m1=temp1%rows()
      n1=temp1%cols()
      m2=temp2%rows()
      n2=temp2%cols()

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp1%data(1,1) * temp2
        call remove_scratch_memory_lock(temp1%data)
      elseif ((m2.EQ.1).AND.(n2.EQ.1)) then 
        mat = temp1 * temp2%data(1,1) 
        call remove_scratch_memory_lock(temp2%data)
      else
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_times_temp: incompatible array sizes for matrix multiplication.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of first array --
        mat%data => temp1%data
        mat%data(1:m1,1:n1) = temp1%data(1:m1,1:n1) * temp2%data(1:m1,1:n1)
        !-- Reclaim memory of second array --
        call remove_scratch_memory_lock(temp2%data)
      endif
    end function temp_times_temp  
    


    !--------------------------------------------------------------------------
    ! function temp_times_array
    !
    ! Multiplication of temporary array temp_array_t and array.
    !
    ! Parameters:
    !   temp......scratch matrix
    !   arr.......matrix/array
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_times_array(temp,arr) result(mat)
      implicit none
      !------------------- Parameters --------------------  
      type(temp_array_t), intent(in), target :: temp
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t)                     :: mat
      !----------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !---------------------------------------------------

      if (DEBUG) write(*,*) 'temp_times_array called...'
      
      m1=temp%rows()
      n1=temp%cols()
      m2=size(arr,1)
      n2=size(arr,2)

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp_array(m2,n2)
        mat%data(1:m2,1:n2) = temp%data(1,1) * arr(1:m2,1:n2)
        !-- Remove locks from temporary arrays --
        call remove_scratch_memory_lock(temp%data)
      else 
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_times_array: incompatible array sizes for matrix multiplication.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = temp%data(1:m1,1:n1) * arr(1:m1,1:n1)
      endif  
    end function temp_times_array
      
      
  

    !--------------------------------------------------------------------------
    ! function array_times_temp
    !
    ! Multiplication of array and temporary array temp_array_t.
    !
    ! Parameters:
    !   arr.......matrix/array
    !   temp......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function array_times_temp(arr,temp) result(mat)
      implicit none
      !----------------- Parameters --------------------  
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t), intent(in), target :: temp
      type(temp_array_t)                     :: mat
      !--------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'array_times_temp called...'
      
      m1=size(arr,1)
      n1=size(arr,2)
      m2=temp%rows()
      n2=temp%cols()

      if ((m2.EQ.1).AND.(n2.EQ.1)) then
        mat = temp_array(m1,n2)
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) * temp%data(1,1)
        !-- Remove locks from temporary array --
        call remove_scratch_memory_lock(temp%data)
      else  
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in array_times_temp: incompatible array sizes for matrix multiplication.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) * temp%data(1:m1,1:n1)
      endif  
    end function array_times_temp
  




    !--------------------------------------------------------------------------
    ! function temp_times_scalar
    !
    ! Multiplication of temp_array_t object and a scalar.
    !
    ! Parameters:
    !   temp1......temporary array object
    !   val2.......scalar
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_times_scalar(temp1,val2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      real(WP),           intent(in), value  :: val2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_times_scalar called...'
      
      m=temp1%rows()
      n=temp1%cols()

      mat%data => temp1%data
      mat%data(1:m,1:n) = temp1%data(1:m,1:n) * val2
    end function temp_times_scalar



    !--------------------------------------------------------------------------
    ! function scalar_times_temp
    !
    ! Multiplication of a scalar and temp_array_t object.
    !
    ! Parameters:
    !   val1.......scalar
    !   temp2......temporary array object
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function scalar_times_temp(val1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      real(WP),           intent(in), value  :: val1
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'scalar_times_temp called...'
      
      m=temp2%rows()
      n=temp2%cols()

      mat%data => temp2%data
      mat%data(1:m,1:n) = val1 * temp2%data(1:m,1:n)
    end function scalar_times_temp
    !********************** END MULTIPLICATION BLOCK **************************

    !************************ BEGIN DIVISION BLOCK ****************************

    !--------------------------------------------------------------------------
    ! function temp_divided_by_temp
    !
    ! Division of two temporary array temp_array_t objects.
    !
    ! Parameters:
    !   temp1......scratch matrix
    !   temp2......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_divided_by_temp(temp1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      type(temp_array_t), intent(in)         :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_divided_by_temp called...'
      
      m1=temp1%rows()
      n1=temp1%cols()
      m2=temp2%rows()
      n2=temp2%cols()

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp1%data(1,1) / temp2
        call remove_scratch_memory_lock(temp1%data)
      elseif ((m2.EQ.1).AND.(n2.EQ.1)) then 
        mat = temp1 / temp2%data(1,1) 
        call remove_scratch_memory_lock(temp2%data)
      else
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_divided_by_temp: incompatible array sizes for matrix dvision.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of first array --
        mat%data => temp1%data
        mat%data(1:m1,1:n1) = temp1%data(1:m1,1:n1) / temp2%data(1:m1,1:n1)
        !-- Reclaim memory of second array --
        call remove_scratch_memory_lock(temp2%data)
      endif
    end function temp_divided_by_temp  
    


    !--------------------------------------------------------------------------
    ! function temp_divided_by_array
    !
    ! Division of temporary array temp_array_t and array.
    !
    ! Parameters:
    !   temp......scratch matrix
    !   arr.......matrix/array
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_divided_by_array(temp,arr) result(mat)
      implicit none
      !------------------- Parameters --------------------  
      type(temp_array_t), intent(in), target :: temp
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t)                     :: mat
      !----------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !---------------------------------------------------

      if (DEBUG) write(*,*) 'temp_divided_by_array called...'
      
      m1=temp%rows()
      n1=temp%cols()
      m2=size(arr,1)
      n2=size(arr,2)

      if ((m1.EQ.1).AND.(n1.EQ.1)) then
        mat = temp_array(m2,n2)
        mat%data(1:m2,1:n2) = temp%data(1,1) / arr(1:m2,1:n2)
        !-- Remove locks from temporary arrays --
        call remove_scratch_memory_lock(temp%data)
      else 
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in temp_divided_by_array: incompatible array sizes for matrix dvision.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = temp%data(1:m1,1:n1) / arr(1:m1,1:n1)
      endif  
    end function temp_divided_by_array
      
      
  

    !--------------------------------------------------------------------------
    ! function array_divided_by_temp
    !
    ! Division of array and temporary array temp_array_t.
    !
    ! Parameters:
    !   arr.......matrix/array
    !   temp......scratch matrix
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function array_divided_by_temp(arr,temp) result(mat)
      implicit none
      !----------------- Parameters --------------------  
      real(WP),           intent(in)         :: arr(:,:)
      type(temp_array_t), intent(in), target :: temp
      type(temp_array_t)                     :: mat
      !--------------- Local Variables -----------------
      integer :: m1,n1,m2,n2
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'array_divided_by_temp called...'
      
      m1=size(arr,1)
      n1=size(arr,2)
      m2=temp%rows()
      n2=temp%cols()

      if ((m2.EQ.1).AND.(n2.EQ.1)) then
        mat = temp_array(m1,n2)
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) / temp%data(1,1)
        !-- Remove locks from temporary array --
        call remove_scratch_memory_lock(temp%data)
      else  
        if (DEBUG) then
          if ((m1.NE.m2).OR.(n1.NE.n2)) then
            write(*,*) 'Error in array_divided_by_temp: incompatible array sizes for matrix dvision.'
            call graceful_exit(0)
          endif  
        endif  

        !-- assume ownership of data --
        mat%data => temp%data
        mat%data(1:m1,1:n1) = arr(1:m1,1:n1) / temp%data(1:m1,1:n1)
      endif  
    end function array_divided_by_temp
  

    !--------------------------------------------------------------------------
    ! function temp_divided_by_scalar
    !
    ! Division of temp_array_t object and a scalar.
    !
    ! Parameters:
    !   temp1......temporary array object
    !   val2.......scalar
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function temp_divided_by_scalar(temp1,val2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      type(temp_array_t), intent(in), target :: temp1
      real(WP),           intent(in), value  :: val2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'temp_divided_by_scalar called...'
      
      m=temp1%rows()
      n=temp1%cols()

      mat%data => temp1%data
      mat%data(1:m,1:n) = temp1%data(1:m,1:n) / val2
    end function temp_divided_by_scalar



    !--------------------------------------------------------------------------
    ! function scalar_divided_by_temp
    !
    ! Division of a scalar and temp_array_t object.
    !
    ! Parameters:
    !   val1.......scalar
    !   temp2......temporary array object
    ! Returns:
    !   mat.....resultant matrix as temp_array_t in scratch memory
    !--------------------------------------------------------------------------
    function scalar_divided_by_temp(val1,temp2) result(mat)
      implicit none
      !----------------- Parameters ------------------  
      real(WP),           intent(in), value  :: val1
      type(temp_array_t), intent(in), target :: temp2
      type(temp_array_t)                     :: mat
      !--------------- Local Variables ---------------
      integer :: m,n
      !-----------------------------------------------

      if (DEBUG) write(*,*) 'scalar_divided_by_temp called...'
      
      m=temp2%rows()
      n=temp2%cols()

      mat%data => temp2%data
      mat%data(1:m,1:n) = val1 / temp2%data(1:m,1:n)
    end function scalar_divided_by_temp

    !************************* END DIVISION BLOCK *****************************



    !--------------------------------------------------------------------------
    ! function temp_array_constructor
    !
    ! Constructs a temp_array_t object that reserves space in scratch memory
    !
    ! Parameters:
    !   rows......rows in temporary array
    !   cols......columns in temporary array
    ! Result:
    !   temp_array_t object pointing to scratch memory
    !--------------------------------------------------------------------------
    function temp_array_constructor(rows,cols) result(this)
      implicit none
      !------------------ Parameters -------------------
      integer,            intent(in) :: rows, cols
      type(temp_array_t)             :: this
      !-------------------------------------------------
      this%data => reserve_temporary_memory(rows,cols)
    end function temp_array_constructor  



    !--------------------------------------------------------------------------
    ! subroutine temp_array_finalize
    !
    ! Destructor for temporary array.
    !
    ! WARNING: this destructor does NOT free (lock removed from scratch array) 
    !          the scratch memory. This is because ownership of the data may 
    !          have been transfered to another object. SCRATCH MEMORY MUST BE 
    !          FREED MANUALLY IF IT IS NO LONGER IN USE. 
    !
    ! Dev Note: I hoped to wrap more memory cleanup into the destructors 
    !           ('final' subroutines) however compilers do not appear to 
    !           reliably call these subroutines when the object is no longer in
    !           use. By contrast, C/C++ will call destructors for temporary
    !           objects as soon as they are out of scope (triggers when they 
    !           are no longer used even in a temporary expression) but gfortran
    !           (and possibly other compilers) doesn't seem to finalize until
    !           the subroutine in which the object appears (even if it is 
    !           temporary) terminates.
    !--------------------------------------------------------------------------
    subroutine temp_array_finalize(this)
      implicit none
      !-----------------------------------------
      type(temp_array_t), intent(inout) :: this
      !-----------------------------------------
      nullify(this%data)
    end subroutine


    !--------------------------------------------------------------------------
    ! subroutine temp_array_assign
    !
    ! Assigns the contents of the temporary matrix to the array then frees
    ! the associated scratch memory (removes lock but doesn't nullify the 
    ! pointer due to 'intent(in)' requirement).
    !
    ! Parameters:
    !   mat......target array for assignment
    !   this.....temp_array_t scratch array
    !--------------------------------------------------------------------------
    subroutine temp_array_assign(mat,this)
      implicit none
      !------------------- Parameters --------------------
      real(WP),            intent(out), target :: mat(:,:)
      class(temp_array_t), intent(in),  target :: this
      !---------------- Local Variables ------------------
      integer :: m1,n1,m2,n2
      !---------------------------------------------------

      m1=size(mat,1)
      n1=size(mat,2)
      if (DEBUG) then
        m2=this%rows()
        n2=this%cols()
        if ((m1.NE.m2).OR.(n1.NE.n2)) then
          write(*,*) 'Error in temp_array_assign: incompatible array sizes.'
          call graceful_exit(0)
        endif  
      endif

      mat(1:m1,1:n1) = this%data(1:m1,1:n1)
      call remove_scratch_memory_lock(this%data)
    end subroutine temp_array_assign    



    !--------------------------------------------------------------------------
    ! function temp_array_rows
    !
    ! Parameters:
    ! Result:
    !   number of rows of the temporary array
    !--------------------------------------------------------------------------
    function temp_array_rows(this) result(rows)
      implicit none
      !----------------------------------------
      class(temp_array_t), intent(in) :: this
      integer                         :: rows
      !----------------------------------------
      rows = size(this%data,1)
    end function temp_array_rows



    !--------------------------------------------------------------------------
    ! function temp_array_cols
    !
    ! Parameters:
    ! Result:
    !   number of columns of the temporary array
    !--------------------------------------------------------------------------
    function temp_array_cols(this) result(cols)
      implicit none
      !----------------------------------------
      class(temp_array_t), intent(in) :: this
      integer                         :: cols
      !----------------------------------------
      cols = size(this%data,2)
    end function temp_array_cols



    !--------------------------------------------------------------------------
    ! function transpose_temp_array
    !
    ! Parameters:
    !   trans......array in scratch memory
    ! Result:
    !   transpose of input array (also located in temporary memory)
    !--------------------------------------------------------------------------
    function transpose_temp_array(temp) result(trans)
      !------------ Parameters ----------------
      type(temp_array_t), intent(in) :: temp
      type(temp_array_t)             :: trans
      !--------- Local Variables --------------
      integer :: m,n,i,j
      !----------------------------------------
      m = temp%rows()
      n = temp%cols()

      trans = temp_array(n,m)
      do j=1,n
        do i=1,m
          trans%data(j,i) = temp%data(i,j)
        enddo
      enddo  

      !-- Remove locks from temporary arrays --
      call remove_scratch_memory_lock(temp%data)
    end function transpose_temp_array
    
    



    !--------------------------------------------------------------------------
    ! subroutine init_scratch_arrays()
    !
    ! Allocate scratch arrays if they are not already allocated or the 
    ! requested number of arrays has been modified.
    !--------------------------------------------------------------------------
    subroutine update_scratch_memory()
      use precision_m
      use linked_list_m
      implicit none
      !-----------------------------------------
      
      !-- Add new scratch arrays if necessary --
      do while (scratch_memory%length .LT. NUM_SCRATCH_ARRAYS) 
        call add_scratch_array()
      enddo  

      !-- Resize dormant arrays --
      call scratch_memory%traverse(update_scratch_array) 

      contains 
      
      subroutine update_scratch_array(node)
        implicit none
        !----------------------- Parameters ---------------------
        type(linked_list_node_t), pointer, intent(inout) :: node
        !--------------------------------------------------------

        !-- type guard --
        select type (scratch => node%value)
          class is (scratch_array_t)
            !-- allocate scratch array if not allocated --
            if (.NOT.allocated(scratch%data)) then
              allocate(scratch%data(SCRATCH_MEMORY_SIZE))
              scratch%in_use = .FALSE.
            endif    

            !-- If array is not currently in use, allow resizing --
            if (.NOT.scratch%in_use) then
              if (size(scratch%data).NE.SCRATCH_MEMORY_SIZE) then
                deallocate(scratch%data)
                allocate(scratch%data(SCRATCH_MEMORY_SIZE))
                scratch%in_use = .FALSE.
              endif
            endif  
        end select
      end subroutine  
    end subroutine update_scratch_memory



    subroutine add_scratch_array()
      use precision_m
      use linked_list_m
      implicit none
      !----------------------------------------
      type(scratch_array_t)         :: scratch
      class(*),             pointer :: ptr
      !----------------------------------------
      allocate(ptr,source=scratch)
      call scratch_memory%append(ptr)
    end subroutine add_scratch_array  



    !--------------------------------------------------------------------------
    ! subroutine matrix_multiply
    !
    ! Performs a matrix multiplication of arrays 'mat1' and 'mat2' and store 
    ! the result in the 'result' array. This procedure assumes that 'result'
    ! is NOT 'mat1' or 'mat2'
    ! Dev Note: For speed, it would be best to replace this routine with a
    !           cache conscious matrix multiplication, or rebind it to a
    !           LAPACK subroutine.
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine matrix_multiply(result,mat1,mat2)
      use precision_m
      implicit none
      !------------------------ Parameters ------------------------
      real(WP), intent(out) :: result(:,:)
      real(WP), intent(in)  :: mat1(:,:)
      real(WP), intent(in)  :: mat2(:,:)
      !--------------------- Local Variables ----------------------
      integer  :: m1,n1,m2,n2,m,n
      !------------------------------------------------------------
      !result= matmul(mat1,mat2)
      !return

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      n2 = size(mat2,2)

      if (DEBUG) then
        m2 = size(mat2,1)
        m = size(result,1)
        n = size(result,1)
        
        if ((m.NE.m1).OR.(n.NE.n2).OR.(n1.NE.m2)) then
          write(*,*) 'Error in matrix_multiply: array shapes incompatible for matrix multiplicaiton...'
          call graceful_exit(0)
        endif  
      endif

      !-- Using intrinsic and BLAS for better performance --
      if ((m1*n2).LE.(BLOCK_SIZE*BLOCK_SIZE)) then
        result = matmul(mat1,mat2)
        !call matrix_multiply_naive(result,mat1,mat2)
      else
        call matrix_multiply_blas(result,mat1,mat2)
      endif
      return


      !-- Using only artisan homebrew code :)   --
      if ((m1*n2).LE.(BLOCK_SIZE*BLOCK_SIZE)) then
        call matrix_multiply_naive(result,mat1,mat2)
      else
        call matrix_multiply_blocked(result,mat1,mat2)
      endif

    end subroutine matrix_multiply  



    !--------------------------------------------------------------------------
    ! subroutine matrix_multiply_naive
    !
    ! Performs a matrix multiplication of arrays 'mat1' and 'mat2' and store 
    ! the result in the 'result' array. This procedure assumes that 'result'
    ! is NOT 'mat1' or 'mat2'
    !   -Does NOT check for compatible dimensions
    !   -Is NOT cache conscious
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine matrix_multiply_naive(result,mat1,mat2)
      use precision_m
      implicit none
      !------------------ Externals -------------------
      external :: dgemm
      !------------------ Parameters ------------------
      real(WP), intent(out) :: result(:,:)
      real(WP), intent(in)  :: mat1(:,:)
      real(WP), intent(in)  :: mat2(:,:)
      !---------------- Local Variables ---------------
      integer  :: i,j,k,m1,n1,n2
      real(WP) :: val
      !------------------------------------------------

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      n2 = size(mat2,2)

      do i=1,m1
        do j=1,n2
          val=0.0_WP
          do k=1,n1
            val = val + mat1(i,k)*mat2(k,j) 
          enddo
          result(i,j) = val
        enddo
      enddo    
    end subroutine matrix_multiply_naive


    !--------------------------------------------------------------------------
    ! subroutine matrix_multiply_blas
    !
    ! Maps to BLAS matrix multiplication subroutine
    !   -Does NOT check for compatible dimensions
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine matrix_multiply_blas(result,mat1,mat2)
      use precision_m
      implicit none
      !------------------ Externals -------------------
      external :: dgemm !Double precision
      external :: sgemm !Single precision
      !------------------ Parameters ------------------
      real(WP), intent(out) :: result(:,:)
      real(WP), intent(in)  :: mat1(:,:)
      real(WP), intent(in)  :: mat2(:,:)
      !---------------- Local Variables ---------------
      integer :: m1,n1,n2
      !------------------------------------------------

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      n2 = size(mat2,2)

      !-- Double Precision --
      call dgemm('n','n', m1, n2, n1, 1.0_WP,            &
                  mat1, m1, mat2, n1, 1.0_WP, result, m1)
      
      !-- Single Precision --
      !call dgemm('n','n', m1, n2, n1, 1.0_WP,            &
      !            mat1, m1, mat2, n1, 1.0_WP, result, m1)

    end subroutine matrix_multiply_blas






    !--------------------------------------------------------------------------
    ! subroutine matrix_multiply_blocked_alt
    !
    ! Performs a matrix multiplication of arrays 'mat1' and 'mat2' using basic 
    ! cache-blocking and stores the result in the 'result' array. This 
    ! procedure assumes that 'result' is NOT 'mat1' or 'mat2'
    !   -Does NOT check for compatible dimensions
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine matrix_multiply_blocked_alt(result,mat1,mat2)
      use precision_m
      implicit none
      !------------------ Parameters ------------------
      real(WP), intent(out) :: result(:,:)
      real(WP), intent(in)  :: mat1(:,:)
      real(WP), intent(in)  :: mat2(:,:)
      !---------------- Local Variables ---------------
      integer  :: m1,n1,n2
      integer  :: m1_blocks, n1_blocks, n2_blocks
      integer  :: mb, nb, nk, ib, jb, kb
      integer  :: i_start, i_stop, j_start, j_stop, &
                  k_start, k_stop
      real(WP) :: inv_bs            
      real(WP) :: partition(BLOCK_SIZE,BLOCK_SIZE)
      !------------------------------------------------

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      n2 = size(mat2,2)

      !m1_blocks = CEILING(real(m1)/real(BLOCK_SIZE))
      !n1_blocks = CEILING(real(n1)/real(BLOCK_SIZE))
      !n2_blocks = CEILING(real(n2)/real(BLOCK_SIZE))

      inv_bs = 1.0_WP/real(BLOCK_SIZE)
      m1_blocks = CEILING(m1*inv_bs)
      n1_blocks = CEILING(n1*inv_bs)
      n2_blocks = CEILING(n2*inv_bs)


      !-- Inner 'i' Blocks --
      mb = BLOCK_SIZE
      do ib=1,(m1_blocks-1)
        i_start = (ib-1)*BLOCK_SIZE + 1
        i_stop = i_start+BLOCK_SIZE-1
        nb = BLOCK_SIZE
        !-- Inner 'j' Blocks --
        do jb=1,(n2_blocks-1)
          j_start = (jb-1)*BLOCK_SIZE + 1
          j_stop = j_start+BLOCK_SIZE-1
          partition(1:mb,1:nb) = 0.0_WP

          !-- Inner 'k' Blocks --
          nk = BLOCK_SIZE
          do kb=1,(n1_blocks-1)
            k_start = (kb-1)*BLOCK_SIZE + 1
            k_stop  = k_start+BLOCK_SIZE-1
            call additive_matmul(partition(1:mb,1:nb),                &
                                  mat1(i_start:i_stop,k_start:k_stop), &
                                  mat2(k_start:k_stop,j_start:j_stop), &
                                  mb,nb,nk)
          enddo

          !-- Last k Block --
          k_start = (n1_blocks-1)*BLOCK_SIZE + 1
          k_stop  = n1 
          nk = 1+k_stop-k_start
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
          result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)                     
        enddo

        !-- Last j Block --
        j_start = (n2_blocks-1)*BLOCK_SIZE + 1
        j_stop  = n2 
        nb = 1+j_stop-j_start
        nk = BLOCK_SIZE
        partition(1:mb,1:nb) = 0.0_WP
        !-- Inner k Blocks --
        do kb=1,(n1_blocks-1)
          k_start = (kb-1)*BLOCK_SIZE + 1
          k_stop  = k_start+BLOCK_SIZE-1
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
        enddo

        !-- Last k Block --
        k_start = (n1_blocks-1)*BLOCK_SIZE + 1
        k_stop  = n1 
        nk = 1+k_stop-k_start
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
        result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)
      enddo !ib=1,(m1_blocks-1)

      !-- Last i Block --
      i_start = (m1_blocks-1)*BLOCK_SIZE + 1
      i_stop  = m1 
      mb = 1+i_stop-i_start
      nb = BLOCK_SIZE
      !-- j Blocks --
      do jb=1,(n2_blocks-1)
        j_start = (jb-1)*BLOCK_SIZE + 1
        j_stop = j_start+BLOCK_SIZE-1
        partition(1:mb,1:nb) = 0.0_WP

        !-- Inner k Blocks --
        nk = BLOCK_SIZE
        do kb=1,(n1_blocks-1)
          k_start = (kb-1)*BLOCK_SIZE + 1
          k_stop  = k_start+BLOCK_SIZE-1
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
        enddo

        !-- Last k Block --
        k_start = (n1_blocks-1)*BLOCK_SIZE + 1
        k_stop  = n1 
        nk = 1+k_stop-k_start
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
        result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)
      enddo

      !-- Last j Block --
      j_start = (n2_blocks-1)*BLOCK_SIZE + 1
      j_stop  = n2 
      nb = 1+j_stop-j_start
      nk = BLOCK_SIZE
      partition(1:mb,1:nb) = 0.0_WP
      !-- Inner k Blocks --
      do kb=1,(n1_blocks-1)
        k_start = (kb-1)*BLOCK_SIZE + 1
        k_stop  = k_start+BLOCK_SIZE-1
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
      enddo

      !-- Last k Block --
      k_start = (n1_blocks-1)*BLOCK_SIZE + 1
      k_stop  = n1 
      nk = 1+k_stop-k_start
      call additive_matmul(partition(1:mb,1:nb),                &
                            mat1(i_start:i_stop,k_start:k_stop), &
                            mat2(k_start:k_stop,j_start:j_stop), &
                            mb,nb,nk)
      result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)     
    end subroutine matrix_multiply_blocked_alt





    !--------------------------------------------------------------------------
    ! subroutine matrix_multiply_blocked
    !
    ! Performs a matrix multiplication of arrays 'mat1' and 'mat2' using 
    ! basic cache-blocking and stores the result in the 'result' array. 
    ! This procedure assumes that 'result' is NOT 'mat1' or 'mat2'
    !   -Does NOT check for compatible dimensions
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine matrix_multiply_blocked(result,mat1,mat2)
      use precision_m
      implicit none
      !------------------ Parameters ------------------
      real(WP), intent(out) :: result(:,:)
      real(WP), intent(in)  :: mat1(:,:)
      real(WP), intent(in)  :: mat2(:,:)
      !---------------- Local Variables ---------------
      integer  :: m1,n1,n2
      integer  :: m1_blocks, n1_blocks, n2_blocks
      integer  :: mb, nb, nk, ib, jb, kb
      integer  :: i_start, i_stop, j_start, j_stop, &
                  k_start, k_stop
      real(WP) :: inv_bs            
      real(WP) :: partition(BLOCK_SIZE,BLOCK_SIZE)
      !------------------------------------------------

      m1 = size(mat1,1)
      n1 = size(mat1,2)
      n2 = size(mat2,2)

      !m1_blocks = CEILING(real(m1)/real(BLOCK_SIZE))
      !n1_blocks = CEILING(real(n1)/real(BLOCK_SIZE))
      !n2_blocks = CEILING(real(n2)/real(BLOCK_SIZE))

      inv_bs = 1.0_WP/real(BLOCK_SIZE)
      m1_blocks = CEILING(m1*inv_bs)
      n1_blocks = CEILING(n1*inv_bs)
      n2_blocks = CEILING(n2*inv_bs)


      !-- Inner 'j' Blocks --
      nb = BLOCK_SIZE
      do jb=1,(n2_blocks-1)
        j_start = (jb-1)*BLOCK_SIZE + 1
        j_stop = j_start+BLOCK_SIZE-1
        mb = BLOCK_SIZE
        !-- Inner 'i' Blocks --
        do ib=1,(m1_blocks-1)
          i_start = (ib-1)*BLOCK_SIZE + 1
          i_stop = i_start+BLOCK_SIZE-1
          partition(1:mb,1:nb) = 0.0_WP

          !-- Inner 'k' Blocks --
          nk = BLOCK_SIZE
          do kb=1,(n1_blocks-1)
            k_start = (kb-1)*BLOCK_SIZE + 1
            k_stop  = k_start+BLOCK_SIZE-1
            call additive_matmul(partition(1:mb,1:nb),                &
                                  mat1(i_start:i_stop,k_start:k_stop), &
                                  mat2(k_start:k_stop,j_start:j_stop), &
                                  mb,nb,nk)
          enddo

          !-- Last k Block --
          k_start = (n1_blocks-1)*BLOCK_SIZE + 1
          k_stop  = n1 
          nk = 1+k_stop-k_start
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
          result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)                     
        enddo

        !-- Last i Block --
        i_start = (n2_blocks-1)*BLOCK_SIZE + 1
        i_stop  = n2 
        mb = 1+i_stop-i_start
        nk = BLOCK_SIZE
        partition(1:mb,1:nb) = 0.0_WP
        !-- Inner k Blocks --
        do kb=1,(n1_blocks-1)
          k_start = (kb-1)*BLOCK_SIZE + 1
          k_stop  = k_start+BLOCK_SIZE-1
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
        enddo

        !-- Last k Block --
        k_start = (n1_blocks-1)*BLOCK_SIZE + 1
        k_stop  = n1 
        nk = 1+k_stop-k_start
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
        result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)
      enddo !ib=1,(m1_blocks-1)

      !-- Last j Block --
      j_start = (n2_blocks-1)*BLOCK_SIZE + 1
      j_stop  = n2 
      nb = 1+j_stop-j_start
      mb = BLOCK_SIZE
      !-- i Blocks --
      do ib=1,(m1_blocks-1)
        i_start = (ib-1)*BLOCK_SIZE + 1
        i_stop = i_start+BLOCK_SIZE-1
        partition(1:mb,1:nb) = 0.0_WP

        !-- Inner k Blocks --
        nk = BLOCK_SIZE
        do kb=1,(n1_blocks-1)
          k_start = (kb-1)*BLOCK_SIZE + 1
          k_stop  = k_start+BLOCK_SIZE-1
          call additive_matmul(partition(1:mb,1:nb),                &
                                mat1(i_start:i_stop,k_start:k_stop), &
                                mat2(k_start:k_stop,j_start:j_stop), &
                                mb,nb,nk)
        enddo

        !-- Last k Block --
        k_start = (n1_blocks-1)*BLOCK_SIZE + 1
        k_stop  = n1 
        nk = 1+k_stop-k_start
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
        result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)
      enddo

      !-- Last i Block --
      i_start = (m1_blocks-1)*BLOCK_SIZE + 1
      i_stop  = m1 
      mb = 1+i_stop-i_start
      nk = BLOCK_SIZE
      partition(1:mb,1:nb) = 0.0_WP
      !-- Inner k Blocks --
      do kb=1,(n1_blocks-1)
        k_start = (kb-1)*BLOCK_SIZE + 1
        k_stop  = k_start+BLOCK_SIZE-1
        call additive_matmul(partition(1:mb,1:nb),                &
                              mat1(i_start:i_stop,k_start:k_stop), &
                              mat2(k_start:k_stop,j_start:j_stop), &
                              mb,nb,nk)
      enddo

      !-- Last k Block --
      k_start = (n1_blocks-1)*BLOCK_SIZE + 1
      k_stop  = n1 
      nk = 1+k_stop-k_start
      call additive_matmul(partition(1:mb,1:nb),                &
                            mat1(i_start:i_stop,k_start:k_stop), &
                            mat2(k_start:k_stop,j_start:j_stop), &
                            mb,nb,nk)
      result(i_start:i_stop, j_start:j_stop) = partition(1:mb,1:nb)     
    end subroutine matrix_multiply_blocked






    !--------------------------------------------------------------------------
    ! subroutine additive_matmul
    !
    ! Performs a naive matrix multiplication that ADDS to the 'result' matrix
    !
    ! Parameters:
    !   result.......result of matrix multiplication of mat1 and mat2
    !   mat1.........2D array (treated as a matrix)
    !   mat2.........2D array (treated as a matrix) 
    !--------------------------------------------------------------------------
    subroutine additive_matmul(result,mat1,mat2,m1,n2,n1)
      use precision_m
      implicit none
      !------------------ Parameters ------------------
      real(WP), intent(inout) :: result(:,:)
      real(WP), intent(in)    :: mat1(:,:)
      real(WP), intent(in)    :: mat2(:,:)
      integer,  intent(in), value :: m1,n2,n1
      !---------------- Local Variables ---------------
      integer :: i,j,k
      !------------------------------------------------

      do i=1,m1
        do j=1,n2
          do k=1,n1
            result(i,j) = result(i,j) + mat1(i,k)*mat2(k,j) 
          enddo
        enddo
      enddo    
    end subroutine additive_matmul



  

    subroutine graceful_exit(val)
      implicit none
      integer, intent(in), value :: val
      call abort()
    end subroutine graceful_exit
  
end module mat_mul_m
    
    
    