
module linked_list_m
    implicit none
    
    type, public :: linked_list_node_t
      class(*),                 pointer :: value => null()
      type(linked_list_node_t), pointer :: prev  => null()
      type(linked_list_node_t), pointer :: next  => null()
      contains
      !final :: linked_list_node_finalize
    end type linked_list_node_t
  
  
  
  
    type, public :: linked_list_t
      integer                           :: length = 0
      type(linked_list_node_t), pointer :: head => null()
      type(linked_list_node_t), pointer :: tail => null()
      contains
      procedure, public  :: append => linked_list_append
      procedure, public  :: erase => linked_list_erase_by_ptr
      procedure, public  :: erase_by_index => linked_list_erase_by_index 
      procedure, public  :: traverse => linked_list_traverse
      procedure, private :: cleanup => linked_list_cleanup
      final              :: linked_list_finalize
    end type linked_list_t
  
    interface linked_list_erase
      module procedure :: linked_list_erase_by_index
      module procedure :: linked_list_erase_by_ptr
    end interface
  
  
    contains
  
  
      subroutine linked_list_node_finalize(this)
        implicit none
        !----------------- Parameters --------------------
        type(linked_list_node_t), intent(inout) :: this
        !-------------- Local Variables ------------------
        if (associated(this%value)) then
          deallocate(this%value)
          nullify(this%value)
          nullify(this%prev)
          nullify(this%next) 
        endif   
      end subroutine linked_list_node_finalize
  
  
  
      subroutine linked_list_append(this,ptr)
        implicit none
        !------------------- Parameters --------------------
        class(linked_list_t), intent(inout), target  :: this
        class(*),             intent(in),    pointer :: ptr
        !----------------- Local Variables -----------------
        type(linked_list_node_t), pointer :: new_node_ptr
        !---------------------------------------------------
  
        allocate(new_node_ptr)
        new_node_ptr%value => ptr
        new_node_ptr%prev  => this%tail
        new_node_ptr%next  => null() 
  
        if (.NOT.associated(this%head)) then
          this%head => new_node_ptr
        else 
          this%tail%next => new_node_ptr  
        endif    
        this%tail => new_node_ptr
        this%length = this%length+1
      end subroutine linked_list_append
  
  
      subroutine linked_list_erase_by_index(this,index)
        implicit none
        !--------------------- Parameters ----------------------
        class(linked_list_t), intent(inout), target  :: this
        integer,              intent(in),    value   :: index
        !------------------- Local Variables -------------------
        type(linked_list_node_t), pointer :: node_ptr
        integer                           :: count
        !-------------------------------------------------------
  
        if ((index.LT.1).OR.(index.GT.this%length)) return
  
        count=1
        node_ptr => this%head
        do while (associated(node_ptr))
          count=count+1
          if (index.EQ.count) then
            !call this%erase(node_ptr)
            call linked_list_erase_by_ptr(this,node_ptr)
          endif
          node_ptr => node_ptr%next   
        enddo   
  
      end subroutine linked_list_erase_by_index  
  
  
  
  
      subroutine linked_list_traverse(this, iterator_func)
        !----------------- Parameters -----------------
        class(linked_list_t), intent(inout) :: this
        !--------------- Local Variables --------------
        type(linked_list_node_t), pointer :: node_ptr
        type(linked_list_node_t), pointer :: temp_ptr
        !----------------------------------------------
        interface
          subroutine iterator_func(node)
            import linked_list_node_t
            type(linked_list_node_t), pointer, intent(inout) :: node
          end subroutine iterator_func
        end interface
  
        if (associated(this%head)) then
          node_ptr => this%head
          do while (associated(node_ptr))
            !-- save pointer (in case iterator procedure modifies it) --
            nullify(temp_ptr)
            temp_ptr => node_ptr%next    
            call iterator_func(node_ptr)
            node_ptr => temp_ptr 
          enddo
        endif  
      end subroutine  
  
  
  
      subroutine linked_list_erase_by_ptr(this, node_ptr)
        !------------------------- Parameters ------------------------
        class(linked_list_t),     intent(inout)          :: this
        type(linked_list_node_t), intent(inout), pointer :: node_ptr
        !----------------------- Local Variables ---------------------
        type(linked_list_node_t), pointer :: current_node_ptr
        type(linked_list_node_t), pointer :: next_node_ptr
        type(linked_list_node_t), pointer :: prev_node_ptr
        !-------------------------------------------------------------
  
        current_node_ptr => this%head
        do while(associated(current_node_ptr))
          
          if (associated(current_node_ptr,node_ptr)) then
            prev_node_ptr => current_node_ptr%prev
            next_node_ptr => current_node_ptr%next
  
            if (associated(current_node_ptr,this%head)) then
              this%head => next_node_ptr
            endif 
            
            if (associated(current_node_ptr,this%tail)) then
              this%tail => prev_node_ptr
            endif 
  
            if (associated(prev_node_ptr)) then
              prev_node_ptr%next => next_node_ptr
            endif
            
            if (associated(next_node_ptr)) then
              next_node_ptr%prev => prev_node_ptr
            endif
  
            deallocate(current_node_ptr)
            nullify(current_node_ptr)
  
            exit
          endif
        enddo  
      end subroutine  linked_list_erase_by_ptr   
  
  
      subroutine linked_list_cleanup(this)
        implicit none  
        !----------------- Parameters -----------------
        class(linked_list_t), intent(inout) :: this
        !----------------------------------------------
  
        call this%traverse(destroy_all)
        nullify(this%head)
        nullify(this%tail)
  
        contains
  
          !-- wrapper procedure for traverse interface --
          subroutine destroy_all(node)
            type(linked_list_node_t), pointer, intent(inout) :: node 
            
            this%head => node%next
            deallocate(node)
            nullify(node)
            this%length = this%length - 1
          end subroutine destroy_all   
      end subroutine linked_list_cleanup
  
  
  
  
      subroutine linked_list_finalize(this)
        implicit none
        !----------------- Parameters --------------------
        type(linked_list_t), intent(inout) :: this
        !-------------- Local Variables ------------------
        call this%cleanup()
      end subroutine linked_list_finalize
  
end module linked_list_m
  
  
  
  
  