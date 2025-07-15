program example_state2
  !! This example shows how to set a `type(linalg_state_type)` variable to process output conditions 
  !! out of a simple division routine. The example is meant to highlight: 
  !! 1) the different mechanisms that can be used to initialize the `linalg_state` variable providing 
  !!    strings, scalars, or arrays, on input to it; 
  !! 2) `pure` setup of the error control
  use stdlib_linalg_state, only: linalg_state_type, LINALG_VALUE_ERROR, LINALG_SUCCESS, &
          linalg_error_handling
  implicit none
  type(linalg_state_type) :: err
  real :: a_div_b
  
  
  ! OK
  call very_simple_division(0.0,2.0,a_div_b,err)
  print *, err%print()
  
  ! Division by zero
  call very_simple_division(1.0,0.0,a_div_b,err)
  print *, err%print()  

  ! Out of bounds
  call very_simple_division(huge(0.0),0.001,a_div_b,err)
  print *, err%print()  
  
  contains
  
     !> Simple division returning an integer flag (LAPACK style)
     subroutine very_simple_division(a,b,a_div_b,err)
        real, intent(in) :: a,b
        real, intent(out) :: a_div_b
        type(linalg_state_type), optional, intent(out) :: err
        
        type(linalg_state_type) :: err0
        real, parameter :: MAXABS = huge(0.0)
        character(*), parameter :: this = 'simple division'
        character(len=512) :: msg
        
        !> Check a
        if (b==0.0) then 
            ! Division by zero
            write(msg, *) 'Division by zero trying ',a,'/',b
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,msg)
        elseif (.not.abs(b)<MAXABS) then 
            ! B is out of bounds
            write(msg, *) 'B is infinity in a/b: [',a,b,']'
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,msg) ! use an array
        elseif (.not.abs(a)<MAXABS) then 
            ! A is out of bounds
            write(msg, * ) 'A is infinity in a/b: a=',a,' b=',b
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR, msg)
        else
            a_div_b = a/b
            if (.not.abs(a_div_b)<MAXABS) then 
                ! Result is out of bounds
                write(msg, * ) 'A/B is infinity in a/b: a=',a,' b=',b
                err0 = linalg_state_type(this,LINALG_VALUE_ERROR,msg)
            else
                err0%state = LINALG_SUCCESS
            end if
        end if
        
        ! Return error flag, or hard stop on failure
        call linalg_error_handling(err0,err)
                
     end subroutine very_simple_division


end program example_state2
