program example_state1
  use stdlib_linalg_state, only: linalg_state_type, LINALG_SUCCESS, LINALG_VALUE_ERROR, &
          operator(/=)
  implicit none
  type(linalg_state_type) :: err
  character(len=512) :: msg

  ! To create a state variable, we enter its integer state flag, followed by a list of variables 
  ! that will be automatically assembled into a formatted error message. No need to provide string formats
  write(msg, *) 'just an example with scalar ',&
                          'integer=',1,'real=',2.0,'complex=',(3.0,1.0),'and array [',1,2,3,'] inputs'
  err = linalg_state_type(LINALG_VALUE_ERROR, msg)

  ! Print flag
  print *, err%print()

  ! Check success
  print *, 'Check error: ',err%error()
  print *, 'Check flag : ',err /= LINALG_SUCCESS


end program example_state1
