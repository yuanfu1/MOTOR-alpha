MODULE errorHandler_m

  !>
  !!=================================================================
  !!  This module handles all error occurances in this model.
  !!
  !!  \author Yuanfu Xie
  !!  \b History
  !!    2018-9 created by Yuanfu Xie
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double

  IMPLICIT NONE

  TYPE :: errorHandler_t
    CHARACTER(LEN=64) :: message
  CONTAINS
    PROCEDURE :: halt
  END TYPE errorHandler_t

CONTAINS
  SUBROUTINE halt(this, ierror)
    CLASS(errorHandler_t), INTENT(IN) :: this
    INTEGER(i_kind), INTENT(IN) :: ierror

    WRITE (*, 1)
    WRITE (*, 2) this%message, ierror
    WRITE (*, 1)
1   FORMAT('#######################################')
2   FORMAT('Error: ', A32, ' with error code: ', I4)
    STOP
  END SUBROUTINE halt

END MODULE errorHandler_m
