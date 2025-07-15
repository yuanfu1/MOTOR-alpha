!!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
! All kins of check's flag decimal position  set up and flag refresh
! HISTORY           : Origionally from RCNMP
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief

MODULE qc_flag_set_m

  USE Meteoro_constants

  IMPLICIT NONE

CONTAINS

  SUBROUTINE flag_replace(a, n, m)

! This subroutine is to replace the check's flag by use of 0,1,2 at the corresponding check position

    IMPLICIT NONE
    !  CLASS(qc_flag_set_t) :: this
    INTEGER :: b, c
    INTEGER, INTENT(inout) :: a   ! input and out flag
    INTEGER :: n                   ! kinds of check flag order position
    INTEGER :: m                   ! sigle check flag
    b = MOD(a, 10**n)
    c = INT(b / 10**(n - 1))
    a = a - c * 10**(n - 1) + m * 10**(n - 1)
  END SUBROUTINE flag_replace
!-------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE flag_decompse(a, n, m)
! This subroutine is to replace the check's flag by use of 0,1,2 at the corresponding check position

    IMPLICIT NONE
    ! CLASS(qc_flag_set_t) :: this
    INTEGER :: b
    INTEGER, INTENT(inout) :: a   ! input and out flag
    INTEGER :: n                   ! kinds of check flag order position
    INTEGER :: m                   ! sigle check flag
    b = MOD(a, 10**n)
    m = INT(b / 10**(n - 1))
!      print *,'a,n,m ', a,n,m
    RETURN
  END SUBROUTINE flag_decompse

END MODULE qc_flag_set_m
