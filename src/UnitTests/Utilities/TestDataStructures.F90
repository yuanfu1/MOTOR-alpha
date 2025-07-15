!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR.Utilities.Utility
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by YUANFU XIE  (yuanfu_xie@yahoo.com), 2024/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a test program for all data structure modules:
PROGRAM main
  USE kinds_m, ONLY: i_kind
  USE linkedList_m, ONLY: linkedList_t
  IMPLICIT NONE
  TYPE(linkedList_t), POINTER :: head => NULL()

  ! INTEGER(i_kind) numElements,i
  ! INTEGER(i_kind), ALLOCATABLE :: values(:)
  LOGICAL :: passed = .TRUE.

  ! call insert(head, 10)
  ! call insert(head, 20)
  ! call insert(head, 30)

  ! print *, "Linked list values:"
  ! call to1DArray(head,numElements,values)
  ! DO i=1,numElements
  !   PRINT*,'Values: ',i,values(i)
  !   IF (values(i) .NE. 10*i) passed = .FALSE.
  ! END DO

  IF (passed) THEN
    WRITE (*, 1)
  ELSE
    WRITE (*, 2)
  END IF
1 FORMAT('Test passed')
2 FORMAT('Test failed')

  ! call free_list(head)
END PROGRAM main
