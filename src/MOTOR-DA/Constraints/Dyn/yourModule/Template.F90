!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE Template_m

  TYPE Template_t
  CONTAINS
    FINAL :: destructor
  END TYPE Template_t

  INTERFACE Template_t
    PROCEDURE :: constructor
  END INTERFACE Template_t

CONTAINS

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(Template_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Template_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

END MODULE Template_m
