!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.AuxTypeObs
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/25, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This type enclosed the essential parameters in constructing types.
MODULE AuxTypeObs_m

  USE MPObs_m, ONLY: MPObs_t
  TYPE AuxTypeObs_t

    TYPE(MPObs_t), POINTER :: mpObs
  CONTAINS

  END TYPE AuxTypeObs_t

  INTERFACE AuxTypeObs_t
    PROCEDURE :: constructor
  END INTERFACE AuxTypeObs_t

CONTAINS

  FUNCTION constructor(mpObs) RESULT(this)
    TYPE(AuxTypeObs_t) :: this
    TYPE(MPObs_t), TARGET :: mpObs

    this%mpObs => mpObs

  END FUNCTION constructor

END MODULE AuxTypeObs_m
