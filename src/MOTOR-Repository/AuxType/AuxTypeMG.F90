!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.AuxType
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/25, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This type enclosed the essential parameters in constructing types.
MODULE AuxTypeMG_m
  USE MultiGrid_m, ONLY: MultiGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t

  TYPE AuxTypeMG_t
    TYPE(mpddGlob_t), POINTER :: mpddGlob
    TYPE(MultiGrid_t), POINTER :: mg
    CHARACTER(LEN=1024) :: m_configFile
  CONTAINS
    FINAL :: destructor
  END TYPE AuxTypeMG_t

  INTERFACE AuxTypeMG_t
    PROCEDURE :: constructor
  END INTERFACE AuxTypeMG_t

CONTAINS

  FUNCTION constructor(configFile, mg) RESULT(this)
    IMPLICIT NONE
    TYPE(AuxTypeMG_t) :: this
    TYPE(MultiGrid_t), TARGET, INTENT(IN) :: mg
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    this%mpddGlob => this%mg%mpddGlob
    this%mg => mg
    this%m_configFile = configFile
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(AuxTypeMG_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor
END MODULE AuxTypeMG_m
