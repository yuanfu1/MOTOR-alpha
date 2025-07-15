!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE ObsAttr_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double

!> @brief
!! Base type of observation Attrmeters
! @see
! @note
! @warning
! @attention
  TYPE ObsAttr_t

  CONTAINS

  END TYPE ObsAttr_t

!!> @brief
!! A data struct containing the model_states
! @see
! @note
! @warning
! @attention
  TYPE, EXTENDS(ObsAttr_t) :: ObsAttrSAT_t
    REAL(r_kind), ALLOCATABLE :: zenangle(:)
    REAL(r_kind), ALLOCATABLE :: azangle(:)
    REAL(r_kind), ALLOCATABLE :: sunzenangle(:)
    REAL(r_kind), ALLOCATABLE :: sunazangle(:)
    REAL(r_kind), ALLOCATABLE :: latitude(:)
    REAL(r_kind), ALLOCATABLE :: longitude(:)
  CONTAINS
    FINAL :: des_ObsAttrSAT
  END TYPE ObsAttrSAT_t

  ! Constructor
  INTERFACE ObsAttr_t
    PROCEDURE :: cons_ObsAttr
  END INTERFACE ObsAttr_t

  INTERFACE ObsAttrSAT_t
    PROCEDURE :: cons_ObsAttrSAT
  END INTERFACE ObsAttrSAT_t

CONTAINS

  FUNCTION cons_ObsAttr(configFile) RESULT(this)
    TYPE(ObsAttr_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

  END FUNCTION

  FUNCTION cons_ObsAttrSAT(configFile) RESULT(this)
    TYPE(ObsAttrSAT_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    this%ObsAttr_t = ObsAttr_t(configFile)
  END FUNCTION

  IMPURE ELEMENTAL SUBROUTINE des_ObsAttrSAT(this)
    IMPLICIT NONE
    TYPE(ObsAttrSAT_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%zenangle)) DEALLOCATE (this%zenangle)
    IF (ALLOCATED(this%azangle)) DEALLOCATE (this%azangle)
    IF (ALLOCATED(this%sunzenangle)) DEALLOCATE (this%sunzenangle)
    IF (ALLOCATED(this%sunazangle)) DEALLOCATE (this%sunazangle)
    IF (ALLOCATED(this%latitude)) DEALLOCATE (this%latitude)
    IF (ALLOCATED(this%longitude)) DEALLOCATE (this%longitude)
  END SUBROUTINE des_ObsAttrSAT

END MODULE ObsAttr_m
