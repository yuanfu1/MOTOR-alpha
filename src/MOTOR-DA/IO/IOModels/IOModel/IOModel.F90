!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/3/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE IOModel_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE State_m, ONLY: State_t

  TYPE IOModel_t
    TYPE(geometry_t), POINTER :: geometry
    CHARACTER(LEN=1024) :: m_configFile

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initializeIOModel
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm_Ens
    PROCEDURE, PUBLIC, PASS(this) :: m_write_Xm_into_bcg
    FINAL :: destructor
  END TYPE IOModel_t

CONTAINS
  SUBROUTINE initializeIOModel(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOModel_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    this%geometry => geometry
    this%m_configFile = configFile
  END SUBROUTINE

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    IMPLICIT NONE
    CLASS(IOModel_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
  END SUBROUTINE m_read_bcg_into_Xm

  SUBROUTINE m_read_bcg_into_Xm_Ens(this, Xm, sg, ensIndx)
    IMPLICIT NONE
    CLASS(IOModel_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: ensIndx
  END SUBROUTINE m_read_bcg_into_Xm_Ens

  SUBROUTINE m_write_Xm_into_bcg(this, Xm, sg, Xb)
    IMPLICIT NONE
    CLASS(IOModel_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(State_t), INTENT(IN) :: Xb
  END SUBROUTINE m_write_Xm_into_bcg

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(IOModel_t), INTENT(INOUT) :: this
  END SUBROUTINE destructor

END MODULE IOModel_m
