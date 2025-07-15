!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/4/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! Under this subdirectory, users of MOTOR can define their own control variables and provide their conversion functions from the control to model states. We start all modules with letter c in the module names List of choices:
!! 1. Vorticity and Divergence control: $\zeta$, $\delta$, $\lnP$, $\rh$, $temp$, the module is named as cVortDive.F90, This abstract module defines the structure of all control variable choices.
MODULE C2MBase_m

  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  TYPE, ABSTRACT :: C2MBase_t
    TYPE(State_t) :: analysis
    CLASS(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE(forward), DEFERRED :: forwardOpr_s
    PROCEDURE(bckward), DEFERRED :: bckwardOpr_s
    PROCEDURE(tangent), DEFERRED :: tangentOpr_s
    PROCEDURE(adjoint), DEFERRED :: adjointOpr_s
  END TYPE C2MBase_t

  ABSTRACT INTERFACE
    FUNCTION forward(this, X) RESULT(Y)
      IMPORT :: i_kind, C2MBase_t, State_t
      CLASS(C2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      TYPE(State_t) :: Y
    END FUNCTION forward
    FUNCTION bckward(this, X) RESULT(Y)
      IMPORT :: i_kind, C2MBase_t, State_t
      CLASS(C2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      TYPE(State_t) :: Y
    END FUNCTION bckward

    FUNCTION tangent(this, X) RESULT(Y)
      IMPORT :: i_kind, C2MBase_t, State_t
      CLASS(C2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      TYPE(State_t) :: Y
    END FUNCTION tangent

    FUNCTION adjoint(this, X, dt) RESULT(Y)
      IMPORT :: i_kind, r_kind, C2MBase_t, State_t
      CLASS(C2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      REAL(r_kind), INTENT(IN), OPTIONAL :: dt
      TYPE(State_t) :: Y
    END FUNCTION adjoint
    FUNCTION bckAdjoint(this, X) RESULT(Y)
      IMPORT :: i_kind, C2MBase_t, State_t
      CLASS(C2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      TYPE(State_t) :: Y
    END FUNCTION bckAdjoint
  END INTERFACE

CONTAINS
  SUBROUTINE initialize_s(this, sg, X, dyc)

    IMPLICIT NONE

    CLASS(C2MBase_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), INTENT(IN) :: X
    CLASS(dyCoresBase_t), TARGET, OPTIONAL :: dyc

    ! The base module provides a default setting for control variables over all space and time:
    ! (see docs: Object-oriented Design.md under etc/doc of MOTOR)

    this%analysis = X
    this%sg => sg

  END SUBROUTINE initialize_s

  SUBROUTINE forwardOpr_s(this, glvl, X, Y)
    CLASS(C2MBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Y

    Y = X ! Identical forward operator
  END SUBROUTINE forwardOpr_s

  SUBROUTINE bckwardOpr_s(this, glvl, X, Y)
    CLASS(C2MBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Y

    Y = X ! Identical forward operator
  END SUBROUTINE bckwardOpr_s

  SUBROUTINE tangentOpr_s(this, glvl, X, Y)
    CLASS(C2MBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Y

    Y = X ! Identical forward operator
  END SUBROUTINE tangentOpr_s

  SUBROUTINE adjointOpr_s(this, glvl, X, Y)
    CLASS(C2MBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Y

    Y = X ! Identical forward operator
  END SUBROUTINE adjointOpr_s

END MODULE C2MBase_m
