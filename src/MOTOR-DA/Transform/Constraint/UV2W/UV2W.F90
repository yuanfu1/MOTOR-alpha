!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/20, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the UV2W Object which is tranform the variables from control variable space to
!! model space.
MODULE UV2W_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: UV2W_t
    ! TYPE(State_t), pointer :: X

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    PROCEDURE, PUBLIC, NOPASS :: transElemForward
    PROCEDURE, PUBLIC, NOPASS :: transElemAdjointMultiply

    FINAL :: destructor
  END TYPE UV2W_t

  ! INTERFACE UV2W_t
  !   PROCEDURE :: constructor
  ! END INTERFACE UV2W_t

CONTAINS

  ! FUNCTION constructor(configFile, X) RESULT(this)
  !   IMPLICIT NONE
  !   TYPE(UV2W_Poisson_t) :: this
  !   CHARACTER(LEN=1024), INTENT(IN) :: configFile
  !   TYPE(State_t), TARGET, INTENT(IN) :: X

  !   this%X => X
  !   IF (X%sg%vLevel < 3) return
  !   CALL this%calParameters
  ! END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(UV2W_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

  END SUBROUTINE transAdjMultiply

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

  END SUBROUTINE

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

  END FUNCTION

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1

  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(UV2W_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t) :: dX1

  END FUNCTION

  FUNCTION transElemForward(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION

  FUNCTION transElemAdjointMultiply(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION
END MODULE UV2W_m
