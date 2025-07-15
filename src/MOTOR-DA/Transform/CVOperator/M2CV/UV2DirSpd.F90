!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/8/2, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
! This module converts u and v control variables into COS(wdir) and SIN(wdir) and wspd for
! use of wind direction and speed as observations.
MODULE UV2DirSpd_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE M2OBase_m, ONLY: M2OBase_t
  USE ObsSet_m, ONLY: ObsSet_t

  TYPE, EXTENDS(M2OBase_t) :: UV2DirSpd_t
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(State_t), POINTER :: X
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

    PROCEDURE, PUBLIC :: transBackward
    PROCEDURE, PUBLIC, NOPASS :: transElemForward
    PROCEDURE, PUBLIC, NOPASS :: transElemAdjointMultiply

    FINAL :: destructor
  END TYPE UV2DirSpd_t

CONTAINS
  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(UV2DirSpd_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), TARGET, INTENT(IN) :: Y
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(UV2DirSpd_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(ObsSet_t), INTENT(in) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = X%zeroCopy()
    CALL this%transAdjMultiply(D, dX, X)
  END FUNCTION

  SUBROUTINE transFwdNonLinear(this, X, Y)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
  END SUBROUTINE transFwdNonLinear

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    IF (X%getVarIdx('rho') .NE. 0) THEN
      X1 = X
      CALL this%transFwdNonLinear(X1)
    ELSE
      PRINT *, "ERROR: No rho in X fields when converting rho to pres"
      STOP
    END IF
  END FUNCTION transFwdNonLinear_opr

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION transFwdTanLinear_opr

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE transFwdTanLinear

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
  END SUBROUTINE transBackward

  FUNCTION transElemAdjointMultiply(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

    valOut = valIn * s1
  END FUNCTION transElemAdjointMultiply

  FUNCTION transElemForward(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

    valOut = valIn * s1
  END FUNCTION transElemForward
END MODULE UV2DirSpd_m
