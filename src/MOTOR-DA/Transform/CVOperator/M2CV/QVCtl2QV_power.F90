!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2022/7/27, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the QVCtl2QV_power Object which is tranform the variables from control variable space to
!! model space.
! This module performs the following conversions:
!
! NL: x = e^[ln(p * x_ctl + 1)/p]
! dx/dx_ctl = e^[ln(p * x_ctl + 1)/p] / ( 1 + p * x_ctl)
! TL: dx= (dx/dx_ctl) * dx_ctl = e^[ln(p * x_ctl + 1)/p] / ( 1 + p * x_ctl) * d(x_ctl)
! AD: grad(x_ctl)_o = grad(x) * dx/d(x_ctl) =  grad(x) * e^[ln(p * x_ctl + 1)/p] / ( 1 + p * x_ctl)
! Backward: x_ctl = (x^p-1)/p
!
! grad(x_ctl)_b = dJ/d(x_ctl) = dJ/dx * dx/d(x_ctl)
! grad(x_ctl)_o = dJ/dH * dH/d(x_ctl) = dJ/dH * dH/dx * dx/d(x_ctl)
! with dJ/dH = 1, then, grad(x_ctl)_o = grad(x) * dx/d(x_ctl)
MODULE QVCtl2QV_power_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: QVCtl2QV_power_t
    TYPE(State_t), POINTER :: X
  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: transBackward

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

    ! GENERIC :: OPERATOR(.yaml.) => transFwdNonLinear_opr
    ! GENERIC :: OPERATOR(.ADJ.) => transAdjMultiply_opr
    FINAL :: destructor
  END TYPE QVCtl2QV_power_t

  INTERFACE QVCtl2QV_power_t
    PROCEDURE :: constructor
  END INTERFACE QVCtl2QV_power_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(QVCtl2QV_power_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(QVCtl2QV_power_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: p = 0.1D0
    REAL(r_kind) :: tmp1, tmp2

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)

      IF ('qvapor' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('qvapor_ctl') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('qvapor_ctl')
        END IF

        X%fields(j)%DATA = (X%fields(j)%DATA**p - 1.0D0) / p
        PRINT *, 'check qvapor_ctl ', MAXVAL(X%fields(j)%DATA), MINVAL(X%fields(j)%DATA)

      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k
    REAL(r_kind) :: p = 0.1D0

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)

      IF ('qvapor_ctl' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('qvapor') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('qvapor')
        END IF

        X%fields(j)%DATA = EXP(LOG(p * X%fields(j)%DATA + 1.0D0) / p)
        PRINT *, 'check NL: qvapor ', MAXVAL(X%fields(j)%DATA), MINVAL(X%fields(j)%DATA)
      END IF
    END DO

  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X%zeroCopy()
    IF (X%getVarIdx('qvapor_ctl') .NE. 0) THEN
      X1%fields(X1%getVarIdx('qvapor_ctl')) = X%fields(X%getVarIdx('qvapor_ctl'))
      CALL this%transFwdNonLinear(X1)
    END IF
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: p = 0.1D0

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)

      IF ('qvapor' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('qvapor_ctl') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('qvapor_ctl')
        END IF

        dX%fields(j)%DATA = dX%fields(j)%DATA * EXP(LOG(p * X%fields(X%getVarIdx('qvapor_ctl'))%DATA + 1.0D0) / p) &
                            / (1.0D0 + p * X%fields(X%getVarIdx('qvapor_ctl'))%DATA)
        PRINT *, 'check AD grad(x_ctl): ', MAXVAL(dX%fields(j)%DATA), MINVAL(dX%fields(j)%DATA)

      END IF
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX%zeroCopy()
    IF (dX%getVarIdx('qvapor') .NE. 0) THEN
      dX1%fields(dX1%getVarIdx('qvapor')) = dX%fields(dX%getVarIdx('qvapor'))
      CALL this%transAdjMultiply(dX1)
    END IF
  END FUNCTION

  FUNCTION transElemForward(valIn, sv) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, sv

    valOut = valIn * sv
  END FUNCTION

  FUNCTION transElemAdjointMultiply(valIn, sv) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, sv

    valOut = valIn * sv
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER :: i, j, k
    REAL(r_kind) :: p = 0.1D0

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)

      IF ('qvapor_ctl' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('qvapor') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('qvapor')
        END IF

        dX%fields(j)%DATA = dX%fields(j)%DATA * EXP(LOG(p * X%fields(X%getVarIdx('qvapor_ctl'))%DATA + 1.0D0) / p) &
                            / (1.0D0 + p * X%fields(X%getVarIdx('qvapor_ctl'))%DATA)
        PRINT *, 'check TL: qvapor ', MAXVAL(dX%fields(j)%DATA), MINVAL(dX%fields(j)%DATA)

      END IF
    END DO

  END SUBROUTINE transFwdTanLinear

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(QVCtl2QV_power_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = dX%zeroCopy()
    IF (dX%getVarIdx('qvapor_ctl') .NE. 0) THEN
      dX1%fields(dX1%getVarIdx('qvapor_ctl')) = dX%fields(dX%getVarIdx('qvapor_ctl'))
      CALL this%transFwdTanLinear(dX1, X)
    END IF

  END FUNCTION

END MODULE QVCtl2QV_power_m
