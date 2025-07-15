!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2022/7/05, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the RhoVCtl2RhoV_sv Object which is tranform the variables from control variable space to
!! model space.
! This module performs the following conversions:
!
! NL: x = x_ctl * p
! dx/dx_ctl = p
! TL: dx = (dx/dx_ctl)*dx_ctl = p * d(x_ctl)
! AD: grad(x_ctl)_o = grad(x) * dx/d(x_ctl) = grad(x) * p
! Backward: x_ctl = x/p
!
! grad(x_ctl)_b = dJ/d(x_ctl) = dJ/dx * dx/d(x_ctl)
! grad(x_ctl)_o = dJ/dH * dH/d(x_ctl) = dJ/dH * dH/dx * dx/d(x_ctl)
! with dJ/dH = 1, then, grad(x_ctl)_o = grad(x) * dx/d(x_ctl)
MODULE RhoVCtl2RhoV_sv_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: RhoVCtl2RhoV_sv_t
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
  END TYPE RhoVCtl2RhoV_sv_t

  INTERFACE RhoVCtl2RhoV_sv_t
    PROCEDURE :: constructor
  END INTERFACE RhoVCtl2RhoV_sv_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(RhoVCtl2RhoV_sv_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(RhoVCtl2RhoV_sv_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i, j, k

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)

      IF ('rhov' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('rhov_ctl') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('rhov_ctl')
        END IF

        DO k = 1, X%sg%tSlots
          X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) / X%sg%SVapor
        END DO
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)

    IMPLICIT NONE
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)

      IF ('rhov_ctl' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('rhov') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('rhov')
        END IF

        DO k = 1, X%sg%tSlots
          X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) * X%sg%SVapor
        END DO

      END IF

    END DO
    ! IF (MAXVAL(X%sg%SVapor) .GT. 1.0)  &
    !     PRINT *, 'WARNING: scaling profile exceed max check ', X%sg%mpddInfo_sg%myrank, MAXVAL(X%sg%SVapor), MINVAL(X%sg%SVapor)
    ! PRINT *, 'check rhov: ', MAXVAL(X%fields(X%getVarIdx('rhov'))%data), MINVAL(X%fields(X%getVarIdx('rhov'))%data)
    ! PRINT *, 'check rhov_ctl profile ', X%fields(X%getVarIdx('rhov_ctl'))%data(:, 1, 1)
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X%zeroCopy()
    IF (X%getVarIdx('rhov_ctl') .NE. 0) THEN
      X1%fields(X1%getVarIdx('rhov_ctl')) = X%fields(X%getVarIdx('rhov_ctl'))
      CALL this%transFwdNonLinear(X1)
    END IF
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    INTEGER(i_kind) :: i, j, k

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)

      IF ('rhov' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('rhov_ctl') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('rhov_ctl')
        END IF

        DO k = 1, dX%sg%tSlots
          dX%fields(j)%DATA(:, :, k) = dX%fields(j)%DATA(:, :, k) * dX%sg%SVapor
        END DO

      END IF

    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX%zeroCopy()
    IF (dX%getVarIdx('rhov') .NE. 0) THEN
      dX1%fields(dX1%getVarIdx('rhov')) = dX%fields(dX%getVarIdx('rhov'))
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
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(RhoVCtl2RhoV_sv_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

END MODULE RhoVCtl2RhoV_sv_m
