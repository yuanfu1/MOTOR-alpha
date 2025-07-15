!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/2/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the QrCtl2Qr Object which is tranform the variables from control variable space to
!! model space.
MODULE QrCtl2Qr_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: QrCtl2Qr_t
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
  END TYPE QrCtl2Qr_t

  INTERFACE QrCtl2Qr_t
    PROCEDURE :: constructor
  END INTERFACE QrCtl2Qr_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(QrCtl2Qr_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(QrCtl2Qr_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i, j, k

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('qrain' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        CALL X%fields(j)%Set_Name('qrain_ctl')
        FORALL (k=1:X%sg%tSlots)
          X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) / X%sg%s1
        END FORALL
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('qrain_ctl' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        CALL X%fields(j)%Set_Name('qrain')
        FORALL (k=1:X%sg%tSlots)
          X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) * X%sg%s1
        END FORALL
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X%zeroCopy()
    IF (X%getVarIdx('qrain_ctl') .NE. 0) THEN
      X1%fields(X1%getVarIdx('qrain_ctl')) = X%fields(X%getVarIdx('qrain_ctl'))
      CALL this%transFwdNonLinear(X1)
    END IF
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    INTEGER(i_kind) :: i, j, k

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
      IF ('qrain' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN
        CALL dX%fields(j)%Set_Name('qrain_ctl')
        FORALL (k=1:dX%sg%tSlots)
          dX%fields(j)%DATA(:, :, k) = dX%fields(j)%DATA(:, :, k) * dX%sg%s1
        END FORALL
      END IF
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(QrCtl2Qr_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX%zeroCopy()
    IF (dX%getVarIdx('qrain') .NE. 0) THEN
      dX1%fields(dX1%getVarIdx('qrain')) = dX%fields(dX%getVarIdx('qrain'))
      CALL this%transAdjMultiply(dX1)
    END IF
  END FUNCTION

  FUNCTION transElemForward(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

    valOut = valIn * s1
  END FUNCTION

  FUNCTION transElemAdjointMultiply(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

    valOut = valIn * s1
  END FUNCTION
END MODULE QrCtl2Qr_m
