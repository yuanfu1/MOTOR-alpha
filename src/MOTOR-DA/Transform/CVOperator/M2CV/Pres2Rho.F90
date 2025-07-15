!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (zilong.qin@gmail.com), 2022/2/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief

MODULE Pres2Rho_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE M2MBase_m, ONLY: M2MBase_t
  USE parameters_m, ONLY: r_d

  TYPE, EXTENDS(M2MBase_t) :: Pres2Rho_t
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
  END TYPE Pres2Rho_t

  INTERFACE Pres2Rho_t
    PROCEDURE :: constructor
  END INTERFACE Pres2Rho_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(Pres2Rho_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Pres2Rho_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: j, sizeDim1, sizeDim2, sizeDim3
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: qRatio

    IF (X%getVarIdx('pres') .EQ. 0) THEN
      PRINT *, "ERROR: No pres in X fields when converting pres to rho"
      STOP
    END IF

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('pres' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        sizeDim1 = SIZE(X%fields(j)%DATA, dim=1)
        sizeDim2 = SIZE(X%fields(j)%DATA, dim=2)
        sizeDim3 = SIZE(X%fields(j)%DATA, dim=3)
        ALLOCATE (qRatio(sizeDim1, sizeDim2, sizeDim3))
        qRatio = 1.0D0
        IF (X%getVarIdx('qvapor') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qvapor'))%DATA
        ELSE
          PRINT *, "ERROR: No qvapor in X fields when converting pres to rho"
          STOP
        END IF
        IF (X%getVarIdx('qrain') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qrain'))%DATA
        END IF
        IF (X%getVarIdx('qice') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('qcloud') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qcloud'))%DATA
        END IF
        IF (X%getVarIdx('qsnow') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qsnow'))%DATA
        END IF
        IF (X%getVarIdx('qgraup') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('temp') .NE. 0) THEN

          IF (X%getVarIdx('rho') .EQ. 0) THEN
            CALL X%fields(j)%Set_Name('rho')
          END IF

          X%fields(j)%DATA = X%fields(j)%DATA / r_d / &
                             X%fields(X%getVarIdx('temp'))%DATA / qRatio
        ELSE
          PRINT *, "ERROR: No temp in X fields when converting pres to rho"
          STOP
        END IF
      END IF
    END DO
    DEALLOCATE (qRatio)
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: j, sizeDim1, sizeDim2, sizeDim3
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: qRatio

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('rho' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        IF (X%getVarIdx('pres') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('pres')
        END IF
        sizeDim1 = SIZE(X%fields(j)%DATA, dim=1)
        sizeDim2 = SIZE(X%fields(j)%DATA, dim=2)
        sizeDim3 = SIZE(X%fields(j)%DATA, dim=3)
        ALLOCATE (qRatio(sizeDim1, sizeDim2, sizeDim3))
        qRatio = 1.0D0
        IF (X%getVarIdx('qvapor') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qvapor'))%DATA
        ELSE
          PRINT *, "ERROR: No qvapor in X fields when converting pres to rho"
          STOP
        END IF
        IF (X%getVarIdx('qrain') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qrain'))%DATA
        END IF
        IF (X%getVarIdx('qice') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('qcloud') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qcloud'))%DATA
        END IF
        IF (X%getVarIdx('qsnow') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qsnow'))%DATA
        END IF
        IF (X%getVarIdx('qgraup') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('temp') .NE. 0) THEN
          X%fields(j)%DATA = X%fields(j)%DATA * r_d * &
                             X%fields(X%getVarIdx('temp'))%DATA * qRatio
        ELSE
          PRINT *, "ERROR: No temp in X fields when converting pres to rho"
          STOP
        END IF
      END IF
    END DO
    IF (ALLOCATED(qRatio)) DEALLOCATE (qRatio)
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    IF (X%getVarIdx('rho') .NE. 0) THEN
      X1 = X
      CALL this%transFwdNonLinear(X1)
    ELSE
      PRINT *, "ERROR: No rho in X fields when converting rho to pres"
      STOP
    END IF
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind) :: j, sizeDim1, sizeDim2, sizeDim3
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: qRatio

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
      IF ('pres' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('rho') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('rho')
        END IF

        sizeDim1 = SIZE(dX%fields(j)%DATA, dim=1)
        sizeDim2 = SIZE(dX%fields(j)%DATA, dim=2)
        sizeDim3 = SIZE(dX%fields(j)%DATA, dim=3)
        ALLOCATE (qRatio(sizeDim1, sizeDim2, sizeDim3))
        qRatio = 1.0D0
        IF (X%getVarIdx('qvapor') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qvapor'))%DATA
        ELSE
          PRINT *, "ERROR: No qvapor in X fields when converting pres to rho"
          STOP
        END IF
        IF (X%getVarIdx('qrain') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qrain'))%DATA
        END IF
        IF (X%getVarIdx('qice') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('qcloud') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qcloud'))%DATA
        END IF
        IF (X%getVarIdx('qsnow') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qsnow'))%DATA
        END IF
        IF (X%getVarIdx('qgraup') .NE. 0) THEN
          qRatio = qRatio + X%fields(X%getVarIdx('qice'))%DATA
        END IF
        IF (X%getVarIdx('temp') .NE. 0) THEN
          dX%fields(j)%DATA = dX%fields(j)%DATA * r_d * &
                              X%fields(X%getVarIdx('temp'))%DATA * qRatio
        ELSE
          PRINT *, "ERROR: No temp in X fields when converting pres to rho in dX"
          STOP
        END IF

      END IF
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(Pres2Rho_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (dX%getVarIdx('pres') .NE. 0) THEN
      dX1 = dX
      CALL this%transAdjMultiply(dX1)
    ELSE
      PRINT *, "ERROR: No pres in dX fields when converting adjoint pres to rho in dX"
      STOP
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

END MODULE Pres2Rho_m
