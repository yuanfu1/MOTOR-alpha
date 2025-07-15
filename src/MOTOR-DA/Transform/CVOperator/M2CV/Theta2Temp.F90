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

MODULE Theta2Temp_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE M2MBase_m, ONLY: M2MBase_t
  USE parameters_m, ONLY: k_d

  TYPE, EXTENDS(M2MBase_t) :: Theta2Temp_t
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
  END TYPE Theta2Temp_t

  INTERFACE Theta2Temp_t
    PROCEDURE :: constructor
  END INTERFACE Theta2Temp_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(Theta2Temp_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Theta2Temp_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: j

    IF (X%getVarIdx('theta') .EQ. 0) THEN
      PRINT *, "ERROR: No theta in X fields when converting theta to temp"
      STOP
    END IF

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('theta' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        IF (X%getVarIdx('pres') .NE. 0) THEN
          X%fields(X%getVarIdx('temp'))%DATA = X%fields(j)%DATA * &
                                               (X%fields(X%getVarIdx('pres'))%DATA / 1.0D05)**k_d
        ELSE
          PRINT *, "ERROR: No pres in X fields when converting theta to temp"
          STOP
        END IF
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: j

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('temp' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
        IF (X%getVarIdx('theta') .EQ. 0) THEN
          CALL X%addVar('theta')
        END IF

        IF (X%getVarIdx('pres') .NE. 0) THEN
          X%fields(X%getVarIdx('theta'))%DATA = X%fields(j)%DATA * &
                                                (1.0D05 / X%fields(X%getVarIdx('pres'))%DATA)**k_d
        ELSE
          PRINT *, "ERROR: No pres in X fields when converting temp to theta"
          STOP
        END IF
      END IF
    END DO
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    IF (X%getVarIdx('temp') .NE. 0) THEN
      X1 = X
      CALL this%transFwdNonLinear(X1)
    ELSE
      PRINT *, "ERROR: No temp in X fields when converting temp to theta"
      STOP
    END IF
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind) :: j

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
      IF ('theta' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('temp') .EQ. 0) THEN
          CALL dX%addVar('temp')
        END IF

        IF (X%getVarIdx('pres') .NE. 0) THEN
          dX%fields(dX%getVarIdx('temp'))%DATA = dX%fields(j)%DATA * &
                                                 (1.0D05 / X%fields(X%getVarIdx('pres'))%DATA)**k_d
        ELSE
          PRINT *, "ERROR: No pres in X fields when converting theta to temp in dX"
          STOP
        END IF

      END IF
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(Theta2Temp_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (dX%getVarIdx('theta') .NE. 0) THEN
      dX1 = dX
      CALL this%transAdjMultiply(dX1, X)
    ELSE
      PRINT *, "ERROR: No theta in dX fields when converting adjoint theta to temp in dX"
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

END MODULE Theta2Temp_m
