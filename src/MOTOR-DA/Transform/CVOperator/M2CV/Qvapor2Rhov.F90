!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2022/2/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief

MODULE Qvapor2Rhov_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE M2MBase_m, ONLY: M2MBase_t
  USE parameters_m, ONLY: k_d

  TYPE, EXTENDS(M2MBase_t) :: Qvapor2Rhov_t
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
  END TYPE Qvapor2Rhov_t

  INTERFACE Qvapor2Rhov_t
    PROCEDURE :: constructor
  END INTERFACE Qvapor2Rhov_t

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(Qvapor2Rhov_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Qvapor2Rhov_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: j

    IF (X%getVarIdx('qvapor') .EQ. 0) THEN
      PRINT *, "ERROR: No qvapor in X fields when converting qvapor to rhov"
      STOP
    END IF

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('qvapor' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('rhov') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('rhov')
        END IF

        IF (X%getVarIdx('rho') .NE. 0) THEN
          X%fields(j)%DATA = X%fields(j)%DATA * X%fields(j)%DATA
        ELSE
          PRINT *, "ERROR: No rho in X fields when converting qvapor to rhov"
          STOP
        END IF
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: j

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF ('rhov' .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

        IF (X%getVarIdx('qvapor') .EQ. 0) THEN
          CALL X%fields(j)%Set_Name('qvapor')
        END IF

        IF (X%getVarIdx('rho') .NE. 0) THEN
          X%fields(j)%DATA = X%fields(j)%DATA / X%fields(X%getVarIdx('rho'))%DATA
        ELSE
          PRINT *, "ERROR: No rho in X fields when converting rhov to qvapor"
          STOP
        END IF
      END IF
    END DO
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    IF (X%getVarIdx('rhov') .NE. 0) THEN
      X1 = X
      CALL this%transFwdNonLinear(X1)
    ELSE
      PRINT *, "ERROR: No rhov in X fields when converting rhov to qvapor"
      STOP
    END IF
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER :: j

    ! q_v = rho_v / rho
    ! d(q_v) = d(rho_v)/rho - d(rho)* rho_v / rho^2
    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)

      IF ('rhov' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('qvapor') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('qvapor')
        END IF

        IF (X%getVarIdx('rho') .NE. 0) THEN
          dX%fields(j)%DATA = dX%fields(j)%DATA/X%fields(X%getVarIdx('rho'))%DATA !&
          ! - dX%fields(dX%getVarIdx('rho'))%data * X%fields(j)%data / (X%fields(X%getVarIdx('rho'))%data)**2
        ELSE
          PRINT *, "ERROR: No rhov in dX fields when converting drhov to dqvapor"
          STOP
        END IF
      END IF
    END DO

  END SUBROUTINE transFwdTanLinear

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    IF (dX%getVarIdx('rhov') .NE. 0) THEN
      dX1 = dX
      CALL this%transFwdTanLinear(dX1, X)
    ELSE
      PRINT *, "ERROR: No rhov in X fields when converting rhov to qvapor"
      STOP
    END IF

  END FUNCTION

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind) :: j

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
      IF ('qvapor' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

        IF (dX%getVarIdx('rhov') .EQ. 0) THEN
          CALL dX%fields(j)%Set_Name('rhov')
        END IF

        IF (X%getVarIdx('rho') .NE. 0) THEN
          dX%fields(j)%DATA = dX%fields(j)%DATA / &
                              X%fields(X%getVarIdx('rho'))%DATA
        ELSE
          PRINT *, "ERROR: No rho in X fields when converting qvapor to rhov in dX"
          STOP
        END IF

      END IF
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(Qvapor2Rhov_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (dX%getVarIdx('qvapor') .NE. 0) THEN
      dX1 = dX
      CALL this%transAdjMultiply(dX1, X)
    ELSE
      PRINT *, "ERROR: No qvapor in dX fields when converting adjoint qvapor to rhov in dX"
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

END MODULE Qvapor2Rhov_m
