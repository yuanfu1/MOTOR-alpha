!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2022/3/7, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the M2O Object which is tranform the variables from model space to
!! Observation space.
MODULE M2ODirect_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2OBase_m, ONLY: M2OBase_t

  TYPE, EXTENDS(M2OBase_t) :: M2ODirect_t
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(State_t), POINTER :: X
  CONTAINS
    FINAL :: destructor
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

    ! GENERIC :: OPERATOR(.TL.) => transFwdNonLinear_opr
    ! GENERIC :: OPERATOR(.ADJ.) => transAdjMultiply_opr
  END TYPE M2ODirect_t

  INTERFACE M2ODirect_t
    PROCEDURE :: constructor
  END INTERFACE M2ODirect_t

CONTAINS

  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(M2ODirect_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), TARGET, INTENT(IN) :: Y
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%Y => Y
    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(M2ODirect_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X, Y)
    CLASS(M2ODirect_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    ! DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
    DO i = 1, UBOUND(Y%ObsFields, 1) ! exclude zero
      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)

        IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
            Y%ObsFields(i)%values(k) = X%fields(j)%Get_Value(Y%ObsFields(i)%idx(k))
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE transFwdTanLinear(this, dX, dY, X)
    CLASS(M2ODirect_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(ObsSet_t), INTENT(INOUT) :: dY
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX, dY)
  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    CLASS(M2ODirect_t) :: this
    TYPE(ObsSet_t), INTENT(in) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
      IF (D%getObsIdx(dX%fields(j)%Get_Name()) .EQ. 0) THEN
        CYCLE
      ELSE
        dX%fields(j)%DATA = 0.0D0
      END IF
      DO i = LBOUND(D%ObsFields, 1), UBOUND(D%ObsFields, 1)
        IF ((TRIM(D%ObsFields(i)%Get_Name())) .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(D%ObsFields(i)%idx, 1), UBOUND(D%ObsFields(i)%idx, 1)
            CALL dX%fields(j)%Add_Value(D%ObsFields(i)%idx(k), D%ObsFields(i)%values(k))
          END DO
        END IF
      END DO
    END DO
    CALL dX%exHalo()
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(M2ODirect_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y
    INTEGER :: i, j, k

    Y = this%Y%zeroCopy()
    CALL this%transFwdNonLinear(X, Y)
  END FUNCTION transFwdNonLinear_opr

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    CLASS(M2ODirect_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(ObsSet_t) :: dY

    dY = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(M2ODirect_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = this%X%zeroCopy()
    CALL this%transAdjMultiply(D, dX)
  END FUNCTION transAdjMultiply_opr

END MODULE M2ODirect_m
