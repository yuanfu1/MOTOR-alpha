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
MODULE LnP2Pres_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: LnP2Pres_t
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

    PROCEDURE, PUBLIC :: transBackward
    FINAL :: destructor
  END TYPE LnP2Pres_t

  ! INTERFACE LnP2Pres_t
  !   PROCEDURE :: constructor
  ! END INTERFACE LnP2Pres_t

CONTAINS

  ! FUNCTION constructor(configFile, X) RESULT(this)
  !   IMPLICIT NONE
  !   TYPE(LnP2Pres_t) :: this
  !   CHARACTER(LEN=1024), INTENT(IN) :: configFile
  !   TYPE(State_t), TARGET, INTENT(IN) :: X

  !   this%X => X
  ! END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(LnP2Pres_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    IF (X%getVarIdx('pres') .EQ. 0) THEN
      PRINT *, "ERROR: No pres in X fields when converting pres to lnp"
      STOP
    ELSE
      ASSOCIATE (pres => X%fields(X%getVarIdx('pres')))
        CALL pres%Set_Name('lnp')
        pres%DATA = LOG(pres%DATA)
      END ASSOCIATE
    END IF
  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    IF (X%getVarIdx('lnp') .EQ. 0) THEN
      PRINT *, "ERROR: No lnp in X fields when converting lnp to pres"
      STOP
    ELSE
      ASSOCIATE (lnp => X%fields(X%getVarIdx('lnp')))
        CALL lnp%Set_Name('pres')
        lnp%DATA = EXP(lnp%DATA)
      END ASSOCIATE
    END IF
  END SUBROUTINE

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (dX%getVarIdx('lnp') .EQ. 0) THEN
      PRINT *, "ERROR: No lnp in X fields when converting lnp to pres"
      STOP
    ELSE
      ASSOCIATE (lnp => dX%fields(dX%getVarIdx('lnp')))
        CALL lnp%Set_Name('pres')
        lnp%DATA = EXP(lnp%DATA) * lnp%DATA
      END ASSOCIATE
    END IF
  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (dX%getVarIdx('lnp') .EQ. 0) THEN
      PRINT *, "ERROR: No lnp in X fields when converting lnp to pres"
      STOP
    ELSE
      ASSOCIATE (lnp => dX%fields(dX%getVarIdx('lnp')))
        CALL lnp%Set_Name('pres')
        lnp%DATA = EXP(lnp%DATA) * lnp%DATA
      END ASSOCIATE
    END IF
  END SUBROUTINE transAdjMultiply

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1

    IF (X%getVarIdx('lnp') .NE. 0) THEN
      X1 = X
      CALL this%transFwdNonLinear(X1)
    ELSE
      PRINT *, "ERROR: No lnp in X fields when converting lnp to pres"
      STOP
    END IF

  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(LnP2Pres_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    IF (X%getVarIdx('lnp') .NE. 0) THEN
      dX1 = dX
      CALL this%transFwdNonLinear(dX1)
    ELSE
      PRINT *, "ERROR: No lnp in X fields when converting lnp to pres"
      STOP
    END IF
  END FUNCTION

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(LnP2Pres_t) :: this
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

END MODULE LnP2Pres_m
