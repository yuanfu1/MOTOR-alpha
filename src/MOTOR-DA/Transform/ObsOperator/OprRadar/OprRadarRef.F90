!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE OprRadarRef_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2OBase_m, ONLY: M2OBase_t

  TYPE, EXTENDS(M2OBase_t) :: OprRadarRef_t
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(State_t), POINTER :: X
    REAL(r_kind) :: a1, a2

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

    PROCEDURE, PUBLIC :: transElemForward
    PROCEDURE, PUBLIC :: transElemAdjointMultiply
    PROCEDURE, PRIVATE :: calRefFactors
    PROCEDURE, PUBLIC :: calRhoRFromRefRho

    FINAL :: destructor
  END TYPE OprRadarRef_t

  INTERFACE OprRadarRef_t
    PROCEDURE :: constructor
  END INTERFACE OprRadarRef_t

CONTAINS

  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(OprRadarRef_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(ObsSet_t), OPTIONAL, TARGET, INTENT(IN) :: Y

    this%X => X
    IF (PRESENT(Y)) this%Y => Y
    CALL this%calRefFactors(43.1D0, 17.5D0)
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(OprRadarRef_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE calRefFactors(this, c1, c2)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    REAL(r_kind) :: c1, c2

    this%a1 = 10**(c1 / 10)
    this%a2 = c2 / 10.0

    this%a1 = c1
    this%a2 = c2
  END SUBROUTINE

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i, j, k

  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X, Y)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y

    INTEGER(i_kind) :: i, j, k

    IF ((X%getVarIdx(TRIM('rhor')) .EQ. 0) .OR. &
        (Y%getObsIdx(TRIM('ref')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in OprRadarRef_t, STOP! '
      STOP
    END IF

    ASSOCIATE (rhor => X%Fields(X%getVarIdx(TRIM('rhor'))), &
               ref => Y%ObsFields(Y%getObsIdx(TRIM('ref'))))

      DO k = LBOUND(ref%idx, 1), UBOUND(ref%idx, 1)
        ref%values(k) = this%transElemForward(rhor%Get_Value(ref%idx(k)))
      END DO

    END ASSOCIATE
  END SUBROUTINE

  SUBROUTINE transFwdTanLinear(this, dX, dY, X)
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(ObsSet_t), INTENT(INOUT) :: dY
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX, dY)
  END SUBROUTINE

  SUBROUTINE calRhoRFromRefRho(this, X, Y)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: val

    IF ((X%getVarIdx(TRIM('rhor')) .EQ. 0)) THEN
      PRINT *, 'Error use the OprRadarRef_t in calRhoRFromRefRho, STOP! '
      STOP
    END IF

    ASSOCIATE (rhor => X%Fields(X%getVarIdx('rhor')), &
               ref => Y%ObsFields(Y%getObsIdx('ref')))

      DO k = LBOUND(ref%idx, 1), UBOUND(ref%idx, 1)
        ! CALL rhor%Set_Value(ref%idx(k), (ref%values(k)/this%a1)**(1/this%a2))
        CALL rhor%Set_Value(ref%idx(k), 10**((ref%values(k) - this%a1) / this%a2))
      END DO
    END ASSOCIATE

  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y
    INTEGER(i_kind) :: i, j, k

    Y = this%Y%zeroCopy()
    CALL this%transFwdNonLinear(X, Y)
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    CLASS(OprRadarRef_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(ObsSet_t) :: dY

    dY = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    INTEGER(i_kind) :: i, j, k

    IF ((X%getVarIdx(TRIM('rhor')) .EQ. 0) .OR. &
        (D%getObsIdx(TRIM('ref')) .EQ. 0)) THEN
      PRINT *, 'Error use the OprRadarRef_t in transAdjMultiply, STOP! '
      STOP
    END IF

    IF (.NOT. PRESENT(X)) THEN
      PRINT *, 'Error use the OprRadarRef_t in transAdjMultiply, STOP! '
      STOP
    END IF

    ASSOCIATE (rhor => X%Fields(X%getVarIdx(TRIM('rhor'))), &
               ref => D%ObsFields(D%getObsIdx(TRIM('ref'))))

      DO k = LBOUND(ref%idx, 1), UBOUND(ref%idx, 1)
        BLOCK
          REAL(r_kind) :: out
          out = this%transElemAdjointMultiply(ref%values(k), rhor%Get_Value(ref%idx(k)))
          CALL dX%Fields(dX%getVarIdx(TRIM('rhor')))%Set_Value(ref%idx(k), out)
        END BLOCK
      END DO

    END ASSOCIATE
    CALL dX%exHalo()
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(OprRadarRef_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = X%zeroCopy()
    CALL this%transAdjMultiply(D, dX, X)
  END FUNCTION

  FUNCTION transElemForward(this, rhor) RESULT(Ref)
    CLASS(OprRadarRef_t), INTENT(IN) :: this

    REAL(r_kind) :: rhor, Ref

    ! ref = this%a1*(rhor)**this%a2
    ref = this%a1 + this%a2 * LOG10(rhor)

    IF (rhor .EQ. 0) ref = 0

  END FUNCTION

  FUNCTION transElemAdjointMultiply(this, in, rhor) RESULT(out)
    CLASS(OprRadarRef_t), INTENT(IN) :: this

    REAL(r_kind) :: in, rhor
    REAL(r_kind) :: out

    ! out = this%a1*(rhor)**(this%a2 - 1)*in      ! d/d(rho) * ()
    out = this%a2 / rhor / LOG(10.0D0) * in      ! d/d(rho) * ()
    ! PRINT *, 'out: ', out
  END FUNCTION

END MODULE OprRadarRef_m
