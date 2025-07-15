!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/21, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BMatrix_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE BFieldBase_m, ONLY: BFieldBase_t
  USE NML_BMatrix_m, ONLY: NML_BMatrix
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE BFieldLaplace_m, ONLY: BFieldLaplace_t
  USE BFieldEnLoc_m, ONLY: BFieldEnLoc_t
  USE BFieldHybInc_m, ONLY: BFieldHybInc_t
  USE Ctl2State_m, ONLY: Ctl2State_t
  USE ObsSet_m, ONLY: ObsSet_t

!#define TRACE_PRESSURE

  TYPE, EXTENDS(AuxTypeSG_t):: BMatrix_t
    CLASS(BFieldBase_t), ALLOCATABLE :: BFields(:)
    TYPE(Ctl2State_t) :: Ctl2State

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC, PASS(this) :: inverse_multiply
    PROCEDURE, PUBLIC, PASS(this) :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC, PASS(this) :: sqrt_inverse_multiply_tl
    PROCEDURE, PUBLIC, PASS(this) :: sqrt_inverse_multiply_adjoint
    PROCEDURE, PUBLIC, PASS(this) :: selectBMat
    PROCEDURE, PUBLIC, PASS(this) :: initialize

    GENERIC :: OPERATOR(.SQRTINVMUL.) => sqrt_inverse_multiply
    ! GENERIC :: OPERATOR(.SQRTINVMULTL.) => sqrt_inverse_multiply_tl
    ! GENERIC :: OPERATOR(.SQRTINVMULADJ.) => sqrt_inverse_multiply_adjoint

  END TYPE BMatrix_t
CONTAINS

  SUBROUTINE initialize(this, configFile, sg, Bsolver, Y)
    CLASS(BMatrix_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(len=*), OPTIONAL :: Bsolver
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    !this%AuxTypeSG_t = AuxTypeSG_t(sg)
    CALL this%AuxTypeSG_t%aux_initialize(sg)
    CALL this%selectBMat(configFile, sg, Bsolver, Y)
    CALL this%Ctl2State%initialize(configFile)

  END SUBROUTINE

  SUBROUTINE selectBMat(this, configFile, sg, Bsolver, Y)
    CLASS(BMatrix_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(len=*), INTENT(IN) :: Bsolver
    CHARACTER(LEN=10), ALLOCATABLE :: varList(:)
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    CALL NML_BMatrix(configFile, varList)
    PRINT *, 'varList of B Mat: ', TRIM(Bsolver),varList

    BLOCK
      INTEGER(i_kind) :: i

      IF (TRIM(Bsolver) == 'EnLoc') THEN
        BLOCK
          ALLOCATE (BFieldEnLoc_t::this%BFields(UBOUND(varList, 1)))
        END BLOCK
      ELSE IF (TRIM(Bsolver) == 'HybInc') THEN
        BLOCK
          ALLOCATE (BFieldHybInc_t::this%BFields(UBOUND(varList, 1)))
        END BLOCK
      ELSE ! IF (TRIM(Bsolver) == 'Laplace') THEN ! Yuanfu Xie set the Laplace operator as default 2025-02-13
        ! The original setting of the above IF makes unspecified Bsolver option.
        BLOCK
          ALLOCATE (BFieldLaplace_t::this%BFields(UBOUND(varList, 1)))
        END BLOCK
      END IF

      ! The Bfield types input from the config files
      DO i = LBOUND(varList, 1), UBOUND(varList, 1)
        ASSOCIATE (BField => this%BFields(i))
          CALL BField%initialize(configFile, sg, varList(i), Y)
        END ASSOCIATE
      END DO
    END BLOCK
    IF (ALLOCATED(varList)) DEALLOCATE (varList)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(BMatrix_t), INTENT(INOUT) :: this
    PRINT *, 'Start of destructor of BMatrix.'

    IF (ALLOCATED(this%BFields)) DEALLOCATE (this%BFields)

    PRINT *, 'End of destructor of BMatrix.'

  END SUBROUTINE destructor

  SUBROUTINE inverse_multiply(this, Xm)
    IMPLICIT NONE
    CLASS(BMatrix_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm

  END SUBROUTINE inverse_multiply

  FUNCTION sqrt_inverse_multiply(this, Xm) RESULT(XX)
    IMPLICIT NONE
    CLASS(BMatrix_t), INTENT(IN) :: this
    TYPE(State_t), INTENT(IN) :: Xm
    INTEGER :: i, j
    TYPE(State_t) :: XX

    XX = Xm

    CALL this%Ctl2State%transFwdNonLinear(XX)

#ifdef TRACE_PRESSURE
    ! Debugging: Yuanfu Xie 2022-10-06:
    WRITE (*, 1) (XX%Fields(i)%Get_Name(), i=LBOUND(XX%Fields, 1), UBOUND(XX%Fields, 1))
    WRITE (*, 2) (this%BFields(i)%Get_Name(), i=LBOUND(this%BFields, 1), UBOUND(this%BFields, 1))
1   FORMAT('sqrt_inverse_multiply XX: ', 10(1X, A10))
2   FORMAT('sqrt_inverse_multiply BF: ', 10(1X, A10))
    PRINT *, 'XX pressure:', MINVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :)), &
      MAXVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :))
#endif

    DO j = LBOUND(XX%Fields, 1), UBOUND(XX%Fields, 1)
      DO i = LBOUND(this%BFields, 1), UBOUND(this%BFields, 1)
        ! Ensure the types of calculating B mat and X are matched.
        IF (TRIM(this%BFields(i)%Get_Name()) .EQ. TRIM(XX%fields(j)%Get_Name())) THEN
          CALL this%BFields(i)%sqrtInvMul(XX%fields(j))
        END IF
      END DO
    END DO
#ifdef TRACE_PRESSURE
    PRINT *, 'XX pressure After:', MINVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :)), &
      MAXVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :))
#endif

  END FUNCTION sqrt_inverse_multiply

  FUNCTION sqrt_inverse_multiply_tl(this, Xm, X) RESULT(XX)
    IMPLICIT NONE
    CLASS(BMatrix_t), INTENT(IN) :: this
    TYPE(State_t), INTENT(IN) :: Xm
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER :: i, j
    TYPE(State_t) :: XX

    XX = Xm

    CALL this%Ctl2State%transFwdTanLinear(XX,X)

#ifdef TRACE_PRESSURE
    ! Debugging: Yuanfu Xie 2022-10-06:
    WRITE (*, 1) (XX%Fields(i)%Get_Name(), i=LBOUND(XX%Fields, 1), UBOUND(XX%Fields, 1))
    WRITE (*, 2) (this%BFields(i)%Get_Name(), i=LBOUND(this%BFields, 1), UBOUND(this%BFields, 1))
1   FORMAT('sqrt_inverse_multiply_tl XX: ', 10(1X, A10))
2   FORMAT('sqrt_inverse_multiply_tl BF: ', 10(1X, A10))
    PRINT *, 'XX pressure:', MINVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :)), &
      MAXVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :))
#endif

    DO j = LBOUND(XX%Fields, 1), UBOUND(XX%Fields, 1)
      DO i = LBOUND(this%BFields, 1), UBOUND(this%BFields, 1)
        ! Ensure the types of calculating B mat and X are matched.
        IF (TRIM(this%BFields(i)%Get_Name()) .EQ. TRIM(XX%fields(j)%Get_Name())) THEN
          CALL this%BFields(i)%sqrtInvMul(XX%fields(j))
        END IF
      END DO
    END DO
#ifdef TRACE_PRESSURE
    PRINT *, 'XX pressure After:', MINVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :)), &
      MAXVAL(XX%fields(XX%getVarIdx('pres'))%DATA(1, :, :))
#endif

  END FUNCTION sqrt_inverse_multiply_tl

  FUNCTION sqrt_inverse_multiply_adjoint(this, Xm, X) RESULT(XX)
    IMPLICIT NONE
    CLASS(BMatrix_t), INTENT(IN) :: this
    TYPE(State_t), INTENT(IN) :: Xm
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    INTEGER :: i, j
    TYPE(State_t) :: XX

    XX = Xm
    DO j = LBOUND(XX%Fields, 1), UBOUND(XX%Fields, 1)
      DO i = LBOUND(this%BFields, 1), UBOUND(this%BFields, 1)
        ! Ensure the types of calculating B mat and X are matched.
        IF (TRIM(this%BFields(i)%Get_Name()) .EQ. TRIM(XX%fields(j)%Get_Name())) THEN
          IF (PRESENT(X)) THEN
            CALL this%BFields(i)%sqrtInvMulAdj(XX%fields(j), X%fields(j))
          ELSE
            CALL this%BFields(i)%sqrtInvMulAdj(XX%fields(j))
          END IF
          ! PRINT *, 'Matched! ',TRIM(XX%fields(j)%Get_Name()),' and ',TRIM(this%BFields(j)%Get_Name())
        END IF
      END DO
    END DO

    CALL this%Ctl2State%transAdjMultiply(XX, X)
  END FUNCTION sqrt_inverse_multiply_adjoint

END MODULE BMatrix_m
