!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.RMatrix
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/10, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE RMatrix_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsField_m, ONLY: ObsField_t
  USE AuxTypeObs_m, ONLY: AuxTypeObs_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE RField_m, ONLY: RField_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  IMPLICIT NONE

  TYPE, EXTENDS(AuxTypeObs_t):: RMatrix_t
    TYPE(RField_t), ALLOCATABLE :: RFields(:)
    CHARACTER(LEN=1024) :: configFile

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    PROCEDURE, PUBLIC, PASS(this) :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC, PASS(this) :: sqrt_inverse_multiply_adjoint

    GENERIC :: OPERATOR(.SQRTINVMUL.) => sqrt_inverse_multiply
    GENERIC :: OPERATOR(.SQRTINVMULADJ.) => sqrt_inverse_multiply_adjoint

  END TYPE RMatrix_t

CONTAINS

  SUBROUTINE initialize(this, configFile, Y, sg, Nocoarest)
    CLASS(RMatrix_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), INTENT(IN) :: Y
    INTEGER(i_kind), OPTIONAL, INTENT(IN) :: Nocoarest(:)
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    ! Yuanfu Xie 2024-05-08: deleted unused varList

    this%AuxTypeObs_t = AuxTypeObs_t(Y%mpObs)
    this%configFile = configFile

    BLOCK
      INTEGER(i_kind) :: i

      ! Clear the RFields if they are used before:
      IF (ALLOCATED(this%RFields)) DEALLOCATE (this%RFields)
      ALLOCATE (this%RFields(UBOUND(Y%ObsFields, 1)))

      ! The RField types input from the config files
      DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
        IF (PRESENT(Nocoarest)) THEN
          CALL this%RFields(i)%initialize(configFile, Y%mpObs, Y%ObsFields(i), sg, Nocoarest(i))
        ELSE
          CALL this%RFields(i)%initialize(configFile, Y%mpObs, Y%ObsFields(i), sg)
        END IF
      END DO
    END BLOCK
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(RMatrix_t), INTENT(INOUT) :: this
    PRINT *, 'Start of destructor of RMatrix.'

    IF (ALLOCATED(this%RFields)) DEALLOCATE (this%RFields)

    PRINT *, 'End of destructor of RMatrix.'

  END SUBROUTINE destructor

  SUBROUTINE inverse_multiply(this, Yin)
    IMPLICIT NONE
    CLASS(RMatrix_t) :: this
    TYPE(ObsSet_t), INTENT(INOUT) :: Yin

  END SUBROUTINE inverse_multiply

  FUNCTION sqrt_inverse_multiply(this, Yin) RESULT(Yout)
    IMPLICIT NONE
    CLASS(RMatrix_t), INTENT(IN) :: this
    TYPE(ObsSet_t), INTENT(IN) :: Yin
    TYPE(ObsSet_t) :: Yout
    INTEGER :: i, j

    Yout = Yin
    DO j = LBOUND(Yout%ObsFields, 1), UBOUND(Yout%ObsFields, 1)
      DO i = LBOUND(this%RFields, 1), UBOUND(this%RFields, 1)
        IF (TRIM(this%RFields(i)%Get_Name()) .EQ. TRIM(Yout%ObsFields(j)%Get_Id_Name())) THEN
          CALL this%RFields(i)%sqrt_inverse_multiply(Yout%ObsFields(j))
        END IF
      END DO
    END DO

  END FUNCTION sqrt_inverse_multiply

  FUNCTION sqrt_inverse_multiply_adjoint(this, Yin) RESULT(Yout)
    IMPLICIT NONE
    CLASS(RMatrix_t), INTENT(IN) :: this
    TYPE(ObsSet_t), INTENT(IN) :: Yin
    TYPE(ObsSet_t) :: Yout
    INTEGER :: i, j

    Yout = Yin
    DO j = LBOUND(Yout%ObsFields, 1), UBOUND(Yout%ObsFields, 1)
      DO i = LBOUND(this%RFields, 1), UBOUND(this%RFields, 1)
        IF (TRIM(this%RFields(i)%Get_Name()) .EQ. TRIM(Yout%ObsFields(j)%Get_Id_Name())) THEN
          CALL this%RFields(i)%sqrt_inverse_multiply_adjoint(Yout%ObsFields(j))
        END IF
      END DO
    END DO

  END FUNCTION sqrt_inverse_multiply_adjoint

END MODULE RMatrix_m
