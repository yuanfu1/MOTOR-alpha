!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA/ObsMG
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-04-26   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2024/04/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides arithmatics operations for an array of obsSet data.
MODULE obsSetArray_m
  USE ObsSet_m, ONLY: ObsSet_t
  USE AuxTypeObs_m, ONLY: AuxTypeObs_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MPObs_m, ONLY: MPObs_t

  TYPE, EXTENDS(AuxTypeObs_t) :: obsSetArray_t
    PRIVATE

    INTEGER(i_kind), PUBLIC :: numTypes, glvl
    TYPE(ObsSet_t), PUBLIC, ALLOCATABLE :: obsArray(:)
    CHARACTER(LEN=20), ALLOCATABLE :: obsTypeName(:)
    CHARACTER(LEN=1024) :: yamlFile
    ! TYPE(MPObs_t), PUBLIC :: mpObs

  CONTAINS

    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: add, dot_multiply, subtract
    !PROCEDURE, PUBLIC :: getObsName_f

    FINAL :: destructor
    GENERIC :: OPERATOR(.DOT.) => dot_multiply
    GENERIC :: OPERATOR(+) => add
    GENERIC :: OPERATOR(-) => subtract
  END TYPE obsSetArray_t

CONTAINS

  FUNCTION GetObsName(this, i) RESULT(name)
    IMPLICIT NONE
    CLASS(obsSetArray_t) :: this
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i

    name = this%obsTypeName(i)
  END FUNCTION

  SUBROUTINE initialize_s(this, numTypes, glvl, sg, yamlFile)
    CLASS(obsSetArray_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numTypes, glvl
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile
    ! TYPE(MPObs_t), INTENT(IN) :: mpObs

    this%numTypes = numTypes
    this%yamlFile = yamlFile
    this%glvl = glvl

    ! Initialize the aux:
    this%AuxTypeObs_t = AuxTypeObs_t(this%mpObs)

    IF (.NOT. ALLOCATED(this%obsArray)) ALLOCATE (this%obsArray(this%numTypes))
    DO i = 1, this%numTypes
      this%obsArray(i) = ObsSet_t(this%yamlFile, this%mpobs)
    END DO
  END SUBROUTINE initialize_s

  FUNCTION add(X1, X2) RESULT(X3)
    CLASS(obsSetArray_t), INTENT(IN) :: X1
    TYPE(obsSetArray_t), INTENT(IN) :: X2
    TYPE(obsSetArray_t) :: X3

    ! Local variables:
    INTEGER(i_kind) :: i

    X3 = X1
    DO i = 1, X3%numTypes
      X3%obsArray(i) = X1%obsArray(i) + X2%obsArray(i)
    END DO

  END FUNCTION add

  FUNCTION dot_multiply(X1, X2)
    CLASS(obsSetArray_t), INTENT(IN) :: X1
    TYPE(obsSetArray_t), INTENT(IN) :: X2

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: dot_multiply

    dot_multiply = 0.0D0
    DO i = 1, X1%numTypes
      dot_multiply = dot_multiply + (X1%obsArray(i) .DOT.X2%obsArray(i))
    END DO
  END FUNCTION

  FUNCTION subtract(X1, X2) RESULT(X3)
    CLASS(obsSetArray_t), INTENT(IN) :: X1
    TYPE(obsSetArray_t), INTENT(IN) :: X2
    TYPE(obsSetArray_t) ::X3

    X3 = X1
    DO i = 1, X3%numTypes
      ! Note: even though substract is implemented in reversed order,
      !       this is to call an abstract data type obsSet_t type to
      !       subtract and the order is kept as a nature subtract:
      X3%obsArray(i) = X1%obsArray(i) - X2%obsArray(i)
    END DO
  END FUNCTION subtract

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(obsSetArray_t), INTENT(inout) :: this

    IF (ALLOCATED(this%obsArray)) DEALLOCATE (this%obsArray)
  END SUBROUTINE destructor
END MODULE obsSetArray_m
