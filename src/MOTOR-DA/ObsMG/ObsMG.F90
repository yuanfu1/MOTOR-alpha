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
!! This module provides an interface of accessing background fields for multigrid.
MODULE obsMG_m

  ! Observation types:
  USE ObsBase_m
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t
  USE ObsRadarRef_m, ONLY: ObsRadarRef_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  !USE ObsAGRI_m, ONLY: ObsAGRI_t
  !USE ObsGIIRS_m, ONLY: ObsGIIRS_t
  USE Satob_m, ONLY: Satob_t
  USE ObsGWST_m, ONLY: ObsGWST_t
  USE ObsLBUOY_m, ONLY: ObsLBUOY_t
  USE ObsSING_m, ONLY: ObsSING_t

  USE ObsSetArray_m, ONLY: ObsSetArray_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE YAMLRead_m
  USE ObsUtilities_m

  USE geometry_m, ONLY: geometry_t

  TYPE :: ObsDataList_t
    CLASS(ObsBase_t), POINTER :: object
  END TYPE ObsDataList_t

  TYPE :: obsMG_t
    CHARACTER(LEN=8), ALLOCATABLE :: obsList(:)
    CHARACTER(LEN=1024) :: yamlFile
    INTEGER(i_kind) :: numObsTyps
    ! Data and matrices for each observation types:
    CLASS(ObsDataList_t), POINTER :: observations(:)
    TYPE(RMatrix_t), ALLOCATABLE :: rMatrices(:)

    ! Thinned obs:
    TYPE(ObsSetArray_t) :: thinnedObs
    ! Obs on multi-processors:
    TYPE(MPObs_t), ALLOCATABLE :: mpObs(:)
    TYPE(geometry_t), POINTER :: geo
  CONTAINS
    PROCEDURE :: initialize_s
    PROCEDURE :: ingestions_s
    PROCEDURE :: thinningSG_s
    PROCEDURE :: R_Matrices_s
    PROCEDURE :: forwardOpr_s
    PROCEDURE :: tangentOpr_s
    PROCEDURE :: adjointOpr_s

    ! R matrices and their adjoints:
    PROCEDURE :: sqrt_inverse_multiply
    PROCEDURE :: sqrt_inverse_multiply_adjoint

    ! Subroutine version of the above operators:
    PROCEDURE :: rSqrtInverseMultiply
    PROCEDURE :: rSqrtInverseMultiplyAdjoint

    PROCEDURE, PUBLIC :: addition
    PROCEDURE, PUBLIC :: subtract
    PROCEDURE, PUBLIC :: multiply
    ! PROCEDURE, PUBLIC :: dotPrdct

    !GENERIC :: OPERATOR(.DOT.) => dotPrdct
    GENERIC :: OPERATOR(*) => multiply
    GENERIC :: OPERATOR(+) => addition
    GENERIC :: OPERATOR(-) => subtract

    GENERIC :: OPERATOR(.SQRTINVMUL.) => sqrt_inverse_multiply
    GENERIC :: OPERATOR(.SQRTINVMULADJ.) => sqrt_inverse_multiply_adjoint

    FINAL :: destroy_s
  END TYPE obsMG_t

CONTAINS
  ! Destructor:
  SUBROUTINE destroy_s(this)

    IMPLICIT NONE

    TYPE(obsMG_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i

    ! Deallocate memory of observation types:
    IF (ASSOCIATED(this%observations)) THEN
      DO i = 1, this%numObsTyps
        DEALLOCATE (this%observations(i)%object)
      END DO
      DEALLOCATE (this%observations)
    END IF

    IF (ALLOCATED(this%rMatrices)) DEALLOCATE (this%rMatrices)

    ! Deallocate memory of mpObs:
    IF (ALLOCATED(this%mpObs)) DEALLOCATE (this%mpObs)

    PRINT *, 'ObsMG destroyed...'
  END SUBROUTINE destroy_s

  SUBROUTINE initialize_s(this, geo, yamlFile)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geo
    CHARACTER(LEN=1024) :: yamlFile

    ! Local variables:
    INTEGER(i_kind) :: i, istatus, mgStart, mgEnd

    this%yamlFile = yamlFile

    this%geo => geo

    mgStart = this%geo%mg%mg_coarsest
    mgEnd = this%geo%mg%mg_finest

    istatus = yaml_get_var(yamlFile, 'IO', 'obsList', this%obsList)

    this%numObsTyps = UBOUND(this%obsList, 1)

    ! Allocate variable arrays:
    ALLOCATE (this%observations(this%numObsTyps), & !this%thinnedObs(this%numObsTyps), &
              this%rMatrices(this%numObsTyps))
    ALLOCATE (this%mpObs(mgStart:mgEnd))

    ! Initialize observation parallel processes:
    DO i = mgStart, mgEnd
      CALL this%mpObs(i)%initializeMPObs(this%geo%mg%sg(i))
    END DO

  END SUBROUTINE initialize_s

  ! This routine reads in all observations at a given grid:
  SUBROUTINE ingestions_s(this, state)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(State_t), INTENT(IN) :: state

    ! Local variables:
    INTEGER(i_kind) :: i

    PRINT *, 'Entering ingestion_s...', this%numObsTyps
    DO i = 1, this%numObsTyps
      SELECT CASE (this%obsList(i))
      CASE ('surfaces')
        ALLOCATE (ObsSurface_t:: this%observations(i)%object)
      CASE ('sounding')
        ALLOCATE (ObsSound_t:: this%observations(i)%object)
      CASE ('profiler')
        ALLOCATE (ObsVwpw_t:: this%observations(i)%object)
      CASE ('radarRef')
        ALLOCATE (ObsRadarRef_t:: this%observations(i)%object)
      CASE ('radarVel')
        ALLOCATE (ObsRadarVel_t:: this%observations(i)%object)
      ! CASE ('fy4_agri')
      !   ALLOCATE (ObsAGRI_t:: this%observations(i)%object)
      CASE ('cloudWnd')
        ALLOCATE (Satob_t:: this%observations(i)%object)
      CASE ('shipRpts')
        ALLOCATE (ObsGWST_t:: this%observations(i)%object)
      CASE ('bouyance')
        ALLOCATE (ObsLBUOY_t:: this%observations(i)%object)
      CASE ('pilotRpt')
        ALLOCATE (ObsSING_t:: this%observations(i)%object)
      END SELECT

      ! Unified process of data types:
      CALL this%observations(i)%object%ObsInitial(this%yamlFile)
      CALL this%observations(i)%object%ObsIngest(state)
    END DO
  END SUBROUTINE ingestions_s

  ! Apply thinning schemes on a given grid level of state input:
  SUBROUTINE thinningSG_s(this, state)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(State_t), INTENT(IN) :: state

    ! Local variables:
    INTEGER(i_kind) :: i

    ! Initialize thinnedObs on a given sg:
    CALL this%thinnedObs%initialize_s(this%numObsTyps, state%sg%gLevel, state%sg, this%yamlFile)

    DO i = 1, this%numObsTyps
      CALL this%observations(i)%object%ObsThinning( &
        state, this%thinnedObs%obsArray(i), this%mpObs(state%sg%gLevel), .TRUE., .FALSE.)
    END DO
  END SUBROUTINE thinningSG_s

  ! Apply thinning schemes on a given grid level of state input:
  SUBROUTINE R_Matrices_s(this, sg)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    ! Local variables:
    INTEGER(i_kind) :: i

    DO i = 1, this%numObsTyps
      CALL this%rMatrices(i)%initialize(this%yamlFile, this%thinnedObs%obsArray(i), sg)
    END DO
  END SUBROUTINE R_Matrices_s

  ! Apply observation forward operator on a given model state to ObsSet:
  FUNCTION forwardOpr_s(this, state) RESULT(Y)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(State_t), INTENT(IN) :: state
    TYPE(ObsSetArray_t) :: Y

    ! Local variables:
    CHARACTER(LEN=1024) :: configfile ! Dummy as the current ObsSet_t does not use it
    INTEGER(i_kind) :: i

    ! Initialize output Y:
    Y = this%thinnedObs

    DO i = 1, this%numObsTyps
      ! Use thinned obs structure but not use its values
      Y%obsArray(i) = this%observations(i)%object%ObsForward(state, Y%obsArray(i))
    END DO
  END FUNCTION forwardOpr_s

  ! Apply observation forward operator on a given model state to ObsSet:
  FUNCTION tangentOpr_s(this, dstate, state) RESULT(Y)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(State_t), INTENT(IN) :: dstate, state
    TYPE(ObsSet_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i
    TYPE(ObsSet_t) :: yi(this%numObsTyps)

    DO i = 1, this%numObsTyps
      yi(i) = this%observations(i)%object%ObsForward(state, this%thinnedObs%obsArray(i)) &
              - this%thinnedObs%obsArray(i)
    END DO
    CALL ObsConcat_s(yi, Y)
  END FUNCTION tangentOpr_s

  ! Apply observation forward operator on a given model state to ObsSet:
  FUNCTION adjointOpr_s(this, D, X) RESULT(dX)

    IMPLICIT NONE

    CLASS(obsMG_t) :: this
    TYPE(ObsSetArray_t) :: D
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: dX

    ! Local variables:
    INTEGER(i_kind) :: i, j, k

    dX = X%zeroCopy()

    DO i = 1, this%numObsTyps
      dX = dX + this%observations(i)%object%ObsAdjoint(D%obsArray(i), X)
    END DO
    CALL dX%exHalo()
  END FUNCTION adjointOpr_s

  ! R matrices:

  FUNCTION sqrt_inverse_multiply(this, Yin) RESULT(Yout)
    IMPLICIT NONE
    CLASS(obsMG_t), INTENT(IN) :: this
    TYPE(ObsSet_t), INTENT(IN) :: Yin(this%numObsTyps)
    TYPE(ObsSet_t), ALLOCATABLE :: Yout(:)

    ! Local variable:
    INTEGER(i_kind) :: i
    TYPE(ObsSet_t), ALLOCATABLE :: temp(:)

    ! Check the memory:
    ! IF (ALLOCATED(Yout)) DEALLOCATE(Yout)
    ! ALLOCATE(Yout(this%numObsTyps),temp(this%numObsTyps))

    ! temp = yin

    ! ! For each obs type:
    ! DO i=1,this%numObsTyps
    !   Yout(i) = this%rMatrices(i).SQRTINVMUL.temp(i)
    ! END DO

  END FUNCTION sqrt_inverse_multiply

  FUNCTION sqrt_inverse_multiply_adjoint(this, Yin) RESULT(Yout)
    IMPLICIT NONE
    CLASS(obsMG_t), INTENT(IN) :: this
    TYPE(ObsSetArray_t), INTENT(IN) :: Yin
    TYPE(ObsSetArray_t) :: Yout
    INTEGER(i_kind) :: i

    ! Initialize Yout:
    Yout = Yin

    ! For each obs type:
    DO i = 1, this%numObsTyps
      Yout%obsArray(i) = this%RMatrices(i) .SQRTINVMULADJ.Yin%obsArray(i)
    END DO
  END FUNCTION sqrt_inverse_multiply_adjoint

  SUBROUTINE rSqrtInverseMultiply(this, Yin, Yout)
    IMPLICIT NONE
    CLASS(obsMG_t), INTENT(IN) :: this
    TYPE(ObsSetArray_t), INTENT(IN) :: Yin
    TYPE(ObsSetArray_t), INTENT(OUT) :: Yout

    ! Local variable:
    INTEGER(i_kind) :: i

    ! Initialize Yout:
    Yout = Yin

    ! For each obs type:
    DO i = 1, this%numObsTyps
      Yout%obsArray(i) = this%rMatrices(i) .SQRTINVMUL.Yin%obsArray(i)
    END DO

  END SUBROUTINE rSqrtInverseMultiply

  SUBROUTINE rSqrtInverseMultiplyAdjoint(this, Yin, Yout)
    IMPLICIT NONE
    CLASS(obsMG_t), INTENT(IN) :: this
    TYPE(ObsSetArray_t), INTENT(IN) :: Yin
    TYPE(ObsSetArray_t), ALLOCATABLE :: Yout

    ! Local variables:
    INTEGER(i_kind) :: i

    ! Initialize Yout:
    Yout = Yin

    ! For each obs type:
    DO i = 1, this%numObsTyps
      Yout%obsArray(i) = this%RMatrices(i) .SQRTINVMULADJ.Yin%obsArray(i)
    END DO
  END SUBROUTINE rSqrtInverseMultiplyAdjoint

  ! Arithematics operators:

  FUNCTION addition(X1, X2) RESULT(X3)

    IMPLICIT NONE

    CLASS(obsMG_t), INTENT(IN) :: X1
    TYPE(obsMG_t), INTENT(IN) :: X2
    TYPE(obsMG_t) ::X3

    ! Local variables:
    INTEGER(i_kind) :: i, j

    X3 = X1
    DO j = 1, UBOUND(X3%thinnedObs%obsArray, 1)
      DO i = 1, UBOUND(X3%thinnedObs%obsArray(j)%ObsFields, 1)
        X3%thinnedObs%obsArray(j)%ObsFields(i) = &
          X2%thinnedObs%obsArray(j)%ObsFields(i) + X1%thinnedObs%obsArray(j)%ObsFields(i)
      END DO
    END DO

  END FUNCTION addition

  FUNCTION subtract(X1, X2) RESULT(X3)

    IMPLICIT NONE

    CLASS(obsMG_t), INTENT(IN) :: X1
    TYPE(obsMG_t), INTENT(IN) :: X2
    TYPE(obsMG_t) ::X3

    ! Local variables:
    INTEGER(i_kind) :: i, j

    X3 = X1
    DO j = 1, UBOUND(X3%thinnedObs%obsArray, 1)
      DO i = 1, UBOUND(X3%thinnedObs%obsArray(j)%ObsFields, 1)
        X3%thinnedObs%obsArray(j)%ObsFields(i) = &
          X2%thinnedObs%obsArray(j)%ObsFields(i) - X1%thinnedObs%obsArray(j)%ObsFields(i)
      END DO
    END DO

  END FUNCTION subtract

  FUNCTION multiply(X1, c) RESULT(X3)

    IMPLICIT NONE

    CLASS(obsMG_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: c
    TYPE(obsMG_t) ::X3

    ! Local variables:
    INTEGER(i_kind) :: i, j

    X3 = X1
    DO j = 1, UBOUND(X3%thinnedObs%obsArray, 1)
      DO i = 1, UBOUND(X3%thinnedObs%obsArray(j)%ObsFields, 1)
        X3%thinnedObs%obsArray(j)%ObsFields(i) = &
          X1%thinnedObs%obsArray(j)%ObsFields(i) * c
      END DO
    END DO

  END FUNCTION multiply
END MODULE obsMG_m
