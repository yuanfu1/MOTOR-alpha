! A new module for observation parent abstract class
! Created by Ting Shu @20230202
MODULE ObsBase_m2
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE YAMLRead_m
  USE mpObs_m, ONLY: mpObs_t
  IMPLICIT NONE

  TYPE, ABSTRACT :: ObsBase_t2
    CHARACTER(LEN=20) :: obsType
    INTEGER(i_kind) :: numFiles, numVars
    CHARACTER(LEN=1024), ALLOCATABLE :: fileNames(:)
    CHARACTER(LEN=20), ALLOCATABLE :: varNames(:)
    ! REAL(r_kind) :: ztop
    ! black site yaml file
    ! CHARACTER(LEN=1024) :: blackList = ''
    CHARACTER(LEN=1024) :: configFile
    REAL(r_kind) :: correlation_threshold

  CONTAINS
    ! PROCEDURE(oneAug), DEFERRED :: obsInitial
    PROCEDURE, PUBLIC :: obsInitial
    PROCEDURE(ingest), DEFERRED :: obsIngest
    ! PROCEDURE(noAugs), DEFERRED :: obsQC
    PROCEDURE(thinning), DEFERRED :: obsThinning
    PROCEDURE(release), DEFERRED :: obsDeallocate

    ! PROCEDURE, PUBLIC :: obsThinning
    ! PROCEDURE, PUBLIC :: obsDeallocate
    ! PROCEDURE, PUBLIC :: obsProcess

    GENERIC :: obsQC => obsQC_conv, obsQC_radar
    PROCEDURE :: obsQC_conv, obsQC_radar
  END TYPE ObsBase_t2

  !TYPE :: Valided_Data
  !    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: valided_obsData, valided_obsErrs, valided_pos
  !    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: valided_time
  !    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: valided_ID
  !    INTEGER(i_kind), ALLOCATABLE :: idx_grd_valid(:, :), idx_hgt_valid(:, :, :), idx_t_valid(:, :)
  !    REAL(r_kind), ALLOCATABLE ::  coe_grd_valid(:, :), coe_hgt_valid(:, :, :), coe_t_valid(:, :)
  !END TYPE Valided_Data

  ABSTRACT INTERFACE
    SUBROUTINE ingest(this, X)
      IMPORT
      CLASS(ObsBase_t2) :: this
      TYPE(State_t) :: X
    END SUBROUTINE ingest

    ! SUBROUTINE noAugs(this)
    !     IMPORT
    !     CLASS(ObsBase_t2) :: this
    ! END SUBROUTINE noAugs

    SUBROUTINE thinning(this, state, thinObs, mpObs, useBkg)
      IMPORT
      CLASS(ObsBase_t2) :: this
      TYPE(State_t), INTENT(IN) :: state
      TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
      TYPE(mpObs_t), INTENT(IN) :: mpObs
      LOGICAL, INTENT(IN) :: useBkg
    END SUBROUTINE thinning

    SUBROUTINE release(this)
      IMPORT
      CLASS(ObsBase_t2) :: this
    END SUBROUTINE release
  END INTERFACE
CONTAINS
  SUBROUTINE obsQC_conv(this, state)
    IMPLICIT NONE
    CLASS(ObsBase_t2) :: this
    TYPE(State_t), INTENT(IN) :: state
  END SUBROUTINE obsQC_conv

  SUBROUTINE obsQC_radar(this, state, radar_types)
    IMPLICIT NONE
    CLASS(ObsBase_t2) :: this
    TYPE(State_t), INTENT(IN) :: state
    CHARACTER(LEN=20), INTENT(IN) :: radar_types
  END SUBROUTINE obsQC_radar

  SUBROUTINE obsInitial(this, configFile, obsType)
    IMPLICIT NONE
    CLASS(ObsBase_t2) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(LEN=*), INTENT(IN) :: obsType

    ! local variables
    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=20) :: obsInfo
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileNames(:)
    INTEGER(i_kind) :: i, status
    LOGICAL :: blackExist

    this%obsType = obsType
    obsInfo = TRIM(obsType)//"_info"

    this%configFile = configFile

    status = yaml_get_var(configFile, obsInfo, 'inputDir', obsFileDir)
    status = yaml_get_var(configFile, obsInfo, 'list', obsFileNames)
    this%numFiles = UBOUND(obsFileNames, 1)
    ALLOCATE (this%fileNames(this%numFiles), STAT=status)
    DO i = 1, this%numFiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileNames(i))
    END DO
    status = yaml_get_var(configFile, obsInfo, 'varNames', this%varNames)
    this%numVars = UBOUND(this%varNames, 1)

    ! get DA ztop
    ! status = yaml_get_var(configFile, 'DASpace', 'ztop', this%ztop)

    ! get black site yaml file address, each type of observation can add black list
    ! status = yaml_get_var(configFile, obsInfo, 'blackExist', blackExist)
    ! IF (blackExist) status = yaml_get_var(configFile, 'Verify', TRIM(this%obsType)//'_BlackListAdd', this%blackList)
    this%correlation_threshold = 0.25D0
  END SUBROUTINE obsInitial

  ! SUBROUTINE obsThinning(this, thinObs)
  !     IMPLICIT NONE
  !     CLASS(ObsBase_t2) :: this
  !     TYPE(ObsSet_t), INTENT(OUT) :: thinObs
  ! END SUBROUTINE obsThinning

  ! SUBROUTINE obsDeallocate(this)
  !     IMPLICIT NONE
  !     CLASS(ObsBase_t2) :: this
  ! END SUBROUTINE obsDeallocate

  ! SUBROUTINE obsProcess(this, configFile, obsType, X, thinObs, mpObs)
  !     IMPLICIT NONE
  !     CLASS(ObsBase_t2) :: this
  !     CHARACTER(LEN=1024), INTENT(IN) :: configFile
  !     CHARACTER(LEN=10), INTENT(IN) :: obsType
  !     TYPE(State_t), INTENT(IN) :: X
  !     TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
  !     TYPE(mpObs_t), INTENT(IN) :: mpObs

  !     CALL this%obsInitial(configFile, obsType)
  !     CALL this%obsIngest(X)
  !     CALL this%obsQC(X)
  !     CALL this%obsThinning(X, thinObs, mpObs, .TRUE.)
  !     CALL this%obsDeallocate()
  ! END SUBROUTINE obsProcess
END MODULE ObsBase_m2
