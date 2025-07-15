MODULE ECMIO_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE AdvanceTime_m
  USE parameters_m, ONLY: degree2radian

  TYPE ECMIO_t
    CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: inFileNames
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: varNames
    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: time_unix
    INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE :: time_gmt
    INTEGER(i_kind)                                      :: numFiles, numVars, nx, ny, nt
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: lon2D, lat2D
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: precipitation, & ! X x Y x T
                                                       t2m, &
                                                       q2m, &
                                                       u10m, &
                                                       v10m, &
                                                       u200m, &
                                                       v200m, &
                                                       u100m, &
                                                       v100m, &
                                                       pressure
  CONTAINS
    PROCEDURE, PRIVATE, PASS(this) :: read_value
    FINAL :: destructor
  END TYPE ECMIO_t

  INTERFACE ECMIO_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  FUNCTION constructor(ECMFileNames, ECMVarNames, mdate) RESULT(this)
    IMPLICIT NONE
    TYPE(ECMIO_t) :: this
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: ECMFileNames
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: ECMVarNames
    INTEGER(i_kind), DIMENSION(:, :), INTENT(IN) :: mdate
    INTEGER :: i

    this%numFiles = UBOUND(ECMFileNames, 1)
    this%numVars = UBOUND(ECMVarNames, 1)
    this%nt = this%numFiles

    ALLOCATE (this%inFileNames(this%numFiles))
    ALLOCATE (this%time_gmt(6, this%numFiles))
    ALLOCATE (this%varNames(this%numVars))

    this%inFileNames = ECMFileNames
    this%time_gmt = mdate
    this%varNames = ECMVarNames

    ALLOCATE (this%time_unix(this%nt))
    DO i = 1, this%nt
      PRINT *, "GMT Time: ", this%time_gmt(:, i)
      CALL Time_GMT_to_Unix(this%time_gmt(:, i), this%time_unix(i))
      PRINT *, "Unix Time: ", this%time_unix(i)
    END DO

    CALL this%read_value()

  END FUNCTION constructor

  SUBROUTINE read_value(this)
    IMPLICIT NONE
    CLASS(ECMIO_t)                                  :: this
    TYPE(NcDataset)                                 :: dataset
    TYPE(NcVariable)                                :: var
    TYPE(NcDimension)                               :: dim
    INTEGER(i_kind)                                 :: t, i, tmpIdx
    CHARACTER(LEN=1024)                             :: inFileName
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: tmpVar

    dataset = NcDataset(TRIM(this%inFileNames(1)), "r")
    ! Get dimensions
    dim = dataset%getDimension("lon"); this%nx = dim%getLength()
    dim = dataset%getDimension("lat"); this%ny = dim%getLength()

    PRINT *, "nx, ny, nt:", this%nx, this%ny, this%nt

    ! Allocate necessary arrays
    ALLOCATE (this%lon2D(this%nx, this%ny))
    ALLOCATE (this%lat2D(this%nx, this%ny))

    ! Read longitude and latitude
    var = dataset%getVariable("lon")
    CALL var%getData(this%lon2D)
    this%lon2D = this%lon2D * degree2radian
    var = dataset%getVariable("lat")
    CALL var%getData(this%lat2D)
    this%lat2D = this%lat2D * degree2radian

    ! Allocate and read other variables
    ALLOCATE (this%t2m(this%nx, this%ny, this%nt))
    ALLOCATE (this%q2m(this%nx, this%ny, this%nt))
    ALLOCATE (this%precipitation(this%nx, this%ny, this%nt))
    ALLOCATE (this%u10m(this%nx, this%ny, this%nt))
    ALLOCATE (this%v10m(this%nx, this%ny, this%nt))
    ALLOCATE (this%u200m(this%nx, this%ny, this%nt))
    ALLOCATE (this%v200m(this%nx, this%ny, this%nt))
    ALLOCATE (this%u100m(this%nx, this%ny, this%nt))
    ALLOCATE (this%v100m(this%nx, this%ny, this%nt))
    ALLOCATE (this%pressure(this%nx, this%ny, this%nt))

    DO t = 1, this%nt
      PRINT *, "Reading file hours: ", t, " of file ", this%nt
      dataset = NcDataset(TRIM(this%inFileNames(t)), "r")
      tmpIdx = this%time_gmt(4, t) + 1
      PRINT *, "Reading time index: ", tmpIdx

      ! Read all data at once
      IF (ALLOCATED(this%t2m)) THEN
        var = dataset%getVariable("t2m")
        CALL var%getData(tmpVar)
        PRINT *, "SHAPE", SHAPE(tmpVar)
        this%t2m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading t2m"
      END IF

      IF (ALLOCATED(this%q2m)) THEN
        var = dataset%getVariable("q2m")
        CALL var%getData(tmpVar)
        this%q2m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading q2m"
      END IF

      IF (ALLOCATED(this%precipitation)) THEN
        var = dataset%getVariable("precipitation")
        CALL var%getData(tmpVar)
        this%precipitation(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading precipitation"
      END IF

      IF (ALLOCATED(this%u10m)) THEN
        var = dataset%getVariable("u10m")
        CALL var%getData(tmpVar)
        this%u10m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading u10m"
      END IF

      IF (ALLOCATED(this%v10m)) THEN
        var = dataset%getVariable("v10m")
        CALL var%getData(tmpVar)
        this%v10m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading v10m"
      END IF

      IF (ALLOCATED(this%u200m)) THEN
        var = dataset%getVariable("u200m")
        CALL var%getData(tmpVar)
        this%u200m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading u200m"
      END IF

      IF (ALLOCATED(this%v200m)) THEN
        var = dataset%getVariable("v200m")
        CALL var%getData(tmpVar)
        this%v200m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading v200m"
      END IF

      IF (ALLOCATED(this%u100m)) THEN
        var = dataset%getVariable("u100m")
        CALL var%getData(tmpVar)
        this%u100m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading u100m"
      END IF

      IF (ALLOCATED(this%v100m)) THEN
        var = dataset%getVariable("v100m")
        CALL var%getData(tmpVar)
        this%v100m(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading v100m"
      END IF

      IF (ALLOCATED(this%pressure)) THEN
        var = dataset%getVariable("pressure")
        CALL var%getData(tmpVar)
        this%pressure(:, :, t) = tmpVar(:, :, tmpIdx)
        PRINT *, "Reading pressure"
      END IF
    END DO

    PRINT *, "DONE ECMIO_t"

    DEALLOCATE (tmpVar)
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(ECMIO_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%t2m)) DEALLOCATE (this%t2m)
    IF (ALLOCATED(this%q2m)) DEALLOCATE (this%q2m)
    IF (ALLOCATED(this%precipitation)) DEALLOCATE (this%precipitation)
    IF (ALLOCATED(this%u10m)) DEALLOCATE (this%u10m)
    IF (ALLOCATED(this%v10m)) DEALLOCATE (this%v10m)
    IF (ALLOCATED(this%u200m)) DEALLOCATE (this%u200m)
    IF (ALLOCATED(this%v200m)) DEALLOCATE (this%v200m)
    IF (ALLOCATED(this%u100m)) DEALLOCATE (this%u100m)
    IF (ALLOCATED(this%v100m)) DEALLOCATE (this%v100m)
    IF (ALLOCATED(this%pressure)) DEALLOCATE (this%pressure)
    IF (ALLOCATED(this%lon2D)) DEALLOCATE (this%lon2D)
    IF (ALLOCATED(this%lat2D)) DEALLOCATE (this%lat2D)
    IF (ALLOCATED(this%time_unix)) DEALLOCATE (this%time_unix)

    PRINT *, 'destructor ECMIO works'

  END SUBROUTINE destructor

END MODULE ECMIO_m
