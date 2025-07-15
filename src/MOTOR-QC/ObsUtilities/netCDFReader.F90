!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-Utilities/Utility.netCDFReader
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
! This is a universal netCDF file reader. It uses a yaml file as a netCDF table to extract data from
! the given file.
MODULE netCDFReader_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE netCDFTable_m, ONLY: netCDFTable_t

  TYPE :: netCDFReader_t
  CONTAINS
    PROCEDURE :: ncReader
  END TYPE netCDFReader_t

  ! INTERFACE netCDFReader_t
  !     PROCEDURE :: constructor
  ! END INTERFACE netCDFReader_t

CONTAINS

  SUBROUTINE ncReader(this, table, ncFile, domain, misVal, ncData, llht)
    CLASS(netCDFReader_t) :: this
    TYPE(netCDFTable_t), INTENT(IN) :: table
    CHARACTER(LEN=*), INTENT(IN) :: ncFile
    REAL(r_kind), INTENT(IN) :: domain(2, 2), misVal
    REAL(r_kind), INTENT(OUT), ALLOCATABLE :: ncData(:, :), llht(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iv, numVars, numValid
    LOGICAL :: Missing
    TYPE(NcVariable) :: var
    TYPE(NcDataset) :: nc
    REAL(r_kind), ALLOCATABLE :: temp(:, :), nc2Data(:, :, :), latlon(:, :, :), time(:, :), height(:, :)

    nc = NcDataset(ncFile, "r")

    numVars = UBOUND(table%VarNames, 1)

    ! Get the latlon:
    var = nc%getVariable(TRIM(table%latlonNames(1))); CALL var%getData(temp)
    ALLOCATE (latlon(2, UBOUND(temp, 1), UBOUND(temp, 2)), &
              time(UBOUND(temp, 1), UBOUND(temp, 2)), &
              height(UBOUND(temp, 1), UBOUND(temp, 2)))
    latlon(1, :, :) = temp
    var = nc%getVariable(TRIM(table%latlonNames(2))); CALL var%getData(temp)
    latlon(2, :, :) = temp

    ! Time and height:
    time = misVal; height = misVal
    IF (TRIM(table%timeName(1)) .NE. 'NONE') THEN
      var = nc%getVariable(TRIM(table%timeName(1))); CALL var%getData(temp)
      time = temp
    END IF
    IF (TRIM(table%heightName(1)) .NE. 'NONE') THEN
      var = nc%getVariable(TRIM(table%heightName(1))); CALL var%getData(temp)
      height = temp
    END IF

    ! Get the data:
    DO iv = 1, numVars
      var = nc%getVariable(TRIM(table%VarNames(iv))); CALL var%getData(temp)
      IF (iv .EQ. 1) ALLOCATE (nc2Data(UBOUND(temp, 1), UBOUND(temp, 2), numVars))
      nc2Data(:, :, iv) = temp
    END DO

    ! Remove data with all values missing:
    numValid = 0
    DO j = 1, UBOUND(temp, 2)
      DO i = 1, UBOUND(temp, 1)
        Missing = .FALSE.
        IF (latlon(1, i, j) .EQ. table%missingValue .OR. &
            latlon(2, i, j) .EQ. table%missingValue .OR. &
            latlon(1, i, j) .LT. domain(1, 1) .OR. &
            latlon(1, i, j) .GT. domain(2, 1) .OR. &
            latlon(2, i, j) .LT. domain(1, 2) .OR. &
            latlon(2, i, j) .GT. domain(2, 2)) THEN
          Missing = .TRUE.
          CYCLE
        END IF
        DO iv = 1, numVars
          IF (table%required(iv) .EQ. 1 .AND. &
              nc2Data(i, j, iv) .EQ. table%missingValue) THEN
            Missing = .TRUE.
            CYCLE
          END IF
        END DO
        numValid = numValid + 1
      END DO
    END DO
    PRINT *, 'Valid: ', numValid

    ! Pass on the data as output and the latlonhgttim:
    ALLOCATE (ncData(numValid, numVars), llht(numValid, 4))
    numValid = 0
    DO j = 1, UBOUND(temp, 2)
      DO i = 1, UBOUND(temp, 1)
        IF (latlon(1, i, j) .EQ. table%missingValue .OR. &
            latlon(2, i, j) .EQ. table%missingValue .OR. &
            latlon(1, i, j) .LT. domain(1, 1) .OR. &
            latlon(1, i, j) .GT. domain(2, 1) .OR. &
            latlon(2, i, j) .LT. domain(1, 2) .OR. &
            latlon(2, i, j) .GT. domain(2, 2)) THEN
          CYCLE
        END IF
        DO iv = 1, numVars
          IF (table%required(iv) .EQ. 1 .AND. &
              nc2Data(i, j, iv) .EQ. table%missingValue) THEN
            CYCLE
          END IF
        END DO
        numValid = numValid + 1
        ncData(numValid, :) = nc2Data(i, j, :)
        llht(numValid, 1:2) = latlon(1:2, i, j)
        llht(numValid, 3) = time(i, j)
        llht(numValid, 4) = height(i, j)
      END DO
    END DO
    DO iv = 1, numVars
      PRINT *, 'Maxmin obs: ', MAXVAL(ncData(:, iv)), MINVAL(ncData(:, iv)), table%varNames(iv), iv
    END DO
    PRINT *, 'Number of valid data: ', numValid
  END SUBROUTINE ncReader
END MODULE netCDFReader_m
