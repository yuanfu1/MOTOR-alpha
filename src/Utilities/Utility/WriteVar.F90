MODULE WriteVar_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE netcdf

CONTAINS
  SUBROUTINE writeVar2Dnc(outFile, NX, NY, LON2D, LAT2D, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outFile
    INTEGER(i_kind), INTENT(IN) :: NX, NY
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: LON2D, LAT2D
    REAL(r_kind), DIMENSION(:, :, :), INTENT(IN) :: varData
    CHARACTER*(*), INTENT(IN) :: varName

    INTEGER(i_kind) :: ncid, x_dimid, y_dimid, t_dimid
    INTEGER(i_kind) :: lon_varid, lat_varid, varid
    INTEGER(i_kind) :: dimidsvar(3), dimidsloc(2)

    ! Create the netCDF file.
    CALL check(nf90_create(TRIM(outFile), NF90_CLOBBER, ncid))

    ! Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
    CALL check(nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_dimid))
    dimidsvar = (/x_dimid, y_dimid, t_dimid/)
    dimidsloc = (/x_dimid, y_dimid/)

    ! Define variables
    CALL check(nf90_def_var(ncid, "lon", NF90_REAL, dimidsloc, lon_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_REAL, dimidsloc, lat_varid))
    CALL check(nf90_def_var(ncid, TRIM(varName), NF90_REAL, dimidsvar, varid))

    ! End define mode.
    CALL check(nf90_enddef(ncid))

    ! Write Data
    CALL check(nf90_put_var(ncid, lon_varid, LON2D))
    CALL check(nf90_put_var(ncid, lat_varid, LAT2D))
    CALL check(nf90_put_var(ncid, varid, varData))
    CALL check(nf90_close(ncid))

  END SUBROUTINE writeVar2Dnc

  SUBROUTINE writeVar3Dnc(outFile, NX, NY, NZ, LON2D, LAT2D, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outFile
    INTEGER(i_kind), INTENT(IN) :: NX, NY, NZ
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: LON2D, LAT2D
    REAL(r_kind), DIMENSION(:, :, :, :), INTENT(IN) :: varData
    CHARACTER*(*), INTENT(IN) :: varName

    INTEGER(i_kind) :: ncid, x_dimid, y_dimid, z_dimid, t_dimid
    INTEGER(i_kind) :: lon_varid, lat_varid, varid
    INTEGER(i_kind) :: dimidsvar(4), dimidsloc(2)

    ! Create the netCDF file.
    CALL check(nf90_create(TRIM(outFile), NF90_CLOBBER, ncid))

    ! Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
    CALL check(nf90_def_dim(ncid, "lev", NZ, z_dimid))
    CALL check(nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_dimid))
    dimidsvar = (/x_dimid, y_dimid, z_dimid, t_dimid/)
    dimidsloc = (/x_dimid, y_dimid/)

    ! Define variables
    CALL check(nf90_def_var(ncid, "lon", NF90_REAL, dimidsloc, lon_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_REAL, dimidsloc, lat_varid))
    CALL check(nf90_def_var(ncid, TRIM(varName), NF90_REAL, dimidsvar, varid))

    ! End define mode.
    CALL check(nf90_enddef(ncid))

    ! Write Data
    CALL check(nf90_put_var(ncid, lon_varid, LON2D))
    CALL check(nf90_put_var(ncid, lat_varid, LAT2D))
    CALL check(nf90_put_var(ncid, varid, varData))

    CALL check(nf90_close(ncid))

  END SUBROUTINE writeVar3Dnc

  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check

  SUBROUTINE writeVar1Dtxt(outFile, varData, outFMT, outStatus)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outFile
    REAL(r_kind), INTENT(IN) :: varData(:)
    CHARACTER(LEN=*), INTENT(IN) :: outFMT
    CHARACTER(LEN=*), INTENT(IN) :: outStatus

    INTEGER(i_kind) :: i

    OPEN (1, FILE=TRIM(outFile), STATUS=TRIM(outStatus))
    ! WRITE(UNIT=1, FMT='(I10,I10)') shape_var(1), shape_var(2)

    DO i = LBOUND(varData, 1), UBOUND(varData, 1)
      WRITE (UNIT=1, FMT=outFMT) varData(i)
    END DO

    CLOSE (1)
  END SUBROUTINE writeVar1Dtxt

  SUBROUTINE writeVar2Dtxt(outFile, varData, outFMT, outStatus)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outFile
    REAL(r_kind), INTENT(IN) :: varData(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: outFMT
    CHARACTER(LEN=*), INTENT(IN) :: outStatus

    INTEGER(i_kind) :: shape_var(2), i

    OPEN (1, FILE=TRIM(outFile), STATUS=TRIM(outStatus))
    ! WRITE(UNIT=1, FMT='(I10,I10)') shape_var(1), shape_var(2)

    ! shape_var = SHAPE(varData)
    ! DO i = 1, shape_var(1)
    DO i = LBOUND(varData, 1), UBOUND(varData, 1)
      WRITE (UNIT=1, FMT=outFMT) varData(i, 1), varData(i, 2)
    END DO

    CLOSE (1)
  END SUBROUTINE writeVar2Dtxt

  SUBROUTINE writeStr1Dtxt(outFile, varData, outFMT, outStatus)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outFile
    CHARACTER(*), INTENT(IN) :: varData(:)
    CHARACTER(LEN=*), INTENT(IN) :: outFMT
    CHARACTER(LEN=*), INTENT(IN) :: outStatus

    INTEGER(i_kind) :: i

    OPEN (1, FILE=TRIM(outFile), STATUS=TRIM(outStatus))
    ! WRITE(UNIT=1, FMT='(I10,I10)') shape_var(1), shape_var(2)

    DO i = LBOUND(varData, 1), UBOUND(varData, 1)
      WRITE (UNIT=1, FMT=outFMT) varData(i)
    END DO

    CLOSE (1)
  END SUBROUTINE writeStr1Dtxt

END MODULE WriteVar_m
