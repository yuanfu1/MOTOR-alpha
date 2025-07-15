MODULE ncWriteGrid_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE netcdf

CONTAINS
  SUBROUTINE writegrid(outfile, NX, NY, xpos, ypos, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outfile
    INTEGER(i_kind), INTENT(IN) :: NX, NY
    REAL(r_kind), DIMENSION(NX), INTENT(IN) :: xpos
    REAL(r_kind), DIMENSION(NY), INTENT(IN) :: ypos
    CHARACTER*(*), INTENT(IN) :: varName
    REAL(r_kind), DIMENSION(NX, NY), INTENT(IN) :: varData
    INTEGER(i_kind) :: ncid, x_dimid, y_dimid
    INTEGER(i_kind) :: x_varid, y_varid, varid
    INTEGER(i_kind), DIMENSION(2) :: dimids

!Create the netCDF file.
    CALL check(nf90_create(TRIM(outfile), NF90_CLOBBER, ncid))

!Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
!Define coordinate variables
    CALL check(nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, x_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, y_varid))
    dimids = (/x_dimid, y_dimid/)
!Define variable
    CALL check(nf90_def_var(ncid, varName, NF90_DOUBLE, dimids, varid))
    CALL check(nf90_enddef(ncid)) !End Definitions
!Write Data
    CALL check(nf90_put_var(ncid, x_varid, xpos))
    CALL check(nf90_put_var(ncid, y_varid, ypos))
    CALL check(nf90_put_var(ncid, varid, varData))
    CALL check(nf90_close(ncid))
  END SUBROUTINE writegrid

  SUBROUTINE writegrid3D(outfile, NX, NY, NZ, xpos, ypos, zpos, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outfile
    INTEGER(i_kind), INTENT(IN) :: NX, NY, NZ
    REAL(r_kind), DIMENSION(NX), INTENT(IN) :: xpos
    REAL(r_kind), DIMENSION(NY), INTENT(IN) :: ypos
    REAL(r_kind), DIMENSION(NZ), INTENT(IN) :: zpos
    CHARACTER(*), INTENT(IN) :: varName
    REAL(r_kind), DIMENSION(NX, NZ, NY), INTENT(IN) :: varData
    INTEGER(i_kind) :: ncid, x_dimid, y_dimid, z_dimid
    INTEGER(i_kind) :: x_varid, y_varid, z_varid, varid
    INTEGER(i_kind), DIMENSION(3) :: dimids

    !Create the netCDF file.
    CALL check(nf90_create(TRIM(outfile), NF90_CLOBBER, ncid))

    !Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
    CALL check(nf90_def_dim(ncid, "vertC", NZ, z_dimid))

    !Define coordinate variables
    CALL check(nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, x_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, y_varid))
    CALL check(nf90_def_var(ncid, "vertC", NF90_DOUBLE, z_dimid, z_varid))

    dimids = (/x_dimid, z_dimid, y_dimid/)
    !Define variable
    CALL check(nf90_def_var(ncid, varName, NF90_DOUBLE, dimids, varid))
    CALL check(nf90_enddef(ncid)) !End Definitions
    !Write Data

    CALL check(nf90_put_var(ncid, x_varid, xpos))
    CALL check(nf90_put_var(ncid, y_varid, ypos))
    CALL check(nf90_put_var(ncid, z_varid, zpos))
    CALL check(nf90_put_var(ncid, varid, varData))
    CALL check(nf90_close(ncid))
  END SUBROUTINE writegrid3D

  SUBROUTINE writegrid3D_XYZ(outfile, NX, NY, NZ, xpos, ypos, zpos, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outfile
    INTEGER(i_kind), INTENT(IN) :: NX, NY, NZ
    REAL(r_kind), DIMENSION(NX), INTENT(IN) :: xpos
    REAL(r_kind), DIMENSION(NY), INTENT(IN) :: ypos
    REAL(r_kind), DIMENSION(NZ), INTENT(IN) :: zpos
    CHARACTER(*), INTENT(IN) :: varName
    REAL(r_kind), DIMENSION(NX, NY, NZ), INTENT(IN) :: varData
    INTEGER(i_kind) :: ncid, x_dimid, y_dimid, z_dimid
    INTEGER(i_kind) :: x_varid, y_varid, z_varid, varid
    INTEGER(i_kind), DIMENSION(3) :: dimids

    !Create the netCDF file.
    CALL check(nf90_create(TRIM(outfile), NF90_CLOBBER, ncid))

    !Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
    CALL check(nf90_def_dim(ncid, "vertC", NZ, z_dimid))

    !Define coordinate variables
    CALL check(nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, x_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, y_varid))
    CALL check(nf90_def_var(ncid, "vertC", NF90_DOUBLE, z_dimid, z_varid))

    dimids = (/x_dimid, y_dimid, z_dimid/)
    !Define variable
    CALL check(nf90_def_var(ncid, varName, NF90_DOUBLE, dimids, varid))
    CALL check(nf90_enddef(ncid)) !End Definitions
    !Write Data

    CALL check(nf90_put_var(ncid, x_varid, xpos))
    CALL check(nf90_put_var(ncid, y_varid, ypos))
    CALL check(nf90_put_var(ncid, z_varid, zpos))
    CALL check(nf90_put_var(ncid, varid, varData))
    CALL check(nf90_close(ncid))
  END SUBROUTINE writegrid3D_XYZ

  SUBROUTINE writegrid3D_s(outfile, NX, NY, NZ, xpos, ypos, zpos, varName, varData)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: outfile
    INTEGER(i_kind), INTENT(IN) :: NX, NY, NZ
    REAL(r_single), DIMENSION(NX), INTENT(IN) :: xpos
    REAL(r_single), DIMENSION(NY), INTENT(IN) :: ypos
    REAL(r_single), DIMENSION(NZ), INTENT(IN) :: zpos
    CHARACTER(*), INTENT(IN) :: varName
    REAL(r_single), DIMENSION(NX, NZ, NY), INTENT(IN) :: varData
    INTEGER(i_kind) :: ncid, x_dimid, y_dimid, z_dimid
    INTEGER(i_kind) :: x_varid, y_varid, z_varid, varid
    INTEGER(i_kind), DIMENSION(3) :: dimids

    !Create the netCDF file.
    CALL check(nf90_create(TRIM(outfile), NF90_CLOBBER, ncid))

    !Define the dimensions.
    CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
    CALL check(nf90_def_dim(ncid, "sigma", NZ, z_dimid))

    !Define coordinate variables
    CALL check(nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, x_varid))
    CALL check(nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, y_varid))
    CALL check(nf90_def_var(ncid, "sigma", NF90_DOUBLE, z_dimid, z_varid))

    dimids = (/x_dimid, z_dimid, y_dimid/)
    !Define variable
    CALL check(nf90_def_var(ncid, varName, NF90_FLOAT, dimids, varid))
    CALL check(nf90_enddef(ncid)) !End Definitions
    !Write Data

    CALL check(nf90_put_var(ncid, x_varid, xpos))
    CALL check(nf90_put_var(ncid, y_varid, ypos))
    CALL check(nf90_put_var(ncid, z_varid, zpos))
    CALL check(nf90_put_var(ncid, varid, varData))
    CALL check(nf90_close(ncid))
  END SUBROUTINE writegrid3D_s

  SUBROUTINE check(status)
    INTEGER, INTENT(in) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check

END MODULE ncWriteGrid_m
