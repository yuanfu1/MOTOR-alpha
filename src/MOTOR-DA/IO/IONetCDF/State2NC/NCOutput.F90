!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE NCOutput_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup

  TYPE :: NCOutput_t
    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dimx, dimy, dimz, dimt
    TYPE(NcVariable)  :: varx, vary, varz, vart

  CONTAINS
    PROCEDURE :: CLOSE
    PROCEDURE :: addVar

  END TYPE NCOutput_t

  INTERFACE NCOutput_t
    PROCEDURE :: constructor
  END INTERFACE NCOutput_t

CONTAINS

  FUNCTION constructor(outfile, xpos, ypos, zpos, tpos, NX, NY, NZ, NT) RESULT(this)
    TYPE(NCOutput_t) :: this
    INTEGER(i_kind) :: NX, NY, NZ, NT
    REAL(r_kind) :: xpos(NX), ypos(NY), zpos(NZ), tpos(NT)
    CHARACTER(*), INTENT(IN) :: outfile

    ASSOCIATE (nc => this%nc, dimx => this%dimx, &
               dimy => this%dimy, dimz => this%dimz, dimt => this%dimt, &
               varx => this%varx, vary => this%vary, varz => this%varz, vart => this%vart)

      nc = NcDataset(TRIM(outfile), "w")

      dimx = nc%setDimension("lon", NX)
      dimy = nc%setDimension("lat", NY)
      dimz = nc%setDimension("sigma", NZ)
      dimt = nc%setDimension("time")

      varx = nc%setVariable("lon", "f32", (/dimx/)); CALL varx%setData(xpos)
      vary = nc%setVariable("lat", "f32", (/dimy/)); CALL vary%setData(ypos)
      varz = nc%setVariable("sigma", "f32", (/dimz/)); CALL varz%setData(zpos)
      vart = nc%setVariable("time", "f32", (/dimt/)); CALL vart%setData(tpos - tpos(1))

      CALL nc%setAttribute("TITLE", "OUTPUT FROM MOTOR V0.1 DA")

      CALL varx%setAttribute("description", "Longitude, west is negative")
      CALL varx%setAttribute("units", "degree_east")
      CALL varx%setAttribute("long_name", "longitude")
      CALL varx%setAttribute("axis", "X")
      CALL varx%setAttribute("_CoordinateAxisType", "Lon")
      CALL varx%setAttribute("_CoordinateAliasForDimension", "west_east")

      CALL vary%setAttribute("description", "Latitude, south is negative")
      CALL vary%setAttribute("units", "degree_north")
      CALL vary%setAttribute("long_name", "latitude")
      CALL vary%setAttribute("axis", "Y")
      CALL vary%setAttribute("_CoordinateAxisType", "Lat")
      CALL vary%setAttribute("_CoordinateAliasForDimension", "south_north")

      CALL varz%setAttribute("description", "Sigma coordinate, from 0 to ztop")
      CALL varz%setAttribute("long_name", "Sigma")
      CALL varz%setAttribute("units", "meters")
      CALL varz%setAttribute("axis", "Z")
      CALL varz%setAttribute("_CoordinateAxisType", "Z")
      CALL varz%setAttribute("_CoordinateAliasForDimension", "bottom_top")

      BLOCK
        USE AdvanceTime_m
        CHARACTER(len=19) :: timeStr(1)
        INTEGER(i_kind) :: gTime(6), uTime_sec

1010    FORMAT(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2)
        uTime_sec = tpos(1)
        CALL Time_Unix_to_GMT(uTime_sec, gTime)
        WRITE (timeStr(1), 1010) gTime(1), '-', gTime(2), '-', gTime(3), ' ', gTime(4), ':', gTime(5), ':', gTime(6)
        CALL vart%setAttribute("description", "seconds since "//timeStr(1))
        CALL vart%setAttribute("units", "seconds since "//timeStr(1))
        CALL vart%setAttribute("axis", "T")
        CALL vart%setAttribute("_CoordinateAxisType", "Time")
        CALL vart%setAttribute("_CoordinateAliasForDimension", "time")
      END BLOCK
    END ASSOCIATE
  END FUNCTION

  SUBROUTINE addVar(this, varName, varData, nz, nx, ny, nt)
    CLASS(NCOutput_t) :: this
    CHARACTER(*), INTENT(IN) :: varName
    REAL(r_kind) :: varData(nz, nx, ny, nt)
    TYPE(NcVariable)  :: var
    REAL(r_kind), ALLOCATABLE :: valueState(:, :, :, :)
    INTEGER(i_kind) :: i

    ALLOCATE (valueState(nx, ny, nz, nt))
    FORALL (i=1:nz) valueState(:, :, i, :) = varData(i, :, :, :)

    var = this%nc%setVariable(varName, "f64", &
                              (/this%dimx, this%dimy, this%dimz, this%dimt/)); CALL var%setData(valueState)
    CALL var%setAttribute("_CoordinateAxes", "lon lat sigma time")

    DEALLOCATE (valueState)
  END SUBROUTINE addVar

  SUBROUTINE CLOSE (this)
    CLASS(NCOutput_t) :: this

    CALL this%nc%CLOSE()
  END SUBROUTINE

END MODULE NCOutput_m
