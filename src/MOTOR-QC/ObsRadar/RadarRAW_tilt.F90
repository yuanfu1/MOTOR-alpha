!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.RadarRAW
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/12/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
MODULE RadarRAW_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double, i_short, i_long, i_llong, i_byte
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE parameters_m, ONLY: degree2radian

  TYPE :: RadarRAW_t

    CHARACTER(LEN=1024) :: configFile
    REAL(r_kind), ALLOCATABLE :: latRef(:), &      !< in Radius
                                 lonRef(:), &      !< in Radius
                                 timRef(:), &
                                 altRef(:), &
                                 valRef(:)         !< dBZ

    REAL(r_kind), ALLOCATABLE :: latVel(:), &      !< in Radius
                                 lonVel(:), &      !< in Radius
                                 timVel(:), &
                                 altVel(:), &
                                 valVel(:)         !< meter/second, in radial direction

  CONTAINS
    FINAL :: destructor
    PROCEDURE :: readRadarData
    PROCEDURE :: outputNCForTest
    PROCEDURE, NOPASS :: countsTodBZ
    PROCEDURE, NOPASS :: countsToVel
    PROCEDURE, NOPASS :: extractVarsFromData
  END TYPE RadarRAW_t

  INTERFACE RadarRAW_t
    PROCEDURE :: constructor
  END INTERFACE RadarRAW_t

CONTAINS

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(RadarRAW_t) :: this
    CHARACTER(LEN=1024) :: configFile

    this%configFile = configFile
  END FUNCTION constructor

  SUBROUTINE readRadarData(this, filename)
    USE parameters_m, ONLY: degree2radian
    IMPLICIT NONE

    CLASS(RadarRAW_t) :: this
    CHARACTER(LEN=*) :: filename
    INTEGER(i_kind) :: i, j, k, l

    TYPE(NcDataset)   :: nc
    TYPE(NcVariable)  :: var
    TYPE(NcDimension)  :: dim
    TYPE(NcGroup)     :: grp

    REAL(r_kind), ALLOCATABLE :: radialAzim(:, :), &
                                 radialElev(:, :), &
                                 radialTime(:, :), &
                                 resolutionV(:)

    INTEGER(i_byte), ALLOCATABLE :: Vel(:, :, :), &
                                    Ref(:, :, :), &
                                    SpmW(:, :, :)

    REAL(r_kind), ALLOCATABLE :: distV(:, :), &
                                 distZ(:, :)

    INTEGER(i_kind) :: elevation, radial, V_bin, Z_bin
    REAL(r_kind) :: latSite, lonSite, altSite
    INTEGER(i_byte) :: V_fillValue, Z_fillValue

    ! Read data from nc files.
    nc = NcDataset(filename, "r")

    var = nc%getVariable("V"); CALL var%getData(Vel)
    CALL var%getAttribute('_FillValue', V_fillValue)

    var = nc%getVariable("Z"); CALL var%getData(Ref)
    CALL var%getAttribute('_FillValue', Z_fillValue)

    var = nc%getVariable("W"); CALL var%getData(SpmW)

    var = nc%getVariable("radialAzim"); CALL var%getData(radialAzim)
    var = nc%getVariable("radialElev"); CALL var%getData(radialElev)
    var = nc%getVariable("radialTime"); CALL var%getData(radialTime)
    var = nc%getVariable("resolutionV"); CALL var%getData(resolutionV)

    var = nc%getVariable("distanceV"); CALL var%getData(distV)
    distV = distV * 1000.
    var = nc%getVariable("distanceZ"); CALL var%getData(distZ)
    distZ = distZ * 1000.

    var = nc%getVariable("siteLat"); CALL var%getData(latSite)
    var = nc%getVariable("siteLon"); CALL var%getData(lonSite)
    var = nc%getVariable("siteAlt"); CALL var%getData(altSite)

    dim = nc%getDimension("elevation"); elevation = dim%getLength()
    dim = nc%getDimension("radial"); radial = dim%getLength()
    dim = nc%getDimension("V_bin"); V_bin = dim%getLength()
    dim = nc%getDimension("Z_bin"); Z_bin = dim%getLength()

    PRINT *, 'elevation: ', elevation
    PRINT *, 'radial: ', radial
    PRINT *, 'V_bin: ', V_bin
    PRINT *, 'Z_bin: ', Z_bin
    PRINT *, 'V_fillValue', V_fillValue
    PRINT *, 'Shape of Ref', SHAPE(Ref)
    PRINT *, 'Shape of radialAzim', SHAPE(radialAzim)
    PRINT *, 'Shape of radialElev', SHAPE(radialElev)
    PRINT *, 'Shape of distZ', SHAPE(distZ)
    PRINT *, 'siteLat', latSite
    PRINT *, 'siteLon', lonSite
    PRINT *, 'siteAlt', altSite
    PRINT *, 'Z: ', Ref(2, 1, :)

    CALL nc%CLOSE()
    ! End of the file reading.

    BLOCK
      REAL(r_kind), PARAMETER :: Z_offset = 66.0
      REAL(r_kind), PARAMETER :: Z_scale = 2.0
      REAL(r_kind), PARAMETER :: V_offset = 129.0
      REAL(r_kind), PARAMETER :: V_scale = 2.0

      CALL extractVarsFromData(radialAzim, radialElev, radialTime, distZ, resolutionV, Ref, &
                               latSite, lonSite, altSite, &
                               elevation, radial, Z_bin, &
                               this%latRef, this%lonRef, this%altRef, this%timRef, this%valRef, &
                               Z_offset, Z_scale, Z_fillValue, 'Ref')

      CALL extractVarsFromData(radialAzim, radialElev, radialTime, distV, resolutionV, Vel, &
                               latSite, lonSite, altSite, &
                               elevation, radial, V_bin, &
                               this%latVel, this%lonVel, this%altVel, this%timVel, this%valVel, &
                               V_offset, V_scale, V_fillValue, 'Vel')

    END BLOCK

    IF (ALLOCATED(radialAzim)) DEALLOCATE (radialAzim)
    IF (ALLOCATED(radialElev)) DEALLOCATE (radialElev)
    IF (ALLOCATED(radialTime)) DEALLOCATE (radialTime)
    IF (ALLOCATED(resolutionV)) DEALLOCATE (resolutionV)
    IF (ALLOCATED(Vel)) DEALLOCATE (Vel)
    IF (ALLOCATED(Ref)) DEALLOCATE (Ref)
    IF (ALLOCATED(SpmW)) DEALLOCATE (SpmW)
    IF (ALLOCATED(distV)) DEALLOCATE (distV)
    IF (ALLOCATED(distZ)) DEALLOCATE (distZ)

  END SUBROUTINE readRadarData

!> @brief
!! Extract data from raw format to lat/lon/alt/time/values
  SUBROUTINE extractVarsFromData(radialAzim, radialElev, radialTime, dist, resolutionV, values, &
                                 latSite, lonSite, altSite, &
                                 elevation, radial, R_bin, &
                                 latVar, lonVar, altVar, timVar, valVar, &
                                 offset, scale, fillvalue, varType)

    REAL(r_kind), INTENT(IN) :: offset, scale
    INTEGER(i_byte), INTENT(IN) :: fillvalue
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: radialAzim(:, :), &
                                             radialElev(:, :), &
                                             radialTime(:, :), &
                                             resolutionV(:), &
                                             dist(:, :)
    INTEGER(i_byte), ALLOCATABLE, INTENT(IN) :: values(:, :, :)

    REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
    INTEGER(i_kind), INTENT(IN) :: elevation, radial, R_bin
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: latVar(:), lonVar(:), altVar(:), timVar(:), &
                                                valVar(:)
    CHARACTER(*), INTENT(IN) :: varType

    INTEGER(i_long) :: i, j, k

    REAL(r_kind), ALLOCATABLE :: latSwap(:, :, :), &
                                 lonSwap(:, :, :), &
                                 altSwap(:, :, :), &
                                 timSwap(:, :, :), &
                                 ValSwap(:, :, :)
    INTEGER(i_long), ALLOCATABLE :: idxKJI(:, :)
    INTEGER(i_long) :: curIdx

    ALLOCATE (latSwap(R_bin, radial, elevation), &
              lonSwap(R_bin, radial, elevation), &
              altSwap(R_bin, radial, elevation), &
              timSwap(R_bin, radial, elevation), &
              ValSwap(R_bin, radial, elevation), &
              idxKJI(R_bin * radial * elevation, 3))

    curIdx = 1
    DO i = 1, elevation
      DO j = 1, radial
        DO k = 1, R_bin
          IF (TRIM(varType) .EQ. 'Ref') THEN ! For Reflectivity
            ValSwap(k, j, i) = countsTodBZ(values(k, j, i), fillValue, offset, scale)
          ELSE IF (TRIM(varType) .EQ. 'Vel') THEN ! For Velocity
            ValSwap(k, j, i) = countsToVel(values(k, j, i), fillValue, offset, scale, resolutionV(i))
          END IF

          IF (ValSwap(k, j, i) .NE. fillValue) THEN
            CALL calLatLon(radialAzim(j, i) * degree2radian, radialElev(j, i) * degree2radian, dist(k, i), &
                           latSite * degree2radian, lonSite * degree2radian, altSite, &
                           latSwap(k, j, i), lonSwap(k, j, i), altSwap(k, j, i))
            timSwap(k, j, i) = radialTime(j, i)
            idxKJI(curIdx, :) = [k, j, i]
            curIdx = curIdx + 1
          END IF
        END DO
      END DO
    END DO

    BLOCK
      INTEGER(i_kind) :: countPts
      countPts = COUNT(ValSwap .NE. fillValue)
      ! countPts = size(ValSwap)
      PRINT *, 'Count of points: ', countPts

      ALLOCATE (latVar(countPts), lonVar(countPts), &
                timVar(countPts), altVar(countPts), valVar(countPts))

      FORALL (i=1:countPts)
        latVar(i) = latSwap(idxKJI(i, 1), idxKJI(i, 2), idxKJI(i, 3))
        lonVar(i) = lonSwap(idxKJI(i, 1), idxKJI(i, 2), idxKJI(i, 3))
        altVar(i) = altSwap(idxKJI(i, 1), idxKJI(i, 2), idxKJI(i, 3))
        timVar(i) = timSwap(idxKJI(i, 1), idxKJI(i, 2), idxKJI(i, 3))
        valVar(i) = ValSwap(idxKJI(i, 1), idxKJI(i, 2), idxKJI(i, 3))

      END FORALL
    END BLOCK

    DEALLOCATE (latSwap, lonSwap, altSwap, timSwap, ValSwap)
  END SUBROUTINE extractVarsFromData

!> @brief
!! Convert from polar coordinates to lat/lon/alt
  SUBROUTINE calLatLon(azim, elev, rbin, latSite, lonSite, altSite, lat, lon, alt)
    USE parameters_m, ONLY: ae

    ! Azimuth clockwise
    REAL(r_kind), INTENT(IN) :: azim, elev, rbin
    REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
    REAL(r_kind), INTENT(OUT) :: lat, lon, alt

    REAL(r_kind) :: gcr
    REAL(r_kind) :: thetaLat
    REAL(r_kind) :: thetaLon
    REAL(r_kind) :: thetaAll
    REAL(r_kind) :: alphaLat, alphaLon, alphaAll
    REAL(r_kind) :: rLat
    REAL(r_kind) :: rLon

    ! Real great circle distance
    gcr = ae + altSite

    ! Theta at lat and lon direction
    thetaLat = ATAN((rbin * DSIN(elev)), (rbin * DCOS(elev) * DCOS(azim)))
    thetaLon = ATAN((rbin * DSIN(elev)), (rbin * DCOS(elev) * DSIN(azim)))
    ! PRINT *, 'thetaLat', thetaLat, elev, azim, rbin

    ! r distance at lat and lon direction
    rLat = (rbin * DSIN(elev)) / (DSIN(thetaLat))
    rLon = (rbin * DSIN(elev)) / (DSIN(thetaLon))

    !
    alphaLat = ATAN((rLat * DCOS(thetaLat)), (gcr + rLat * DSIN(thetaLat)))
    alphaLon = ATAN((rLon * DCOS(thetaLon)), (gcr + rLon * DSIN(thetaLon)))
    ! PRINT *, 'thetaLat', rbin

    !
    lat = latSite + alphaLat
    lon = lonSite + alphaLon

    !
    alphaAll = ATAN((rbin * DCOS(elev)), (gcr + rbin * DSIN(elev)))
    alt = (gcr / (DCOS(alphaAll))) + (rbin * DSIN(elev) / DCOS(alphaAll)) - gcr
  END SUBROUTINE calLatLon

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(RadarRAW_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%latRef)) DEALLOCATE (this%latRef)
    IF (ALLOCATED(this%lonRef)) DEALLOCATE (this%lonRef)
    IF (ALLOCATED(this%timRef)) DEALLOCATE (this%timRef)
    IF (ALLOCATED(this%AltRef)) DEALLOCATE (this%AltRef)
    IF (ALLOCATED(this%ValRef)) DEALLOCATE (this%ValRef)

    IF (ALLOCATED(this%latVel)) DEALLOCATE (this%latVel)
    IF (ALLOCATED(this%lonVel)) DEALLOCATE (this%lonVel)
    IF (ALLOCATED(this%timVel)) DEALLOCATE (this%timVel)
    IF (ALLOCATED(this%AltVel)) DEALLOCATE (this%AltVel)
    IF (ALLOCATED(this%ValVel)) DEALLOCATE (this%ValVel)
  END SUBROUTINE destructor

!> @brief
!! Convert integer Z count value to dbz
  PURE FUNCTION countsTodBZ(zcounts, fillValue, Z_offset, Z_scale) RESULT(dBZ)
    INTEGER(i_byte), INTENT(in) :: zcounts, fillValue
    REAL(r_kind), INTENT(in) :: Z_offset, Z_scale
    REAL(r_kind) :: dBZ
!      From the NetCDF header
!      Z:valid_range = 2b, -2b ;     (2 through 254)
!      Z:below_threshold = 0b ;      (0)
!      Z:range_ambiguous = 1b ;      (1)
!      Z:_FillValue = -1b ;          (255)

    dBZ = zcounts

    IF (zcounts .EQ. fillValue) THEN
      dBZ = fillValue
      RETURN
    END IF

    ! IF (dBZ .gt. 127.) THEN
    !   PRINT *, 'error in Reflectivity: ', dbz_hold, zcounts
    !   dBZ = fillValue
    ! END IF
    ! Convert from signed to unsigned
    IF (dBZ .LT. 0.) dBZ = 256 + dBZ

    IF (dBZ .EQ. 1. .OR. dbZ .EQ. 0.) THEN
      dBZ = fillValue  ! Range Ambiguous
      RETURN
    END IF

    dBZ = (dBZ - Z_offset) / Z_scale
  END FUNCTION countsTodBZ

!> @brief
!! Convert integer V count value to radial velocity (meters/sec)
  PURE FUNCTION countsToVel(vcounts, fillValue, V_offset, V_scale, V_resolution) RESULT(Vel)
    INTEGER(i_byte), INTENT(in) :: vcounts, fillValue
    REAL(r_kind), INTENT(in) :: V_offset, V_scale, V_resolution
    REAL(r_kind) :: Vel

    Vel = vcounts

    IF (vcounts .EQ. fillValue) THEN
      Vel = fillValue
      RETURN
    END IF

    ! IF (Vel .gt. 127.) THEN
    !   PRINT *, 'error in Velocity: ', Vel
    !   Vel = 0.
    ! ENDIF

    !  Convert from signed to unsigned
    IF (Vel .LT. 0.) Vel = 256.+Vel

    IF (Vel .EQ. 1. .OR. Vel .EQ. 0.) THEN
      Vel = fillValue  ! Invalid Measurement
      RETURN
    END IF

    IF (V_resolution .EQ. 0.) THEN ! QC Check
      Vel = fillValue
      RETURN
    END IF

    Vel = (Vel - V_offset) / V_scale
  END FUNCTION

  SUBROUTINE outputNCForTest(this, filename)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    IMPLICIT NONE
    CLASS(RadarRAW_t) :: this
    CHARACTER(*) :: filename

    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1
    TYPE(NcVariable)  :: var, scalar
    TYPE(NcGroup)     :: grp

    INTEGER(i_kind) :: recNum
! open a dataset in write mode
! args:
!     filename
!     mode ("w": write, "r": read-only, "a": read-write)
    nc = NcDataset(filename, "w")

    recNum = SIZE(this%latRef)
    dim1 = nc%setDimension("recNum", recNum)

    var = nc%setVariable("latRef", "f64", (/dim1/)); CALL var%setData(this%latRef / degree2radian)
    var = nc%setVariable("lonRef", "f64", (/dim1/)); CALL var%setData(this%lonRef / degree2radian)
    var = nc%setVariable("timRef", "f64", (/dim1/)); CALL var%setData(this%timRef)
    var = nc%setVariable("altRef", "f64", (/dim1/)); CALL var%setData(this%altRef)
    var = nc%setVariable("valRef", "f64", (/dim1/)); CALL var%setData(this%valRef)

    ! close the dataset
    CALL nc%CLOSE()
  END SUBROUTINE

END MODULE RadarRAW_m
