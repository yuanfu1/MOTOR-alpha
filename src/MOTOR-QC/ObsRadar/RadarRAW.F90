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
  IMPLICIT NONE

  TYPE :: RadarRAW_t

    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind) :: sizeRef, sizeVel
    REAL(r_kind), ALLOCATABLE :: latRef(:), &      !< in Radius
                                 lonRef(:), &      !< in Radius
                                 timRef(:), &
                                 altRef(:), &
                                 valRef(:)         !< dBZ

    REAL(r_kind), ALLOCATABLE :: latVel(:), &      !< in Radius
                                 lonVel(:), &      !< in Radius
                                 timVel(:), &
                                 altVel(:), &
                                 valVel(:), &      !< meter/second, in radial direction
                                 valSpmW(:)

    REAL(r_kind) :: locRadar(3)
    CHARACTER(len=20) :: station
    REAL(r_kind), ALLOCATABLE :: unAmbigRangeRef(:)
    REAL(r_kind), ALLOCATABLE :: unAmbigRangeVel(:)

  CONTAINS
    FINAL :: destructor
    PROCEDURE :: readRadarData
    PROCEDURE :: read_SAD_Data
    PROCEDURE :: read_DXK_Data
    PROCEDURE :: outputNCForTest
    PROCEDURE, NOPASS :: countsTodBZ, countsTodBZ_DXK
    PROCEDURE, NOPASS :: countsToVel, countsToVel_DXK
    PROCEDURE, NOPASS :: countsToSpmW
    PROCEDURE, NOPASS :: extractVarsFromData, extractVarsFromDataForVel
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

  FUNCTION getStaName(filename) RESULT(staName)
    CHARACTER(LEN=*) :: filename
    CHARACTER(LEN=5) :: staName
    TYPE(NcDataset)   :: nc

    nc = NcDataset(filename, "r")
    CALL nc%getAttribute("Station", staName)
    CALL nc%CLOSE()

  END FUNCTION

  SUBROUTINE readRadarData(this, filename, varname)
    USE parameters_m, ONLY: degree2radian
    IMPLICIT NONE
    CLASS(RadarRAW_t) :: this
    CHARACTER(LEN=*) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: varname
    TYPE(NcDataset) :: nc
    CHARACTER(len=20) :: summary

    nc = NcDataset(filename, "r")
    CALL nc%getAttribute("Summary", summary); 
    PRINT *, 'summary: ', summary(13:15)
    PRINT *, 'The radar type is ', summary(13:15), '.'

    SELECT CASE (summary(13:15))
    CASE ("SAD")
      PRINT *, 'Going here.'
      CALL this%read_SAD_Data(varname, nc)
    CASE ("DXK")
      CALL this%read_DXK_Data(varname, nc)
    CASE DEFAULT
      CALL this%read_SAD_Data(varname, nc)
    END SELECT

    CALL nc%CLOSE()

  END SUBROUTINE

  SUBROUTINE read_DXK_Data(this, varname, nc)
    CLASS(RadarRAW_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: varname
    TYPE(NcDataset), INTENT(INOUT) :: nc

    INTEGER(i_kind) :: i, j, k, l, posixOnDate
    TYPE(NcVariable)  :: var
    TYPE(NcDimension)  :: dim
    TYPE(NcGroup)     :: grp

    REAL(r_kind), ALLOCATABLE :: radialAzimR(:, :), &
                                 radialElevR(:, :), &
                                 radialTimeR(:, :), &
                                 radialAzimV(:, :), &
                                 radialElevV(:, :), &
                                 radialTimeV(:, :)
    !  resolutionV(:)

    INTEGER(i_short), ALLOCATABLE :: Vel(:, :, :), &
                                     Ref(:, :, :), &
                                     SpmW(:, :, :)

    REAL(r_kind), ALLOCATABLE :: distV(:), &
                                 distR(:)

    INTEGER(i_kind) :: elevNumV, elevNumR, &
                       radialNumV, radialNumR, &
                       gateNumV, gateNumR
    REAL(r_kind) :: latSite, lonSite, altSite
    INTEGER(i_byte), ALLOCATABLE :: V_fillValue(:), R_fillValue(:)

    REAL(r_kind) :: Z_offset
    REAL(r_kind) :: Z_scale
    REAL(r_kind) :: V_offset
    REAL(r_kind) :: V_scale

    ! PRINT *, 'V_scale: ', V_scale
    ! PRINT *, 'V_offset: ', V_offset
    ! PRINT *, 'Z_scale: ', Z_scale
    ! PRINT *, 'Z_offset: ', Z_offset

    ! Some QC data has not SpectrumWidth variable
    ! var = nc%getVariable("SpectrumWidth"); CALL var%getData(SpmW)
    ! var = nc%getVariable("resolutionV"); CALL var%getData(resolutionV)

    CALL nc%getAttribute("Station", this%Station); 
    CALL nc%getAttribute("StationLatitude", latSite); 
    CALL nc%getAttribute("StationLongitude", lonSite); 
    CALL nc%getAttribute("StationElevationInMeters", altSite); 
    this%locRadar = (/latSite * degree2radian, lonSite * degree2radian, altSite/)
    PRINT *, 'siteLat', latSite
    PRINT *, 'siteLon', lonSite
    PRINT *, 'siteAlt', altSite
    PRINT *, 'Station: ', this%Station

    ! End of the file reading.
    BLOCK
      USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
      CHARACTER(len=50) :: dateStr
      INTEGER(i_kind) :: YYYY, MM, DD

      CALL nc%getAttribute("base_date", dateStr); 
      PRINT *, 'dateStr: ', dateStr
      READ (dateStr, "(I4,1x,I2,1x,I2)") YYYY, MM, DD
      PRINT *, 'date: ', YYYY, MM, DD
      CALL Time_GMT_to_Unix((/YYYY, MM, DD, 0, 0, 0/), posixOnDate)
      PRINT *, 'posixOnDate: ', posixOnDate
    END BLOCK

    IF (TRIM(varname) == 'ref') THEN

      var = nc%getVariable("Reflectivity"); CALL var%getData(Ref)
      ! CALL var%getAttribute('missing_value', R_fillValue)
      CALL var%getAttribute("scale_factor", Z_scale); 
      CALL var%getAttribute("add_offset", Z_offset); 
      ! Z_scale = 1/Z_scale
      ! Z_offset = Z_offset*2*(-1)
      ALLOCATE (R_fillValue(1))
      R_fillValue = (/0/)

      var = nc%getVariable("azimuth"); CALL var%getData(radialAzimR)
      var = nc%getVariable("elevation"); CALL var%getData(radialElevR)
      var = nc%getVariable("time"); CALL var%getData(radialTimeR)
      radialTimeR = radialTimeR / 1000.0 + posixOnDate

      var = nc%getVariable("unambiguousRange"); CALL var%getData(this%unAmbigRangeRef)
      this%unAmbigRangeRef = this%unAmbigRangeRef * 1000.0D0

      dim = nc%getDimension("scan"); elevNumR = dim%getLength()
      dim = nc%getDimension("radial"); radialNumR = dim%getLength()
      dim = nc%getDimension("gate"); gateNumR = dim%getLength()

      var = nc%getVariable("distance"); CALL var%getData(distR)

      PRINT *, 'elevNumR: ', elevNumR
      PRINT *, 'radialNumR: ', radialNumR
      PRINT *, 'gateNumR: ', gateNumR
      PRINT *, 'R_fillValue', R_fillValue
      PRINT *, 'radialTimeR: ', radialTimeR(1, 1)
      PRINT *, 'Z_offset: ', Z_offset
      PRINT *, 'Z_scale: ', Z_scale

      CALL extractVarsFromData(radialAzimR, radialElevR, radialTimeR, distR, Ref, &
                               latSite, lonSite, altSite, &
                               elevNumR, radialNumR, gateNumR, &
                               this%latRef, this%lonRef, this%altRef, this%timRef, this%valRef, &
                               Z_offset, Z_scale, R_fillValue, 'Ref', this%unAmbigRangeRef, "DXK")
      this%sizeRef = SIZE(this%latRef)

    ELSE IF (TRIM(varname) == 'vel') THEN
      var = nc%getVariable("RadialVelocity"); CALL var%getData(Vel)
      ! CALL var%getAttribute('missing_value', V_fillValue)
      CALL var%getAttribute("scale_factor", V_scale); 
      CALL var%getAttribute("add_offset", V_offset); 
      ! V_scale = 1/V_scale
      ! V_offset = V_offset*2*(-1)
      ALLOCATE (V_fillValue(5))
      V_fillValue = (/0, 1, 2, 3, 4/)

      var = nc%getVariable("azimuth"); CALL var%getData(radialAzimV)
      var = nc%getVariable("elevation"); CALL var%getData(radialElevV)
      var = nc%getVariable("time"); CALL var%getData(radialTimeV)
      var = nc%getVariable("unambiguousRange"); CALL var%getData(this%unAmbigRangeVel)
      this%unAmbigRangeVel = this%unAmbigRangeVel * 1000.0D0
      var = nc%getVariable("distance"); CALL var%getData(distV)
      radialTimeV = radialTimeV / 1000.0 + posixOnDate

      dim = nc%getDimension("scan"); elevNumV = dim%getLength()
      dim = nc%getDimension("radial"); radialNumV = dim%getLength()
      dim = nc%getDimension("gate"); gateNumV = dim%getLength()
      PRINT *, 'elevNumV: ', elevNumV
      PRINT *, 'radialNumV: ', radialNumV
      PRINT *, 'gateNumV: ', gateNumV
      PRINT *, 'V_fillValue', V_fillValue
      PRINT *, 'radialTimeV: ', radialTimeV(1, 1)
      PRINT *, 'V_offset: ', V_offset
      PRINT *, 'V_scale: ', V_scale

      CALL extractVarsFromData(radialAzimV, radialElevV, radialTimeV, distV, Vel, &
                               latSite, lonSite, altSite, &
                               elevNumV, radialNumV, gateNumV, &
                               this%latVel, this%lonVel, this%altVel, this%timVel, this%valVel, &
                               V_offset, V_scale, V_fillValue, 'Vel', this%unAmbigRangeVel, "DXK")
      this%sizeVel = SIZE(this%latVel)
    END IF

    IF (ALLOCATED(radialAzimR)) DEALLOCATE (radialAzimR)
    IF (ALLOCATED(radialElevR)) DEALLOCATE (radialElevR)
    IF (ALLOCATED(radialTimeR)) DEALLOCATE (radialTimeR)

    IF (ALLOCATED(radialAzimV)) DEALLOCATE (radialAzimV)
    IF (ALLOCATED(radialElevV)) DEALLOCATE (radialElevV)
    IF (ALLOCATED(radialTimeV)) DEALLOCATE (radialTimeV)

    ! IF (ALLOCATED(resolutionV)) DEALLOCATE (resolutionV)
    IF (ALLOCATED(Vel)) DEALLOCATE (Vel)
    IF (ALLOCATED(Ref)) DEALLOCATE (Ref)
    IF (ALLOCATED(SpmW)) DEALLOCATE (SpmW)
    IF (ALLOCATED(distV)) DEALLOCATE (distV)
    IF (ALLOCATED(distR)) DEALLOCATE (distR)
  END SUBROUTINE

  SUBROUTINE read_SAD_Data(this, varname, nc)
    CLASS(RadarRAW_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: varname
    TYPE(NcDataset), INTENT(INOUT) :: nc
    INTEGER(i_kind) :: i, j, k, l, posixOnDate
    TYPE(NcVariable)  :: var
    TYPE(NcDimension)  :: dim
    TYPE(NcGroup)     :: grp

    REAL(r_kind), ALLOCATABLE :: radialAzimR(:, :), &
                                 radialElevR(:, :), &
                                 radialTimeR(:, :), &
                                 radialAzimV(:, :), &
                                 radialElevV(:, :), &
                                 radialTimeV(:, :)
    !  resolutionV(:)

    INTEGER(i_short), ALLOCATABLE :: Vel(:, :, :), &
                                     Ref(:, :, :), &
                                     SpmW(:, :, :)

    REAL(r_kind), ALLOCATABLE :: distV(:), &
                                 distR(:)

    INTEGER(i_kind) :: elevNumV, elevNumR, &
                       radialNumV, radialNumR, &
                       gateNumV, gateNumR
    REAL(r_kind) :: latSite, lonSite, altSite
    INTEGER(i_byte), ALLOCATABLE :: V_fillValue(:), R_fillValue(:)

    REAL(r_kind) :: Z_offset
    REAL(r_kind) :: Z_scale
    REAL(r_kind) :: V_offset
    REAL(r_kind) :: V_scale

    ! PRINT *, 'V_scale: ', V_scale
    ! PRINT *, 'V_offset: ', V_offset
    ! PRINT *, 'Z_scale: ', Z_scale
    ! PRINT *, 'Z_offset: ', Z_offset

    ! Some QC data has not SpectrumWidth variable
    ! var = nc%getVariable("SpectrumWidth"); CALL var%getData(SpmW)
    ! var = nc%getVariable("resolutionV"); CALL var%getData(resolutionV)

    CALL nc%getAttribute("Station", this%Station); 
    CALL nc%getAttribute("StationLatitude", latSite); 
    CALL nc%getAttribute("StationLongitude", lonSite); 
    CALL nc%getAttribute("StationElevationInMeters", altSite); 
    this%locRadar = (/latSite * degree2radian, lonSite * degree2radian, altSite/)
    PRINT *, 'siteLat', latSite
    PRINT *, 'siteLon', lonSite
    PRINT *, 'siteAlt', altSite

    ! End of the file reading.
    BLOCK
      USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
      CHARACTER(len=10) :: dateStr
      INTEGER(i_kind) :: YYYY, MM, DD

      CALL nc%getAttribute("base_date", dateStr); 
      PRINT *, 'dateStr: ', dateStr
      READ (dateStr, "(I4,1x,I2,1x,I2)") YYYY, MM, DD
      PRINT *, 'date: ', YYYY, MM, DD
      CALL Time_GMT_to_Unix((/YYYY, MM, DD, 0, 0, 0/), posixOnDate)
      PRINT *, 'posixOnDate: ', posixOnDate
    END BLOCK

    IF (TRIM(varname) == 'ref') THEN

      var = nc%getVariable("Reflectivity"); CALL var%getData(Ref)
      ! CALL var%getAttribute('missing_value', R_fillValue)
      CALL var%getAttribute("scale_factor", Z_scale); 
      CALL var%getAttribute("add_offset", Z_offset); 
      ALLOCATE (R_fillValue(5))
      R_fillValue = (/0, 1, 2, 3, 4/)

      var = nc%getVariable("azimuthR"); CALL var%getData(radialAzimR)
      var = nc%getVariable("elevationR"); CALL var%getData(radialElevR)
      var = nc%getVariable("timeR"); CALL var%getData(radialTimeR)
      var = nc%getVariable("unambiguousRangeR"); CALL var%getData(this%unAmbigRangeRef)
      this%unAmbigRangeRef = this%unAmbigRangeRef * 1000.0D0
      var = nc%getVariable("distanceR"); CALL var%getData(distR)
      radialTimeR = radialTimeR / 1000.0 + posixOnDate

      dim = nc%getDimension("scanR"); elevNumR = dim%getLength()
      dim = nc%getDimension("radialR"); radialNumR = dim%getLength()
      dim = nc%getDimension("gateR"); gateNumR = dim%getLength()
      PRINT *, 'elevNumR: ', elevNumR
      PRINT *, 'radialNumR: ', radialNumR
      PRINT *, 'gateNumR: ', gateNumR
      PRINT *, 'R_fillValue', R_fillValue
      PRINT *, 'radialTimeR: ', radialTimeR(1, 1)

      CALL extractVarsFromData(radialAzimR, radialElevR, radialTimeR, distR, Ref, &
                               latSite, lonSite, altSite, &
                               elevNumR, radialNumR, gateNumR, &
                               this%latRef, this%lonRef, this%altRef, this%timRef, this%valRef, &
                               Z_offset, Z_scale, R_fillValue, 'Ref', this%unAmbigRangeRef, "SAD")
      this%sizeRef = SIZE(this%latRef)

    ELSE IF (TRIM(varname) == 'vel') THEN
      var = nc%getVariable("RadialVelocity"); CALL var%getData(Vel)
      var = nc%getVariable("SpectrumWidth"); CALL var%getData(SpmW)
      ! CALL var%getAttribute('missing_value', V_fillValue)
      CALL var%getAttribute("scale_factor", V_scale); 
      CALL var%getAttribute("add_offset", V_offset); 
      ALLOCATE (V_fillValue(5))
      V_fillValue = (/0, 1, 2, 3, 4/)

      var = nc%getVariable("azimuthV"); CALL var%getData(radialAzimV)
      var = nc%getVariable("elevationV"); CALL var%getData(radialElevV)
      var = nc%getVariable("timeV"); CALL var%getData(radialTimeV)
      ! PRINT *, 'Reading velocity data. 1',
      BLOCK
        REAL(r_kind) :: ccRV
        IF (UBOUND(radialTimeV, 2) == 1) THEN
          var = nc%getVariable("unambiguousRangeV"); CALL var%getData(ccRV)
          ALLOCATE (this%unAmbigRangeVel(1))
          this%unAmbigRangeVel(1) = ccRV
        ELSE
          var = nc%getVariable("unambiguousRangeV"); CALL var%getData(this%unAmbigRangeVel)
        END IF

      END BLOCK
      ! PRINT *, 'Reading velocity data. 2'

      this%unAmbigRangeVel = this%unAmbigRangeVel * 1000.0D0
      var = nc%getVariable("distanceV"); CALL var%getData(distV)
      radialTimeV = radialTimeV / 1000.0 + posixOnDate

      dim = nc%getDimension("scanV"); elevNumV = dim%getLength()
      dim = nc%getDimension("radialV"); radialNumV = dim%getLength()
      dim = nc%getDimension("gateV"); gateNumV = dim%getLength()
      PRINT *, 'elevNumV: ', elevNumV
      PRINT *, 'radialNumV: ', radialNumV
      PRINT *, 'gateNumV: ', gateNumV
      PRINT *, 'V_fillValue', V_fillValue
      PRINT *, 'radialTimeV: ', radialTimeV(1, 1)

      CALL extractVarsFromDataForVel(radialAzimV, radialElevV, radialTimeV, distV, Vel, SpmW, &
                                     latSite, lonSite, altSite, &
                                     elevNumV, radialNumV, gateNumV, &
                                     this%latVel, this%lonVel, this%altVel, this%timVel, this%valVel, &
                                     V_offset, V_scale, V_fillValue, 'Vel', this%unAmbigRangeVel, "SAD")
      this%sizeVel = SIZE(this%latVel)
      ! PRINT *, 'sizeVel: ', this%sizeVel, this%unAmbigRangeVel
      ! STOP
    END IF

    IF (ALLOCATED(radialAzimR)) DEALLOCATE (radialAzimR)
    IF (ALLOCATED(radialElevR)) DEALLOCATE (radialElevR)
    IF (ALLOCATED(radialTimeR)) DEALLOCATE (radialTimeR)

    IF (ALLOCATED(radialAzimV)) DEALLOCATE (radialAzimV)
    IF (ALLOCATED(radialElevV)) DEALLOCATE (radialElevV)
    IF (ALLOCATED(radialTimeV)) DEALLOCATE (radialTimeV)

    ! IF (ALLOCATED(resolutionV)) DEALLOCATE (resolutionV)
    IF (ALLOCATED(Vel)) DEALLOCATE (Vel)
    IF (ALLOCATED(Ref)) DEALLOCATE (Ref)
    IF (ALLOCATED(SpmW)) DEALLOCATE (SpmW)
    IF (ALLOCATED(distV)) DEALLOCATE (distV)
    IF (ALLOCATED(distR)) DEALLOCATE (distR)

  END SUBROUTINE

!> @brief
!! Extract data from raw format to lat/lon/alt/time/values
  SUBROUTINE extractVarsFromDataForVel(radialAzim, radialElev, radialTime, dist, velRaw, SpmWRaw, &
                                       latSite, lonSite, altSite, &
                                       elevNum, radialNum, gateNum, &
                                       latVar, lonVar, altVar, timVar, valVar, &
                                       offset, scale, fillvalue, varType, unAmbigRange, radarType)

    REAL(r_kind), INTENT(IN) :: offset, scale
    INTEGER(i_byte), INTENT(IN) :: fillvalue(:)
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: radialAzim(:, :), &
                                             radialElev(:, :), &
                                             radialTime(:, :), &
                                             !  resolutionV(:), &
                                             dist(:)
    INTEGER(i_short), ALLOCATABLE, INTENT(IN) :: velRaw(:, :, :), SpmWRaw(:, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: unAmbigRange(:)

    REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
    INTEGER(i_kind), INTENT(IN) :: elevNum, radialNum, gateNum
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: latVar(:), lonVar(:), altVar(:), timVar(:), &
                                                valVar(:)
    CHARACTER(*), INTENT(IN) :: varType
    CHARACTER(*), INTENT(IN) :: radarType

    INTEGER(i_long) :: i, j, k

    REAL(r_kind), ALLOCATABLE :: latSwap(:, :, :), &
                                 lonSwap(:, :, :), &
                                 altSwap(:, :, :), &
                                 timSwap(:, :, :), &
                                 ValSwap(:, :, :)
    REAL(r_kind), ALLOCATABLE :: val(:, :, :), SpmW(:, :, :)

    ALLOCATE (val(gateNum, radialNum, elevNum))
    ALLOCATE (SpmW(gateNum, radialNum, elevNum))

    ! print *, 'velRaw', velRaw
    val = 0.0D0
    SpmW = 0.0D0
    IF (radarType == "SAD") THEN
      IF (TRIM(varType) .EQ. 'Ref') THEN ! For Reflectivity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsTodBZ(velRaw(k, j, i), offset, scale)
      ELSE IF (TRIM(varType) .EQ. 'Vel') THEN ! For Velocity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i))
          val(k, j, i) = countsToVel(velRaw(k, j, i), offset, scale)
          SpmW(k, j, i) = countsToSpmW(SpmWRaw(k, j, i), offset, scale)
        END FORALL
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i) &
                .AND. SpmW(k, j, i) > 3.5D0)
          val(k, j, i) = 0.0D0
        END FORALL
        ! print*, SpmW
      END IF
    ELSEIF (radarType == "DXK") THEN
      IF (TRIM(varType) .EQ. 'Ref') THEN ! For Reflectivity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsTodBZ_DXK(velRaw(k, j, i), offset, scale)
      ELSE IF (TRIM(varType) .EQ. 'Vel') THEN ! For Velocity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsToVel_DXK(velRaw(k, j, i), offset, scale)
      END IF
    END IF

    BLOCK ! Interpolation
      REAL(r_kind), ALLOCATABLE :: radialAzimInterp(:, :), &
                                   radialElevInterp(:, :), &
                                   radialTimeInterp(:, :), &
                                   !  resolutionV(:), &
                                   distInterp(:)
      ! INTEGER(i_byte), ALLOCATABLE :: velRawInterp(:, :, :)
      INTEGER(i_kind) :: elevNumInterp, radialNumInterp, gateNumInterp

      INTEGER :: elevScaleFactor

      elevScaleFactor = 0

      gateNumInterp = gateNum
      radialNumInterp = radialNum
      elevNumInterp = elevNum + (elevNum - 1) * elevScaleFactor

      ALLOCATE (radialAzimInterp(radialNumInterp, elevNumInterp), &
                radialElevInterp(radialNumInterp, elevNumInterp), &
                radialTimeInterp(radialNumInterp, elevNumInterp), &
                ValSwap(gateNumInterp, radialNumInterp, elevNumInterp))

      DO i = 0, elevNum - 2
        DO j = 1, elevScaleFactor + 1
          radialAzimInterp(:, i * (elevScaleFactor + 1) + j) = radialAzim(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialAzim(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          radialElevInterp(:, i * (elevScaleFactor + 1) + j) = radialElev(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialElev(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          radialTimeInterp(:, i * (elevScaleFactor + 1) + j) = radialTime(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialTime(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          ValSwap(:, :, i * (elevScaleFactor + 1) + j) = val(:, :, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + val(:, :, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
        END DO
      END DO
! ValSwap = val
      radialAzimInterp(:, elevNumInterp) = radialAzim(:, elevNum)
      radialElevInterp(:, elevNumInterp) = radialElev(:, elevNum)
      radialTimeInterp(:, elevNumInterp) = radialTime(:, elevNum)
      ValSwap(:, :, elevNumInterp) = val(:, :, elevNum)

      distInterp = dist
      !

      ALLOCATE (latSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                lonSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                altSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                timSwap(gateNumInterp, radialNumInterp, elevNumInterp))

      altSwap = 0.0D0

      DO i = 1, elevNumInterp
        DO j = 1, radialNumInterp
          DO k = 1, gateNumInterp
            IF (ValSwap(k, j, i) .NE. 0) THEN
              CALL calLatLon(radialAzimInterp(j, i) * degree2radian, radialElevInterp(j, i) * degree2radian, distInterp(k), &
                             latSite * degree2radian, lonSite * degree2radian, altSite, &
                             latSwap(k, j, i), lonSwap(k, j, i), altSwap(k, j, i))
              timSwap(k, j, i) = radialTimeInterp(j, i)
              IF (altSwap(k, j, i) > 9000D0) ValSwap(k, j, i) = 0
              ! IF (i >= 6)   ValSwap(k, j, i) = 0
              IF (altSwap(k, j, i) < 500D0) ValSwap(k, j, i) = 0
              ! IF (radialElevInterp(j, i)>4) ValSwap(k, j, i) = 0
            END IF
          END DO
        END DO
      END DO

!       WRITE (*, 1) MAXVAL(altSwap(:, :, 1)), MAXLOC(altSwap(:, :, 1)), MAXVAL(radialElevInterp(:, 1)), &
!         MAXVAL(ValSwap), MINVAL(ValSwap)
! 1     FORMAT('ExtractRdr - alt max: ', D10.3, ' at: ', 2I4, ' Elev max: ', D10.3, ' Ref max/min: ', 2D10.3)

      ! Pack all non-zero values
      latVar = PACK(latSwap, ValSwap /= 0)
      lonVar = PACK(lonSwap, ValSwap /= 0)
      altVar = PACK(altSwap, ValSwap /= 0)
      timVar = PACK(timSwap, ValSwap /= 0)
      valVar = PACK(ValSwap, ValSwap /= 0)

      PRINT *, 'AAA latVar: ', SIZE(latVar), MAXVAL(valVar), MINVAL(valVar)

      ! Yuanfu Xie changed the radar data to all non-negative values:
      !latVar = PACK(latSwap, ValSwap >= 0)
      !lonVar = PACK(lonSwap, ValSwap >= 0)
      !altVar = PACK(altSwap, ValSwap >= 0)
      !timVar = PACK(timSwap, ValSwap >= 0)
      !valVar = PACK(ValSwap, ValSwap >= 0)

      DEALLOCATE (radialAzimInterp, radialElevInterp, radialTimeInterp)
    END BLOCK

    DEALLOCATE (latSwap, lonSwap, altSwap, timSwap, ValSwap, val, SpmW)
  END SUBROUTINE
!> @brief
!! Extract data from raw format to lat/lon/alt/time/values
  SUBROUTINE extractVarsFromData(radialAzim, radialElev, radialTime, dist, values, &
                                 latSite, lonSite, altSite, &
                                 elevNum, radialNum, gateNum, &
                                 latVar, lonVar, altVar, timVar, valVar, &
                                 offset, scale, fillvalue, varType, unAmbigRange, radarType)

    REAL(r_kind), INTENT(IN) :: offset, scale
    INTEGER(i_byte), INTENT(IN) :: fillvalue(:)
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: radialAzim(:, :), &
                                             radialElev(:, :), &
                                             radialTime(:, :), &
                                             !  resolutionV(:), &
                                             dist(:)
    INTEGER(i_short), ALLOCATABLE, INTENT(IN) :: values(:, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: unAmbigRange(:)

    REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
    INTEGER(i_kind), INTENT(IN) :: elevNum, radialNum, gateNum
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: latVar(:), lonVar(:), altVar(:), timVar(:), &
                                                valVar(:)
    CHARACTER(*), INTENT(IN) :: varType
    CHARACTER(*), INTENT(IN) :: radarType

    INTEGER(i_long) :: i, j, k

    REAL(r_kind), ALLOCATABLE :: latSwap(:, :, :), &
                                 lonSwap(:, :, :), &
                                 altSwap(:, :, :), &
                                 timSwap(:, :, :), &
                                 ValSwap(:, :, :)
    REAL(r_kind), ALLOCATABLE :: val(:, :, :)

    ALLOCATE (val(gateNum, radialNum, elevNum))

    ! DO i = 1, elevNum
    !   DO j = 1, radialNum
    !     DO k = 1, gateNum
    !       IF (TRIM(varType) .eq. 'Ref') THEN ! For Reflectivity
    !         val(k, j, i) = countsTodBZ(values(k, j, i), offset, scale)
    !       ELSE IF (TRIM(varType) .eq. 'Vel') then ! For Velocity
    !         val(k, j, i) = countsToVel(values(k, j, i), offset, scale)
    !       END IF
    !     END DO
    !   END DO
    ! END DO

    ! print *, 'values', values
    val = 0.0D0
    IF (radarType == "SAD") THEN
      IF (TRIM(varType) .EQ. 'Ref') THEN ! For Reflectivity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsTodBZ(values(k, j, i), offset, scale)
      ELSE IF (TRIM(varType) .EQ. 'Vel') THEN ! For Velocity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i) .AND. dist(k) > 20000) &
          val(k, j, i) = countsToVel(values(k, j, i), offset, scale)
      END IF
    ELSEIF (radarType == "DXK") THEN
      IF (TRIM(varType) .EQ. 'Ref') THEN ! For Reflectivity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsTodBZ_DXK(values(k, j, i), offset, scale)
      ELSE IF (TRIM(varType) .EQ. 'Vel') THEN ! For Velocity
        FORALL (i=1:elevNum, j=1:radialNum, k=1:gateNum, dist(k) < unAmbigRange(i)) &
          val(k, j, i) = countsToVel_DXK(values(k, j, i), offset, scale)
      END IF
    END IF

    BLOCK ! Interpolation
      REAL(r_kind), ALLOCATABLE :: radialAzimInterp(:, :), &
                                   radialElevInterp(:, :), &
                                   radialTimeInterp(:, :), &
                                   !  resolutionV(:), &
                                   distInterp(:)
      ! INTEGER(i_byte), ALLOCATABLE :: valuesInterp(:, :, :)
      INTEGER(i_kind) :: elevNumInterp, radialNumInterp, gateNumInterp

      INTEGER :: elevScaleFactor
      !!
      ! radialAzimInterp = radialAzim
      ! radialElevInterp = radialElev
      ! radialTimeInterp = radialTime
      ! ValSwap = val
      ! distInterp = dist

      ! gateNumInterp = gateNum
      ! radialNumInterp = radialNum
      ! elevNumInterp = elevNum
      !
      elevScaleFactor = 0

      gateNumInterp = gateNum
      radialNumInterp = radialNum
      elevNumInterp = elevNum + (elevNum - 1) * elevScaleFactor

      ALLOCATE (radialAzimInterp(radialNumInterp, elevNumInterp), &
                radialElevInterp(radialNumInterp, elevNumInterp), &
                radialTimeInterp(radialNumInterp, elevNumInterp), &
                ValSwap(gateNumInterp, radialNumInterp, elevNumInterp))

      DO i = 0, elevNum - 2
        DO j = 1, elevScaleFactor + 1
          radialAzimInterp(:, i * (elevScaleFactor + 1) + j) = radialAzim(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialAzim(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          radialElevInterp(:, i * (elevScaleFactor + 1) + j) = radialElev(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialElev(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          radialTimeInterp(:, i * (elevScaleFactor + 1) + j) = radialTime(:, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + radialTime(:, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
          ValSwap(:, :, i * (elevScaleFactor + 1) + j) = val(:, :, i + 1) * ((elevScaleFactor + 2 - j) * 1.0D0 / (elevScaleFactor + 1)) + val(:, :, i + 2) * ((j - 1) * 1.0D0 / (elevScaleFactor + 1))
        END DO
      END DO
! ValSwap = val
      radialAzimInterp(:, elevNumInterp) = radialAzim(:, elevNum)
      radialElevInterp(:, elevNumInterp) = radialElev(:, elevNum)
      radialTimeInterp(:, elevNumInterp) = radialTime(:, elevNum)
      ValSwap(:, :, elevNumInterp) = val(:, :, elevNum)

      distInterp = dist
      !

      ALLOCATE (latSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                lonSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                altSwap(gateNumInterp, radialNumInterp, elevNumInterp), &
                timSwap(gateNumInterp, radialNumInterp, elevNumInterp))

      altSwap = 0.0D0

      ! FORALL (i=1:elevNumInterp, j=1:radialNumInterp, k=1:gateNumInterp, ValSwap(k, j, i) /= 0)
      !   CALL calLatLon(radialAzimInterp(j, i)*degree2radian, radialElevInterp(j, i)*degree2radian, distInterp(k), &
      !                  latSite*degree2radian, lonSite*degree2radian, altSite, &
      !                  latSwap(k, j, i), lonSwap(k, j, i), altSwap(k, j, i))
      !   timSwap(k, j, i) = radialTimeInterp(j, i)
      ! END FORALL

      DO i = 1, elevNumInterp
        DO j = 1, radialNumInterp
          DO k = 1, gateNumInterp
            IF (ValSwap(k, j, i) .NE. 0) THEN
              CALL calLatLon(radialAzimInterp(j, i) * degree2radian, radialElevInterp(j, i) * degree2radian, distInterp(k), &
                             latSite * degree2radian, lonSite * degree2radian, altSite, &
                             latSwap(k, j, i), lonSwap(k, j, i), altSwap(k, j, i))
              timSwap(k, j, i) = radialTimeInterp(j, i)
              IF (altSwap(k, j, i) > 9000D0) ValSwap(k, j, i) = 0
              ! IF (i >= 6)   ValSwap(k, j, i) = 0
              IF (altSwap(k, j, i) < 2500D0) ValSwap(k, j, i) = 0
            END IF
          END DO
        END DO
      END DO

!       WRITE (*, 1) MAXVAL(altSwap(:, :, 1)), MAXLOC(altSwap(:, :, 1)), MAXVAL(radialElevInterp(:, 1)), &
!         MAXVAL(ValSwap), MINVAL(ValSwap)
! 1     FORMAT('ExtractRdr - alt max: ', D10.3, ' at: ', 2I4, ' Elev max: ', D10.3, ' Ref max/min: ', 2D10.3)

      ! Pack all non-zero values
      latVar = PACK(latSwap, ValSwap /= 0)
      lonVar = PACK(lonSwap, ValSwap /= 0)
      altVar = PACK(altSwap, ValSwap /= 0)
      timVar = PACK(timSwap, ValSwap /= 0)
      valVar = PACK(ValSwap, ValSwap /= 0)

      PRINT *, 'AAA latVar: ', SIZE(latVar), MAXVAL(valVar), MINVAL(valVar)

      ! Yuanfu Xie changed the radar data to all non-negative values:
      !latVar = PACK(latSwap, ValSwap >= 0)
      !lonVar = PACK(lonSwap, ValSwap >= 0)
      !altVar = PACK(altSwap, ValSwap >= 0)
      !timVar = PACK(timSwap, ValSwap >= 0)
      !valVar = PACK(ValSwap, ValSwap >= 0)

      DEALLOCATE (radialAzimInterp, radialElevInterp, radialTimeInterp)
    END BLOCK

    DEALLOCATE (latSwap, lonSwap, altSwap, timSwap, ValSwap, val)
  END SUBROUTINE extractVarsFromData

  SUBROUTINE calLatLon(azim, elev, rbin, latSite, lonSite, altSite, lat, lon, alt)
    USE parameters_m, ONLY: ae
    ! Azimuth clockwise
    REAL(r_kind), INTENT(IN) :: azim, elev, rbin
    REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
    REAL(r_kind), INTENT(OUT) :: lat, lon, alt
    REAL(r_kind) :: gcr
    REAL(r_kind) :: SIN_ELEV, COS_ELEV, gcrPrbinXSIN_ELEV

    ! ! Real great circle distance
    ! gcr = ae + altSite
    ! gcr = 8500000 + altSite

    ! SIN_ELEV = DSIN(elev)
    ! COS_ELEV = DCOS(elev)
    ! gcrPrbinXSIN_ELEV = gcr + rbin*SIN_ELEV

    ! !
    ! lat = latSite + ATAN((rbin*(COS_ELEV*DCOS(azim)))/(gcrPrbinXSIN_ELEV))
    ! lon = lonSite + ATAN((rbin*(COS_ELEV*DSIN(azim)))/(gcrPrbinXSIN_ELEV))
    ! alt = gcrPrbinXSIN_ELEV/DCOS(ATAN((rbin*COS_ELEV)/gcrPrbinXSIN_ELEV)) - gcr

    CALL radar_to_latlon(lat, lon, alt, azim / degree2radian, rbin, elev / degree2radian, latSite / degree2radian, lonSite / degree2radian, altSite)

    ! print*, 'aa', azim/degree2radian, elev/degree2radian, rbin
    ! print*, 'bb', lat, lon, alt
    lat = lat * degree2radian
    lon = lon * degree2radian
  END SUBROUTINE calLatLon

! !> @brief
! !! Convert from polar coordinates to lat/lon/alt
!   SUBROUTINE calLatLon(azim, elev, rbin, latSite, lonSite, altSite, lat, lon, alt)
!     USE parameters_m, ONLY: ae

!     ! Azimuth clockwise
!     REAL(r_kind), INTENT(IN) :: azim, elev, rbin
!     REAL(r_kind), INTENT(IN) :: latSite, lonSite, altSite
!     REAL(r_kind), INTENT(OUT) :: lat, lon, alt

!     REAL(r_kind) :: gcr
!     REAL(r_kind) :: thetaLat
!     REAL(r_kind) :: thetaLon
!     REAL(r_kind) :: thetaAll
!     REAL(r_kind) :: alphaLat, alphaLon, alphaAll
!     REAL(r_kind) :: rLat
!     REAL(r_kind) :: rLon

!     ! Real great circle distance
!     gcr = ae + altSite

!     ! Theta at lat and lon direction
!     thetaLat = ATAN((rbin*DSIN(elev)), (rbin*DCOS(elev)*DCOS(azim)))
!     thetaLon = ATAN((rbin*DSIN(elev)), (rbin*DCOS(elev)*DSIN(azim)))
!     ! PRINT *, 'thetaLat', thetaLat, elev, azim, rbin

!     ! r distance at lat and lon direction
!     rLat = (rbin*DSIN(elev))/(DSIN(thetaLat))
!     rLon = (rbin*DSIN(elev))/(DSIN(thetaLon))

!     !
!     alphaLat = ATAN((rLat*DCOS(thetaLat)), (gcr + rLat*DSIN(thetaLat)))
!     alphaLon = ATAN((rLon*DCOS(thetaLon)), (gcr + rLon*DSIN(thetaLon)))
!     ! PRINT *, 'thetaLat', rbin

!     !
!     lat = latSite + alphaLat
!     lon = lonSite + alphaLon

!     !
!     alphaAll = ATAN((rbin*DCOS(elev)), (gcr + rbin*DSIN(elev)))
!     alt = (gcr/(DCOS(alphaAll))) + (rbin*DSIN(elev)/DCOS(alphaAll)) - gcr
!   END SUBROUTINE calLatLon

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(RadarRAW_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%latRef)) DEALLOCATE (this%latRef)
    IF (ALLOCATED(this%lonRef)) DEALLOCATE (this%lonRef)
    IF (ALLOCATED(this%timRef)) DEALLOCATE (this%timRef)
    IF (ALLOCATED(this%altRef)) DEALLOCATE (this%altRef)
    IF (ALLOCATED(this%valRef)) DEALLOCATE (this%valRef)

    IF (ALLOCATED(this%latVel)) DEALLOCATE (this%latVel)
    IF (ALLOCATED(this%lonVel)) DEALLOCATE (this%lonVel)
    IF (ALLOCATED(this%timVel)) DEALLOCATE (this%timVel)
    IF (ALLOCATED(this%AltVel)) DEALLOCATE (this%AltVel)
    IF (ALLOCATED(this%ValVel)) DEALLOCATE (this%ValVel)

    IF (ALLOCATED(this%unAmbigRangeVel)) DEALLOCATE (this%unAmbigRangeVel)
    IF (ALLOCATED(this%unAmbigRangeRef)) DEALLOCATE (this%unAmbigRangeRef)

  END SUBROUTINE destructor

!> @brief
!! Convert integer Z count value to dbz
  PURE ELEMENTAL FUNCTION countsTodBZ_DXK(zcounts, Z_offset, Z_scale) RESULT(dBZ)
    INTEGER(i_short), INTENT(in) :: zcounts
    ! INTEGER(i_byte), INTENT(IN) :: fillValue(:)
    REAL(r_kind), INTENT(in) :: Z_offset, Z_scale
    REAL(r_kind) :: dBZ
    INTEGER :: i

    dBZ = zcounts
    IF (dBZ .EQ. 255) THEN
      dBZ = 0  ! Range Ambiguous
      RETURN
    END IF

    dBZ = dBZ * Z_scale + Z_Offset
    IF (dBZ .LT. 0) dBZ = 0.0D0
  END FUNCTION

!> @brief
!! Convert integer V count value to radial velocity (meters/sec)
  PURE ELEMENTAL FUNCTION countsToVel_DXK(vcounts, V_offset, V_scale) RESULT(Vel)
    INTEGER(i_short), INTENT(in) :: vcounts
    ! INTEGER(i_byte), INTENT(IN) :: fillValue(:)
    REAL(r_kind), INTENT(in) :: V_offset, V_scale
    REAL(r_kind) :: Vel

    Vel = vcounts
    !  Convert from signed to unsigned
    IF (Vel .EQ. 255) THEN
      Vel = 0  ! Invalid Measurement
      RETURN
    END IF

    Vel = Vel * V_scale + V_offset
  END FUNCTION

!> @brief
!! Convert integer Z count value to dbz
  PURE ELEMENTAL FUNCTION countsTodBZ(zcounts, Z_offset, Z_scale) RESULT(dBZ)
    INTEGER(i_short), INTENT(in) :: zcounts
    ! INTEGER(i_byte), INTENT(IN) :: fillValue(:)
    REAL(r_kind), INTENT(in) :: Z_offset, Z_scale
    REAL(r_kind) :: dBZ
    INTEGER :: i
!      From the NetCDF header
!      Z:valid_range = 2b, -2b ;     (2 through 254)
!      Z:below_threshold = 0b ;      (0)
!      Z:range_ambiguous = 1b ;      (1)
!      Z:_FillValue = -1b ;          (255)

    dBZ = zcounts

    ! DO i = LBOUND(fillValue, 1), UBOUND(fillValue, 1)
    !   IF (zcounts .eq. fillValue(i)) THEN
    !     dBZ = 0
    !     RETURN
    !   END IF
    ! END DO

    ! IF (dBZ .gt. 127.) THEN
    !   PRINT *, 'error in Reflectivity: ', dbz_hold, zcounts
    !   dBZ = fillValue
    ! END IF
    ! Convert from signed to unsigned
    IF (dBZ .LT. 0.) dBZ = 256 + dBZ

    IF (dBZ .EQ. 1. .OR. dbZ .EQ. 0.) THEN
      dBZ = 0  ! Range Ambiguous
      RETURN
    END IF

    dBZ = dBZ * Z_scale + Z_Offset

    IF (dBZ .LT. 0) dBZ = 0.0D0
  END FUNCTION countsTodBZ

!> @brief
!! Convert integer V count value to radial velocity (meters/sec)
  PURE ELEMENTAL FUNCTION countsToVel(vcounts, V_offset, V_scale) RESULT(Vel)
    INTEGER(i_short), INTENT(in) :: vcounts
    ! INTEGER(i_byte), INTENT(IN) :: fillValue(:)
    REAL(r_kind), INTENT(in) :: V_offset, V_scale
    REAL(r_kind) :: Vel

    Vel = vcounts
    ! DO i = LBOUND(fillValue, 1), UBOUND(fillValue, 1)
    !   IF (vcounts .eq. fillValue(i)) THEN
    !     Vel = 0
    !     RETURN
    !   END IF
    ! END DO

    ! IF (Vel .gt. 127.) THEN
    !   PRINT *, 'error in Velocity: ', Vel
    !   Vel = 0.
    ! ENDIF

    !  Convert from signed to unsigned
    IF (Vel .LT. 0.) Vel = 256.+Vel

    IF (Vel .EQ. 1. .OR. Vel .EQ. 0.) THEN
      Vel = 0  ! Invalid Measurement
      RETURN
    END IF

    ! IF (V_resolution .eq. 0.) THEN ! QC Check
    !   Vel = fillValue
    !   RETURN
    ! END IF

    Vel = Vel * V_scale + V_offset
    ! IF (Vel .LT. 0) Vel = 0.0D0
  END FUNCTION

  PURE ELEMENTAL FUNCTION countsToSpmW(vcounts, V_offset, V_scale) RESULT(Vel)
    INTEGER(i_short), INTENT(in) :: vcounts
    ! INTEGER(i_byte), INTENT(IN) :: fillValue(:)
    REAL(r_kind), INTENT(in) :: V_offset, V_scale
    REAL(r_kind) :: Vel

    Vel = vcounts
    ! DO i = LBOUND(fillValue, 1), UBOUND(fillValue, 1)
    !   IF (vcounts .eq. fillValue(i)) THEN
    !     Vel = 0
    !     RETURN
    !   END IF
    ! END DO

    ! IF (Vel .gt. 127.) THEN
    !   PRINT *, 'error in Velocity: ', Vel
    !   Vel = 0.
    ! ENDIF

    !  Convert from signed to unsigned
    IF (Vel .LT. 0.) Vel = 256.+Vel

    IF (Vel .EQ. 1. .OR. Vel .EQ. 0.) THEN
      Vel = 0  ! Invalid Measurement
      RETURN
    END IF

    ! IF (V_resolution .eq. 0.) THEN ! QC Check
    !   Vel = fillValue
    !   RETURN
    ! END IF

    Vel = Vel * V_scale + V_offset
    ! IF (Vel .LT. 0) Vel = 0.0D0
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

    var = nc%setVariable("latVel", "f64", (/dim1/)); CALL var%setData(this%latVel / degree2radian)
    var = nc%setVariable("lonVel", "f64", (/dim1/)); CALL var%setData(this%lonVel / degree2radian)
    var = nc%setVariable("timVel", "f64", (/dim1/)); CALL var%setData(this%timVel)
    var = nc%setVariable("altVel", "f64", (/dim1/)); CALL var%setData(this%altVel)
    var = nc%setVariable("valVel", "f64", (/dim1/)); CALL var%setData(this%valVel)

    PRINT *, 'Output radar files.'

    ! close the dataset
    CALL nc%CLOSE()
  END SUBROUTINE

END MODULE RadarRAW_m
