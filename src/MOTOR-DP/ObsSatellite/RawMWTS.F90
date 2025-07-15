!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.RawMWTS
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2023/04/19, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE RawMWTS_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE mpObs_m, ONLY: mpObs_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE Read_HDF5_interfaces_m
  USE Dp_utils_m
  USE Satellite_utils_m

  IMPLICIT NONE
  LOGICAL, PARAMETER :: DEBUG_MODE = .FALSE.

  TYPE :: RawMWTS_t

    CHARACTER(LEN=1024) :: configFile
    TYPE(State_t) :: X
    INTEGER(i_kind) :: nchans = 0, nRawObs = 0
    ! `nchans` will be numVars
    ! Do not deallocate nchans and rttov_chan_lists

    INTEGER(i_kind) :: npixel, nscan, cal_nobs, mwts_nobs
    CHARACTER(len=50) :: inst_name = 'mwts2', platform_name = 'fy3_4' !FY3D
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: chan_lists, ifuse
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: oversea
    ! Defined for ObsMWTS(ObsBase) type
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    REAL(r_single), ALLOCATABLE :: sol_zenith(:), sol_azi(:)
    REAL(r_kind), ALLOCATABLE :: latitude(:), longitude(:)
    REAL(r_kind), ALLOCATABLE :: obsvalue(:, :), obserr(:)
    INTEGER(i_kind), ALLOCATABLE :: obstime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), landseamask(:)
    INTEGER(i_kind), ALLOCATABLE :: iscanpos(:)

  CONTAINS

    PROCEDURE, PUBLIC  :: RawMWTS_setup
    PROCEDURE, PUBLIC  :: RawMWTS_read
    PROCEDURE, PUBLIC  :: Destroy_RawMWTS
    ! PROCEDURE, PRIVATE :: Write_RawObs_fy4_RawMWTS

  END TYPE RawMWTS_t

CONTAINS

  SUBROUTINE RawMWTS_setup(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(RawMWTS_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, INTENT(IN) :: idxFile

    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
    INTEGER(i_kind) :: i, j, k
    CHARACTER(LEN=1024) :: satinfo_file
    INTEGER(i_kind) :: ichan
    LOGICAL :: istat

    IF (yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'satellite', platform_name) == 0) STOP
    IF (yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'instrument', inst_name) == 0) STOP
    this%platform_name = TRIM(platform_name(idxFile))
    this%inst_name = TRIM(inst_name(idxFile))
    DEALLOCATE (platform_name, inst_name)

    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"// &
                   TRIM(this%platform_name)//'-'//TRIM(this%inst_name)//".txt"

    ALLOCATE (this%chan_lists(this%nchans), this%ifuse(this%nchans))
    ALLOCATE (this%obserr(this%nchans))

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//'-'//TRIM(this%inst_name), this%nchans, chan_lists=this%chan_lists, &
                             ifuse=this%ifuse, Err=this%obserr)

    PRINT *, "RawMWTS_setup successfully run"

  END SUBROUTINE RawMWTS_setup

  SUBROUTINE RawMWTS_read(this, radFile)
    USE NETCDF
    USE HDF5
    USE AdvanceTime_m
    IMPLICIT NONE

    CLASS(RawMWTS_t) :: this
    TYPE(Read_HDF5_interfaces_t) :: HDF5_interfaces
    CHARACTER(LEN=*), INTENT(IN) :: radFile
    INTEGER(HID_T) :: rad_fileid
    INTEGER(i_kind) :: num_mwts_obs
    INTEGER(i_kind) :: iret, ncid, cmid, i, j, iobs
    CHARACTER(LEN=40) :: varname
    INTEGER(i_kind), ALLOCATABLE :: obstime(:, :), obstime_2D(:, :)
    INTEGER(i_kind) :: obstimetmp
    REAL(8) :: obstimetmp1
    INTEGER(i_kind), ALLOCATABLE :: gmt_time(:)
    INTEGER(4) :: rank, ichan
    INTEGER(4), ALLOCATABLE :: dims(:)
    CHARACTER(len=100) :: outputDir = "/Users/yaliwu/Downloads/"

    ! Defined for reading/writing raw obs only
    REAL(r_single), ALLOCATABLE :: sat_zenith_2D(:, :), sat_azi_2D(:, :)
    REAL(r_single), ALLOCATABLE :: sol_zenith_2D(:, :), sol_azi_2D(:, :)
    REAL(r_kind), ALLOCATABLE :: latitude_2D(:, :), longitude_2D(:, :), elevation_2D(:, :)
    REAL(r_kind), ALLOCATABLE :: ObsAngle_2D(:, :), LandCover_2D(:, :), LandSeaMask_2D(:, :)
    REAL(r_kind), ALLOCATABLE :: daycnt(:), mscnt(:)
    REAL(r_single), ALLOCATABLE :: obsvalue_2D(:, :, :)
    REAL(r_kind), ALLOCATABLE :: cloud_flag_2D(:, :)
    REAL(4) :: slope, intercept
    INTEGER(8) :: attr_dim
    LOGICAL :: radFile_exists, file_exists
    REAL(4) :: obstime1, obstime2

    IF (DEBUG_MODE) PRINT *, 'Note: the HDF5 read only works for Base processor for now'
    !IF (.not. sg%isBaseProc()) RETURN

    INQUIRE (file=TRIM(radFile), exist=radFile_exists)
    file_exists = radFile_exists
    IF (.NOT. file_exists) RETURN

    !!!!!!!!!!!!!!!!!!!!!!test only!!!
    CALL HDF5_interfaces%open_hdf(TRIM(radFile), rad_fileid)
    CALL HDF5_interfaces%get_hdf_array_rank(rad_fileid, '/Data/Earth_Obs_BT', rank)
    ALLOCATE (dims(rank))
    CALL HDF5_interfaces%get_hdf_array_dims(rad_fileid, '/Data/Earth_Obs_BT', dims)
    PRINT *, 'dims of Earth_Obs_BT: ', dims(1), dims(2), dims(3)
    IF (DEBUG_MODE) PRINT *, 'dims of Earth_Obs_BT: ', dims

    this%nchans = dims(3)
    this%nscan = dims(2)
    this%npixel = dims(1)
    this%mwts_nobs = this%npixel * this%nscan
    this%nRawObs = this%mwts_nobs
    DEALLOCATE (dims)

    ! Defined in this LOCAL SUBROUTINE
    ALLOCATE (obstime(8, this%nscan))
    ALLOCATE (obstime_2D(this%npixel, this%nscan))

    ! Defined and saved for ObsMWTS(ObsBase) type
    ALLOCATE (this%sat_zenith(this%nRawObs))
    ALLOCATE (this%sat_azi(this%nRawObs))
    ALLOCATE (this%sol_zenith(this%nRawObs))
    ALLOCATE (this%sol_azi(this%nRawObs))
    ALLOCATE (this%latitude(this%nRawObs))
    ALLOCATE (this%longitude(this%nRawObs))
    ALLOCATE (this%obstime(this%nRawObs))
    ALLOCATE (this%obsvalue(this%nRawObs, this%nchans))
    ALLOCATE (this%cloud_flag(this%nRawObs))
    ALLOCATE (this%landseamask(this%nRawObs))
    ALLOCATE (this%iscanpos(this%nRawObs))

    ! Defined for raw data dimensions and for Write_RawObs_fy5_mwts
    ! Deallocated before the end of this LOCAL SUBROUTINE
    ALLOCATE (sat_zenith_2D(this%npixel, this%nscan))
    ALLOCATE (sat_azi_2D(this%npixel, this%nscan))
    ALLOCATE (sol_zenith_2D(this%npixel, this%nscan))
    ALLOCATE (sol_azi_2D(this%npixel, this%nscan))
    ALLOCATE (latitude_2D(this%npixel, this%nscan))
    ALLOCATE (longitude_2D(this%npixel, this%nscan))
    ALLOCATE (elevation_2D(this%npixel, this%nscan))
    ALLOCATE (ObsAngle_2D(this%npixel, this%nscan))
    ALLOCATE (LandCover_2D(this%npixel, this%nscan))
    ALLOCATE (LandSeaMask_2D(this%npixel, this%nscan))
    ALLOCATE (daycnt(this%nscan))
    ALLOCATE (mscnt(this%nscan))
    ALLOCATE (cloud_flag_2D(this%npixel, this%nscan))
    ALLOCATE (obsvalue_2D(this%npixel, this%nscan, this%nchans))

    IF (DEBUG_MODE) PRINT *, 'ALLOCATION OVER'

    ! Get TBB data and Time info
    ! TBB unit: K
    varname = '/Data/Earth_Obs_BT'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), obsvalue_2D)
    attr_dim = 1
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Slope', (/attr_dim/), slope)
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Intercept', (/attr_dim/), intercept)
    obsvalue_2D = obsvalue_2D * slope + intercept
    PRINT *, 'slope = ', slope
    PRINT *, 'intercept = ', intercept

    ! varname = '/Data/Time'
    ! CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), obstime)
    ! varname = 'EpochTime'
    ! PRINT *, 'read_hdf_global_attribute of ', TRIM(varname)
    ! CALL HDF5_interfaces%read_hdf_global_attribute(rad_fileid, TRIM(varname), (/attr_dim/), obstime1)
    ! varname = 'Observing Beginning Time'
    ! PRINT *, 'read_hdf_global_attribute of ', TRIM(varname)
    ! CALL HDF5_interfaces%read_hdf_global_attribute(rad_fileid, TRIM(varname), (/attr_dim/), obstime2)
    ! PRINT *, 'obstime = ', obstime1, obstime2

    ! Get latitude/longitude
    ! latitude: -90~90; longitude: -180~180 degree
    varname = '/Geolocation/Latitude'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), latitude_2D)
    varname = '/Geolocation/Longitude'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), longitude_2D)
    ! ! convert (-90,90) to (0,90)
    ! WHERE(latitude_2D < 0.0 ) latitude_2D = latitude_2D + 180.0
    ! convert (-180,180) to (0,360)
    WHERE (longitude_2D < 0.0) longitude_2D = longitude_2D + 360.0

    ! Get evelation etc
    ! elevation Unit: m
    ! ObsAngle Unit: degree
    ! landcover Description: "The type of land cover, 0 Water 1; Evergreen NeedleleafForest 2; Evergreen Broadleaf Forest;
    !                                                 3 Deciduous NeedleleafForest; 4 Deciduous Broadleaf Forest; 5 Mixed Forests;
    !                                                 6 Closed Shrublands; 7 Open Shrublands; 8 Woody Savannas; 9 Savannas; 10 Grasslands;
    !                                                11 Permanent Wetlands; 12 Croplands; 13 Urban and Built-Up; 14 Cropland/Natural Vegetation Mosaic;
    !                                                15 Snow and Ice; 16 Barren or Sparsely Vegetated;
    !                           17 (IGBP Water Bodies, recoded to 0 for MODIS Land Product consistency.) 254 Unclassified 255 Fill Value";
    ! landseamask Description: 1 land, 2 continentalwater, 3 sea, 5 boundary
    ! SensorZenith/SolarZenith Unit: degree, valid range: 0, 18000
    ! SensorAzimuth/SolarAzimuth Unit: degree, valid range: -18000, 18000
    !
    ! varname = 'Geolocation/Altitude'
    ! CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), elevation_2D)
    ! varname = '/Geolocation/LandCover'
    ! CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), LandCover_2D)
    ! varname = '/Geolocation/LandSeaMask'
    ! CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), LandSeaMask_2D)
    varname = '/Geolocation/Scnlin_daycnt'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), daycnt)
    varname = '/Geolocation/Scnlin_mscnt'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), mscnt)
    attr_dim = 1
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Slope', (/attr_dim/), slope)
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Intercept', (/attr_dim/), intercept)
    mscnt = mscnt * slope + intercept
    varname = '/Geolocation/SensorAzimuth'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), sat_azi_2D)
    varname = '/Geolocation/SensorZenith'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), sat_zenith_2D)
    varname = '/Geolocation/SolarAzimuth'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), sol_azi_2D)
    varname = '/Geolocation/SolarZenith'
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), sol_zenith_2D)
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Slope', (/attr_dim/), slope)
    CALL HDF5_interfaces%read_hdf_attribute_real4(rad_fileid, TRIM(varname), 'Intercept', (/attr_dim/), intercept)
    sat_zenith_2D = sat_zenith_2D * slope + intercept
    sat_azi_2D = sat_azi_2D * slope + intercept
    sol_zenith_2D = sol_zenith_2D * slope + intercept
    sol_azi_2D = sol_azi_2D * slope + intercept
    WHERE (sat_azi_2D < 0.0) sat_azi_2D = sat_azi_2D + 360.0
    WHERE (sol_azi_2D < 0.0) sol_azi_2D = sol_azi_2D + 360.0
    WHERE (sat_zenith_2D > 180.0 .OR. sat_zenith_2D < 0.0) sat_zenith_2D = missing
    WHERE (sol_zenith_2D > 180.0 .OR. sol_zenith_2D < 0.0) sol_zenith_2D = missing
    WHERE (sat_azi_2D > 360.0 .OR. sat_azi_2D < 0.0) sat_azi_2D = missing
    WHERE (sol_azi_2D > 360.0 .OR. sol_azi_2D < 0.0) sol_azi_2D = missing

    ! time format conversion
    ALLOCATE (gmt_time(6))
    DO j = 1, this%nscan
      obstimetmp1 = daycnt(j) + DBLE(mscnt(j)) / 86400000.
      CALL mjd2cal(obstimetmp1, gmt_time(1), gmt_time(2), gmt_time(3), gmt_time(4), gmt_time(5), gmt_time(6))
      ! PRINT *, 'GMT time = ', gmt_time
      CALL Time_GMT_to_Unix(gmt_time, obstimetmp)
      obstime_2D(:, j) = obstimetmp
    END DO
    this%obsTime = RESHAPE(obstime_2D, (/this%nRawObs/))
    ! PRINT *, 'this%obsTime = ', this%obsTime(200)
    DEALLOCATE (obstime, obstime_2D, gmt_time)

    CALL HDF5_interfaces%close_hdf()

    ! reshape arrays from 2D to 1D
    this%latitude = RESHAPE(latitude_2D, (/this%nRawObs/)) * degree2radian
    this%longitude = RESHAPE(longitude_2D, (/this%nRawObs/)) * degree2radian
    this%sat_zenith = RESHAPE(sat_zenith_2D, (/this%nRawObs/)) * degree2radian
    this%sol_zenith = RESHAPE(sol_zenith_2D, (/this%nRawObs/)) * degree2radian
    this%sat_azi = RESHAPE(sat_azi_2D, (/this%nRawObs/)) * degree2radian
    this%sol_azi = RESHAPE(sol_azi_2D, (/this%nRawObs/)) * degree2radian
    this%landseamask = RESHAPE(LandSeaMask_2D, (/this%nRawObs/))
    DO ichan = 1, this%nchans
      this%obsvalue(:, ichan) = RESHAPE(obsvalue_2D(:, :, ichan), (/this%nRawObs/))
    END DO

    this%iscanpos = missing !???

    ! IF (DEBUG_MODE) PRINT *, 'CLM=', maxval(this%cloud_flag), minval(this%cloud_flag)
    IF (DEBUG_MODE) PRINT *, 'obsvalue=', MAXVAL(this%obsvalue), MINVAL(this%obsvalue)
    IF (DEBUG_MODE) PRINT *, 'latitude=', MAXVAL(this%latitude) / degree2radian, MINVAL(this%latitude) / degree2radian
    IF (DEBUG_MODE) PRINT *, 'longitude=', MAXVAL(this%longitude) / degree2radian, MINVAL(this%longitude) / degree2radian

    ! PRINT *, "mwts_read successfully run"

    ! PRINT *, 'DEBUG only: Write_RawObs_fy5_mwts3'
    ! IF(yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
    ! CALL Write_RawObs_fy5_mwts3(this, TRIM(outputDir)//'/Raw_FY3_MWTS2.nc')

    DEALLOCATE (sat_zenith_2D, sat_azi_2D, sol_zenith_2D, sol_azi_2D)
    DEALLOCATE (latitude_2D, longitude_2D, obsvalue_2D, cloud_flag_2D)
    DEALLOCATE (elevation_2D, ObsAngle_2D, LandCover_2D, daycnt, mscnt)
    ! PRINT *, 'RawMWTS_read is successfully run'

  END SUBROUTINE RawMWTS_read

  SUBROUTINE Write_RawObs_fy5_mwts(this, filename)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    CLASS(RawMWTS_t) :: this
    ! TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(*), INTENT(IN) :: filename
    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1, dim2
    TYPE(NcVariable)  :: var
    TYPE(NcGroup)     :: grp
    INTEGER :: Npixels, Nchans
    REAL(r_kind), ALLOCATABLE    :: out_array(:)

    ! IF (.not. sg%isBaseProc()) RETURN

    nc = NcDataset(TRIM(filename), "w")
    Npixels = SIZE(this%latitude, 1)
    Nchans = this%nchans
    dim1 = nc%setDimension("Npixels", Npixels)
    dim2 = nc%setDimension("Channel", Nchans)
    PRINT *, 'output NC dims: ', dim1, dim2
    ALLOCATE (out_array(Npixels))

    ! Reverse latitude for plotting
    var = nc%setVariable("lat", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%latitude .LT. missing - 1.0)
      out_array = this%latitude * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    PRINT *, 'lat = ', MAXVAL(out_array), MINVAL(out_array), MAXVAL(this%latitude), MINVAL(this%latitude)
    CALL var%setData(out_array)

    var = nc%setVariable("lon", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%longitude .LT. missing - 1.0)
      out_array = this%longitude * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    PRINT *, 'lon = ', MAXVAL(out_array), MINVAL(out_array), MAXVAL(this%longitude), MINVAL(this%longitude)
    CALL var%setData(out_array)

    var = nc%setVariable("tbb02", "f32", (/dim1, dim2/))
    CALL var%setFillValue(9999.0)
    PRINT *, 'tbb02 = ', MAXVAL(this%ObsValue(:, 2)), MINVAL(this%ObsValue(:, 2))
    CALL var%setData(this%ObsValue(:, 2))

    var = nc%setVariable("tbb03", "f32", (/dim1, dim2/))
    CALL var%setFillValue(9999.0)
    CALL var%setData(this%ObsValue(:, 3))

    var = nc%setVariable("sat_zenith", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%sat_zenith .LT. missing - 1.0)
      out_array = this%sat_zenith * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    PRINT *, 'sat_zenith = ', MAXVAL(out_array), MINVAL(out_array), MAXVAL(this%sat_zenith), MINVAL(this%sat_zenith)
    CALL var%setData(out_array)

    var = nc%setVariable("sat_azimuth", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%sat_azi .LT. missing - 1.0)
      out_array = this%sat_azi * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    CALL var%setData(out_array)

    var = nc%setVariable("cloud_flag", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    CALL var%setData(this%cloud_flag)

    ! close the dataset
    CALL nc%CLOSE()
    DEALLOCATE (out_array)

  END SUBROUTINE Write_RawObs_fy5_mwts

  SUBROUTINE Destroy_RawMWTS(this)
    CLASS(RawMWTS_t) :: this

    ! Deallocate memories:
    IF (ALLOCATED(this%chan_lists)) DEALLOCATE (this%chan_lists)
    IF (ALLOCATED(this%ifuse)) DEALLOCATE (this%ifuse)
    IF (ALLOCATED(this%oversea)) DEALLOCATE (this%oversea)
    IF (ALLOCATED(this%obserr)) DEALLOCATE (this%obserr)
    IF (ALLOCATED(this%sat_zenith)) DEALLOCATE (this%sat_zenith)
    IF (ALLOCATED(this%sat_azi)) DEALLOCATE (this%sat_azi)
    IF (ALLOCATED(this%sol_zenith)) DEALLOCATE (this%sol_zenith)
    IF (ALLOCATED(this%sol_azi)) DEALLOCATE (this%sol_azi)
    IF (ALLOCATED(this%latitude)) DEALLOCATE (this%latitude)
    IF (ALLOCATED(this%longitude)) DEALLOCATE (this%longitude)
    IF (ALLOCATED(this%obsvalue)) DEALLOCATE (this%obsvalue)
    IF (ALLOCATED(this%cloud_flag)) DEALLOCATE (this%cloud_flag)
    IF (ALLOCATED(this%landseamask)) DEALLOCATE (this%landseamask)
    IF (ALLOCATED(this%obstime)) DEALLOCATE (this%obstime)
    IF (ALLOCATED(this%iscanpos)) DEALLOCATE (this%iscanpos)

  END SUBROUTINE Destroy_RawMWTS

END MODULE RawMWTS_m
