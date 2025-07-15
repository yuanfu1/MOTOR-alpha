!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.RawAGRI
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2024/10/30, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE RawAGRI_m
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
  USE DP_utils_m
  USE Satellite_utils_m

  IMPLICIT NONE
  LOGICAL, PARAMETER :: DEBUG_MODE = .FALSE.

  TYPE :: RawAGRI_t

    CHARACTER(LEN=1024) :: configFile
    TYPE(State_t) :: X
    INTEGER(i_kind) :: nchans = 0, nRawObs = 0
    ! `nchans` will be numVars
    ! Do not deallocate nchans and rttov_chan_lists

    ! Satellite information
    ! Read from the static/Satellite/satinfo_fy4_1-agri.txt file
    INTEGER(i_kind) :: reglength, regwidth, cal_nobs, agri_nobs
    REAL(r_kind) :: agri_longitude
    REAL(r_kind) :: agri_reso
    REAL(r_kind) :: col_offset, col_scale, radius_earth_a, radius_earth_b, agri_h
    CHARACTER(len=50) :: inst_name = 'agri', platform_name = 'fy4_1'
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: chan_lists, ifuse

    ! Defined for ObsAGRI(ObsBase) type
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    REAL(r_single), ALLOCATABLE :: sol_zenith(:), sol_azi(:)
    REAL(r_kind), ALLOCATABLE :: latitude(:), longitude(:)
    REAL(r_kind), ALLOCATABLE :: obsvalue(:, :), obserr(:)
    INTEGER(i_kind), ALLOCATABLE :: obstime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:)
    INTEGER(i_kind), ALLOCATABLE :: iscanpos(:)

  CONTAINS

    PROCEDURE, PUBLIC  :: RawAGRI_setup
    PROCEDURE, PUBLIC  :: RawAGRI_read
    PROCEDURE, PUBLIC  :: Destroy_RawAGRI
    ! PROCEDURE, PRIVATE :: Write_RawObs_fy4_RawAGRI
    PROCEDURE, PUBLIC :: Convert2TBB
    PROCEDURE, PUBLIC :: Convert2latlon

  END TYPE RawAGRI_t

CONTAINS

  SUBROUTINE RawAGRI_setup(this, configFile)
    IMPLICIT NONE
    CLASS(RawAGRI_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    INTEGER(i_kind) :: i, j, k, istatus
    CHARACTER(LEN=1024) :: satinfo_file
    INTEGER(i_kind) :: ichan
    LOGICAL :: istat
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: ifuse

    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'satellite', this%platform_name)
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"// &
                   TRIM(this%platform_name)//'-'//TRIM(this%inst_name)//".txt"

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//'-'//TRIM(this%inst_name), this%nchans)

    ALLOCATE (this%chan_lists(this%nchans), this%ifuse(this%nchans))
    ALLOCATE (this%obserr(this%nchans))

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//'-'//TRIM(this%inst_name), this%nchans, chan_lists=this%chan_lists, &
                             ifuse=this%ifuse, Err=this%obserr, lon=this%agri_longitude, reso=this%agri_reso, col_offset=this%col_offset, &
                             col_scale=this%col_scale, radius_earth_a=this%radius_earth_a, radius_earth_b=this%radius_earth_b, satellite_h=this%agri_h)

  IF (yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'ChanList', ifuse) .NE. 0) THEN
    PRINT *, '====== Determine chan usage from the satinfo file ======'
  ELSE
    PRINT *, '====== Determine chan usage from the YAML file ======'
    PRINT *, ifuse
    this%ifuse = ifuse
  END IF
      
  IF (ALLOCATED(ifuse)) DEALLOCATE(ifuse)
  
  PRINT *, "RawAGRI_setup successfully run"

  END SUBROUTINE RawAGRI_setup

  SUBROUTINE RawAGRI_read(this, radFile, angleFile, cmFile)
    USE NETCDF
    USE HDF5
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    IMPLICIT NONE

    CLASS(RawAGRI_t) :: this
    TYPE(Read_HDF5_interfaces_t) :: HDF5_interfaces
    CHARACTER(LEN=*), INTENT(IN) :: radFile, angleFile, cmFile
    INTEGER(HID_T) :: rad_fileid, angle_fileid
    INTEGER(i_kind) :: num_agri_obs
    INTEGER(i_kind) :: iret, ncid, cmid, i, j, iobs
    CHARACTER(LEN=2)  :: chan_str
    CHARACTER(LEN=40) :: varname, GroupName
    INTEGER(4), ALLOCATABLE :: LineNumber(:, :), ColumnNumber(:, :), nomvalue(:, :, :)
    REAL(8), ALLOCATABLE :: obstime(:, :), obstime_2D(:, :), obstime_unix(:)
    INTEGER(i_kind), ALLOCATABLE :: gmt_time(:)
    REAL(r_single), ALLOCATABLE :: calvalue(:, :)
    INTEGER(4) :: rank
    INTEGER(4), ALLOCATABLE :: dims(:)
    CHARACTER(len=100) :: outputDir = "/public/home/wuyl/run/motor/MOTOR02/tests_DA/20220527_GGF_SURF_SOUND/tmp/"

    ! Defined for reading/writing raw obs only
    REAL(r_single), ALLOCATABLE :: sat_zenith_2D(:, :), sat_azi_2D(:, :)
    REAL(r_single), ALLOCATABLE :: sol_zenith_2D(:, :), sol_azi_2D(:, :)
    REAL(r_kind), ALLOCATABLE :: latitude_2D(:, :), longitude_2D(:, :)
    REAL(r_single), ALLOCATABLE :: obsvalue_2D(:, :, :)
    REAL(r_kind), ALLOCATABLE :: cloud_flag_2D(:, :)
    LOGICAL :: radFile_exists, angleFile_exists, cm_exists, file_exists

    IF (DEBUG_MODE) PRINT *, 'Note: the HDF5 read only works for Base processor for now'
    !IF (.not. sg%isBaseProc()) RETURN

    INQUIRE (file=TRIM(radFile), exist=radFile_exists)
    INQUIRE (file=TRIM(angleFile), exist=angleFile_exists)
    INQUIRE (file=TRIM(cmFile), exist=cm_exists)
    file_exists = radFile_exists .AND. angleFile_exists .AND. cm_exists
    IF (.NOT. file_exists) RETURN

    CALL HDF5_interfaces%open_hdf(TRIM(radFile), rad_fileid)
    IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
    IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'Data'
    CALL HDF5_interfaces%get_hdf_array_rank(rad_fileid, '/'//TRIM(GroupName)//'/NOMChannel01', rank)
    ALLOCATE (dims(rank))
    CALL HDF5_interfaces%get_hdf_array_dims(rad_fileid, '/'//TRIM(GroupName)//'/NOMChannel01', dims)
    IF (DEBUG_MODE) PRINT *, 'dims of NOMChannel01: ', dims

    this%regwidth = dims(1)
    this%reglength = dims(2)
    this%agri_nobs = this%reglength * this%regwidth
    this%nRawObs = this%agri_nobs

    CALL HDF5_interfaces%open_hdf(TRIM(radFile), rad_fileid)
    IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
    IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'Calibration'
    CALL HDF5_interfaces%get_hdf_array_rank(rad_fileid, '/'//TRIM(GroupName)//'/CALChannel01', rank)
    DEALLOCATE (dims)
    ALLOCATE (dims(rank))
    CALL HDF5_interfaces%get_hdf_array_dims(rad_fileid, '/'//TRIM(GroupName)//'/CALChannel01', dims)
    IF (DEBUG_MODE) PRINT *, 'dims of CALChannel01: ', dims
    this%cal_nobs = dims(1)
    DEALLOCATE (dims)

    IF (DEBUG_MODE) PRINT *, 'DEBUG: ', this%reglength, this%regwidth, this%cal_nobs, this%agri_nobs

    ! Defined in this LOCAL SUBROUTINE
    ALLOCATE (LineNumber(this%regwidth, this%reglength))
    ALLOCATE (ColumnNumber(this%regwidth, this%reglength))
    ALLOCATE (nomvalue(this%regwidth, this%reglength, this%nchans))
    ALLOCATE (calvalue(this%cal_nobs, this%nchans))
    ALLOCATE (obstime(this%reglength, 2))
    ALLOCATE (obstime_2D(this%regwidth, this%reglength))
    ALLOCATE (obstime_unix(this%nRawObs))

    ! Defined and saved for ObsAGRI(ObsBase) type
    ALLOCATE (this%sat_zenith(this%nRawObs))
    ALLOCATE (this%sat_azi(this%nRawObs))
    ALLOCATE (this%sol_zenith(this%nRawObs))
    ALLOCATE (this%sol_azi(this%nRawObs))
    ALLOCATE (this%latitude(this%nRawObs))
    ALLOCATE (this%longitude(this%nRawObs))
    ALLOCATE (this%obstime(this%nRawObs))
    ALLOCATE (this%obsvalue(this%nRawObs, this%nchans))
    ALLOCATE (this%cloud_flag(this%nRawObs))
    ALLOCATE (this%iscanpos(this%nRawObs))

    ! Defined for raw data dimensions and for Write_RawObs_fy4_agri
    ! Deallocated before the end of this LOCAL SUBROUTINE
    ALLOCATE (sat_zenith_2D(this%regwidth, this%reglength))
    ALLOCATE (sat_azi_2D(this%regwidth, this%reglength))
    ALLOCATE (sol_zenith_2D(this%regwidth, this%reglength))
    ALLOCATE (sol_azi_2D(this%regwidth, this%reglength))
    ALLOCATE (latitude_2D(this%regwidth, this%reglength))
    ALLOCATE (longitude_2D(this%regwidth, this%reglength))
    ALLOCATE (cloud_flag_2D(this%regwidth, this%reglength))
    ALLOCATE (obsvalue_2D(this%regwidth, this%reglength, this%nchans))

    IF (DEBUG_MODE) PRINT *, 'ALLOCATION OVER'

    IF (DEBUG_MODE) PRINT *, "agri_read angleFile=", TRIM(angleFile)

    ! Read NOMChannel
    DO j = 1, this%nchans
      WRITE (chan_str, '(i2.2)') this%chan_lists(j)
      IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
      IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'Data'
      varname = '/'//TRIM(GroupName)//'/NOMChannel'//chan_str
      IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
      CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), nomvalue(:, :, j))
    END DO

    ! Read CALChannel
    DO j = 1, this%nchans
      WRITE (chan_str, '(i2.2)') this%chan_lists(j)
      IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
      IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'Calibration'
      varname = '/'//TRIM(GroupName)//'/CALChannel'//chan_str
      IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
      CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), calvalue(:, j))
    END DO

    this%regwidth = SIZE(nomvalue, 1)
    this%reglength = SIZE(nomvalue, 2)
    this%cal_nobs = SIZE(calvalue, 1)
    this%agri_nobs = this%reglength * this%regwidth
    this%nRawObs = this%agri_nobs

    ! Read NOMObsTime
    IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
    IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'NOMObs'
    varname = '/'//TRIM(GroupName)//'/NOMObsTime'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), obstime)
    IF (DEBUG_MODE) PRINT *, obstime(1450:1470, 1)
    ! Since difference of the starttime and the endtime is less than 1s, we simplify the obstime as the starttime
    DO j = 1, this%reglength
      obstime_2D(:, j) = obstime(j, 1)
    END DO
    ! print*, 'obstime_2D:',maxval(obstime_2D), minval(obstime_2D)
    obstime_unix = RESHAPE(obstime_2D, (/this%nRawObs/))
    ! time format conversion
    ALLOCATE (gmt_time(6))
    DO iobs = 1, this%nRawObs
      IF (obstime_unix(iobs) > missing + 1) THEN
        gmt_time(1) = INT(obstime_unix(iobs) / 1E13)
        gmt_time(2) = INT(MOD(obstime_unix(iobs), 1E13) / 1E11)
        gmt_time(3) = INT(MOD(obstime_unix(iobs), 1E11) / 1E9)
        gmt_time(4) = INT(MOD(obstime_unix(iobs), 1E9) / 1E7)
        gmt_time(5) = INT(MOD(obstime_unix(iobs), 1E7) / 1E5)
        gmt_time(6) = INT(MOD(obstime_unix(iobs), 1E5) / 1E3)
        ! PRINT *, 'obstime_unix = ', obstime_unix(iobs)
        ! PRINT *, 'gmt_time = ', gmt_time
        CALL Time_GMT_to_Unix(gmt_time, this%obsTime(iobs))
      ELSE
        this%obsTime(iobs) = missing
      END IF
    END DO
    DEALLOCATE (obstime, obstime_2D, obstime_unix, gmt_time)

    CALL HDF5_interfaces%close_hdf()

    CALL HDF5_interfaces%open_hdf(TRIM(angleFile), angle_fileid)

    ! Read NOMSatelliteZenith
    IF (TRIM(this%platform_name) .EQ. 'fy4_1') GroupName = ''
    IF (TRIM(this%platform_name) .EQ. 'fy4_2') GroupName = 'Navigation'
    varname = '/'//TRIM(GroupName)//'/NOMSatelliteZenith'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), sat_zenith_2D)
    WHERE (sat_zenith_2D > 400.0) sat_zenith_2D = missing
    this%sat_zenith = RESHAPE(sat_zenith_2D, (/this%nRawObs/)) * degree2radian

    ! Read NOMSatelliteAzimuth
    varname = '/'//TRIM(GroupName)//'/NOMSatelliteAzimuth'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), sat_azi_2D)
    WHERE (sat_azi_2D > 400.0) sat_azi_2D = missing
    this%sat_azi = RESHAPE(sat_azi_2D, (/this%nRawObs/)) * degree2radian

    ! Read NOMSunZenith
    varname = '/'//TRIM(GroupName)//'/NOMSunZenith'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), sol_zenith_2D)
    WHERE (sol_zenith_2D > 400.0) sol_zenith_2D = missing
    this%sol_zenith = RESHAPE(sol_zenith_2D, (/this%nRawObs/)) * degree2radian

    ! Read NOMSunAzimuth
    varname = '/'//TRIM(GroupName)//'/NOMSunAzimuth'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), sol_azi_2D)
    WHERE (sol_azi_2D > 400.0) sol_azi_2D = missing
    this%sol_azi = RESHAPE(sol_azi_2D, (/this%nRawObs/)) * degree2radian

    ! Read LineNumber
    varname = '/'//TRIM(GroupName)//'/LineNumber'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), LineNumber)

    ! Read ColumnNumber
    varname = '/'//TRIM(GroupName)//'/ColumnNumber'
    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
    CALL HDF5_interfaces%read_hdf(angle_fileid, TRIM(varname), ColumnNumber)

    CALL HDF5_interfaces%close_hdf()

    ! Read cloud mask
    iret = nf90_open(TRIM(cmFile), nf90_NOWRITE, ncid)

    IF (iret /= 0) THEN
      PRINT *, "Cannot open NETCDF4 file "//cmFile
      RETURN
    END IF

    iret = nf90_inq_varid(ncid, 'CLM', cmid)
    iret = nf90_get_var(ncid, cmid, cloud_flag_2D)
    WHERE (cloud_flag_2D > 120) cloud_flag_2D = missing
    this%cloud_flag = RESHAPE(cloud_flag_2D, (/this%nRawObs/))
    ! cloud_flag:
    ! 0:cloud,1:probably cloud,2:probably clear,3:clear,126:space,127:fillvalue

    iret = nf90_close(ncid)

    ! convert to BT
    DO j = 1, this%nchans
      CALL Convert2TBB(this, nomvalue(:, :, j), calvalue(:, j), this%obsvalue(:, j), obsvalue_2D(:, :, j))
    END DO
    ! PRINT *, ' maxval/minval nomvalue =', maxval(nomvalue(:, :, 1)), minval(nomvalue(:, :, 1))
    ! PRINT *, ' maxval/minval calvalue =', maxval(calvalue(:, 1)), minval(calvalue(:, 1))

    ! convert LineNumber/ColumnNumber to latitude/longitude
    CALL Convert2latlon(this, ColumnNumber, LineNumber, latitude_2D, longitude_2D)
    this%latitude = RESHAPE(latitude_2D, (/this%nRawObs/))
    this%longitude = RESHAPE(longitude_2D, (/this%nRawObs/))

    this%iscanpos = missing
    ! CALL Calc_scan_angle(this, this%agri_longitude, this%longitude, this%latitude)

    IF (DEBUG_MODE) PRINT *, 'CLM=', MAXVAL(this%cloud_flag), MINVAL(this%cloud_flag)
    IF (DEBUG_MODE) PRINT *, 'obsvalue=', MAXVAL(this%obsvalue), MINVAL(this%obsvalue)
    IF (DEBUG_MODE) PRINT *, 'latitude=', MAXVAL(this%latitude) / degree2radian, MINVAL(this%latitude) / degree2radian
    IF (DEBUG_MODE) PRINT *, 'longitude=', MAXVAL(this%longitude) / degree2radian, MINVAL(this%longitude) / degree2radian

    DEALLOCATE (LineNumber, ColumnNumber, nomvalue, calvalue)

    ! PRINT *, "agri_read successfully run"

    ! PRINT *, 'DEBUG only: Write_RawObs_fy4_agri'
    IF (yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
    ! IF (this%X%sg%isBaseProc()) &
    ! CALL Write_RawObs_fy4_agri(this, TRIM(outputDir)//'/Raw_FY4_AGRI.nc')
    ! END IF
    DEALLOCATE (sat_zenith_2D, sat_azi_2D, sol_zenith_2D, sol_azi_2D)
    DEALLOCATE (latitude_2D, longitude_2D, obsvalue_2D, cloud_flag_2D)
    ! PRINT *, 'RawAGRI_read is successfully run'

  END SUBROUTINE RawAGRI_read

  SUBROUTINE Write_RawObs_fy4_agri(this, filename)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    CLASS(RawAGRI_t) :: this
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
    ALLOCATE (out_array(Npixels))

    ! Reverse latitude for plotting
    var = nc%setVariable("lat", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%latitude .LT. missing - 1.0)
      out_array = this%latitude * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    ! print *, 'lat = ', maxval(out_array), minval(out_array), maxval(this%latitude), minval(this%latitude)
    CALL var%setData(out_array)

    var = nc%setVariable("lon", "f32", (/dim1/))
    CALL var%setFillValue(9999.0)
    WHERE (this%longitude .LT. missing - 1.0)
      out_array = this%longitude * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    ! print *, 'lon = ', maxval(out_array), minval(out_array), maxval(this%longitude), minval(this%longitude)
    CALL var%setData(out_array)

    var = nc%setVariable("tbb02", "f32", (/dim1, dim2/))
    CALL var%setFillValue(9999.0)
    ! print *, 'tbb02 = ', maxval(this%ObsValue(:, 2)), minval(this%ObsValue(:, 2))
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

  END SUBROUTINE Write_RawObs_fy4_agri

  SUBROUTINE Convert2TBB(this, nom_data, cal_data, tbb, tbb2D)
    IMPLICIT NONE
    CLASS(RawAGRI_t) :: this
    INTEGER(4), DIMENSION(:, :), INTENT(IN) :: nom_data
    REAL(r_single), DIMENSION(:), INTENT(IN)      :: cal_data
    REAL(r_kind), DIMENSION(:), INTENT(OUT)     :: tbb
    REAL(r_single), DIMENSION(:, :), INTENT(OUT), OPTIONAL :: tbb2D
    ! *** Local vars ***
    INTEGER(4)                  :: tmp1D
    INTEGER :: i, j, iobs

    !*** Initialization of full disk brightness temperature data field
    tbb2D = 0.

    DO i = 1, this%regwidth
      DO j = 1, this%reglength
        tmp1D = nom_data(i, j)
        IF (tmp1D > 0 .AND. tmp1D < this%cal_nobs) THEN
          tbb2D(i, j) = cal_data(tmp1D)
        ELSE
          tbb2D(i, j) = missing
        END IF
      END DO
    END DO

    tbb = RESHAPE(tbb2D, (/this%nRawObs/))

  END SUBROUTINE Convert2TBB

  SUBROUTINE Convert2latlon(this, column, line, lat, lon)
    IMPLICIT NONE
    CLASS(RawAGRI_t) :: this
    REAL(r_kind), PARAMETER :: para = 2.0**(-16)
    REAL(r_kind) :: COFF ! column offset
    REAL(r_kind) :: CFAC ! column scale factor
    INTEGER(4), DIMENSION(:, :), INTENT(IN) :: column
    INTEGER(4), DIMENSION(:, :), INTENT(IN) :: line
    REAL(r_kind), DIMENSION(:, :), INTENT(OUT) :: lon
    REAL(r_kind), DIMENSION(:, :), INTENT(OUT) :: lat

    REAL(r_kind) :: LOFF, LFAC
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: x, y, cosx, cosy, siny, cos2y, hcosxcosy, cos2y_ea_eb_siny_2, &
                                                  sd, sn, s1, s2, s3, sxy
    REAL(r_kind) :: ea, eb! earch radius (km)
    REAL(r_kind) :: h ! height of satellite
    REAL(r_kind) :: lamdaD  ! sub sallite
    INTEGER(i_kind) :: icol, iline

    ALLOCATE (x(this%regwidth, this%reglength))
    ALLOCATE (y(this%regwidth, this%reglength))
    ALLOCATE (cosx(this%regwidth, this%reglength))
    ALLOCATE (cosy(this%regwidth, this%reglength))
    ALLOCATE (siny(this%regwidth, this%reglength))
    ALLOCATE (cos2y(this%regwidth, this%reglength))
    ALLOCATE (hcosxcosy(this%regwidth, this%reglength))
    ALLOCATE (cos2y_ea_eb_siny_2(this%regwidth, this%reglength))
    ALLOCATE (sd(this%regwidth, this%reglength))
    ALLOCATE (sn(this%regwidth, this%reglength))
    ALLOCATE (s1(this%regwidth, this%reglength))
    ALLOCATE (s2(this%regwidth, this%reglength))
    ALLOCATE (s3(this%regwidth, this%reglength))
    ALLOCATE (sxy(this%regwidth, this%reglength))

    COFF = this%col_offset
    CFAC = this%col_scale
    ea = this%radius_earth_a  ! earth radius (km)
    eb = this%radius_earth_b  ! earch radius (km)
    h = this%agri_h  ! height of satellite
    lamdaD = degree2radian * this%agri_longitude
    LOFF = COFF  ! line offset
    LFAC = CFAC  ! line scale factor

    x = degree2radian * ((column - COFF) / (para * CFAC))
    y = degree2radian * ((line - LOFF) / (para * LFAC))

    cosx = COS(x)
    cosy = COS(y)
    siny = SIN(y)
    cos2y = cosy**2
    hcosxcosy = h * cosx * cosy
    cos2y_ea_eb_siny_2 = cos2y + (ea / eb * siny)**2
    sd = SQRT(hcosxcosy**2.0 - cos2y_ea_eb_siny_2 * (h**2.0 - ea**2.0))
    sn = (hcosxcosy - sd) / cos2y_ea_eb_siny_2
    s1 = h - sn * cosx * cosy
    s2 = sn * SIN(x) * cosy
    s3 = -sn * siny
    sxy = SQRT(s1**2.0 + s2**2.0)

    ! radian already
    lon = ATAN(s2 / s1) + lamdaD
    lat = ATAN(ea**2.0 / eb**2.0 * s3 / sxy)

    ! PRINT *, ' maxval/minval lat =', maxval(lat)/degree2radian, minval(lat)/degree2radian
    ! PRINT *, ' maxval/minval lon =', maxval(lon)/degree2radian, minval(lon)/degree2radian
    WHERE (lon / degree2radian > 360.0 .OR. lon / degree2radian < 0.0 .OR. ISNAN(lon)) lon = missing
    WHERE (lat / degree2radian > 90.0 .OR. lat / degree2radian < -90.0 .OR. ISNAN(lat)) lat = missing

    DEALLOCATE (x, y, cosx, cosy, siny, cos2y, hcosxcosy, cos2y_ea_eb_siny_2)
    DEALLOCATE (sd, sn, s1, s2, s3, sxy)

  END SUBROUTINE Convert2latlon

  ! SUBROUTINE Calc_scan_angle(this)
  !   IMPLICIT NONE
  !   CLASS(RawAGRI_t) :: this
  !   REAL(r_kind), ALLOCATABLE(:) :: tmp

  !   tmp = this%longitude
  !   tmp = tand(this%agri_longitude - this%longitude) / sind(this%latitude)
  !   this%scan_angle = atand(tmp)

  !   DEALLOCATE(tmp)

  ! END SUBROUTINE Calc_scan_angle

  SUBROUTINE Destroy_RawAGRI(this)
    CLASS(RawAGRI_t) :: this

    ! Deallocate memories:
    IF (ALLOCATED(this%chan_lists)) DEALLOCATE (this%chan_lists)
    IF (ALLOCATED(this%ifuse)) DEALLOCATE (this%ifuse)
    IF (ALLOCATED(this%obserr)) DEALLOCATE (this%obserr)
    IF (ALLOCATED(this%sat_zenith)) DEALLOCATE (this%sat_zenith)
    IF (ALLOCATED(this%sat_azi)) DEALLOCATE (this%sat_azi)
    IF (ALLOCATED(this%sol_zenith)) DEALLOCATE (this%sol_zenith)
    IF (ALLOCATED(this%sol_azi)) DEALLOCATE (this%sol_azi)
    IF (ALLOCATED(this%latitude)) DEALLOCATE (this%latitude)
    IF (ALLOCATED(this%longitude)) DEALLOCATE (this%longitude)
    IF (ALLOCATED(this%obsvalue)) DEALLOCATE (this%obsvalue)
    IF (ALLOCATED(this%cloud_flag)) DEALLOCATE (this%cloud_flag)
    IF (ALLOCATED(this%obstime)) DEALLOCATE (this%obstime)
    IF (ALLOCATED(this%iscanpos)) DEALLOCATE (this%iscanpos)

  END SUBROUTINE Destroy_RawAGRI

END MODULE RawAGRI_m
