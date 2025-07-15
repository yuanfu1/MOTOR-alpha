!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.RawAGRI
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang, Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie for debugging the cloud mask with more output info 2022/12/24, @GBA-MWF, Shenzhen
!   Moving the zsland.dat to static directory
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE RawGIIRS_m
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
  LOGICAL, PARAMETER :: DEBUG_MODE = .TRUE.

  TYPE :: RawGIIRS_t

    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind) :: nchans, nchan_used, nRawObs

    ! `nchans` will be numVars
    ! Do not deallocate nchans and rttov_chan_lists

    ! Read from the static/Satellite/satinfo_fy4_1-giirs.txt file
    !INTEGER(i_kind) :: raw_ncol, raw_nline, raw_cy, raw_nobs
    INTEGER(i_kind) :: reglength, regwidth, cal_nobs, giirs_nobs
    REAL(r_kind) :: raw_longitude
    REAL(r_kind) :: raw_reso
    REAL(r_kind) :: col_offset, col_scale, radius_earth_a, radius_earth_b, giirs_h

    CHARACTER(len=50) :: inst_name = 'giirs', platform_name = 'fy4_1'
    INTEGER*4                                    :: yyyy, mn, dd, hh, mm
    ! REAL*8                                       :: obstime
    INTEGER(i_kind), ALLOCATABLE :: chan_lists(:), rttov_chan_lists(:)
    INTEGER(i_kind), ALLOCATABLE :: ifuse(:)
    INTEGER(i_kind), ALLOCATABLE :: iscanlines(:), iscanpos(:)
    ! REAL*8 , ALLOCATABLE :: obstime(:)
    INTEGER(i_kind), ALLOCATABLE :: obstime(:)
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    REAL(r_single), ALLOCATABLE :: sol_zenith(:), sol_azi(:)
    REAL(r_single), ALLOCATABLE :: scan_angle(:)
    REAL(r_kind), ALLOCATABLE :: latitude(:), longitude(:)
    REAL(r_single), ALLOCATABLE :: isurf_height(:)
    REAL(r_single), ALLOCATABLE :: isurf_type(:)
    REAL(r_kind), ALLOCATABLE :: obsvalue(:, :)
    REAL(r_kind), ALLOCATABLE :: obserr(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:)
    REAL(r_kind), ALLOCATABLE :: ihirsflag(:)
    REAL(r_kind), ALLOCATABLE :: iprepro(:, :)
    REAL(r_kind), ALLOCATABLE :: iavhrr(:, :)
    REAL(r_kind), ALLOCATABLE :: ts(:)
    REAL(r_kind), ALLOCATABLE :: tctop(:)

    ! TYPE(SingleGrid_t), POINTER :: sg

  CONTAINS

    PROCEDURE, PUBLIC  :: RawGIIRS_setup
    PROCEDURE, PUBLIC  :: RawGIIRS_read
    PROCEDURE, PUBLIC  :: Destroy_RawGIIRS
    ! PROCEDURE, PRIVATE :: READ_GIIRS_HDF_FILE
    !PROCEDURE, PRIVATE :: get_r2t
    !PROCEDURE, PRIVATE :: READ_GIIRS_ANGLE
    !PROCEDURE, PRIVATE :: landmask
    !PROCEDURE, PRIVATE :: tbb2apo_tbb
    !PROCEDURE, PRIVATE :: readCLMdat
    !PROCEDURE, PRIVATE :: lonlat2lc

  END TYPE RawGIIRS_t

  ! INTERFACE RawGIIRS_t
  !  PROCEDURE :: constructor_RawGIIRS
  ! END INTERFACE RawGIIRS_t

CONTAINS

  !FUNCTION constructor_RawGIIRS(configFile) RESULT(this)
  !  IMPLICIT NONE
  !  TYPE(RawGIIRS_t) :: this
  !  CHARACTER(LEN=1024), INTENT(IN) :: configFile
  ! TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg

  !  this%configFile = configFile

  !END FUNCTION constructor_RawGIIRS

  SUBROUTINE RawGIIRS_setup(this, configFile)
    ! This subroutine sets up:
    ! (1) namelist parameters
    ! (2) observation file names, including radiance, angle, latlon, etc
    ! (3) satellite info files
    ! (4) HDF/NC variable names
    IMPLICIT NONE
    CLASS(RawGIIRS_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    !CHARACTER(LEN=20), ALLOCATABLE :: rad_varname(:, :), angle_varname(:), latlon_varname(:)
    INTEGER(i_kind) :: nprefix, numvar_rad, numvar_angle, numvar_latlon
    INTEGER(i_kind) :: num_insts, num_plats, num_satinfos, num_satpaths, num_obsformats, numVars, obsFilesNum
    INTEGER(i_kind) :: i, j, k
    INTEGER(i_kind) :: II, nobs, ich, NCH, istatus
    REAL            :: surfem_tmp, raderr_tmp
    CHARACTER(LEN=1024) :: satinfo_file
    CHARACTER(LEN=2) :: chan_str
    INTEGER(i_kind) :: ichan, giirs_loc, fy4_1_loc
    INTEGER(i_kind), ALLOCATABLE :: rttov_chan_lists(:)

    LOGICAL :: istat

    ! (3) Prepare satllite info files. This is related to an individual instrument.
    ! Get channel lists, observation errors, and parameters for allsky DA.

    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS', 'satellite', this%platform_name)
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    !satinfo_file='/sources/MOTOR/static'
    satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"// &
                   TRIM(this%platform_name)//'-'//TRIM(this%inst_name)//".txt"
    INQUIRE (file=satinfo_file, exist=istat)
    IF (istat) THEN
      ! IF (DEBUG_MODE) PRINT *, 'RawGIIRS satinfo_file=', trim(satinfo_file)
      OPEN (unit=121, file=satinfo_file, status='old')
      READ (121, *)
      READ (121, *) this%nchans
      ! print *,"NCHAN:",this%nchans
      ALLOCATE (this%chan_lists(this%nchans))
      ALLOCATE (rttov_chan_lists(this%nchans))
      ALLOCATE (this%ifuse(this%nchans))
      ALLOCATE (this%obserr(this%nchans))

      NCH = 0
      DO ich = 1, this%nchans
        READ (121, *) this%chan_lists(ich), II, this%ifuse(ich), surfem_tmp, raderr_tmp
        IF (this%ifuse(ich) .GT. 0) THEN                      !! =1: use, / 0: not use
          NCH = NCH + 1
          this%obserr(ich) = raderr_tmp
          rttov_chan_lists(NCH) = this%chan_lists(ich)
          !print *,ich, NCH,"chan_lists:",this%chan_lists(ich),rttov_chan_lists(NCH),this%ifuse(ich),surfem_tmp,raderr_tmp
        END IF
      END DO
      ! this%ifuse = 1
      this%nchan_used = NCH
      ALLOCATE (this%rttov_chan_lists(this%nchan_used))
      this%rttov_chan_lists(1:this%nchan_used) = rttov_chan_lists(1:this%nchan_used)

      READ (121, *)
      READ (121, *)
      READ (121, *)
      READ (121, *) this%raw_longitude
      READ (121, *)
      READ (121, *) this%raw_reso
      READ (121, *)
      READ (121, *) this%col_offset, this%col_scale, this%radius_earth_a, this%radius_earth_b, this%giirs_h
      CLOSE (121)
      DEALLOCATE (rttov_chan_lists)
    ELSE
      PRINT *, '----------satinfo file was not found in RawGIIRS----------'
      STOP
    END IF

    ! PRINT *, "RawGIIRS_setup successfully run",NCH

  END SUBROUTINE RawGIIRS_setup

  SUBROUTINE RawGIIRS_read(this, nfile, CLMfilename)
    USE NETCDF
    USE HDF5
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    IMPLICIT NONE

    CLASS(RawGIIRS_t) :: this
    !  INTEGER*4,         INTENT(IN) :: num_file
    CHARACTER(LEN=*), INTENT(IN) :: nfile
    CHARACTER(LEN=*), INTENT(IN) :: CLMfilename

    CHARACTER(len=200)  :: filename
    CHARACTER(len=200)  :: txtfilename
    CHARACTER(len=200)  :: outfilename
    CHARACTER(LEN=20)   :: varname
    CHARACTER(LEN=2)    :: chan_str
    CHARACTER(len=200)  :: hourname
    CHARACTER(len=200)  :: minutename
    CHARACTER(len=45)   :: xx
    INTEGER(HID_T)      :: rad_fileid

    CHARACTER(len=9)    :: dsetname1 = "ES_RealLW"
    CHARACTER(len=14)   :: dsetlonlw = "IRLW_Longitude"
    CHARACTER(len=13)   :: dsetlatlw = "IRLW_Latitude"
    CHARACTER(len=23)   :: dsetQFlw = "QF_LWElementExploration"
    CHARACTER(len=10)   :: dsettime = 'NOMObsTime'
    CHARACTER(len=30)   :: anglelwname(4)

    INTEGER*4, PARAMETER :: num_detec = 128, num_lw = 689, num_mw = 961  !!nobs=num_file*num_detec

    REAL*4, DIMENSION(num_detec, num_lw)  :: rad_arr_lw
    REAL*4, DIMENSION(num_detec, num_lw)  :: bt_arr_lw
    REAL*4, DIMENSION(num_detec)         :: lon_arr_lw, lat_arr_lw, QF_lw
    REAL*4, DIMENSION(num_detec, 4)       :: angle_arr_lw

    INTEGER*4                                    :: i, j, k, iobs, nobs
    INTEGER*4                                    :: ii, jj, kk
    INTEGER*4                                    :: yyyy, mn, dd, hh, mm
    REAL*4, DIMENSION(num_detec, num_lw + num_mw)         :: tbb, tbb_apo
    REAL*4, DIMENSION(num_detec)                       :: isurf_height, isurf_type
    REAL*4, DIMENSION(num_detec, 4)                     :: angle_data
    REAL*4, DIMENSION(num_detec)                       :: lon, lat, QF, tlat, tlon
    REAL*4, DIMENSION(num_lw)                     :: wn_lw
    REAL*4                                       :: wn_temp, rad_temp, bt_temp
    ! REAL(8), ALLOCATABLE :: obstime(:)
    INTEGER(i_kind), ALLOCATABLE :: gmt_time(:)
    REAL(8)     :: obstime
    CHARACTER(len=14)   ::ctime
    INTEGER(4), ALLOCATABLE :: nomvalue(:)
    INTEGER*4            ::  errStatus, ioss     ! error status: 0 OK, +ve for fatal error
    TYPE(Read_HDF5_interfaces_t) :: HDF5_interfaces

    ! TYPE(type_rad)     :: giirs(nobs)

    REAL*4               :: cloudmask(num_detec)
    DATA anglelwname/"IRLW_SatelliteZenith", "IRLW_SatelliteAzimuth", "IRLW_SolarZenith", "IRLW_SolarAzimuth"/

    nobs = num_detec
    this%nRawObs = num_detec
    tbb = missing
    QF = missing

    ALLOCATE (this%sat_zenith(this%nRawObs))
    ALLOCATE (this%sat_azi(this%nRawObs))
    ALLOCATE (this%sol_zenith(this%nRawObs))
    ALLOCATE (this%sol_azi(this%nRawObs))
    ALLOCATE (this%latitude(this%nRawObs))
    ALLOCATE (this%longitude(this%nRawObs))
    ALLOCATE (this%iscanlines(this%nRawObs))
    ALLOCATE (this%iscanpos(this%nRawObs))
    ALLOCATE (this%obsvalue(this%nRawObs, this%nchans))
    ALLOCATE (this%obstime(this%nRawObs))
    ! ALLOCATE (obstime(this%nRawObs))
    ALLOCATE (this%cloud_flag(this%nRawObs))
    ALLOCATE (this%iavhrr(this%nRawObs, 13))
    ALLOCATE (this%ihirsflag(this%nRawObs))
    ALLOCATE (this%iprepro(this%nRawObs, 5))
    ALLOCATE (this%ts(this%nRawObs))
    ALLOCATE (this%tctop(this%nRawObs))
    ALLOCATE (this%isurf_height(this%nRawObs))
    ALLOCATE (this%isurf_type(this%nRawObs))

![1] Read parameter from shell
!READ(*,*)CLMfilename
!READ(*,*)outfilename
!READ(*,"(A200)")temp

![2] Get LW wave number
!wn_lw = (/(wn_temp,wn_temp=700.0,1130.0,0.625)/)
    DO i = 1, num_lw
      wn_lw(i) = 700.0 + (i - 1) * 0.625
    END DO

![3] Read GIIRS HDF FILE
!O i = 1,num_file
    ! WRITE(filename,'(A100)')nfile
    filename = TRIM(nfile)
    ! print *,"TBB read filename:",TRIM(filename),TRIM(CLMfilename)

    !CALL READ_GIIRS_HDF_FILE(filename,dsetname1,dsetlonlw,dsetlatlw,dsetQFlw,dsettime,num_detec,num_lw,&
    !                       rad_arr_lw,lon,lat,QF,obstime)
    CALL READ_GIIRS_HDF_FILE(filename, dsetname1, dsetlonlw, dsetlatlw, dsetQFlw, num_detec, num_lw, &
                             rad_arr_lw, lon, lat, QF)
!                         rad_arr_lw(:,:),lon(:),lat(:), QF(:))
!  CALL HDF5_interfaces%open_hdf(TRIM(nfile), rad_fileid)
!
    ! Read NOMChannel
!  DO j = 1, this%nchans
!    WRITE (chan_str, '(i2.2)') this%chan_lists(j)
!    varname = 'NOMChannel'//chan_str
!    IF (DEBUG_MODE) PRINT *, 'agri_read DEBUG varname=', varname
!    CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), nomvalue(:, :, j))
!  END DO

!END DO
!print *,"j,lat(j),lon(j)"
!do j=1,num_detec
!  print *,j,lat(j),lon(j)
!enddo
![4] Convert radio to TBB
    kk = 0
    ii = 0
!DO i = 1,num_file
    DO j = 1, num_detec
      DO k = 1, num_lw
        wn_temp = wn_lw(k)
        rad_temp = rad_arr_lw(j, k)
        IF (rad_temp >= -300.0 .AND. rad_temp <= 300.0) THEN
          CALL get_r2t(bt_temp, rad_temp, wn_temp)
          tbb(j, k) = bt_temp                      ! tbb((i-1)*num_detec+j,k) = bt_temp
        ELSE
          tbb(j, k) = missing
        END IF
      END DO
      IF (lat(j) .LT. -888.0 .OR. lon(j) .LT. -888.0) THEN
        lat(j) = missing
        lon(j) = missing
        kk = kk + 1
        IF (kk == num_detec) PRINT *, "RawGIIRS_read: All lat&lon are missing"
      ELSE
        ii = ii + 1
        tlat(ii) = lat(j)
        tlon(ii) = lon(j)
      END IF
    END DO
    IF (ii .GT. 0) PRINT *, 'Cover lat:', MINVAL(tlat(1:ii)), MAXVAL(tlat(1:ii)), 'lon:', MINVAL(tlon(1:ii)), MAXVAL(tlon(1:ii))
    PRINT *, 'raw giirs: tbb = ', MAXVAL(tbb), MINVAL(tbb)
!END DO
!  print *,"Convert radio to TBB over"
    !do j=1,nobs
    ! print *,j,lat(j),lon(j),(tbb(j,k),k=1,18)
    !enddo

![5] READ GIIRS HDF ANGLE
!DO i = 1,num_file
!  WRITE(filename,'(A100)')nfile
    filename = TRIM(nfile)
    DO j = 1, 4
      CALL READ_GIIRS_ANGLE(filename, anglelwname(j), num_detec, angle_data(:, j))
    END DO
!END DO
    ! print *,"READ GIIRS HDF ANGLE over"
    !do j=1,nobs
    !print *,j,lat(j),lon(j),(angle_data(j,k),k=1,4)
! enddo

![6] Convert 3D 2D to 1D
!k=1
!DO i = 1,num_file
!  DO j = 1,num_detec
!    tbb(k,1:689) = bt_arr_lw(i,j,:)
    !  angle_data(k,:) = angle_arr_lw(i,j,:)
    !  lon(k) = lon_arr_lw(i,j)
    ! lat(k) = lat_arr_lw(i,j)
    !  QF(k) = QF_lw(i,j)
    !  k = k+1
!  END DO
!END DO

![7] apodization function, num_mw no data?
    CALL tbb2apo_tbb(nobs, num_lw, num_mw, tbb, tbb_apo)
! print *,"apodization over"

![8] generate landmask
    CALL landmask(nobs, lat, lon, isurf_height, isurf_type)
! print *,"landmask over"

![9] generate cloudmask

!BLOCK
!  USE netcdf

    !INTEGER :: iret, ncid,  cmid
    !INTEGER(i_byte) :: cloud(2748,2748)

    ! iret = nf90_open('/mnt/d/motor/TestCase/2021050410/input/obs/satellite/L2_GLB_CLM_20210504070000_20210504071459_4000M_V0001.NC'&
    ! , nf90_NOWRITE, ncid)
    ! IF (iret /= 0) THEN
    !   PRINT *, "Cannot open NETCDF4 file "//CLMfilename
    !   RETURN
    ! END IF
    ! print *,"open CLM"
    ! iret = nf90_inq_varid(ncid, 'CLM', cmid)
    ! print *,"inq CLM", iret,cmid,shape(cloud)
    ! iret = nf90_get_var(ncid, cmid, cloud)
    ! print *,"get CLM"
!END BLOCK
! print *,'Inside RawGIIRS_read - readCLMdat: ',TRIM(CLMfilename)
    CALL readCLMdat(CLMfilename, nobs, lon, lat, cloudmask)
! print *,"cloudmask over"

! Read NOMObsTime
!varname = 'NOMObsTime'
!IF (DEBUG_MODE) PRINT *, 'giirs_read DEBUG varname=', varname
!CALL HDF5_interfaces%open_hdf(TRIM(filename), rad_fileid)
!CALL HDF5_interfaces%read_hdf(rad_fileid, TRIM(varname), obstime)
!IF (DEBUG_MODE) PRINT *, 'time=', obstime(1:5)
! print *,"FileTime-Xie debug:",TRIM(filename)

! This is hardcoded for now, one needs to count the length of character of filename to find the date string!!!
read(filename,"(A130,A14)") xx,ctime
PRINT *, 'file: ', xx,"  ctime: ",ctime
read(ctime,"(F14.0)") obstime
! print *,"ObsTime: ",obstime
    ALLOCATE (gmt_time(6))
    DO iobs = 1, this%nRawObs
      IF (obstime > missing + 1) THEN
        gmt_time(1) = INT(obstime / 1E10)
        gmt_time(2) = INT(MOD(obstime, 1E10) / 1E8)
        gmt_time(3) = INT(MOD(obstime, 1E8) / 1E6)
        gmt_time(4) = INT(MOD(obstime, 1E6) / 1E4)
        gmt_time(5) = INT(MOD(obstime, 1E4) / 1E2)
        gmt_time(6) = INT(MOD(obstime, 1E2) / 1E0)
        ! PRINT *, 'obstime_unix = ', obstime_unix(iobs)
        ! PRINT *, iobs,'gmt_time = ', gmt_time
        CALL Time_GMT_to_Unix(gmt_time, this%obsTime(iobs))
      ELSE
        this%obsTime(iobs) = missing
      END IF
    END DO
    DEALLOCATE (gmt_time)

    DO iobs = 1, nobs
      ! this%obstime(iobs)          = obstime
      this%iscanlines(iobs) = INT((iobs - 1) / 128) + 1     ! file number(1-num_file)
      this%iscanpos(iobs) = MOD(iobs - 1, 128) + 1       ! detec number(1-128)
      this%latitude(iobs) = lat(iobs) * degree2radian
      this%longitude(iobs) = lon(iobs) * degree2radian
      this%isurf_height(iobs) = isurf_height(iobs)
      this%isurf_type(iobs) = isurf_type(iobs)
      this%sat_zenith(iobs) = angle_data(iobs, 1) * degree2radian
      this%sat_azi(iobs) = angle_data(iobs, 2) * degree2radian
      this%sol_zenith(iobs) = angle_data(iobs, 3) * degree2radian
      this%sol_azi(iobs) = angle_data(iobs, 4) * degree2radian
      DO j = 1, 127      !num_lw+num_mw
        this%obsvalue(iobs, j) = tbb_apo(iobs, j)
      END DO
      this%obsvalue(iobs, 128) = tbb_apo(iobs, 362)
      this%iavhrr(iobs, 1:13) = missing
      this%ihirsflag(iobs) = QF(iobs)
      this%iprepro(iobs, 1:5) = 0
      this%cloud_flag(iobs) = cloudmask(iobs)                 !clearsky_data(iobs)
      this%ts(iobs) = missing
      this%tctop(iobs) = missing

      !print *,iobs,"RawGIIRS_read:",this%obstime(iobs),this%latitude(iobs)/degree2radian,this%longitude(iobs)/degree2radian,  &
      !        this%isurf_type(iobs),this%cloud_flag(iobs),this%obsvalue(iobs,24)
    END DO
    CALL Write_RawObs_fy4_giirs(this, '/public/home/wuyl/run/motor/intel/MOTOR/output/'//'Raw_FY4_GIIRS.nc')
  END SUBROUTINE RawGIIRS_read

  SUBROUTINE Write_RawObs_fy4_giirs(this, filename)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    CLASS(RawGIIRS_t) :: this
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

    var = nc%setVariable("tbb", "f32", (/dim1, dim2/))
    CALL var%setFillValue(9999.0)
    ! print *, 'tbb02 = ', maxval(this%ObsValue(:, 2)), minval(this%ObsValue(:, 2))
    CALL var%setData(this%ObsValue)

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

  END SUBROUTINE Write_RawObs_fy4_giirs

!SUBROUTINE READ_GIIRS_HDF_FILE(filename,dataname,lonname,latname,QFname,timename,d1,d2,data_arr,lon_arr,lat_arr,QF_arr,time_arr)
  SUBROUTINE READ_GIIRS_HDF_FILE(filename, dataname, lonname, latname, QFname, d1, d2, data_arr, lon_arr, lat_arr, QF_arr)
    USE HDF5
    IMPLICIT NONE

    CHARACTER(len=*)    ::   filename
    CHARACTER(len=*)      ::   dataname
    CHARACTER(len=*)     ::   lonname
    CHARACTER(len=*)     ::   latname
    CHARACTER(len=*)     ::   QFname
    !CHARACTER(len=*)     ::   timename

    INTEGER               ::  d1, d2, i, j
    INTEGER(HID_T)        ::  file_id
    INTEGER(HID_T)        ::  data_id
    INTEGER(HID_T)        ::   lon_id
    INTEGER(HID_T)        ::   lat_id
    INTEGER(HID_T)        ::    QF_id
    ! INTEGER(HID_T)        ::   time_id
    INTEGER               ::    error

    REAL, DIMENSION(d1, d2) :: data_arr
    REAL, DIMENSION(d1)    :: lon_arr, lat_arr, QF_arr
    REAL(8), DIMENSION(d1) :: time_arr

    INTEGER(HSIZE_T), DIMENSION(2)  :: dims
    INTEGER(HSIZE_T), DIMENSION(1)  :: dimss

    dims(1) = d1
    dims(2) = d2
!print *,"Read filename:",filename," data:",dataname," lon:",lonname," lat:",latname," QF:",QFname
    CALL h5open_f(error)
    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)

    CALL h5dopen_f(file_id, dataname, data_id, error)
    CALL h5dread_f(data_id, H5T_NATIVE_REAL, data_arr, dims, error)

    dimss(1) = d1

    CALL h5dopen_f(file_id, lonname, lon_id, error)
    CALL h5dread_f(lon_id, H5T_NATIVE_REAL, lon_arr, dimss, error)
!do j=1,d1
    ! print *,j,lon_id,lon_arr(j),error
!nddo

    CALL h5dopen_f(file_id, latname, lat_id, error)
    CALL h5dread_f(lat_id, H5T_NATIVE_REAL, lat_arr, dimss, error)

    CALL h5dopen_f(file_id, QFname, QF_id, error)
    CALL h5dread_f(QF_id, H5T_NATIVE_REAL, QF_arr, dimss, error)

!CALL h5dopen_f(file_id,timename,time_id,error)
!CALL h5dread_f(time_id,H5T_NATIVE_REAL,time_arr,dimss,error)
!do j=1,d1
!  print *,'QF_arr: ',j,lat_arr(j),lon_arr(j),QF_arr(j)
!enddo
    CALL h5close_f(error)

  END SUBROUTINE READ_GIIRS_HDF_FILE

  SUBROUTINE get_r2t(bt, radiation, wn)
    IMPLICIT NONE
    REAL         ::  bt, radiation, wn
    REAL         ::  H, C, B, C1, C2
    H = 6.6237E-27
    C = 2.99791E+10
    B = 1.38024E-16
    C1 = 2.*H * C * C
    C2 = H * C / B
    bt = C2 * wn / (LOG(C1 * wn * wn * wn / radiation + 1))

  END SUBROUTINE get_r2t

  SUBROUTINE READ_GIIRS_ANGLE(filename, anglelwname, d1, angle_arr)
    USE HDF5
    IMPLICIT NONE
    CHARACTER(len=200)    ::   filename
    CHARACTER(len=30)     ::   anglelwname

    INTEGER               ::  d1
    INTEGER(HID_T)        ::  file_id
    INTEGER(HID_T)        ::  angle_id
    INTEGER               ::  error

    REAL, DIMENSION(d1)    :: angle_arr

    INTEGER(HSIZE_T), DIMENSION(1)  :: dims

    dims(1) = d1

    CALL h5open_f(error)
    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)

    CALL h5dopen_f(file_id, TRIM(anglelwname), angle_id, error)
    CALL h5dread_f(angle_id, H5T_NATIVE_REAL, angle_arr, dims, error)

    CALL h5close_f(error)

  END SUBROUTINE READ_GIIRS_ANGLE

  SUBROUTINE landmask(nobs, lat, lon, isurf_height, isurf_type)

    INTEGER, INTENT(IN) :: nobs
    REAL, DIMENSION(nobs), INTENT(INOUT) :: lat, lon
    REAL, DIMENSION(nobs), INTENT(OUT)   :: isurf_height, isurf_type
    CHARACTER(len=1024)    ::   filename

!------local variables---------------
    REAL              :: xlat, xlon
    INTEGER           :: iy, jx, i, j, irec, Num_obs, Num_bad
    REAL              :: hei, lnd, r
    REAL              :: hh(720, 360), ls(720, 360)
    REAL              :: a1, a2, a3, a4
    REAL              :: alpha, alpha1, beta, beta1

!!open (68,file='zsland.dat',convert='big_endian',form='unformatted',access='direct',recl=720*360)
!!irec=1
!!read(68,rec=irec) hh
!!irec=2
!!read(68,rec=irec) ls
!!close(68)
!filename='/mnt/d/motor/TestGIIRS/rtstatic/'//'zsland.dat'
! Yuanfu Xie changed the zsland.dat to static/Satellite directory:
!filename='/Volumes/MacBackup/caseStudy/GIIRS-DATA/'//'zsland.dat'
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", filename)
    filename = TRIM(filename)//'/Satellite/zsland.dat'

    PRINT *, "zsland:", TRIM(filename)
    OPEN (68, FILE=filename, FORM='formatted')
    DO j = 1, 360
      READ (68, *) hh(:, j)
    END DO
    DO j = 1, 360
      READ (68, *) ls(:, j)
    END DO
    CLOSE (68)

    Num_bad = 0
    DO Num_obs = 1, nobs
      xlat = lat(Num_obs)
      xlon = lon(Num_obs)
      !print *,Num_obs,xlat,xlon,isurf_height(Num_obs),isurf_type(Num_obs)
      IF ((xlat .GE. 9999.0) .OR. (xlon .GE. 9999.0)) THEN
        lat(Num_obs) = missing
        lon(Num_obs) = missing
        !  isurf_height(Num_obs)=missing
        !  isurf_type(Num_obs)=missing
        Num_bad = Num_bad + 1
        CYCLE
      END IF
      !print *,Num_obs,xlat,xlon

      IF (xlon < 0) xlon = xlon + 360
      CALL Intp2DWeight(xlat, xlon, 0.5, 0.5, -89.75, 0.0, 360, 720, jx, iy, alpha, alpha1, beta, beta1)
      a1 = alpha1 * beta1
      a2 = alpha1 * beta
      a3 = alpha * beta
      a4 = alpha * beta1
      hei = a1 * hh(iy, jx) + a2 * hh(iy, jx + 1) + a3 * hh(iy + 1, jx + 1) + a4 * hh(iy + 1, jx)
      lnd = a1 * ls(iy, jx) + a2 * ls(iy, jx + 1) + a3 * ls(iy + 1, jx + 1) + a4 * ls(iy + 1, jx)

      isurf_height(Num_obs) = hei
      IF (isurf_height(Num_obs) .GT. 8888.0 .OR. isurf_height(Num_obs) .LT. 0.0) isurf_height(Num_obs) = 0
      isurf_type(Num_obs) = lnd
      IF (ABS(lnd - 1) .LT. 0.5) isurf_type(Num_obs) = 2
      IF (ABS(lnd - 2) .LT. 0.5) isurf_type(Num_obs) = 0
    END DO
    r = 0.0
    IF (nobs > 0) r = REAL(Num_bad) / REAL(nobs)
    PRINT *, "landmask: used ", r
    !DO Num_obs = 1,nobs
    !  print *,Num_obs,isurf_type(Num_obs),isurf_height(Num_obs)
    !ENDDO

  END SUBROUTINE landmask

  SUBROUTINE tbb2apo_tbb(nobs, num_lw, num_mw, tbb_0, tbb_apo)

    IMPLICIT NONE

    INTEGER                                 :: nobs, num_lw, num_mw
    REAL*4, DIMENSION(nobs, num_lw + num_mw)     :: tbb_0, tbb_apo
    INTEGER              :: iobs, ich
    REAL*4, DIMENSION(num_lw + num_mw)        :: tbb_temp, tbb_apo_temp

    DO iobs = 1, nobs
      tbb_temp(:) = tbb_0(iobs, :)
      ! ich = 1
      tbb_apo_temp(1) = 0.54 * tbb_temp(1) + 0.23 * tbb_temp(2)
      ! ich = num_lw/689
      tbb_apo_temp(num_lw) = 0.54 * tbb_temp(num_lw) + 0.23 * tbb_temp(num_lw - 1)
      ! ich = num_lw+1/690
      tbb_apo_temp(num_lw + 1) = 0.54 * tbb_temp(num_lw + 1) + 0.23 * tbb_temp(num_lw + 2)
      ! ich = num_lw+num_mw/1650
      tbb_apo_temp(num_lw + num_mw) = 0.54 * tbb_temp(num_lw + num_mw) + 0.23 * tbb_temp(num_lw + num_mw - 1)

      DO ich = 2, num_lw - 1
        tbb_apo_temp(ich) = 0.54 * tbb_temp(ich) + 0.23 * tbb_temp(ich - 1) + 0.23 * tbb_temp(ich + 1)
      END DO

      DO ich = num_lw + 2, num_lw + num_mw - 1
        tbb_apo_temp(ich) = 0.54 * tbb_temp(ich) + 0.23 * tbb_temp(ich - 1) + 0.23 * tbb_temp(ich + 1)
      END DO

      tbb_apo(iobs, :) = tbb_apo_temp(:)
    END DO

  END SUBROUTINE tbb2apo_tbb

  SUBROUTINE readCLMdat(CLMfilename, nobs, lon, lat, cloudmask)
    USE NETCDF
    USE HDF5
    IMPLICIT NONE

    CHARACTER(len=*)          :: CLMfilename     !="CLM_20171118_01.dat"
    INTEGER                     :: nobs            !=540*128
    REAL*4, DIMENSION(*)      :: lon, lat

    INTEGER           :: iobs, l, c, k, iline, icol, n
    REAL*4            :: lon_temp, lat_temp

    REAL        :: cloud_temp(16)
    REAL        :: cloudmask(*)
    INTEGER(1), ALLOCATABLE :: cloud_arr(:, :)
    INTEGER(i_kind) :: iret, ncid, cmid, i, j

    !******** read CLM file ******
    ! print *,"Inside readCLMdat: ",TRIM(CLMfilename)
    ! OPEN(98,FILE=TRIM(CLMfilename),FORM="unformatted",ACCESS="direct",RECL=2748*2748,Status="old")
    ! READ(98,REC=1)CLM
    ! CLOSE(98)
    ! print *,"CLM:",(CLM(n),n=1,100)
    !  cloud_arr = transpose(reshape(CLM,(/2748,2748/)))
    ! cloud=reshape(CLM,(/2748,2748/))
    ! print *,"reshape"
    ! print *,shape(cloud)
    ! print*, transpose(cloud)
    ! cloud_arr = transpose(cloud)
    ! print *,"transpose"
    ! Read cloud mask
    ! ALLOCATE (cloud_arr(1092,2748))
    ALLOCATE (cloud_arr(2748, 2748))

    ! PRINT*,'test 1', cloud_arr(1,1)
    cloud_arr = 0
    ! PRINT*, 'test 2'

    iret = nf90_open(TRIM(CLMfilename), nf90_NOWRITE, ncid)

    IF (iret /= 0) THEN
      PRINT *, "Cannot open NETCDF4 file "//CLMfilename
      RETURN
    END IF
    ! print *,"open CLM",iret,ncid
    iret = nf90_inq_varid(ncid, 'CLM', cmid)
    ! print *,"inq CLM", iret,ncid,cmid,shape(cloud_arr)
    iret = nf90_get_var(ncid, cmid, cloud_arr)
    ! print *,"get CLM",iret,ncid,cmid
    ! WHERE(cloud > 120) cloud = missing
    ! this%cloud_flag = reshape(cloud, (/2748*2748/))
    ! do j=1,2740,274
    !    print *,j,"CLM:",(cloud_arr(i,j),i=1,2740,274)
    ! enddo
    ! cloud_flag:
    ! 0:cloud,1:probably cloud,2:probably clear,3:clear,126:space,127:fillvalue
    !print *, 'Going here. 1'

    iret = nf90_close(ncid)

    !print *, 'Going here. 2'
    !*****************************

    n = 0
    DO iobs = 1, nobs
      cloudmask(iobs) = missing
      lon_temp = lon(iobs)
      lat_temp = lat(iobs)
      ! print *,"iobs:",iobs,lon_temp,lat_temp
      IF ((lat_temp .GT. 8888.0) .OR. (lon_temp .GT. 8888.0)) THEN
        n = n + 1
        IF (n .EQ. nobs) PRINT *, "lat/lon are all missing", n
        CYCLE
      END IF
      ! print *,"readCLMdat - start lonlat2lc: ",iobs,lon_temp,lat_temp
      CALL lonlat2lc(lon_temp, lat_temp, l, c)
      ! print *,"lonlat2lc:",lon_temp,lat_temp,l,c

      k = 1
      DO iline = 1, 4
        DO icol = 1, 4

          IF (cloud_arr(l - 2 + iline, c - 2 + icol) .EQ. 0) THEN
            cloud_temp(k) = 1.-0./3.
          ELSE IF (cloud_arr(l - 2 + iline, c - 2 + icol) .EQ. 1) THEN
            cloud_temp(k) = 1.-1./3.
          ELSE IF (cloud_arr(l - 2 + iline, c - 2 + icol) .EQ. 2) THEN
            cloud_temp(k) = 1.-2./3.
          ELSE IF (cloud_arr(l - 2 + iline, c - 2 + icol) .EQ. 3) THEN
            cloud_temp(k) = 1.-3./3.
          END IF
          k = k + 1
        END DO
      END DO

      cloudmask(iobs) = SUM(cloud_temp(1:16)) / 16.
      ! WRITE(*,*) iobs,lon_temp,lat_temp
      ! WRITE(*,*) l,c
      ! WRITE(*,*) cloudmask(iobs),cloud_temp(:)
      ! IF (cloudmask(iobs) .LE. 1.0D-4) &
!       WRITE(*,1) iobs,lon_temp,lat_temp,l,c,cloudmask(iobs) !,cloud_temp(1:16)
! 1     FORMAT('CMread: iobs: ',I3,' latlon: ',2D12.4,' lc: ',2I5,' cloudMask: ',D12.4)
    END DO

    !DEALLOCATE (cloud_arr)

  END SUBROUTINE readCLMdat

  SUBROUTINE lonlat2lc(lon0, lat0, l, c)

    IMPLICIT NONE

    REAL  :: lon0, lat0
    REAL  :: lon_rad, lat_rad
    INTEGER  :: l, c
    REAL  :: lambda_e, phi_e, Re, r1, r2, r3, rn, x, y
    REAL  :: n1, n2, n3
    REAL  :: COE = 2.0**(-16)
    REAL, PARAMETER :: ea = 6378.137                          ! the long half axis of the earth (km)
    REAL, PARAMETER :: eb = 6356.7523                         ! the short half axis of the earth (km)
    REAL, PARAMETER :: h = 42164.                             ! The distance from the center of the earth to the center of the satellite (km)
    REAL, PARAMETER :: lambda_d = 104.7 * pi / 180                ! The longitude of the subpoint of the satellite
    REAL, PARAMETER :: COFF = 1375.5, LOFF = 1375.5          !! values for resolution(4km)
    INTEGER, PARAMETER :: CFAC = 10233137, LFAC = 10233137   !! values for resolution(4km)

    !!!  step 1
!! check geographical lon and lat

!!!  step 2
!! convert unit of lon/lat(o --> rad)

    lon_rad = lon0 * pi / 180.
    lat_rad = lat0 * pi / 180.

    !!!  step 3
    !! convert geographical lon/lat to the earth's lon/lat

    lambda_e = lon_rad
    n1 = TAN(lat_rad)
    phi_e = ATAN(eb * eb * n1 / ea / ea)

    !!! step 4
    !! calculate Re

    n2 = COS(phi_e)
    n3 = (ea * ea - eb * eb) * n2 * n2 / ea / ea
    Re = eb / SQRT(1 - n3)

    !!!  step 5
    !! calculate r1,r2,r3

    r1 = h - (Re * COS(phi_e) * COS(lambda_e - lambda_d))
    r2 = -Re * COS(phi_e) * SIN(lambda_e - lambda_d)
    r3 = Re * SIN(phi_e)

    !!! step 6
    !! calculate rn,x,y

    rn = SQRT(r1 * r1 + r2 * r2 + r3 * r3)
    x = (ATAN(-r2 / r1)) * 180 / pi
    y = (ASIN(-r3 / rn)) * 180 / pi

!!!  step 7
!! calculate c,l

    c = NINT(COFF + x * COE * CFAC)
    l = NINT(LOFF + y * COE * LFAC)

  END SUBROUTINE lonlat2lc

  SUBROUTINE Intp2DWeight(lat_ob, lon_ob, dlat, dlon, lats, lonw, mi, mj, &
                          i, j, alpha, alpha1, beta, beta1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Discription:
!    Transfer obs. x to grid j and calculate its nondimensional
!    distance to grid j and j+1
!
!  Method:
!    obs. coord at G
!    oy = (lat_ob-latS)/dlat
!    ox = (lon_ob-latW)/dlon
!  where:
!    Lon_ob    :  longitude degree of obs at G.
!    Lat_ob    :  latitude degree of obs at G.
!    latS      :  the start lat. degree of grid,  at the southwestern corner;
!    latW      :  the start lon. degree of grid.
!    dlat      :  mesh space in degree in y direction
!    dlon      :  mesh space in degree in x direction
!    mi        :  Total number of grid in y direction
!    mj        :  Total number of grid in x direction
!
!    i+1,j                 i+1,j+1
!    |-----------------------|-^-
!    |           |           | |
!    |           |           |alpha1
!    |           |           | |
!    |           |           | |
!    |-----------O-----------|-v-
!    |           |           | ^
!    |           |           | |
!    |           |           | |
!    |           |           |alpha
!    |           |           | |
!    |           |           | |
!    |-----------|-----------|-v-
!    |<---beta-->|<--beta1-->|
!    i,j                   i,j+1
!
! --------   ---------    --------
! Current code owner: RCNMP
!
! History:
! Version      Date       Comment
! --------   ---------    --------
!   1.0      2001.08.12   Coding           Zhuang Shiyu
!   1.0      2001.10.25   Debug & test     ZHANG HUA
!
! Code Description:
!   Language:           Fortran 90
!   Software Standard:
!
! Parent module:
!   module_setupstructures
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! USE module_constants, ONLY : latS, lonW, dlat, dlon

    IMPLICIT NONE

    !* Subroutine arguments:
    REAL, INTENT(in)  :: lat_ob, lon_ob, lats, lonw
    REAL, INTENT(in)  :: dlat, dlon
    INTEGER, INTENT(in)  :: mi, mj
    REAL, INTENT(out) :: alpha, alpha1, beta, beta1
    INTEGER, INTENT(out) :: i, j

    !* Local parameters:
    REAL                                :: ox, oy
    !- End of header ----------------------------------------------

    oy = (lat_ob - lats) / dlat + 1.0
    ox = (lon_ob - lonw) / dlon + 1.0

!print *,'(lat_ob-latS)/dlat',lat_ob,latS,dlat
!print *,'(lon_ob-lonW)/dlon',lon_ob,lonW,dlon

    i = INT(oy + 1.0E-4)
    j = INT(ox + 1.0E-4)

    alpha = oy - REAL(i)
    alpha1 = 1.0 - alpha

    beta = ox - REAL(j)
    beta1 = 1.0 - beta

    IF (i <= 0) i = 1
    IF (i >= mi) i = mi - 1
    IF (j <= 0) j = 1
    IF (j >= mj) j = mj - 1

  END SUBROUTINE Intp2DWeight

  SUBROUTINE Destroy_RawGIIRS(this)
    CLASS(RawGIIRS_t) :: this
    ! Deallocate memories:

    IF (ALLOCATED(this%chan_lists)) DEALLOCATE (this%chan_lists)
    IF (ALLOCATED(this%rttov_chan_lists)) DEALLOCATE (this%rttov_chan_lists)
    IF (ALLOCATED(this%ifuse)) DEALLOCATE (this%ifuse)
    IF (ALLOCATED(this%obserr)) DEALLOCATE (this%obserr)
    IF (ALLOCATED(this%iscanlines)) DEALLOCATE (this%iscanlines)
    IF (ALLOCATED(this%iscanpos)) DEALLOCATE (this%iscanpos)
    IF (ALLOCATED(this%sat_zenith)) DEALLOCATE (this%sat_zenith)
    IF (ALLOCATED(this%sat_azi)) DEALLOCATE (this%sat_azi)
    IF (ALLOCATED(this%sol_zenith)) DEALLOCATE (this%sol_zenith)
    IF (ALLOCATED(this%sol_azi)) DEALLOCATE (this%sol_azi)
    IF (ALLOCATED(this%latitude)) DEALLOCATE (this%latitude)
    IF (ALLOCATED(this%longitude)) DEALLOCATE (this%longitude)
    IF (ALLOCATED(this%obsvalue)) DEALLOCATE (this%obsvalue)
    IF (ALLOCATED(this%cloud_flag)) DEALLOCATE (this%cloud_flag)
    IF (ALLOCATED(this%ihirsflag)) DEALLOCATE (this%ihirsflag)
    IF (ALLOCATED(this%iprepro)) DEALLOCATE (this%iprepro)
    IF (ALLOCATED(this%iavhrr)) DEALLOCATE (this%iavhrr)
    IF (ALLOCATED(this%ts)) DEALLOCATE (this%ts)
    IF (ALLOCATED(this%tctop)) DEALLOCATE (this%tctop)
    IF (ALLOCATED(this%obstime)) DEALLOCATE (this%obstime)
    IF (ALLOCATED(this%isurf_height)) DEALLOCATE (this%isurf_height)
    IF (ALLOCATED(this%isurf_type)) DEALLOCATE (this%isurf_type)

  END SUBROUTINE Destroy_RawGIIRS

END MODULE RawGIIRS_m
