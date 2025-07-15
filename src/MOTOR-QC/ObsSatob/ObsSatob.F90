!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsSatob
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang
! VERSION           : V 0.0
! HISTORY           :
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a Satob data structure.
MODULE ObsSatob_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m
  USE slint, ONLY: slint_init, src_grid, tgt_grid
  USE mpObs_m, ONLY: mpObs_t
  USE NMLRead_m
  USE ncReadVar_m
  USE WriteVar_m
  USE AdvanceTime_m
  USE conversions_m
  USE parameters_m, ONLY: degree2radian
  USE GtsQC_m

  IMPLICIT NONE

  TYPE type_obs_field
    REAL    :: VALUE       ! Observation
    REAL    :: corret_value       ! corret_value
    REAL    :: corret_mothod      ! corret_mothod
    REAL    :: error       ! Observational error
    INTEGER :: flag        ! Observational quality identifier
    REAL    :: o_b          ! o_b
    REAL    :: o_a         ! o_a
    REAL    :: bias      ! the bias correction
  END TYPE type_obs_field

  TYPE type_header_info
    CHARACTER(LEN=30)  :: name          ! Station name
    CHARACTER(LEN=9)   :: call_name     ! call of ships or plane
    INTEGER               :: itime(5)      ! yyyy-mm-dd-hh-mm
    REAL                  :: stn_id        ! 5 digit station identifer
    REAL                  :: obs_type_id   !
    REAL                  :: lat           ! Latitude in degree
    REAL                  :: lon           ! Longitude in degree
    INTEGER               :: position      ! the lat lon check flag
    TYPE(type_obs_field) :: elv           ! Elevation in m
    TYPE(type_obs_field) :: ps            ! Surface pressure
    INTEGER               :: levels        !(for temp it means number of levels)
    INTEGER               :: instru_code   ! instrument
    INTEGER               :: radia_correct ! the radiation correction method for temp
    REAL                  :: baseline_bias ! the first thick bias
    REAL                  :: QI            ! the qi index for satob
  END TYPE type_header_info

  TYPE type_level_data
    TYPE(type_obs_field) :: p ! Pressure in Pa
    TYPE(type_obs_field) :: u ! Wind x-component in m/s
    TYPE(type_obs_field) :: v ! Wind y-component in m/s
    TYPE(type_obs_field) :: dd ! Wind direction
    TYPE(type_obs_field) :: ff ! Wind speed
    TYPE(type_obs_field) :: rh ! relative humidity in %
    TYPE(type_obs_field) :: phi ! Height in gm
    TYPE(type_obs_field) :: T ! Temperature in K
    TYPE(type_obs_field) :: Td !
    TYPE(type_obs_field) :: Vtan ! Tangent wind for bda
  END TYPE type_level_data

  TYPE type_satob
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
  END TYPE type_satob

  TYPE, EXTENDS(ObsBase_t) :: ObsSatob_t
    INTEGER            :: num_satob
    LOGICAL            :: logic_satob_qc
    INTEGER, PARAMETER :: max_satob = 1200000,   ! ... of satellite wind obs.
    INTEGER, PARAMETER :: iobsmissing = -888888
    REAL, PARAMETER :: robsmissing = -888888.0
    REAL, PARAMETER :: inobsmissing = 999999.0
    REAL, PARAMETER :: elvmissing = 9999.0
    REAL, PARAMETER :: e = 0.000001

  CONTAINS
    PROCEDURE :: ObsInitial => satobInitial
    PROCEDURE :: ObsIngest => satobIngest
    PROCEDURE :: ObsQC => satobQC
  END TYPE ObsSatob_t

CONTAINS

  SUBROUTINE satobInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsSatob_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    INTEGER :: i

    ! read the obs satob files names
    CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", obsFileDir)
    obsFileDir = '/sources/data/satob202206/'
    CALL namelist_read(TRIM(configFile), "obsFileList_Satob", obsFileName)
    this%numfiles = UBOUND(obsFileName, 1)
    PRINT *, "amount of satob obs files: ", this%numfiles

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numfiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of satob obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! Read in an option:
    CALL namelist_read(TRIM(configFile), "logic_satob_qc", this%logic_satob_qc)
    CALL namelist_read(TRIM(configFile), "interpolation", this%interpolation)

    ! read the var names of obs files in configFile
    this%numVars = 4   !1!5
    ALLOCATE (this%varNames(this%numVars))

    this%varNames(1) = "uwnd"
    this%varNames(2) = "vwnd"
    this%varNames(3) = "temp"
    this%varNames(4) = "pres"

    this%obsType = "SATOB"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ! this%missing = 999999.00
    this%missing = 999999999.0D0
    this%invalid = 3.4D38

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius = 5.0D2        ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(3, :) = 5.0D2   ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.2D0 ! No influence under this threshold value

    ! Threshold values for QC:
    this%qcThreshold(4) = 500.0D0
    this%qcThreshold(5) = 500.0D0
    this%qcThreshold(1) = 500.0D0
    this%qcThreshold(2) = 500.0D0
    this%qcThreshold(3) = 1500.0D0 * 100.0D0

  END SUBROUTINE satobInitial

  SUBROUTINE satobIngest0(this, X)
    IMPLICIT NONE

    CLASS(ObsSatob_t) :: this
    TYPE(State_t) :: X

    ! local variables
    INTEGER :: nf_op, nc_id, var_id, i, num_obs, num_dims, nf90_max_va_dims
    INTEGER :: ifile, num_obs_total, num_dataFile, ind_obsin_beg, ind_obsin_end
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_obsin

    ! local variables, for get the unique obs in multi-files
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), obsErrs(:, :), olatlon(:, :)
    REAL(r_kind), ALLOCATABLE :: obsHght(:), land(:), obspres(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:), num_sitein(:), num_obsin(:)
    CHARACTER(LEN=20), ALLOCATABLE :: stNames(:)
    CHARACTER(LEN=10), ALLOCATABLE :: timestr(:)
    CHARACTER(LEN=30), ALLOCATABLE :: Id_time(:), Id_time_unique(:)
    LOGICAL, ALLOCATABLE :: index_unique_judge(:)
    INTEGER, ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER :: num_obs_total_unique, j
    CHARACTER(LEN=1024) :: fntxt_Idtime_a, fntxt_Idtime_u

    CHARACTER(LEN=5)  :: obsTypeh, siteId
    CHARACTER(LEN=20) :: nhead
    INTEGER :: ctime(4), sitenum, ios, ilevels, imm, isite, iobs, tmp
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
    REAL(r_kind) :: missingValue = 999999.0D0

    ! Get the total amount of observations in multi-files. -Jiongming Pang
    num_dataFile = this%numFiles
    ALLOCATE (num_sitein(num_dataFile))
    ALLOCATE (num_obsin(num_dataFile))
    num_obsin(:) = 0
    DO ifile = 1, num_dataFile
      PRINT *, 'filename: (', TRIM(this%fileNames(ifile)), ')'
      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a5,5i8)', iostat=ios) obsTypeh, ctime(1:4), sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        num_sitein(ifile) = 0
        CLOSE (101)
        CYCLE
      END IF
      num_sitein(ifile) = sitenum
      DO isite = 1, num_sitein(ifile)
        READ (101, '(3x,a5,5f10.2,i6)') siteId, rlat, rlon, relv, rday, rhm, ilevels
        num_obsin(ifile) = num_obsin(ifile) + ilevels
        DO j = 1, ilevels
          READ (101, *)
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file ', ifile, ':', num_obsin(ifile)
    END DO
    num_obs_total = SUM(num_obsin)
    PRINT *, 'Total amount of Obs in all files: ', num_obs_total

    ! Allocated memory for arrays for all data from multi-files. -Jiongming Pang
    ALLOCATE (obsData(num_obs_total, this%numVars))
    ALLOCATE (obsErrs(num_obs_total, this%numVars))
    ALLOCATE (olatlon(2, num_obs_total))
    ALLOCATE (obsTime(num_obs_total))
    ALLOCATE (obsHght(num_obs_total))
    ALLOCATE (land(num_obs_total))
    ALLOCATE (stNames(num_obs_total))
    ALLOCATE (timestr(num_obs_total))
    ALLOCATE (Id_time(num_obs_total))
    ALLOCATE (obspres(num_obs_total))

    ! obs file loop
    iobs = 0
    DO ifile = 1, num_dataFile
      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a5,5i8)', iostat=ios) obsTypeh, ctime(1:4), sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        CLOSE (101)
        CYCLE
      END IF
      DO isite = 1, sitenum
        READ (101, '(3x,a5,5f10.2,i6)') siteId, rlat, rlon, relv, rday, rhm, ilevels
        DO j = 1, ilevels
          iobs = iobs + 1
          ! READ(101, '(10f10.2)') p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
          READ (101, *) p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d

          IF (ABS(dd - missingValue) .LT. 0.01D0 .OR. ABS(ff - missingValue) .LT. 0.01D0) THEN
            obsData(iobs, 4) = this%missing   ! uwnd
            obsData(iobs, 5) = this%missing   ! vwnd
          ELSE
            CALL wind_to_uv(ff, dd, obsData(iobs, 4), obsData(iobs, 5))
          END IF
          IF (ABS(t - missingValue) .GT. 0.01D0 .AND. t .LT. 200.0D0) obsData(iobs, 1) = t + 273.15
          IF (ABS(t - missingValue) .LT. 0.01D0 .OR. t .GE. 200.0D0) obsData(iobs, 1) = this%missing

          IF (ABS(td - missingValue) .GT. 0.01D0 .AND. td .LT. 200.0D0) obsData(iobs, 2) = td + 273.15
          IF (ABS(td - missingValue) .LT. 0.01D0 .OR. td .GE. 200.0D0) obsData(iobs, 2) = this%missing

          IF (ABS(p - missingValue) .GT. 0.01D0 .AND. p .LT. 1500.0D0 * 100.0D0) obsData(iobs, 3) = p * 100.0D0
          IF (ABS(p - missingValue) .LT. 0.01D0 .OR. p .GE. 1500.0D0 * 100.0D0) obsData(iobs, 3) = this%missing

          obspres(iobs) = obsData(iobs, 3)

          IF (ABS(hgh - missingValue) .LT. 0.01D0) THEN
            obsHght(iobs) = this%missing
          ELSE
            obsHght(iobs) = hgh
          END IF

          imm = INT(rhm - ctime(4) * 100)
          CALL Time_GMT_to_Unix((/ctime(1:4), imm, 0/), tmp)
          IF (ABS(ttd - missingValue) .LT. 0.01D0) THEN
            obsTime(iobs) = tmp
          ELSE
            obsTime(iobs) = tmp + ttd
          END IF

          IF (ABS(lat_d - missingValue) .LT. 0.01D0 .OR. ABS(lon_d - missingValue) .LT. 0.01D0) THEN
            olatlon(1, iobs) = rlat   ! lat
            olatlon(2, iobs) = rlon   ! lon
          ELSE
            olatlon(1, iobs) = rlat + lat_d   ! lat
            olatlon(2, iobs) = rlon + lon_d   ! lon
          END IF

          land(iobs) = 1.0D0
          obsErrs(iobs, :) = this%missing * (-1.0D0)   !-999999999.0D0
          stNames(iobs) = siteId
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file ', ifile, ':', num_obsin(ifile)
    END DO

    ! get the index of the used elements
    ALLOCATE (index_unique_tmp(num_obs_total))
    DO i = 1, num_obs_total
      index_unique_tmp(i) = i
    END DO
    WHERE (ABS(obspres - this%missing) .LT. 0.01D0) index_unique_tmp = 0
    WHERE (ABS(obsHght - this%missing) .LT. 0.01D0) index_unique_tmp = 0
    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    PRINT *, "the amount of all input obs and the unique: ", num_obs_total, num_obs_total_unique

    ! get the unique elements
    ALLOCATE (this%obsData(num_obs_total_unique, this%numVars))
    ALLOCATE (this%obsErrs(num_obs_total_unique, this%numVars))
    ALLOCATE (this%olatlon(2, num_obs_total_unique))
    ALLOCATE (this%obsTime(num_obs_total_unique))
    ALLOCATE (this%obsHght(num_obs_total_unique))
    ALLOCATE (this%land(num_obs_total_unique))
    ALLOCATE (this%StNames(num_obs_total_unique))

    this%numObs = num_obs_total_unique
    this%obsData = obsData(index_unique, :)
    this%obsErrs = obsErrs(index_unique, :)
    this%olatlon = olatlon(:, index_unique) * degree2radian   ! convert degree to radian. - JMPang
    this%obsTime = obsTime(index_unique)
    this%obsHght = obsHght(index_unique)
    this%land = land(index_unique)
    this%StNames = stNames(index_unique)

    DO i = 1, this%numVars
      PRINT *, "max and min of ", TRIM(this%varNames(i)), " in Uni obs: ", &
        MAXVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. this%missing)), &
        MINVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. this%missing))
    END DO

    ! deallocate the local arrays
    IF (ALLOCATED(obsData)) DEALLOCATE (obsData)
    IF (ALLOCATED(obsErrs)) DEALLOCATE (obsErrs)
    IF (ALLOCATED(olatlon)) DEALLOCATE (olatlon)
    IF (ALLOCATED(obsTime)) DEALLOCATE (obsTime)
    IF (ALLOCATED(obsHght)) DEALLOCATE (obsHght)
    IF (ALLOCATED(land)) DEALLOCATE (land)
    IF (ALLOCATED(stNames)) DEALLOCATE (stNames)
    IF (ALLOCATED(timestr)) DEALLOCATE (timestr)
    IF (ALLOCATED(Id_time)) DEALLOCATE (Id_time)
    IF (ALLOCATED(index_unique_tmp)) DEALLOCATE (index_unique_tmp)
    IF (ALLOCATED(index_unique)) DEALLOCATE (index_unique)

    PRINT *, "DONE satobIngest"

  END SUBROUTINE satobIngest0

  SUBROUTINE satobIngest(this, X, satob, num_satob)

    IMPLICIT NONE

    CLASS(ObsSatob_t) :: this
    TYPE(State_t) :: X
    TYPE(type_satob), INTENT(inout)   :: satob(1:max_satob) ! observation data
    INTEGER, INTENT(inout)   :: num_satob

    ! local variables
    INTEGER :: nf_op, nc_id, var_id, i, num_obs, num_dims, nf90_max_va_dims
    INTEGER :: ifile, num_obs_total, num_dataFile, ind_obsin_beg, ind_obsin_end
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_obsin

    ! local variables, for get the unique obs in multi-files
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), obsErrs(:, :), olatlon(:, :)
    REAL(r_kind), ALLOCATABLE :: obsHght(:), land(:), obspres(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:), num_sitein(:), num_obsin(:)
    CHARACTER(LEN=20), ALLOCATABLE :: stNames(:)
    CHARACTER(LEN=10), ALLOCATABLE :: timestr(:)
    CHARACTER(LEN=30), ALLOCATABLE :: Id_time(:), Id_time_unique(:)
    LOGICAL, ALLOCATABLE :: index_unique_judge(:)
    INTEGER, ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER :: num_obs_total_unique, j
    CHARACTER(LEN=1024) :: fntxt_Idtime_a, fntxt_Idtime_u

    CHARACTER(LEN=5)  :: obsTypeh, siteId
    CHARACTER(LEN=20) :: nhead
    INTEGER :: ctime(4), sitenum, ios, ilevels, imm, isite, iobs, tmp
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
    REAL(r_kind) :: missingValue = 999999.0D0
    !------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len=3)            :: RTYPE
    INTEGER                     :: yyyy, mon, dd1, hh
    INTEGER                     :: num_head, num_dat, num_qcm, num_shift
    INTEGER                     :: nlev ! the number of levels in one station (in)
    INTEGER                     :: k, i
    REAL                        :: p, id, T, dd, ff  !content
    REAL                        :: flag_p, flag_T, flag_dd, flag_ff
    REAL                        :: a, b
    INTEGER                     :: step       ! flag position
    REAL                        :: flag
    INTEGER                     :: i1
    REAL                        :: stn_id, lat, lon
    REAL                        :: instru_code
    INTEGER                     :: status
!----- add whitelist -------------------------
    INTEGER                   :: num_whtlst, iend, unit_wht, whitelist, whtid
    LOGICAL                   :: alive
    CHARACTER(LEN=120)        :: filewht
    INTEGER, ALLOCATABLE      :: idflag(:)
    ALLOCATE (idflag(2000))
!--------------------------------------------------------------------
!   read list
!---------------------------------------------------------------------
    whitelist = 0
    num_whtlst = 0
    unit_wht = 60
    idflag = 0
    INTEGER, ALLOCATABLE      :: idflag(:)
    ALLOCATE (idflag(2000))
!--------------------------------------------------------------------
!   read list
!---------------------------------------------------------------------

!--------------------------------------------------------------------
! 1.0  read the date, actual data length and number of obs. stations
!---------------------------------------------------------------------
    step = element_lib_flag
    unit_in = 61
    PRINT *, 'step1 begin read the satob data', this%fileNames(1)
    OPEN (101, FILE=TRIM(this%fileNames(1)), STATUS='old', ACTION='READ')
    OPEN (unit_in, file=TRIM(this%fileNames(1)), FORM='FORMATTED', iostat=status, status='old', blocksize=1048576, buffered='yes', buffercount=1)
    IF (status /= 0) THEN
      PRINT *, 'file_in_satob is not exist=', TRIM(this%fileNames(1))
      num_satob = 0
      RETURN
    END IF
    i = 0
    READ (unit_in, *, END=333) RTYPE, yyyy, mon, dd1, hh, num_satob, num_head, num_dat, num_qcm
    PRINT *, 'RTYPE: ', RTYPE, yyyy, mon, dd1, hh, num_satob, num_head, num_dat, num_qcm

!--------------------------------------------------------------------
!   check the input file date  ! use obs_size_constants
!--------------------------------------------------------------------
    IF (yyyy /= yyyy0) THEN
      PRINT *, 'wrong input', yyyy, 'please check file=', file_in_satob
      RETURN
    END IF

    IF (mon /= mon0) THEN
      PRINT *, 'wrong input', mon, 'please check file=', file_in_satob
      RETURN
    END IF

    IF (dd1 /= dd0) THEN
      PRINT *, 'wrong input', dd, 'please check file=', file_in_satob
      RETURN
    END IF

    IF (hh /= hh0) THEN
      PRINT *, 'wrong input', hh, 'please check file=', file_in_satob
      RETURN
    END IF
!--------------------------------------------------------------------
!    read out the data
!--------------------------------------------------------------------
    PRINT *, 'original num_satob', num_satob
    DO i1 = 1, num_satob
      READ (unit_in, *, END=333) stn_id, lat, lon, a, b, instru_code
      PRINT *, i1, 'stn_id: ', stn_id, lat, lon, a, b, instru_code
      READ (unit_in, '(8f10.2)', END=333) p, T, dd, ff, flag_p, flag_T, flag_dd, flag_ff
      PRINT *, 'p/T/dd/ff: ', p, T, dd, ff, flag_p, flag_T, flag_dd, flag_ff

      IF (idflag(INT(stn_id)) .EQ. 0) THEN
        PRINT *, 'idflag(int(stn_id)): ', stn_id, idflag(INT(stn_id))
        CYCLE
      END IF
      IF (lat .GT. 89.0 .OR. lat .LT. -89.0) CYCLE
      IF (lon .GT. 180.0 .OR. lon .LT. -180.0) CYCLE

!   read data with QI >80
      IF (flag_ff > 100.) CYCLE
      IF (flag_ff < 80.) CYCLE
      !  if(instru_code==4.0)cycle
      !   if(instru_code==5.0)cycle
      !   if(instru_code==6.0)cycle
      !   if(instru_code>2)cycle

      IF (stn_id == 173.0) THEN
        IF (instru_code > 2 .AND. flag_ff < 85.) CYCLE
      ELSE
        IF (instru_code > 2) CYCLE
      END IF

      IF (p > 950.0) CYCLE
      IF (p < 150.0) CYCLE

      !   if(abs(lat)<20.and.p<125)cycle
      !  if(abs(lat)>20.and.p<175)cycle
!    if(abs(stn_id -471) < e)cycle    !for korea
!    if(abs(stn_id -810) < e)cycle    !for india
      !   remove out the station with no obs data
      IF (ABS(p - inobsmissing) < e) CYCLE
      IF (ABS(dd - inobsmissing) < e .OR. ABS(ff - inobsmissing) < e) CYCLE
      !------------------------------------------------------------------
      ! remove the data acord the EC data selection added by tianwh 20140504
      !-------------------------------------------------------------------
      IF (instru_code == 1.0 .AND. p >= 700.0 .AND. flag_ff < 85) CYCLE
      !  if(instru_code==1.0.and.p>=400.0.and.lat<20.0.and.lat>-20.0.and.flag_ff<85)cycle
      !    if(instru_code==1.0.and.p<700.0.and.p>400..and.flag_ff<90)cycle
      ! if(instru_code==2.0.and.p<=700.0)cycle
      IF (instru_code == 3.0 .AND. p >= 401.0) CYCLE
      IF (instru_code == 5.0 .AND. p >= 401.0) CYCLE
      IF (instru_code == 7.0 .AND. p >= 401.0) CYCLE
      IF (stn_id == 783.0 .AND. lat > -60 .AND. lat < 60) CYCLE
      IF (stn_id == 784.0 .AND. lat > -60 .AND. lat < 60) CYCLE

      !-------------------------------------------------------------------
      IF ((i + 1) > max_satob) THEN
        PRINT *, 'the max_satob is too small num_satob=', max_satob
        i = max_satob
        EXIT
      END IF
      IF ((lat <= tkmaxlat) .AND. (lat >= tkminlat)) THEN
      IF ((lon <= tkmaxlon) .AND. (lon >= tkminlon)) THEN

        i = i + 1
        satob(i)%header_info%STempn_id = stn_id
        satob(i)%header_info%lat = lat
        satob(i)%header_info%lon = lon
        satob(i)%header_info%instru_code = instru_code
        satob(i)%header_info%QI = flag_dd

!-------------------------------------------
!  change  999999.0 to -888888.0
!-------------------------------------------

        IF (ABS(p - inobsmissing) < e) p = robsmissing
        IF (ABS(T - inobsmissing) < e) T = robsmissing
        IF (ABS(dd - inobsmissing) < e) dd = robsmissing
        IF (ABS(ff - inobsmissing) < e) ff = robsmissing

        !-------------------------------------------
        ! set the initial QC flag
        !-------------------------------------------
        satob(i)%level_data%p%flag = flag_missing
        satob(i)%level_data%T%flag = flag_missing
        satob(i)%level_data%dd%flag = flag_missing
        satob(i)%level_data%ff%flag = flag_missing
        !-------------------------------------------
        ! observation time
        !-------------------------------------------
        satob(i)%header_info%itime(1) = yyyy            ! year
        satob(i)%header_info%itime(2) = mon             ! month
        satob(i)%header_info%itime(3) = INT(a)          ! dd
        satob(i)%header_info%itime(4) = INT(b / 100)      ! hh
        IF (dd1 == 1 .AND. a >= 28) THEN
          satob(i)%header_info%itime(2) = mon - 1
          IF (mon == 1) THEN
            satob(i)%header_info%itime(2) = 12
            satob(i)%header_info%itime(1) = yyyy - 1
          END IF
        END IF
      END IF
      satob(i)%header_info%itime(5) = MOD(INT(b), 100) ! mm
      !-----------------------------------
      satob(i)%level_data%p%VALUE = p
      satob(i)%level_data%T%VALUE = T
      satob(i)%level_data%dd%VALUE = dd
      satob(i)%level_data%ff%VALUE = ff

      !      call flag_replace(satob(i)%level_data%p%flag ,step, int(flag_p))
      !      call flag_replace(satob(i)%level_data%T%flag,step , int(flag_T))
      !      call flag_replace(satob(i)%level_data%dd%flag,step , int(flag_dd))
      !      call flag_replace(satob(i)%level_data%ff%flag,step , int(flag_ff))
      !-------------------------------------------------------
      END IF
      END IF

      !-------------------------------------------
      END DO  ! the loop of the station
333   CONTINUE
      num_satob = i
      DEALLOCATE (idflag)
      CLOSE (unit_in)

      PRINT *, 'step1 read the satob data is over ', 'num_satob  = ', num_satob

      RETURN

      END SUBROUTINE satobIngest

      SUBROUTINE satobQC(satob, num_satob)
        USE qc_flag_set_m
        !use   index_satob
        IMPLICIT NONE

        TYPE(type_satob), INTENT(inout)   :: satob(1:max_satob) ! observation data
        INTEGER, INTENT(inout)   :: num_satob
        ! TYPE (type_ob)              :: yo
        ! TYPE (type_xb)              :: xb ! Background structures

        INTEGER                     :: i, num_stn
        REAL :: stn_h, stn_p, c1, stn_id1

        !---------------------------
        !  read out the satob data
        !---------------------------

        IF (num_satob == 0) logic_satob_qc = .FALSE.
        IF (logic_satob_qc) THEN
          !----------------------------------------
          ! begin the QC
          !----------------------------------------
          CALL satob_extreme_check(satob, num_satob)

          !-----------------------------------------------

          CALL satob_inter_consist_check(satob, num_satob)

          !------------------------------------------------
          ! get the u v from dd ff !

          CALL satob_ddff_uv(satob, num_satob)
          !------------------------------------------------

          ! if(logic_first_guess_check)then
          !  call satob_bg_check(yo%satob,yo%num_satob,xb)
          ! print*,'the satob_first_guess_check is over'
          ! endif

          !-------------------------------------

          CALL satob_complex_dma(satob, num_satob)
          PRINT *, 'satob_complex_dma is over'
          !----------------------------------------------
        END IF
        RETURN

      END SUBROUTINE satobQC

      SUBROUTINE satob_extreme_check(satob, num_satob)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Description:
        ! Climate extremet check on temperature,wind
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        USE qc_flag_set_m

        IMPLICIT NONE

        TYPE(type_satob), INTENT(inout)   :: satob(1:max_satob) ! observation data
        INTEGER, INTENT(inout)   :: num_satob

        ! Local arguments:

        INTEGER :: num_stn, iflag
        REAL    :: tmin1, tmin2, tt, tmax1, tmax2, stn_id, fff, dd, p
!---------------------------------------------------------------------------
        PRINT *, 'step2 begin satob  extreme_check'

        !---------------------------------------------------------------------------
        DO num_stn = 1, num_satob
          !  wind  extreme check
          fff = satob(num_stn)%level_data%ff%VALUE
          IF (ABS(fff - robsmissing) > e) THEN
            CALL w_extreme_qc(satob(num_stn)%level_data%p%VALUE, satob(num_stn)%header_info%lat, &
                              fff, satob(num_stn)%level_data%ff%flag)
            satob(num_stn)%level_data%dd%flag = satob(num_stn)%level_data%ff%flag
          END IF

          ! pressure extreme check
          p = satob(num_stn)%level_data%p%VALUE
          IF (ABS(p - robsmissing) > e) THEN
            IF (p < 150) THEN
              iflag = 2
              !print *,'p < 150',stn_id,p
              CALL flag_replace(satob(num_stn)%level_data%p%flag, extreme_flag, iflag)
              CALL flag_replace(satob(num_stn)%level_data%p%flag, DMA_flag, iflag)
            END IF
          END IF
        END DO
        PRINT *, 'step2 satob extreme_check is over'

        RETURN
      END SUBROUTINE satob_extreme_check

      SUBROUTINE satob_inter_consist_check(satob, num_satob)
        !+++++++++++++++++++++++++++++++++++++++++++++++++
        !Description:
        !Internal consistency check on T, Td, wind
!+++++++++++++++++++++++++++++++++++++++++++++++++++

        IMPLICIT NONE
        TYPE(type_satob), INTENT(inout) :: satob(1:max_satob) !observation data

!  Local arguments:
        INTEGER :: num_stn, num_satob
!---------------------------------------
        PRINT *, 'step3 begin satob inter_consist_check'
        DO num_stn = 1, num_satob
          CALL satob_wind_consist(num_stn, satob)
        END DO
        PRINT *, 'step3 satob_consist_check is over'
        RETURN
      END SUBROUTINE satob_inter_consist_check

      SUBROUTINE satob_wind_consist(num_stn, satob)
        !********************************************************************
        ! this subroutine is for checking consistency about wind on direction and speed, and flag is set.
        !
        ! modified by tianwh 2012.07.13
        !********************************************************************

        ! USE inter_consist_check

        IMPLICIT NONE

        !  Subroutine arguments:
        TYPE(type_satob), INTENT(inout) :: satob(1:max_satob) ! observation data
        INTEGER, INTENT(in)    :: num_stn
        !  Locale arguments:
        INTEGER                        :: iol, iflag
        REAL                           :: dd, ff

        dd = satob(num_stn)%level_data%dd%VALUE
        ff = satob(num_stn)%level_data%ff%VALUE

        IF (MOD(satob(num_stn)%level_data%dd%flag, 10) == 2 .OR. &
            MOD(satob(num_stn)%level_data%ff%flag, 10) == 2) RETURN
        IF (MOD(satob(num_stn)%level_data%dd%flag, 10) == 3 .OR. &
            MOD(satob(num_stn)%level_data%ff%flag, 10) == 3) RETURN
        IF (ABS(dd - robsmissing) < e .AND. ABS(ff - robsmissing) < e) RETURN

        CALL wind_consist_check(dd, ff, iflag)

        CALL flag_replace(satob(num_stn)%level_data%dd%flag, inter_consist_flag, iflag)

        CALL flag_replace(satob(num_stn)%level_data%ff%flag, inter_consist_flag, iflag)

        IF (iflag == 2) THEN
          CALL flag_replace(satob(num_stn)%level_data%dd%flag, DMA_flag, iflag)

          CALL flag_replace(satob(num_stn)%level_data%ff%flag, DMA_flag, iflag)

          !print *,'dd and fff inconsist',satob(num_stn)%header_info%STempn_id,dd,ff
          !    print*,satob(num_stn)%header_info%STempn_id,satob(num_stn)%level_data%p%value,u,v
        END IF

        RETURN
      END SUBROUTINE satob_wind_consist

      SUBROUTINE satob_ddff_uv(satob, num_satob)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Description:
        !  purpose: computer the u v from dd ff
        !
        ! method :
        !          ANG=450-(DD+180)
        !          u=FF*cos(ANG)
        !          v=FF*sin(ANG)
        !
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IMPLICIT NONE
        TYPE(type_satob)            :: satob(max_satob)
        INTEGER                      :: num_satob
        REAL                         :: ANG
        INTEGER                      :: i, k, flag
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        PRINT *, 'begin compute the u v from dd ff'
        DO i = 1, num_satob    !the loop of the station numbers
          flag = MAX(satob(i)%level_data%dd%flag, satob(i)%level_data%ff%flag)
          CALL ddff_uv(satob(i)%level_data%dd%VALUE, satob(i)%level_data%ff%VALUE, &
                       satob(i)%level_data%u%VALUE, satob(i)%level_data%v%VALUE, &
                       flag, robsmissing)
          satob(i)%level_data%u%flag = flag
          satob(i)%level_data%v%flag = flag

        END DO

        PRINT *, 'get the u v from dd ff'
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        RETURN
      END SUBROUTINE satob_ddff_uv

    END MODULE ObsSatob_m
