!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Satob
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang
! VERSION           : V 0.0
! HISTORY           : Origionally from RCNMP
!                   : Modified by Yuanfu Xie for turning off debugging info 2022-11-27
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, @GBA-MWF, Shenzhen for adding obs
!     forward, tangent and adjoint operators.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a Satob data structure.
MODULE Satob_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE NMLRead_m, ONLY: namelist_read
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
  USE parameters_m, ONLY: degree2radian, missing, KelvinsT0
  USE Meteoro_constants, ONLY: threshold
  USE Meteoro_type_define
  USE extreme_value_m
  USE GtsQC_m
  USE qc_flag_set_m
  USE interp_m
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: Satob_t
    TYPE(type_satob), ALLOCATABLE  :: satob(:)
    !INTEGER            :: num_satob
    INTEGER(i_kind)    :: vLevel
    LOGICAL            :: logic_satob_qc
    REAL(r_kind), ALLOCATABLE  :: domain_latlon(:)
    REAL(r_kind), ALLOCATABLE  :: stn_id(:)
    REAL(r_kind), ALLOCATABLE  :: QI(:)

  CONTAINS
    PROCEDURE :: ObsInitial => satob_setup ! satobInitial
    PROCEDURE :: ObsIngest => satobIngest
    PROCEDURE :: ObsForward => satobForward
    PROCEDURE :: ObsTangent => satobTangent
    PROCEDURE :: ObsAdjoint => satobAdjoint
    PROCEDURE :: ObsQC => satob_qc    ! satobQC
    PROCEDURE :: satob_setup
    PROCEDURE :: satob_read
    PROCEDURE :: satob_qc
    PROCEDURE :: satob_assign
    PROCEDURE :: satob_err

    ! Yuanfu Xie adds a routine filling cloud drift obs heights with the background heights:
    PROCEDURE :: backHeights
  END TYPE Satob_t

CONTAINS

  SUBROUTINE satob_setup(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(Satob_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: fileNames(:)
    CHARACTER(len=200)  :: obsFileList_Satob, file

    INTEGER :: i, datetime, istatus
    INTEGER, ALLOCATABLE :: timeArray(:)

    ! read the obs satob files names
    ! CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", obsFileDir)
    ! CALL namelist_read(TRIM(configFile), "obsFileList_Satob", obsFileList_Satob)
    ! obsFileDir='/sources/data/satob202206/'
    ! istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Satob', obsFileList_Satob)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'input_dir_CloudWind', obsFileDir)
    IF (istatus .NE. 0) THEN
      PRINT *, 'It seems you are running cloud wind but no providing the data file directory. Check and rerun!'
      STOP
    END IF
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Satob', this%fileNames)
    IF (istatus .NE. 0) THEN
      PRINT *, 'It seems you are running cloud wind but no providing the data files. Check and rerun!'
      STOP
    END IF

    ! file=TRIM(obsFileDir)//'/'//TRIM(obsFileList_Satob)
    ! print *,"obsFileList_Satob: ",TRIM(file)
    ! OPEN(10,file=trim(file))
    ! CALL LineNumber(this%numfiles,10)
    ! REWIND(10)
    ! ALLOCATE(fileNames(this%numfiles))
    ! ALLOCATE(this%fileNames(this%numfiles))
    this%numfiles = UBOUND(this%fileNames, 1)
    PRINT *, "numfiles=", this%numfiles

    DO i = 1, this%numfiles
      ! READ(10,*) fileNames(i)
      this%fileNames(i) = TRIM(obsFileDir)//'/'//TRIM(this%fileNames(i))
      PRINT *, i, " this%fileNames(i)=", TRIM(this%fileNames(i))
    END DO
    CLOSE (10)

    ! Read in an interpolation option:
    !CALL namelist_read(TRIM(configFile), "interpolation",this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)
    IF (istatus .EQ. 0) THEN
      IF (this%interpolation .EQ. 1) THEN
        WRITE (*, 21)
21      FORMAT('Slint interpolation package is used in this analysis')
      ELSE
        WRITE (*, 22)
22      FORMAT('Bilinear interpolation scheme is used in the analysis')
      END IF
    END IF

    ! Read in an option:
    !CALL namelist_read(TRIM(configFile), "logic_satob_qc",this%logic_satob_qc)
    istatus = yaml_get_var(TRIM(configFile), 'ObsFlag_satob', 'logic_satob_qc', this%logic_satob_qc)
    ! CALL namelist_read(TRIM(configFile), "interpolation",this%interpolation)
    !----------------------------------------------------------------------------
    !CALL namelist_read(TRIM(configFile), "datetime",this%datetime)!----------------------------------------------------------------------------
    istatus = yaml_get_var(TRIM(configFile), 'analysis_para', 'end_time', timeArray)
    WRITE (*, 1) timeArray(1:UBOUND(timeArray, 1))
1   FORMAT('END TIME from yaml: ', 10I6)

    ! datetime = 2022060100
    ! this%yyyy0 = datetime/1000000
    ! this%mon0 = mod(datetime,1000000)/10000
    ! this%dd0 = mod(datetime,10000)/100
    ! this%hh0 = mod(datetime,100)

    ! read the var names of obs files in configFile
    ALLOCATE (this%domain_latlon(4))
    ! CALL namelist_read(TRIM(configFile), 'vLevel', this%vLevel)
    istatus = yaml_get_var(TRIM(configFile), 'modelState', 'vLevel', this%vLevel)

    ! Covered this hard coded info by Yuanfu Xie 2022-11-28
    this%domain_latlon(1) = 18.0
    this%domain_latlon(2) = 27.0
    this%domain_latlon(3) = 107.0
    this%domain_latlon(4) = 109.0

    !CALL namelist_read(TRIM(configFile), 'domain_latlon', this%domain_latlon)
    ! print *,'domain:',this%domain_latlon
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
    this%qcThreshold(1) = 20.0D0
    this%qcThreshold(2) = 20.0D0
    this%qcThreshold(3) = 10.0D0

  END SUBROUTINE satob_setup

  SUBROUTINE satobIngest(this, X)

    IMPLICIT NONE

    CLASS(Satob_t) :: this
    TYPE(State_t) :: X

    ! Wrapper: containing Dr. Zhang Hua's data read, QC and assign plus Yuanfu Xie's background height fill:
    PRINT *, 'Cloud Drift Wind at G: ', X%sg%glevel, X%sg%mpddInfo_sg%myrank
    CALL this%satob_read(X)
    PRINT *, 'step 1: satob data is read in with numObs: ', this%numObs
    CALL this%satob_qc
    CALL this%satob_assign

    IF (this%numObs .GT. 0) THEN
      CALL this%backHeights(X)
      WRITE (*, 2) this%numObs, X%sg%glevel, X%sg%mpddInfo_sg%myrank
2     FORMAT(I6, ' total Cloud Drift Wind processed at G: ', I2, ' pc: ', I2)
    ELSE
      WRITE (*, 1) X%mpddGlob%myrank
1     FORMAT('satobIngest: No cloud drift wind is ingested at proc: ', I3, /, ' check your files and rerun if needed.')
    END IF
  END SUBROUTINE satobIngest

  SUBROUTINE satob_read(this, X)

    IMPLICIT NONE

    CLASS(Satob_t) :: this
    TYPE(State_t) :: X
    ! TYPE(type_satob),  INTENT(inout)   :: satob(1:max_satob) ! observation data
    ! INTEGER ,          INTENT(inout)   :: num_satob

    ! local variables
    ! Yuanfu Xie adds an time array for removing data not in analysis window:
    INTEGER(i_kind) :: itime(6), i4time

    INTEGER :: nf_op, nc_id, var_id, i, num_obs, num_dims, nf90_max_va_dims
    INTEGER :: ifile, num_obs_total, num_dataFile, ind_obsin_beg, ind_obsin_end, num_satob
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
    INTEGER :: sitenum, ios, ilevels, imm, isite, iobs, tmp
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, p, rid, hgh, td, dd, ff, ttd, lat_d, lon_d
    !REAL(r_kind) :: missingValue = 999999.0D0
    !------------------------------------------------

    CHARACTER(len=3)            :: RTYPE
    INTEGER                     :: yyyy, mon, dd1, hh
    INTEGER                     :: num_head, num_dat, num_qcm, num_shift
    INTEGER                     :: nlev ! the number of levels in one station (in)
    INTEGER                     :: k
    REAL                        :: id, T  !content
    REAL                        :: flag_p, flag_T, flag_dd, flag_ff
    REAL                        :: a, b
    INTEGER                     :: step       ! flag position
    REAL                        :: flag
    INTEGER                     :: i1
    REAL                        :: stn_id, lat, lon
    REAL                        :: tkminlat, tkmaxlat, tkminlon, tkmaxlon
    REAL                        :: instru_code
    INTEGER                     :: status
!----- add whitelist -------------------------
    INTEGER                   :: num_whtlst, iend, unit_wht, unit_in, whitelist, whtid
    LOGICAL                   :: alive
    CHARACTER(LEN=120)        :: filewht, file_in_satob, satinfo_file
    INTEGER, ALLOCATABLE      :: idflag(:)
    INTEGER                   :: iflag, ilatlon, icode, ipres

    ! Yuanfu Xie added a satob debugging flag for turning off some printing info:
!#define SATOB_CHECK

    ! Yuanfu Xie assigns the last element of itime to 0:
    itime(6) = 0

    this%domain_latlon(1) = 18.0
    this%domain_latlon(2) = 27.0
    this%domain_latlon(3) = 107.0
    this%domain_latlon(4) = 109.0
    PRINT *, 'Inside read satob Xlatlon: ', X%sg%maxLatGlob / degree2radian, &
      X%sg%minLatGlob / degree2radian, X%sg%maxLonGlob / degree2radian, X%sg%minlonGlob / degree2radian

    ALLOCATE (idflag(2000))
!--------------------------------------------------------------------
!   read list
!---------------------------------------------------------------------
    whitelist = 0
    num_whtlst = 0
    unit_wht = 60
    idflag = 0

    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    !filewht='/sources/data/static/satob_list'
    filewht = TRIM(satinfo_file)//'/'//'satob_list'
    PRINT *, 'White file: ', TRIM(filewht)
    INQUIRE (file=TRIM(filewht), exist=alive)
    IF (.NOT. alive) THEN
      WRITE (0, *) ' Error: no satob whitelist file!'
      num_satob = 0
      RETURN
    ELSE
      whitelist = 1
      OPEN (unit_wht, FILE=TRIM(filewht))
      DO WHILE (.TRUE.)
        READ (unit_wht, *, iostat=iend) whtid
        IF (iend /= 0) EXIT
        idflag(whtid) = 1
        num_whtlst = num_whtlst + 1
      END DO
      CLOSE (unit_wht)
    END IF
    PRINT *, 'read satob whitelist over', num_whtlst

!--------------------------------------------------------------------
!   domain
!---------------------------------------------------------------------
    tkminlat = this%domain_latlon(1)
    tkmaxlat = this%domain_latlon(2)
    tkminlon = this%domain_latlon(3)
    tkmaxlon = this%domain_latlon(4)
    PRINT *, 'domain:', tkminlat, tkmaxlat, tkminlon, tkmaxlon
!--------------------------------------------------------------------
! 1.0  read the date, actual data length and number of obs. stations
!---------------------------------------------------------------------
    unit_in = 61
    num_dataFile = this%numFiles
    ALLOCATE (this%satob(max_satob)) ! max_satob is defined in meteoro_constants.F90
    ALLOCATE (num_sitein(num_dataFile))
    ALLOCATE (num_obsin(num_dataFile))
    num_obsin(:) = 0
    i = 0
    DO ifile = 1, num_dataFile
!  DO ifile=1, 1
      file_in_satob = this%fileNames(ifile)
      PRINT *, 'step1 begin read the satob data', ifile, TRIM(this%fileNames(ifile))
      !OPEN(unit_in,file=trim(this%fileNames(ifile)),FORM='FORMATTED',iostat=status,status='old',blocksize=1048576,buffered='yes',buffercount=1)
      OPEN (unit_in, file=TRIM(this%fileNames(ifile)), FORM='FORMATTED', iostat=status, status='old')
      IF (status /= 0) THEN
        PRINT *, 'file_in_satob is not exist=', TRIM(this%fileNames(ifile))
        num_satob = 0
        ! return  ! Yuanfu Xie does not think this is good for giving up other data file when one is missing
        CYCLE   ! Yuanfu Xie tries to read this next file
      END IF
      READ (unit_in, *, END=333) RTYPE, yyyy, mon, dd1, hh, num_satob, num_head, num_dat, num_qcm
      PRINT *, 'RTYPE: ', RTYPE, yyyy, mon, dd1, hh, num_satob, num_head, num_dat, num_qcm

!--------------------------------------------------------------------
!   check the input file date  ! use obs_size_constants
!--------------------------------------------------------------------
      !if(yyyy/=this%yyyy0)then
      !    print*,'wrong input', yyyy, 'please check file=',file_in_satob
      !     return
      !endif

      !if(mon/=this%mon0)then
      !    print*,'wrong input', mon, 'please check file=',file_in_satob
      !     return
      !endif

      ! if(dd1/=this%dd0)then
      !   print*,'wrong input', dd , 'please check file=',file_in_satob
      !   return
      ! endif

      !if(hh/=this%hh0)then
      !      print*,'wrong input', hh, 'please check file=',file_in_satob
      !     return
      ! endif
!--------------------------------------------------------------------
!    read out the data
!--------------------------------------------------------------------
      PRINT *, 'original num_satob', num_satob
      iflag = 0
      ilatlon = 0
      icode = 0
      ipres = 0
      DO i1 = 1, num_satob
        READ (unit_in, *, END=333) stn_id, lat, lon, a, b, instru_code
        READ (unit_in, '(8f10.2)', END=333) p, T, dd, ff, flag_p, flag_T, flag_dd, flag_ff
#ifdef SATOB_CHECK
        PRINT *, i1, 'stn_id: ', stn_id, lat, lon, a, b, instru_code
        PRINT *, 'p/T/dd/ff: ', p, T, dd, ff, flag_p, flag_T, flag_dd, flag_ff
#endif

        IF (idflag(INT(stn_id)) .EQ. 0) THEN
#ifdef SATOB_CHECK
          ! print *,'idflag(int(stn_id)): ',stn_id,idflag(int(stn_id))
          PRINT *, 'Station is blocked or white out: ', INT(stn_id)
#endif
          CYCLE
        END IF
        IF (lat .GT. 89.0 .OR. lat .LT. -89.0) CYCLE
        IF (lon .GT. 180.0 .OR. lon .LT. -180.0) CYCLE

!   read data with QI >80 /domain
        IF (flag_ff > 100.) CYCLE
        IF (flag_ff < 80.) CYCLE
        iflag = iflag + 1
        !print *,ifile,i1,'iflag= ',iflag
        !  if(instru_code==4.0)cycle
        !   if(instru_code==5.0)cycle
        !   if(instru_code==6.0)cycle
        !   if(instru_code>2)cycle

        IF (stn_id == 173.0) THEN
          IF (instru_code > 2 .AND. flag_ff < 85.) CYCLE
        ELSE
          IF (instru_code > 2) CYCLE
        END IF
        icode = icode + 1
        ! print *,ifile,i1,'icode= ',icode
        IF (p > 950.0) CYCLE
        IF (p < 150.0) CYCLE
        ipres = ipres + 1
#ifdef SATOB_CHECK
        PRINT *, ifile, i1, 'iflag/icode/ipres= ', iflag, icode, ipres, p, t
#endif
        !   remove out the station with no obs data
        IF (ABS(p - inobsmissing) < threshold) CYCLE
        IF (ABS(dd - inobsmissing) < threshold .OR. ABS(ff - inobsmissing) < threshold) CYCLE
        !------------------------------------------------------------------
        ! remove the data acord the EC data selection added by tianwh 20140504
        ! Why? questioned by Yuanfu Xie 2022-11-27
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
          PRINT *, 'the max_satob is too small num_satob=', max_satob, ' Stopping read data any further...'
          i = max_satob
          EXIT
        END IF
#ifdef SATOB_CHECK
        PRINT *, i, 'domain check lat:', lat, tkminlat, tkmaxlat, 'lon:', lon, tkminlon, tkmaxlon
#endif
        !  IF( (lat<=tkmaxlat).and.(lat>=tkminlat) ) THEN
        !  IF( (lon<=tkmaxlon).and.(lon>=tkminlon) ) THEN
        !IF ((lat<=tkmaxlat).and.(lat>=tkminlat) .AND. &
        !    (lon<=tkmaxlon).and.(lon>=tkminlon)) THEN
        IF ((lat .LE. X%sg%maxLatGlob / degree2radian) .AND. &
            (lat .GE. X%sg%minLatGlob / degree2radian) .AND. &
            (lon .LE. X%sg%maxLonGlob / degree2radian) .AND. &
            (lon .GE. X%sg%minLonGlob / degree2radian)) THEN

          !-------------------------------------------
          ! observation time: checking if it is inside the analysis window: Yuanfu Xie 2022-11-30
          !-------------------------------------------
          ! this%satob(i)%header_info%itime(1)=yyyy            ! year
          ! this%satob(i)%header_info%itime(2)=mon             ! month
          ! this%satob(i)%header_info%itime(3)=int(a)          ! dd
          ! this%satob(i)%header_info%itime(4)=int(b/100)      ! hh
          !  if(dd1 == 1 .and.a >= 28) then
          !   this%satob(i)%header_info%itime(2)=mon -1
          !      if(mon == 1) then
          !       this%satob(i)%header_info%itime(2)= 12
          !       this%satob(i)%header_info%itime(1)=yyyy -1
          !       end if
          !  end if

          !  this%satob(i)%header_info%itime(5)=mod(int(b),100) ! mm

          itime(1) = yyyy            ! year
          itime(2) = mon             ! month
          itime(3) = INT(a)          ! dd
          itime(4) = INT(b / 100)      ! hh
          IF (dd1 == 1 .AND. a >= 28) THEN
            itime(2) = mon - 1
            IF (mon == 1) THEN
              itime(2) = 12
              itime(1) = yyyy - 1
            END IF
          END IF

          itime(5) = MOD(INT(b), 100) ! mm
          CALL Time_GMT_to_Unix(itime, i4time)
#ifdef SATOB_CHECK
          !print *,i,'itime:',itime,i4time-X%sg%tt(1),X%sg%tt(X%sg%tSlots)-i4time
          WRITE (*, 1) i, ifile, itime, i4time - INT(X%sg%tt(1)), INT(X%sg%tt(X%sg%tSlots)) - i4time, i4time, INT(X%sg%tt(1)), &
            X%sg%mpddInfo_sg%myrank
1         FORMAT('itime: ', I5, I3, I5, 5I3, ' t-a: ', I6, ' b-t: ', I6, ' i4:', 2I12, ' pc', I2)
#endif

          ! Check if the obs time is inside the analysis window:
          IF (i4time - INT(X%sg%tt(1)) .LT. 0 .OR. INT(X%sg%tt(X%sg%tSlots)) - i4time .LT. 0) CYCLE

          i = i + 1
          this%satob(i)%header_info%stn_id = stn_id
          this%satob(i)%header_info%lat = lat
          this%satob(i)%header_info%lon = lon
          this%satob(i)%header_info%instru_code = instru_code
          this%satob(i)%header_info%QI = flag_dd

          ! Save obs time in Unix form: setting all elements to i4time
          this%satob(i)%header_info%itime = i4time

#ifdef SATOB_CHECK
          PRINT *, i, 'In domain:', stn_id, lat, lon, instru_code, instru_code
#endif
!-------------------------------------------
!  change  999999.0 to -888888.0
!-------------------------------------------
          ! Yuanfu Xie added +1.0D1 avoiding the round off effect as p, T, dd, ff are single precision variables.
          IF (ABS(p - inobsmissing) < threshold + 1.0D1) p = missing
          IF (ABS(T - inobsmissing) < threshold + 1.0D1) THEN
            T = missing
          ELSE
            T = T + KelvinsT0   ! Change temperature from celsius to kelvin
          END IF
          IF (ABS(dd - inobsmissing) < threshold + 1.0D1) dd = missing
          IF (ABS(ff - inobsmissing) < threshold + 1.0D1) ff = missing

          !-------------------------------------------
          ! set the initial QC flag
          !-------------------------------------------
          this%satob(i)%level_data%p%flag = flag_missing
          this%satob(i)%level_data%T%flag = flag_missing
          this%satob(i)%level_data%dd%flag = flag_missing
          this%satob(i)%level_data%ff%flag = flag_missing
          !-----------------------------------
          this%satob(i)%level_data%p%VALUE = p
          this%satob(i)%level_data%T%VALUE = T
          this%satob(i)%level_data%dd%VALUE = dd
          this%satob(i)%level_data%ff%VALUE = ff

          !      call flag_replace(satob(i)%level_data%p%flag ,step, int(flag_p))
          !      call flag_replace(satob(i)%level_data%T%flag,step , int(flag_T))
          !      call flag_replace(satob(i)%level_data%dd%flag,step , int(flag_dd))
          !      call flag_replace(satob(i)%level_data%ff%flag,step , int(flag_ff))
          !-------------------------------------------------------
          ! ENDIF
        END IF

        !-------------------------------------------
      END DO  ! the loop of the station
    END DO
333 CONTINUE
    this%numObs = i

    DEALLOCATE (idflag)
    CLOSE (unit_in)

    RETURN

  END SUBROUTINE satob_read

  FUNCTION satobForward(this, X, O) RESULT(Y)

    IMPLICIT NONE

    CLASS(Satob_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j, k

    Y = O
    DO i = 1, UBOUND(Y%ObsFields, 1)
      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
            Y%ObsFields(i)%values(k) = X%fields(j)%Get_Value(Y%ObsFields(i)%idx(k))
          END DO
        END IF
      END DO
    END DO
  END FUNCTION satobForward

  FUNCTION satobTangent(this, dX, X) RESULT(Y)

    IMPLICIT NONE

    CLASS(Satob_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'satob: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION satobTangent

  FUNCTION satobAdjoint(this, dY, X) RESULT(dX)

    IMPLICIT NONE

    CLASS(Satob_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    ! Local variables:
    INTEGER(i_kind) :: i, j, k

    dX = X%zeroCopy()

    DO i = 1, UBOUND(dY%ObsFields, 1)
      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF ((TRIM(dY%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(dY%ObsFields(i)%idx, 1), UBOUND(dY%ObsFields(i)%idx, 1)
            CALL dX%fields(j)%Add_Value(dY%ObsFields(i)%idx(k), dY%ObsFields(i)%values(k))
          END DO
        END IF
      END DO
    END DO
  END FUNCTION satobAdjoint

  SUBROUTINE satob_qc(this)
    ! USE extreme_value_m
    !use   index_satob
    IMPLICIT NONE
    CLASS(Satob_t) :: this
    !TYPE(type_satob),  INTENT(inout)   :: satob(1:max_satob) ! observation data
    !INTEGER ,          INTENT(inout)   :: num_satob
    !TYPE(extreme_value_t)         :: extreme
    ! TYPE (type_ob)              :: yo
    ! TYPE (type_xb)              :: xb ! Background structures

    INTEGER                     :: i, num_stn
    REAL :: stn_h, stn_p, c1, stn_id1

    !---------------------------
    !  read out the satob data
    !---------------------------

    IF (this%numObs == 0) this%logic_satob_qc = .FALSE.

    PRINT *, 'Satob QC ... ', this%numObs
    IF (this%logic_satob_qc) THEN
      !----------------------------------------
      ! begin the QC
      !----------------------------------------
      CALL read_temp_extreme
      PRINT *, 'read_temp_extreme', this%numObs
      CALL satob_extreme_check(this%satob, this%numObs)
      PRINT *, 'satob_extreme_check'
      !-----------------------------------------------

      CALL satob_inter_consist_check(this%satob, this%numObs)
      PRINT *, 'satob_inter_consist_check'
      !------------------------------------------------
      ! get the u v from dd ff !

      CALL satob_ddff_uv(this%satob, this%numObs)
      PRINT *, 'satob_ddff_uv'
      !------------------------------------------------

      ! if(logic_first_guess_check)then
      !  call satob_bg_check(yo%satob,yo%num_satob,xb)
      ! print*,'the satob_first_guess_check is over'
      ! endif

      !-------------------------------------

      CALL satob_complex_dma(this%satob, this%numObs)
      PRINT *, 'satob_complex_dma is over'
      !----------------------------------------------
    ELSE
      PRINT *, 'logic_satob_qc is FALSE, so QC is skipped'
    END IF
    RETURN

  END SUBROUTINE satob_qc

  SUBROUTINE satob_assign(this)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Description:
    !
    !  purpose:  write out the input observation of satob
    !
    !  method :  write out the data after QC
    !
    ! modified : by Yuanfu Xie for cleaning up output messages: 2022-12-01
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! use   obs_size_constants,ONLY : max_satob,max_mut_levels
    ! use   Meteoro_constants
    ! use   obs_type_define
    ! use   namelist_define
    ! use   qc_flag_set
    !------------------------------------------------
    IMPLICIT NONE

    CLASS(Satob_t) :: this
    !TYPE (type_satob)           :: satob(max_satob)
    !integer                     :: num_satob
    CHARACTER(len=5)            :: RTYPE
    INTEGER                     :: i, k
    INTEGER                     :: yyyy, mon, dd, hh
    REAL                        :: u, v
    INTEGER                     :: unit_out_satob
    INTEGER                     :: m, num, num0
    INTEGER                     :: kk, ii, ll, mm, nn, gg
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), obsErrs(:, :), olatlon(:, :)
    REAL(r_kind), ALLOCATABLE :: obsHght(:), land(:), obspres(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    INTEGER(i_kind), ALLOCATABLE :: num_sitein(:), num_obsin(:), stn_id(:), QI(:)
    CHARACTER(LEN=20), ALLOCATABLE :: stNames(:)
    CHARACTER(LEN=10), ALLOCATABLE :: timestr(:)
    CHARACTER(LEN=30), ALLOCATABLE :: Id_time(:), Id_time_unique(:)
    !---------------------------------------------------------
    !-1.0 write out the data in satob GRAPES var can read direct
    !---------------------------------------------------------
    PRINT *, 'begin write out the satob data', this%numObs, this%numVars

    IF (this%numObs .EQ. 0) THEN
      WRITE (*, 1) this%numObs, this%numVars
1     FORMAT('No cloud drift wind data...', 2I3)
      RETURN
    END IF

    ! Found cloud drift wind data:
    this%obsType = 'CDW'

    ALLOCATE (obsData(this%numObs, this%numVars))
    !  ALLOCATE(obsErrs(this%numVars))
    ALLOCATE (olatlon(2, this%numObs))
    ALLOCATE (obsTime(this%numObs))
    ALLOCATE (obsHght(this%numObs))
    ALLOCATE (stn_id(this%numObs))
    ALLOCATE (QI(this%numObs))
    ALLOCATE (timestr(this%numObs))
    ALLOCATE (Id_time(this%numObs))
    ALLOCATE (obspres(this%numObs))

    num = 0
    num0 = 0
    DO i = 1, this%numObs

      IF (ABS(this%satob(i)%level_data%p%VALUE - missing) < threshold .OR. MOD(this%satob(i)%level_data%p%flag, 10) == 2) CYCLE
      IF (MOD(this%satob(i)%level_data%p%flag, 10) == 1) CYCLE
      IF (this%satob(i)%header_info%position == 2) CYCLE
      CALL flag_decompse(this%satob(i)%level_data%u%flag, 9, m)
      IF (m == 1) num0 = num0 + 1
      IF (MOD(this%satob(i)%level_data%u%flag, 10) == 2 .OR. m == 1) CYCLE
      IF (MOD(this%satob(i)%level_data%v%flag, 10) == 2 .OR. m == 1) CYCLE

      IF (MOD(this%satob(i)%level_data%u%flag, 10) == 1) CYCLE
      IF (MOD(this%satob(i)%level_data%v%flag, 10) == 1) CYCLE
      IF (ABS(this%satob(i)%level_data%u%VALUE - missing) < threshold) CYCLE
      IF (ABS(this%satob(i)%level_data%v%VALUE - missing) < threshold) CYCLE
      !satob(i)%level_data%elevation%value = 900000.000
      !satob(i)%level_data%ps%value = missing
      this%satob(i)%level_data%phi%VALUE = missing
      this%satob(i)%level_data%rh%VALUE = missing
      this%satob(i)%header_info%levels = 1
      this%satob(i)%header_info%stn_id = this%satob(i)%header_info%stn_id * 100 + this%satob(i)%header_info%instru_code

      num = num + 1
      IF (this%satob(i)%header_info%lon < 0.) THEN
        this%satob(i)%header_info%lon = this%satob(i)%header_info%lon + 360.
      END IF

      obsData(num, 1) = this%satob(i)%level_data%u%VALUE
      obsData(num, 2) = this%satob(i)%level_data%v%VALUE
      obsData(num, 3) = this%satob(i)%level_data%T%VALUE
      obsData(num, 4) = this%satob(i)%level_data%p%VALUE
      ! obsErrs(num, :)
      olatlon(1, num) = this%satob(i)%header_info%lat * degree2radian   ! convert degree to radian. - JMPang
      olatlon(2, num) = this%satob(i)%header_info%lon * degree2radian   ! convert degree to radian. - JMPang

      ! Yuanfu Xie converts the obs time to Unix time:
      obsTime(num) = this%satob(i)%header_info%itime(1) ! Note itime(1) is already set to i4time in satob_read
      !obsTime(num) = this%satob(i)%header_info%itime(1)*1000000 + this%satob(i)%header_info%itime(2)*10000 +      &
      !                this%satob(i)%header_info%itime(3)*100 + this%satob(i)%header_info%itime(4)
#ifdef SATOB_CHECK
      PRINT *, i, num, 'Satob_time:', obstime(num), this%satob(i)%header_info%itime(1) * 1000000, this%satob(i)%header_info%itime(2) * 10000, &
        this%satob(i)%header_info%itime(3) * 100, this%satob(i)%header_info%itime(4)
      PRINT *, 'Temperature assigned: ', num, this%satob(i)%level_data%T%VALUE, this%satob(i)%level_data%u%VALUE
#endif
      obsHght(num) = missing  !900000.000 ! Change it to missing value later
      stn_id(num) = this%satob(i)%header_info%stn_id
      QI(num) = this%satob(i)%header_info%QI
      CALL Satob_err(this, num)

    END DO ! the loop of the station number
    CLOSE (unit_out_satob)
    !------------------------------------------------------

    ! Yuanfu Xie added a check for output obs: safe guard for avoiding allocate zero element arrays
    IF (num .GT. 0) THEN
      WRITE (*, 2) num, num0
2     FORMAT('The number of cloud wind obs =', I6, ' number of flag out :', I5)

      ALLOCATE (this%obsData(num, this%numVars))
      ALLOCATE (this%obsErrs(num, this%numVars))
      ALLOCATE (this%olatlon(2, num))
      ALLOCATE (this%obsTime(num))
      ALLOCATE (this%obsHght(num))
      ALLOCATE (this%stn_id(num))
      ALLOCATE (this%QI(num), this%stNames(num))

      this%numObs = num
      !CALL Satob_err(this)

      DO i = 1, this%numObs
        this%obsData(i, 1) = obsData(i, 1)
        this%obsData(i, 2) = obsData(i, 2)
        this%obsData(i, 3) = obsData(i, 3)
        this%obsData(i, 4) = obsData(i, 4) * 100.0D0 ! Yuanfu changes the pressure ob into Pascals
        this%obsErrs(i, 1) = this%satob(i)%level_data%u%error
        this%obsErrs(i, 2) = this%satob(i)%level_data%v%error
        this%olatlon(1, i) = olatlon(1, i)   ! convert degree to radian. - JMPang
        this%olatlon(2, i) = olatlon(2, i)   ! convert degree to radian. - JMPang

        ! Yuanfu Xie corrected the obsTime and
        this%obsTime(i) = obsTime(i)
        this%obsHght(i) = obsHght(i)

        ! Yuanfu Xie temporarily gave a station name for cloud drift wind:
        WRITE (this%stNames(i) (1:10), 11) i
11      FORMAT('CDW', I7.7)
#ifdef SATOB_CHECK
        PRINT *, i, this%olatlon(1, i), this%olatlon(2, i), 'p:', this%obsData(i, 4), 'u:', this%obsData(i, 1), this%obsErrs(i, 1), 'v:', this%obsData(i, 2), this%obsErrs(i, 2)
#endif
      END DO

      DO i = 1, this%numVars
        WRITE (*, 3) TRIM(this%varNames(i)), &
          MAXVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. missing)), &
          MINVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. missing)), &
          COUNT(this%obsData(:, i) .NE. missing)
3       FORMAT("Cloud Wind: max and min of ", A, " in Uni obs: ", 2D12.4, ' NumOb used: ', I5)
      END DO

    ELSE
      WRITE (*, 4) this%numObs
4     FORMAT('satob_assign: the number of Cloud drifted wind is', I2)
    END IF

    ! Deallocate the temporary memory:
    DEALLOCATE (obsData, olatlon, obsTime, obsHght, stn_id, QI, &
                timestr, Id_time, obspres)

  END SUBROUTINE satob_assign

  SUBROUTINE Satob_err(this, istn)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Description:
    ! Allign observation error according the pressure.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IMPLICIT NONE
    !*Subroutine arguments:
    CLASS(Satob_t) :: this
    ! INTEGER, INTENT(in) :: num_sound
    ! TYPE (type_satob), INTENT(inout) :: sound(num_sound)
    !*Local arguments:
    INTEGER, INTENT(in) :: istn ! loop counter: loop of every obs station
    !-----------------------------------------------------------------------------
    ! [1.0] Set up satob obs. error structure:
    !-----------------------------------------------------------------------------
    ! DO istn=1, this%numObs
    !1-- for P, U, V
    IF (this%satob(istn)%level_data%p%VALUE >= 925.0) THEN ! for 1000 hPa
      this%satob(istn)%level_data%u%error = U_err_satob(1)
      this%satob(istn)%level_data%v%error = V_err_satob(1)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 925.0 .AND. & ! for 850 hPa
             this%satob(istn)%level_data%p%VALUE >= 775.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(2)
      this%satob(istn)%level_data%v%error = V_err_satob(2)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 775.0 .AND. & ! for 700 hPa
             this%satob(istn)%level_data%p%VALUE >= 600.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(3)
      this%satob(istn)%level_data%v%error = V_err_satob(3)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 600.0 .AND. & ! for 500 hPa
             this%satob(istn)%level_data%p%VALUE >= 450.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(4)
      this%satob(istn)%level_data%v%error = V_err_satob(4)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 450.0 .AND. & ! for 400 hPa
             this%satob(istn)%level_data%p%VALUE >= 350.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(5)
      this%satob(istn)%level_data%v%error = V_err_satob(5)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 350.0 .AND. & ! for 300 hPa
             this%satob(istn)%level_data%p%VALUE >= 275.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(6)
      this%satob(istn)%level_data%v%error = V_err_satob(6)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 275.0 .AND. & ! for 250 hPa
             this%satob(istn)%level_data%p%VALUE >= 225.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(7)
      this%satob(istn)%level_data%v%error = V_err_satob(7)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 225.0 .AND. & ! for 200 hPa
             this%satob(istn)%level_data%p%VALUE >= 175.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(8)
      this%satob(istn)%level_data%v%error = V_err_satob(8)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 175.0 .AND. & ! for 150 hPa
             this%satob(istn)%level_data%p%VALUE >= 125.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(9)
      this%satob(istn)%level_data%v%error = V_err_satob(9)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 125.0 .AND. & ! for 100 hPa
             this%satob(istn)%level_data%p%VALUE >= 85.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(10)
      this%satob(istn)%level_data%v%error = V_err_satob(10)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 85.0 .AND. & ! for 70 hPa
             this%satob(istn)%level_data%p%VALUE >= 60.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(11)
      this%satob(istn)%level_data%v%error = V_err_satob(11)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 60.0 .AND. & ! for 50 hPa
             this%satob(istn)%level_data%p%VALUE >= 40.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(12)
      this%satob(istn)%level_data%v%error = V_err_satob(12)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 40.0 .AND. & ! for 30 hPa
             this%satob(istn)%level_data%p%VALUE >= 25.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(13)
      this%satob(istn)%level_data%v%error = V_err_satob(13)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 25.0 .AND. & ! for 20 hPa
             this%satob(istn)%level_data%p%VALUE >= 15.0) THEN
      this%satob(istn)%level_data%u%error = U_err_satob(14)
      this%satob(istn)%level_data%v%error = V_err_satob(14)
    ELSE IF (this%satob(istn)%level_data%p%VALUE < 15.0) THEN ! for 10 hPa
      this%satob(istn)%level_data%u%error = U_err_satob(15)
      this%satob(istn)%level_data%v%error = V_err_satob(15)
    END IF
    !END DO
  END SUBROUTINE Satob_err

  SUBROUTINE satob_extreme_check(satob, num_satob)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Description:
    ! Climate extremet check on temperature,wind
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! USE qc_flag_set_m

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
      IF (ABS(fff - missing) > threshold) THEN
        CALL w_extreme_qc(satob(num_stn)%level_data%p%VALUE, satob(num_stn)%header_info%lat, &
                          fff, satob(num_stn)%level_data%ff%flag)
        satob(num_stn)%level_data%dd%flag = satob(num_stn)%level_data%ff%flag
      END IF

      ! pressure extreme check
      p = satob(num_stn)%level_data%p%VALUE
      IF (ABS(p - missing) > threshold) THEN
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
    IF (ABS(dd - missing) < threshold .AND. ABS(ff - missing) < threshold) RETURN

    CALL wind_consist_check(dd, ff, iflag)

    CALL flag_replace(satob(num_stn)%level_data%dd%flag, inter_consist_flag, iflag)

    CALL flag_replace(satob(num_stn)%level_data%ff%flag, inter_consist_flag, iflag)

    IF (iflag == 2) THEN
      CALL flag_replace(satob(num_stn)%level_data%dd%flag, DMA_flag, iflag)

      CALL flag_replace(satob(num_stn)%level_data%ff%flag, DMA_flag, iflag)

      PRINT *, 'dd and fff inconsist', satob(num_stn)%header_info%stn_id, dd, ff
      ! print*,satob(num_stn)%header_info%stn_id,satob(num_stn)%level_data%p%value,u,v
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
                   flag, missing)
      satob(i)%level_data%u%flag = flag
      satob(i)%level_data%v%flag = flag

    END DO

    PRINT *, 'get the u v from dd ff'
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    RETURN
  END SUBROUTINE satob_ddff_uv

  SUBROUTINE satob_complex_dma(satob, num_satob)
    !    complex desition makeing     for satob observation data

    IMPLICIT NONE
    TYPE(type_satob), INTENT(inout) :: satob(1:max_satob) !observation data
    INTEGER                         :: num_satob
    INTEGER                         :: is, k

    DO is = 1, num_satob
      k = 1
      IF (ABS(satob(is)%level_data%p%VALUE - missing) > threshold) &
        CALL complex_dma(satob(is)%level_data%p%flag)

      IF (ABS(satob(is)%level_data%t%VALUE - missing) > threshold) &
        CALL complex_dma(satob(is)%level_data%t%flag)

      IF (ABS(satob(is)%level_data%u%VALUE - missing) > threshold) &
        CALL complex_dma(satob(is)%level_data%u%flag)
      satob(is)%level_data%v%flag = satob(is)%level_data%u%flag

    END DO
    RETURN
  END SUBROUTINE satob_complex_dma

  SUBROUTINE LineNumber(row, unint)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unint
    INTEGER            :: row

    row = 0
10  READ (unint, *, END=100)
    row = row + 1
    GOTO 10
100 CONTINUE
  END SUBROUTINE LineNumber

  ! Yuanfu Xie adds a routine for filling obs heights with the background:
  SUBROUTINE backHeights(this, X)
    USE domainCheck_m, ONLY: domainCheck_t
    IMPLICIT NONE
    CLASS(Satob_t) :: this
    TYPE(State_t) :: X

    ! Local variables:
    INTEGER(i_kind) :: nstencil, io, i, j, ip_bk, iu_bk, iv_bk, it_bk, ip_ob, iu_ob, iv_ob, it_ob, ih, it(2), ihgt, istatus
    INTEGER(i_kind), ALLOCATABLE :: idx(:, :)
    REAL(r_kind) :: tim(2), fit_min, tempH, tempP, tempU, tempV, tempT, distance, Pressure_weight
    REAL(r_kind), ALLOCATABLE :: coe(:, :), hgt(:), backPres(:, :), backHght(:, :), backUwnd(:, :), backVwnd(:, :), backTemp(:, :)
    TYPE(domainCheck_t) :: domain

    INCLUDE 'mpif.h'

    ! Debugging:
    INTEGER(i_kind) :: counter

    ! Read in the weight for pressure in height adjustment:
    istatus = yaml_get_var(TRIM(this%configFile), 'ObsFlag_satob', 'satobPresWeight', Pressure_weight)
    IF (istatus .NE. 0) Pressure_weight = 1.0D2 ! Default value

    WRITE (*, 1) this%numObs, Pressure_weight, X%sg%mpddInfo_sg%myrank, this%varNames(1:UBOUND(this%varNames, 1)) (1:5)
1   FORMAT('Entering CDW height adjustment - numObs: ', I6, ' PressureWeight = ', E10.2, ' pc: ', I2, ' varNames: ', 100A)

    ! Background pressures at obs location, a column array:
    ip_bk = MAX(X%getVarIdx('pres'), X%getVarIdx('psl'))
    iu_bk = X%getVarIdx('uwnd')
    iv_bk = X%getVarIdx('vwnd')
    it_bk = X%getVarIdx('temp')
    ip_ob = this%getVarIdx('pres')
    iu_ob = this%getVarIdx('uwnd')
    iv_ob = this%getVarIdx('vwnd')
    it_ob = this%getVarIdx('temp')

    ! If no uwnd,vwnd,temp,pres in the background or obs, we do not allow cloud drift wind into our DA!
    IF (ip_bk .EQ. 0 .OR. iu_bk .EQ. 0 .OR. iv_bk .EQ. 0 .OR. it_bk .EQ. 0 .OR. &
        ip_ob .EQ. 0 .OR. iu_ob .EQ. 0 .OR. iv_ob .EQ. 0 .OR. it_ob .EQ. 0) THEN
      this%obsHght = missing
      WRITE (*, 2) ip_bk, iu_bk, iv_bk, it_bk, ip_ob, iu_ob, iv_ob, it_ob
2     FORMAT('CDW data is dismissed due to missing variables: BK: ', 4I3, ' OB: ', 4I3)
      RETURN
    ELSE
      PRINT *, 'Starting CDW height adjustment...'
    END IF

    ALLOCATE (hgt(this%numObs), backPres(X%sg%vLevel, this%numObs), backHght(X%sg%vLevel, this%numObs), &
              backUwnd(X%sg%vLevel, this%numObs), backVwnd(X%sg%vLevel, this%numObs), &
              backTemp(X%sg%vLevel, this%numObs))
    backPres = 0.0D0
    backHght = 0.0D0
    backUWND = 0.0D0
    backVWND = 0.0D0
    backTemp = 0.0D0

    ! Spacial interpolation:
    CALL domain%horizontalIntp(X%sg, this%numObs, this%olatlon, idx, coe, this%interpolation, nstencil)

    ! Count valided CDW:
    counter = 0
    DO io = 1, this%numObs
      ! Find the time interval containing the obs time:
      tim = 0.0D0
      it = 0
      DO j = 1, X%sg%tSlots - 1
        tim(1) = X%sg%tt(j + 1) - this%obsTime(io)
        tim(2) = this%obsTime(io) - X%sg%tt(j)
        IF (tim(1) .GE. 0.0D0 .AND. tim(2) .GE. 0.0D0) THEN
          it(1) = j; it(2) = j + 1; 
          tim = tim / (tim(1) + tim(2))
          EXIT
        END IF
      END DO

      IF (MINVAL(idx(:, io)) .GT. 0 .AND. &
          MAXVAL(idx(:, io)) .LE. UBOUND(X%fields(ip_bk)%DATA, 2) .AND. it(1) .GT. 0) THEN
        counter = counter + 1
        DO j = 1, 2  ! Time interpolation
          DO i = 1, nstencil
            backPres(:, io) = backPres(:, io) + &
                              X%fields(ip_bk)%DATA(:, idx(i, io), it(j)) * coe(i, io) * tim(j)
            backUwnd(:, io) = backUwnd(:, io) + &
                              X%fields(iu_bk)%DATA(:, idx(i, io), it(j)) * coe(i, io) * tim(j)
            backVwnd(:, io) = backVwnd(:, io) + &
                              X%fields(iv_bk)%DATA(:, idx(i, io), it(j)) * coe(i, io) * tim(j)
            backTemp(:, io) = backTemp(:, io) + &
                              X%fields(it_bk)%DATA(:, idx(i, io), it(j)) * coe(i, io) * tim(j)
          END DO
        END DO

        DO i = 1, nstencil
          backHght(:, io) = backHght(:, io) + &
                            X%sg%zHght(:, idx(i, io)) * coe(i, io)
        END DO
      ELSE
        backPres(:, io) = missing
        backHght(:, io) = missing
        backUwnd(:, io) = missing
        backVwnd(:, io) = missing
        backTemp(:, io) = missing
      END IF
    END DO
    PRINT *, 'Total Satob inside subdomain: ', counter, X%sg%mpddInfo_sg%myrank

    ! Calculate the background height and pressure at obs locations:
    DO io = 1, this%numObs
      IF (MINVAL(idx(:, io)) .LE. 0 .OR. backPres(1, io) .EQ. missing .OR. ABS(this%obsData(io, it_ob) - missing) .LT. 10.0) CYCLE
      ! Find the pressure levels containing the obs:
      fit_min = 1.0D10 ! default fit value
      this%obsHght(io) = missing
      ihgt = 0

      DO i = 1, X%sg%vLevel - 1
        distance = DSQRT(Pressure_weight * (LOG(backPres(i, io)) - LOG(this%obsData(io, ip_ob)))**2 + &
                         (backUwnd(i, io) - this%obsData(io, iu_ob))**2 + &
                         (backVwnd(i, io) - this%obsData(io, iv_ob))**2 + &
                         (backTemp(i, io) - this%obsData(io, it_ob))**2)
        IF (fit_min .GT. distance) THEN
          fit_min = distance
          this%obsHght(io) = backHght(i, io)
          ihgt = i
#ifdef SATOB_CHECK
!#ifdef DEBUG
          WRITE (*, 3) i, io, DSQRT(Pressure_weight * (LOG(backPres(i, io)) - LOG(this%obsData(io, ip_ob)))**2), &
            DSQRT((backUwnd(i, io) - this%obsData(io, iu_ob))**2), &
            DSQRT((backVwnd(i, io) - this%obsData(io, iv_ob))**2), &
            DSQRT((backTemp(i, io) - this%obsData(io, it_ob))**2), &
            backPres(i, io), this%obsData(io, ip_ob), &
            fit_min, this%obsHght(io), X%sg%mpddInfo_sg%myrank
3         FORMAT('Height intp: ', I3, ' io: ', I4, ' bk: ', 4D10.2, ' PP', 2E11.4, ' fitmin:', E11.4, ' hght:', E11.4, ' pc', I2)
          PRINT *, 'PRESsure: ', i, io, LOG(backPres(i, io)), LOG(this%obsData(io, ip_ob)), X%sg%vLevel, ' pc: ', X%sg%mpddInfo_sg%myrank
          PRINT *, 'UWIND: ', i, io, backUwnd(i, io), this%obsData(io, iu_ob), ' pc: ', X%sg%mpddInfo_sg%myrank
          PRINT *, 'VWIND: ', i, io, backVwnd(i, io), this%obsData(io, iv_ob), ' pc: ', X%sg%mpddInfo_sg%myrank
          PRINT *, 'TEMPerature: ', i, io, backTemp(i, io), this%obsData(io, it_ob), ' pc: ', X%sg%mpddInfo_sg%myrank
#endif
        END IF
      END DO

#ifdef SATOB_CHECK
      IF (ihgt .GT. 0) WRITE (*, 4) io, ihgt, this%obsHght(io), backPres(ihgt, io), this%obsData(io, ip_ob), &
        X%sg%gLevel, X%sg%mpddInfo_sg%myrank
4     FORMAT('SatobHeight: ', I6, ' Hght set: ', I3, E12.4, ' Pres: ', E12.4, ' ObsPres: ', E12.4, ' G', I2, ' pc:', I2)
#endif

    END DO
    ! Yuanfu Xie added MPI_allreduce to broadcast the available heights from any process to all processors to share: 2024-08-17
    BLOCK
      REAL(r_kind) :: broadCastHeight(this%numObs)
      broadCastHeight = this%obsHght

      CALL X%sg%mpddInfo_sg%barrier()
      ! Note: MPI_allreduce and MPI_reduce have different argument list.
      ! MPI_allreduce does not need a root info. But it went through on my laptop.
      CALL MPI_allreduce(broadCastHeight, this%obsHght, this%numObs, MPI_DOUBLE, MPI_MIN, &
                         X%sg%mpddGlob%comm, X%sg%mpddGlob%ierr)
    END BLOCK

    ! Deallocate memory:
    DEALLOCATE (idx, coe, hgt, backPres, backHght, backUwnd, backVwnd, backTemp)
  END SUBROUTINE backHeights

END MODULE Satob_m
