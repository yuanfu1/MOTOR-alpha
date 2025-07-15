!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ReadWriteGIIRS
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang, Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE ReadWriteGIIRS_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE RawGIIRS_m, ONLY: RawGIIRS_t
  USE Dp_utils_m
  USE Satellite_utils_m
  USE RawSatellite_BC_m
  USE RadBase_m, ONLY: RadBase_t
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t

  IMPLICIT NONE

  TYPE, EXTENDS(RadBase_t) :: ReadWriteGIIRS_t
    CHARACTER(LEN=1024) :: inst_name='giirs', platform_name='fy4_1'
    INTEGER(i_kind) :: nchan_used
    INTEGER(i_kind), ALLOCATABLE :: chan_idx(:)
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: radFiles ! radiance(tbb) observation files
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: angleFiles ! satellite angle observation files
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: cmFiles   ! L2 cloud mask files

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => giirs_setup
    PROCEDURE, PUBLIC  :: ObsIngest => giirs_read
    PROCEDURE, PUBLIC  :: ObsQC => giirs_qc

    PROCEDURE, PUBLIC :: ObsPrepareForSg 
    PROCEDURE, PUBLIC :: OutputForThinning
    PROCEDURE, PUBLIC :: destructor 
  END TYPE ReadWriteGIIRS_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    
    CLASS(ReadWriteGIIRS_t) :: this
    TYPE(State_t) :: X
    
    CALL this%initialize(X, this%bias_scheme, this%bias_predictors)
    CALL this%ObsPrepareForSgBase(X, this%configFile, this%inst_name, this%platform_name)

  END SUBROUTINE ObsPrepareForSg

  SUBROUTINE OutputForThinning(this, state, mpObs)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    USE mo_netcdf, only: NcDataset, NcDimension
    IMPLICIT NONE

    CLASS(ReadWriteGIIRS_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs

    CHARACTER(len=200) :: filename
    CHARACTER(len=1024) :: outputDir
    INTEGER(i_kind) :: istatus
        
    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    filename = TRIM(outputDir) //'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    CALL this%OutputForThinningBase(state, mpObs, filename, this%inst_name, this%platform_name)

  END SUBROUTINE OutputForThinning

  SUBROUTINE giirs_qc(this)

    IMPLICIT NONE
    CLASS(ReadWriteGIIRS_t) :: this

  END SUBROUTINE giirs_qc

  SUBROUTINE giirs_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ReadWriteGIIRS_t) :: this
    INTEGER(i_kind), PARAMETER :: nfile_max = 1000
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    CHARACTER(LEN=1024), ALLOCATABLE :: inst_name(:), platform_name(:), satinfos(:), satpaths(:)
    INTEGER(i_kind) :: nprefix, numvar_rad, numvar_angle, numvar_latlon
    INTEGER(i_kind) :: num_insts, num_plats, num_satinfos, num_satpaths, num_obsformats, numVars, obsFilesNum
    INTEGER(i_kind) :: i, j, k
    CHARACTER(LEN=1024) :: satinfo_file
    CHARACTER(LEN=2) :: chan_str
    INTEGER(i_kind) :: ichan, giirs_loc, fy4_1_loc
    INTEGER(i_kind) :: nfiles_fdi, nfiles_clm, ifiles_fdi, ifiles_clm, ntxts
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: fdi_file, clm_file
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: fdi_filenames, clm_filenames, tempClm
    LOGICAL :: istat
    INTEGER(i_kind) :: file_id = 18
    INTEGER(i_kind) :: istatus

    ! (1) Read the namelist file and get necessary parameters
    this%configFile = configFile

    istatus = yaml_get_var(TRIM(this%configFile), 'FY4-GIIRS', 'satellite', this%platform_name)
    CALL Get_rttov_giirs_chan_info(this%nchans, this%nchan_used, this%ifuse, this%chan_idx, TRIM(this%platform_name)//"-giirs", this%chan_lists, this%rttov_chan_lists)
    this%numVars = this%nchan_used + 4

    istatus = yaml_get_var(TRIM(configFile), 'RTTOV',     'inst_name',        inst_name)
    istatus = yaml_get_var(TRIM(configFile), 'RTTOV', 'platform_name',    platform_name)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS',       'satinfo',         satinfos)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS',       'satpath',         satpaths)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS',     'fdi_files',         fdi_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS',     'clm_files',         clm_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS',   'bias_scheme', this%bias_scheme)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-GIIRS', 'bias_predictors', this%bias_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)

    num_insts = UBOUND(inst_name, 1)
    num_plats = UBOUND(platform_name, 1)
    num_satinfos = UBOUND(satinfos, 1)
    num_satpaths = UBOUND(satpaths, 1)
    ntxts = UBOUND(fdi_file, 1)
    WRITE (*, 11) num_insts, num_plats, num_satinfos, ntxts, TRIM(satinfos(1)), &
      TRIM(satpaths(1)), TRIM(fdi_file(1)), TRIM(clm_file(1))
11  FORMAT('giirs_setup:', 4I2, ' satinfos: ', A, ' path:', A, ' fdi:', A, ' clm: ', A)

    ! IF (num_insts .NE. num_satinfos .or. num_insts .NE. num_satpaths .or. num_insts .NE. ntxts) THEN
    !   PRINT *, 'Numbers of inst_name/platform_name and satinfos/satpaths should be identical'
    !   STOP
    ! END IF

    ! Find location/index of the giirs instrument in the namelist file
    ! giirs_loc = FINDLOC([inst_name], 'giirs', dim=1)
    ! fy4_1_loc = FINDLOC([platform_name], 'fy4_1', dim=1)

    BLOCK
      INTEGER IDX
      DO idx = LBOUND(inst_name, 1), UBOUND(inst_name, 1)
        IF (inst_name(idx) == 'giirs') THEN
          giirs_loc = idx
          EXIT
        END IF
      END DO

        DO idx = LBOUND(platform_name, 1), UBOUND(platform_name, 1)
          IF( platform_name(idx) == 'fy4_1' .OR. platform_name(idx) == 'fy4_2') THEN
            fy4_1_loc = idx
            EXIT
          END IF
        END DO
    END BLOCK  

    IF (giirs_loc /= fy4_1_loc) THEN
      PRINT *, 'giirs_loc = ', giirs_loc
      PRINT *, 'fy4_1_loc = ', fy4_1_loc
      PRINT *, '====== Warning: Instrument and platform names should be matched ======'
      PRINT *, '====== BUT agri and giirs can share the same platform name ======'
      ! STOP
    END IF

    ! Read in an interpolation option:
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! (2) Prepare those obs files for later reading
    ALLOCATE (fdi_filenames(nfile_max), clm_filenames(nfile_max), tempClm(nfile_max))

    ifiles_fdi = 1
    DO WHILE (.TRUE.)
      ! OPEN (unit=file_id, file=TRIM(satpaths(giirs_loc))//'/'//TRIM(fdi_file(giirs_loc)))
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(fdi_file(1)))
      READ (file_id, *, END=200) fdi_filenames(ifiles_fdi)
      !  PRINT *,'radFiles: ', TRIM(fdi_filenames(ifiles_fdi))
      ifiles_fdi = ifiles_fdi + 1
    END DO
200 CONTINUE
    CLOSE (file_id)

    ifiles_clm = 1
    DO WHILE (.TRUE.)
      ! OPEN (unit=file_id, file=TRIM(satpaths(giirs_loc))//'/'//TRIM(clm_file(giirs_loc)))
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(clm_file(1)))
      READ (file_id, *, END=202) clm_filenames(ifiles_clm)
      !  PRINT *,'cmFiles: ', TRIM(clm_filenames(ifiles_clm))
      ifiles_clm = ifiles_clm + 1
    END DO
202 CONTINUE
    CLOSE (file_id)

    nfiles_fdi = ifiles_fdi - 1
    nfiles_clm = ifiles_clm - 1
    PRINT *, 'Numbers of FDI/CLM files:', nfiles_fdi, nfiles_clm
    ! TODO: remove those unmatched files in SHELL
    !IF ( nfiles_fdi .NE. nfiles_clm) THEN
    !  PRINT *, 'Numbers of FDI/CLM files should be identical'
    !  STOP
    !END IF
    this%numfiles = nfiles_fdi
    IF (.NOT. ALLOCATED(this%radFiles)) &
      ALLOCATE (this%radFiles(this%numFiles), this%angleFiles(this%numFiles), this%cmFiles(this%numFiles))

    ! Matching ClmFile with FdiFile
    CALL MatchClm(fdi_filenames, clm_filenames, tempClm, nfiles_fdi, nfiles_clm)

    DO i = 1, this%numfiles
      ! this%radFiles(i) = TRIM(satpaths(giirs_loc))//fdi_filenames(i)
      ! this%cmFiles(i) = TRIM(satpaths(giirs_loc))//tempClm(i)
      this%radFiles(i) = TRIM(satpaths(1))//'/'//fdi_filenames(i)
      this%cmFiles(i) = TRIM(satpaths(1))//'/'//tempClm(i)
!       WRITE(*,1) i,TRIM(this%radFiles(i)), TRIM(this%cmFiles(i))
! 1     FORMAT('Giirs setup file: ',I3,' Radfile: ',A,' CMfile: ',A)
    END DO

    IF (ALLOCATED(fdi_file)) DEALLOCATE(fdi_file)
    IF (ALLOCATED(clm_file)) DEALLOCATE(clm_file)
    IF (ALLOCATED(fdi_filenames)) DEALLOCATE(fdi_filenames)
    IF (ALLOCATED(clm_filenames)) DEALLOCATE(clm_filenames)
    DEALLOCATE(inst_name, platform_name, satinfos, satpaths)    
    ! PRINT *, "giirs_setup successfully run"

  END SUBROUTINE giirs_setup

  SUBROUTINE MatchClm(fdi_file, clm_file, cmFiles, numfiles, nfiles_clm)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN)   ::fdi_file(:)
    CHARACTER(len=*), INTENT(IN)   ::clm_file(:)
    INTEGER(i_kind), INTENT(IN)   ::numfiles, nfiles_clm
    CHARACTER(len=*), INTENT(INOUT)   ::cmFiles(:)

    INTEGER(i_kind), ALLOCATABLE :: gmt_fdi(:), gmt_clm(:)
    REAL(8), ALLOCATABLE :: obstime_fdi(:), obstime_clm(:)
    CHARACTER(len=14)   :: xx
    CHARACTER(len=14)   :: ctime
    INTEGER(i_kind)     :: i, j, n
    REAL(r_kind)        :: dt0, dt
    INTEGER(i_kind), ALLOCATABLE :: hhmm_fdi(:), hhmm_clm(:)

    ALLOCATE (obstime_fdi(numfiles))
    ALLOCATE (obstime_clm(numfiles))
    ALLOCATE (hhmm_fdi(numfiles))
    ALLOCATE (hhmm_clm(numfiles))
    ALLOCATE (gmt_fdi(6), gmt_clm(6))
    hhmm_fdi = 0
    hhmm_clm = 0
    DO i = 1, numfiles
      READ (fdi_file(i), "(A14,A14)") xx, ctime
      ! PRINT *, i,"fditime:",ctime
      READ (ctime, "(F14.0)") obstime_fdi(i)
      ! print *,"ObsTime_fdi: ",obstime_fdi(i)

      gmt_fdi(1) = INT(obstime_fdi(i) / 1E10)
      gmt_fdi(2) = INT(MOD(obstime_fdi(i), 1E10) / 1E8)
      gmt_fdi(3) = INT(MOD(obstime_fdi(i), 1E8) / 1E6)
      gmt_fdi(4) = INT(MOD(obstime_fdi(i), 1E6) / 1E4)
      gmt_fdi(5) = INT(MOD(obstime_fdi(i), 1E4) / 1E2)
      gmt_fdi(6) = INT(MOD(obstime_fdi(i), 1E2) / 1E0)
      hhmm_fdi(i) = gmt_fdi(4) * 100 + gmt_fdi(5)
      ! print *,'hhmm_fdi(i)',hhmm_fdi(i),TRIM(fdi_file(i))
    END DO

    DO i = 1, nfiles_clm
      ! read(clm_file(i),"(A11,A14)") xx,ctime  ! GIIRS naming format
      read(clm_file(i),"(A14,A14)") xx,ctime  ! original naming format
      ! PRINT *, i,"clmtime:",ctime
      READ (ctime, "(F14.0)") obstime_clm(i)
      ! print *,"ObsTime_clm: ",obstime_clm(i)

      gmt_clm(1) = INT(obstime_clm(i) / 1E10)
      gmt_clm(2) = INT(MOD(obstime_clm(i), 1E10) / 1E8)
      gmt_clm(3) = INT(MOD(obstime_clm(i), 1E8) / 1E6)
      gmt_clm(4) = INT(MOD(obstime_clm(i), 1E6) / 1E4)
      gmt_clm(5) = INT(MOD(obstime_clm(i), 1E4) / 1E2)
      gmt_clm(6) = INT(MOD(obstime_clm(i), 1E2) / 1E0)
      hhmm_clm(i) = gmt_clm(4) * 100 + gmt_clm(5)
      ! print *,'hhmm_clm(i)',hhmm_clm(i),TRIM(clm_file(i))
    END DO
    PRINT *, 'clm_file(n) size:', SIZE(clm_file), numfiles, nfiles_clm

    n = 1
    DO i = 1, numfiles
      dt0 = 10000.0
      DO j = 1, nfiles_clm
        dt = ABS(REAL(hhmm_fdi(i)) - REAL(hhmm_clm(j)))
        !print *,i,j,dt,dt0,hhmm_fdi(i),hhmm_clm(j)
        IF (dt .LT. dt0) THEN
          !print *,'if<dt0:',i,j,dt,dt0,hhmm_fdi(i),hhmm_clm(j)
          dt0 = dt
          n = j
        END IF
      END DO
      ! print *,'do i-loop',i,n,hhmm_fdi(i),hhmm_clm(n)
      cmFiles(i) = clm_file(n)
    END DO
    DEALLOCATE (obstime_fdi, obstime_clm, hhmm_fdi, hhmm_clm, gmt_fdi, gmt_clm)

  END SUBROUTINE MatchClm

  SUBROUTINE giirs_read(this, X)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    IMPLICIT NONE
    CLASS(ReadWriteGIIRS_t) :: this
    TYPE(State_t) :: X
    TYPE(RawGIIRS_t) :: RawGIIRS
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)

    INTEGER(i_kind) :: num_obs_total, nobs, ichan, i, j, iobs, ich
    ! CHARACTER(LEN=18) :: obstime_str
    CHARACTER(LEN=2)  :: chan_name
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: nRawObs_total, ivar
    INTEGER(i_kind) :: nv, nc, numgrd(2), fileunit
    CHARACTER(len=1024) :: filename_bias_coefs
    REAL(r_kind) :: domain(2,2)
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :)         ! observation errors
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), iscanpos(:)
    REAL(r_kind), ALLOCATABLE :: Angles(:, :)

    ! Get the latlon grid information:
    CALL getDomain_latlon(this%configFile, nv, nc, numgrd, domain)
    WRITE (*, 1) domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2)
1   FORMAT('giirs_read - Analysis domain lats in degrees', 2D12.4, ' lons', 2D12.4)
    domain = domain * degree2radian
    
    nRawObs_total = 0
    IF (this%numFiles > 0) THEN
      CALL RawGIIRS%RawGIIRS_setup(this%configFile)
      CALL RawGIIRS%RawGIIRS_read(TRIM(this%radFiles(1)), TRIM(this%cmFiles(1)))

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawGIIRS%obsvalue,RawGIIRS%nchans,RawGIIRS%ifuse)
      ! Debugging: test no cloud by Yuanfu Xie 2022-12-24
      ! RawGIIRS%cloud_flag = 0.0D0

      DO ichan = 1, SIZE(RawGIIRS%obsvalue, 2)
        CALL Select_Clear('giirs', RawGIIRS%obsvalue(:, ichan), RawGIIRS%cloud_flag, ichan)
      END DO
      ! CALL Select_Clear('giirs',RawGIIRS%obsvalue,RawGIIRS%nRawObs,RawGIIRS%cloud_flag)
!       WRITE(*,15) RawGIIRS%nRawObs,MAXVAL(RawGIIRS%obsvalue(:, :)),MINVAL(RawGIIRS%obsvalue(:,:))
! 15    FORMAT('After Select_Clear nRawObs',I4,' max-min val: ',2D12.4)
      CALL QC_BoundaryCheck(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%iscanpos)
!       WRITE(*,16) RawGIIRS%nRawObs,MAXVAL(RawGIIRS%obsvalue(:, :)),MINVAL(RawGIIRS%obsvalue(:,:))
! 16    FORMAT('After QC_BoundaryCheck nRawObs',I4,' max-min val: ',2D12.4)
      CALL Select_InDomain(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%latitude, RawGIIRS%longitude, &
                           domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))
!       WRITE(*,17) RawGIIRS%nRawObs,MAXVAL(RawGIIRS%obsvalue(:, :)),MINVAL(RawGIIRS%obsvalue(:,:))
! 17    FORMAT('After Select_InDomain nRawObs',I4,' max-min val: ',2D12.4)
      CALL Select_InTimeWin(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%obstime, X%sg%tt(1), X%sg%tt(X%sg%tSlots))
!       WRITE(*,4) RawGIIRS%nRawObs,MAXVAL(RawGIIRS%obsvalue(:, :)),MINVAL(RawGIIRS%obsvalue(:,:))
! 4     FORMAT('END Select OBS QC etc nRawObs',I4,' max-min val: ',2D12.4)

      ! TODO: Need to set num_obs_total for each channel
      DO iobs = 1, RawGIIRS%nRawObs
        IF (MAXVAL(RawGIIRS%obsvalue(iobs, :)) .LT. missing) &
          WRITE (*, 2) iobs, i, RawGIIRS%latitude(iobs), RawGIIRS%longitude(iobs), &
          RawGIIRS%obstime(iobs), missing, UBOUND(RawGIIRS%obsvalue, 2), &
          MAXVAL(RawGIIRS%obsvalue(iobs, :)), MINVAL(RawGIIRS%obsvalue(iobs, :))
2       FORMAT('giirs_read iobs', I4, ' from file: ', I4, ' lat/lon/obstime:', 2D12.4, I16, /, ' missingValue: ', D12.4, &
               ' ObSize:', I3, ' max-min: ', 2D12.4)
        IF (ANY(RawGIIRS%obsvalue(iobs, :) .LT. missing) .AND. &
            RawGIIRS%latitude(iobs) .LT. missing .AND. &
            RawGIIRS%longitude(iobs) .LT. missing .AND. &
            RawGIIRS%obstime(iobs) .GT. missing) THEN
          nRawObs_total = nRawObs_total + 1
          ! print *,'Used obs:',i,iobs,nRawObs_total,TRIM(this%radFiles(i)),RawGIIRS%latitude(iobs),RawGIIRS%longitude(iobs)
        END IF
      END DO

      CALL RawGIIRS%Destroy_RawGIIRS()
      nRawObs_total = nRawObs_total * this%numFiles
    END IF
      
    IF(X%sg%isBaseProc()) PRINT *, 'nRawObs_total = ', nRawObs_total
    this%numObs = nRawObs_total

    ! For satellite observations, number of variables is : nchans + four angles
    !this%numVars = RawGIIRS%nchans + 4
    this%numVars = this%nchan_used + 4
    PRINT *, 'this%numVars:', this%numVars
    ! ObsType: 'fy4_1_giirs'
    ! VarType: 'tbb'

    this%obsType = TRIM(RawGIIRS%platform_name)//'_'//TRIM(RawGIIRS%inst_name)
    this%numObs = nRawObs_total

    ALLOCATE(this%varNames(this%numVars))
    ALLOCATE(obsData(nRawObs_total, this%numVars))
    ALLOCATE(obsErrs(nRawObs_total, this%numVars))
    ALLOCATE(Angles(4, nRawObs_total))
    ALLOCATE(olatlon(2, nRawObs_total))
    ALLOCATE(obsTime(nRawObs_total))
    ALLOCATE(cloud_flag(nRawObs_total))
    ALLOCATE(iscanpos(nRawObs_total))

    obsTime = missing
    obsData = missing
    obsErrs = missing
    Angles = missing
    olatlon = missing 
    cloud_flag = missing
    iscanpos = missing 

    nobs = 0
    DO i = 1, this%numFiles

      !!! Read twice to reduce memory requirement !!!
      ! e.g. 10 obs files, each file has 14 2D arrays, each 2D array has a size of 2074*2074
      ! At the same time ,there are 11 obsbase types
      CALL RawGIIRS%RawGIIRS_setup(this%configFile)

      CALL RawGIIRS%RawGIIRS_read(TRIM(this%radFiles(i)), TRIM(this%cmFiles(i)))
      ! PRINT *, i,'RawGIIRS_read: ', TRIM(this%radFiles(i))
      PRINT *, 'START OBS QC'
      ! CALL Select_active_channels(RawGIIRS%obsvalue,RawGIIRS%nchans,RawGIIRS%ifuse)

      ! Debugging: test no cloud by Yuanfu Xie 2022-12-24
      ! RawGIIRS%cloud_flag = 0.0D0

      DO ichan = 1, SIZE(RawGIIRS%obsvalue, 2)
        CALL Select_Clear('giirs', RawGIIRS%obsvalue(:, ichan), RawGIIRS%cloud_flag, ichan)
      END DO
      ! CALL Select_Clear('giirs',RawGIIRS%obsvalue,RawGIIRS%nRawObs,RawGIIRS%cloud_flag)
      CALL QC_BoundaryCheck(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%iscanpos)
      CALL Select_InDomain(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%latitude, RawGIIRS%longitude, &
                           domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))
      CALL Select_InTimeWin(RawGIIRS%obsvalue, RawGIIRS%nRawObs, RawGIIRS%obstime, X%sg%tt(1), X%sg%tt(X%sg%tSlots))

      DO iobs = 1, RawGIIRS%nRawObs
        IF ( ANY(RawGIIRS%obsvalue(iobs, :) .LT. missing ) .AND. &
        RawGIIRS%latitude(iobs) .LT. missing .AND. &
        RawGIIRS%longitude(iobs) .LT. missing .AND.  &
        RawGIIRS%obstime(iobs) .GT. missing ) THEN 

        nobs = nobs + 1
        obsTime(nobs) = RawGIIRS%obstime(iobs)
        olatlon(1,nobs) = RawGIIRS%latitude(iobs)
        olatlon(2,nobs) = RawGIIRS%longitude(iobs)

        DO ichan = 5, this%numVars
          obsData(nobs,ichan) = RawGIIRS%obsvalue(iobs, ichan-4)
          obsErrs(nobs,ichan) = RawGIIRS%obserr(ichan-4)
        END DO

        cloud_flag(nobs) = RawGIIRS%cloud_flag(iobs)
        obsTime(nobs)    = RawGIIRS%obstime(iobs)
        olatlon(1,nobs) = RawGIIRS%latitude(iobs)
        olatlon(2,nobs) = RawGIIRS%longitude(iobs)

        ! Set for four satellite angles
        Angles(1, nobs) = RawGIIRS%sat_zenith(iobs)
        Angles(2, nobs) = RawGIIRS%sat_azi(iobs)
        Angles(3, nobs) = RawGIIRS%sol_zenith(iobs)
        Angles(4, nobs) = RawGIIRS%sol_azi(iobs)

        iscanpos(nobs) = RawGIIRS%iscanpos(iobs)

        END IF
      END DO

      CALL RawGIIRS%Destroy_RawGIIRS()
    END DO

    ! Set for angles
    DO ichan = 1, 4
      obsData(:,ichan) = Angles(ichan, :)
    END DO
    
    ! Remove missing values
    this%numObs = nobs

    ALLOCATE(this%obsHght(this%numObs))
    this%obsHght = -1.0E8
    this%obsTime = obsTime(1:this%numObs)
    this%obsData = obsData(1:this%numObs, :)
    this%obsData_org = this%obsData
    this%obsErrs = obsErrs(1:this%numObs,:)
    this%Angles = Angles(:, 1:this%numObs)
    this%olatlon = olatlon(:, 1:this%numObs)
    this%cloud_flag = cloud_flag(1:this%numObs)
    this%iscanpos = iscanpos(1:this%numObs)

    this%varNames(1) = 'sat_zenith'
    this%varNames(2) = 'sat_azi'
    this%varNames(3) = 'sol_zenith'
    this%varNames(4) = 'sol_azi'
    this%varNames(5:this%numVars) = 'tbb'

    DEALLOCATE(obsTime, obsData, obsErrs, Angles, olatlon, cloud_flag, iscanpos)
    ! END OF reading 

    this%nchans = this%nchan_used

    ! Read coef files
    fileunit = 101
    this%npred = 3
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", filename_bias_coefs)
    filename_bias_coefs = TRIM(filename_bias_coefs)//'/Satellite/'//'bcor_fy4a_hps'
    PRINT *, 'Read_bias_coefs_WRFDA: ', TRIM(filename_bias_coefs)

    IF (TRIM(this%bias_scheme) .EQ. 'ScanAirs') THEN
      CALL Read_bias_coefs_WRFDA_CMAGFS(fileunit, filename_bias_coefs, this%npred, this%BCOEF, this%SCOEF)
      PRINT *, 'Check coefs values: ', SHAPE(this%bcoef), MAXVAL(this%bcoef), MINVAL(this%bcoef)
    END IF

    ! ! Initialize auxiliary variables
    ! ALLOCATE(this%tb_bakAtmodel(state%sg%num_cell, state%sg%tSlots, this%numVars-4))
    ! this%tb_bakAtmodel = ZERO
    PRINT *, 'giirs_read has been done successfully'
  END SUBROUTINE giirs_read

  SUBROUTINE destructor(this)
    CLASS(ReadWriteGIIRS_t) :: this

    CALL this%destroy()
    IF (ALLOCATED(this%cmFiles)) DEALLOCATE (this%cmFiles)
    IF (ALLOCATED(this%angleFiles)) DEALLOCATE (this%angleFiles)
    IF (ALLOCATED(this%radFiles)) DEALLOCATE (this%radFiles)
    IF (ALLOCATED(this%chan_idx)) DEALLOCATE (this%chan_idx)

  END SUBROUTINE destructor

END MODULE ReadWriteGIIRS_m
