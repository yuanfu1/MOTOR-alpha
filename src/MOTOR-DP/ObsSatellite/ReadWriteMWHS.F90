!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ReadWriteMWHS
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @briefs
!! This module implements a satellite radiance data structure.
MODULE ReadWriteMWHS_m
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
  USE RawMWHS_m, ONLY: RawMWHS_t
  USE Dp_utils_m
  USE Satellite_utils_m
  USE RawSatellite_BC_m
  USE RadBase_m, ONLY: RadBase_t
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t

  IMPLICIT NONE

  TYPE, EXTENDS(RadBase_t) :: ReadWriteMWHS_t

    CHARACTER(LEN=1024) :: inst_name = 'mwhs2', platform_name = 'fy3_4'
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: radFiles ! radiance(tbb) observation files

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => mwhs_setup
    PROCEDURE, PUBLIC  :: ObsIngest => mwhs_read
    PROCEDURE, PUBLIC  :: ObsQC => mwhs_qc

    PROCEDURE, PUBLIC :: ObsPrepareForSg
    PROCEDURE, PUBLIC :: OutputForThinning
    PROCEDURE, PUBLIC :: destructor
  END TYPE ReadWriteMWHS_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)

    CLASS(ReadWriteMWHS_t) :: this
    TYPE(State_t) :: X

    CALL this%initialize(X, this%bias_scheme, this%bias_predictors)
    CALL this%ObsPrepareForSgBase(X, this%configFile, this%inst_name, this%platform_name)

  END SUBROUTINE ObsPrepareForSg

  SUBROUTINE OutputForThinning(this, state, mpObs)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    USE mo_netcdf, ONLY: NcDataset, NcDimension
    IMPLICIT NONE

    CLASS(ReadWriteMWHS_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs

    CHARACTER(len=200) :: filename
    CHARACTER(len=1024) :: outputDir
    INTEGER :: istatus

    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    filename = TRIM(outputDir)//'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    CALL this%OutputForThinningBase(state, mpObs, filename, this%inst_name, this%platform_name)

  END SUBROUTINE OutputForThinning

  SUBROUTINE mwhs_qc(this)
    CLASS(ReadWriteMWHS_t) :: this

  END SUBROUTINE mwhs_qc

  SUBROUTINE mwhs_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ReadWriteMWHS_t) :: this
    INTEGER(i_kind), PARAMETER :: nfile_max = 1000
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    CHARACTER(LEN=1024), ALLOCATABLE :: satpaths(:)
    INTEGER(i_kind) :: numVars
    INTEGER(i_kind) :: i, j, k
    CHARACTER(LEN=1024) :: satinfo_file
    CHARACTER(LEN=2) :: chan_str
    INTEGER(i_kind) :: ichan
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: rad_file, rad_filenames
    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
    INTEGER(i_kind) :: nfiles_rad, ifiles_rad
    INTEGER :: istat
    INTEGER(i_kind) :: file_id = 18
    INTEGER(i_kind) :: istatus
    INTEGER :: nchans

    ! (1) Read the namelist file and get necessary parameters
    this%configFile = configFile

    IF (yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'satellite', platform_name) .NE. 0) STOP
    IF (yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'instrument', inst_name) .NE. 0) STOP
    IF (PRESENT(idxFile)) THEN
      this%idxFile = idxFile
      this%platform_name = TRIM(platform_name(idxFile))
      this%inst_name = TRIM(inst_name(idxFile))
    ELSE
      PRINT *, 'NO instrument is used'
      STOP
    END IF
    DEALLOCATE (platform_name, inst_name)
    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans)
    nchans = this%nchans

    ALLOCATE (this%chan_lists(nchans), this%rttov_chan_lists(nchans), this%ifuse(nchans), this%oversea(nchans))

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans, chan_lists=this%chan_lists, &
                             rttov_chan_lists=this%rttov_chan_lists, ifuse=this%ifuse, only_over_sea=this%oversea)

    this%numVars = this%nchans + 4

    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'satpath', satpaths)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'rad_files', rad_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'bias_scheme', this%bias_scheme)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'bias_predictors', this%bias_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)

    ! Read in an interpolation option:
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! (2) Prepare those obs files for later reading
    ALLOCATE (rad_filenames(nfile_max))

    ifiles_rad = 1
    DO WHILE (.TRUE.)
      OPEN (unit=file_id, file=TRIM(satpaths(idxFile))//'/'//TRIM(rad_file(idxFile)), iostat=istat)
      READ (file_id, *, END=200) rad_filenames(ifiles_rad)
      ! PRINT *,'radFiles: ', TRIM(rad_filenames(ifiles_rad))
      ifiles_rad = ifiles_rad + 1
      IF (istat /= 0) THEN
        PRINT *, 'ifiles_rad not exit'
        CYCLE
      END IF
    END DO
200 CONTINUE
    CLOSE (file_id)

    nfiles_rad = ifiles_rad - 1

    this%numfiles = nfiles_rad
    ! this%numfiles = 1 ! Test only
    IF (.NOT. ALLOCATED(this%radFiles)) &
      ALLOCATE (this%radFiles(this%numFiles))
    DO i = 1, this%numfiles
      this%radFiles(i) = TRIM(satpaths(idxFile))//'/'//rad_filenames(i)
      ! PRINT *, this%radFiles(i)
    END DO

    IF (ALLOCATED(rad_file)) DEALLOCATE (rad_file)
    IF (ALLOCATED(rad_filenames)) DEALLOCATE (rad_filenames)
    DEALLOCATE (satpaths)
    ! PRINT *, "mwhs_setup successfully run"

  END SUBROUTINE mwhs_setup

  SUBROUTINE mwhs_read(this, X)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    USE domainCheck_m, ONLY: domainCheck_t
    IMPLICIT NONE
    CLASS(ReadWriteMWHS_t) :: this
    TYPE(State_t) :: X
    TYPE(RawMWHS_t) :: RawMWHS
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    INTEGER(i_kind), ALLOCATABLE :: pred_use(:, :)

    INTEGER(i_kind) :: num_obs_total, nobs, ichan, i, j, iobs
    CHARACTER(LEN=2)  :: chan_name
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: nRawObs_total, nRawObs, ivar
    INTEGER(i_kind) :: nv, nc, numgrd(2), fileunit
    CHARACTER(len=1024) :: filename_bias_coefs
    REAL(r_kind) :: domain(2, 2)
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :)         ! observation errors
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), iscanpos(:)
    REAL(r_kind), ALLOCATABLE :: Angles(:, :)
    REAL(r_kind), ALLOCATABLE :: org_olatlon(:, :), org_obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: Hidxtgt(:, :), Tidxtgt(:, :), Vidxtgt(:, :, :)
    REAL(r_kind), ALLOCATABLE :: Hcoetgt(:, :), Tcoetgt(:, :), Vcoetgt(:, :, :)
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: maskValided(:)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domainCheck

    ! Get the latlon grid information:
    CALL getDomain_latlon(this%configFile, nv, nc, numgrd, domain)
    domain = domain * degree2radian

    ! domain(1,1) = MINVAL(X%sg%cell_cntr(1,:))
    ! domain(2,1) = MAXVAL(X%sg%cell_cntr(1,:))
    ! domain(1,2) = MINVAL(X%sg%cell_cntr(2,:))
    ! domain(2,2) = MAXVAL(X%sg%cell_cntr(2,:))
    ! PRINT *, 'get domain ', domain(1,1), domain(2,1), domain(1,2), domain(2,2)

    nRawObs_total = 0
    IF (this%numFiles > 0) THEN
      DO i = 1, this%numFiles
        PRINT *, 'read for the first time'
        CALL RawMWHS%RawMWHS_setup(this%configFile, this%idxFile)
        PRINT *, 'RawMWHS_read: ', TRIM(this%radFiles(i))
        CALL RawMWHS%RawMWHS_read(TRIM(this%radFiles(i)))

        PRINT *, 'START Select OBS'
        ! CALL Select_active_channels(RawMWHS%obsvalue,RawMWHS%nchans,RawMWHS%ifuse)
        CALL Select_InDomain(RawMWHS%obsvalue, RawMWHS%nRawObs, RawMWHS%latitude, RawMWHS%longitude, &
                             domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))
        ! CALL Select_InTimeWin(RawMWHS%obsvalue,RawMWHS%nRawObs,RawMWHS%obstime,X%sg%tt(1), X%sg%tt(X%sg%tSlots))
        PRINT *, 'END Select OBS'

        ! TODO: Need to set num_obs_total for each channel
        DO iobs = 1, RawMWHS%nRawObs
          IF (ANY(RawMWHS%obsvalue(iobs, :) .LT. missing) .AND. &
              RawMWHS%latitude(iobs) .LT. missing .AND. &
              RawMWHS%longitude(iobs) .LT. missing) THEN
            nRawObs_total = nRawObs_total + 1
          END IF
        END DO

        CALL RawMWHS%Destroy_RawMWHS()
      END DO
    END IF

    IF (X%sg%isBaseProc()) PRINT *, 'nRawObs_total = ', nRawObs_total
    this%numObs = nRawObs_total

    ! For satellite observations, number of variables is : nchans + four angles
    this%numVars = RawMWHS%nchans + 4
    ! ObsType: 'fy3_5_mwhs3'
    ! VarType: 'tbb'

    this%obsType = TRIM(RawMWHS%platform_name)//'_'//TRIM(RawMWHS%inst_name)
    this%numObs = nRawObs_total

    ALLOCATE (this%varNames(this%numVars))
    ALLOCATE (obsData(nRawObs_total, this%numVars))
    ALLOCATE (obsErrs(nRawObs_total, this%numVars))
    ALLOCATE (Angles(4, nRawObs_total))
    ALLOCATE (olatlon(2, nRawObs_total))
    ALLOCATE (obsTime(nRawObs_total))
    ALLOCATE (cloud_flag(nRawObs_total))
    ALLOCATE (iscanpos(nRawObs_total))

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
      PRINT *, 'read for the second time'
      CALL RawMWHS%RawMWHS_setup(this%configFile, this%idxFile)
      ! PRINT *, 'RawMWHS_read: ', TRIM(this%radFiles(i))
      CALL RawMWHS%RawMWHS_read(TRIM(this%radFiles(i)))

      PRINT *, i, ' file has been read in '

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawMWHS%obsvalue,RawMWHS%nchans,RawMWHS%ifuse)
      nRawObs = SIZE(RawMWHS%latitude, 1)
      ALLOCATE (org_olatlon(2, nRawObs))
      ALLOCATE (org_obsHght(nRawObs))
      ALLOCATE (maskValided(nRawObs))
      org_olatlon(1, :) = RawMWHS%latitude(:)
      org_olatlon(2, :) = RawMWHS%longitude(:)
      org_obsHght(:) = -1.0E8
      CALL domainCheck%validation(X%sg, nRawObs, org_olatlon, org_obsHght, &
                                  RawMWHS%obstime, numValided, maskValided, &
                                  Hidxtgt, Vidxtgt, Tidxtgt, &
                                  Hcoetgt, Vcoetgt, Tcoetgt, &
                                  1, nhstencil, 'Satellite', &
                                  domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))

      DO iobs = 1, nRawObs

        IF (maskValided(iobs) .GT. 0) THEN

          nobs = nobs + 1
          ! print *, 'size of obsTime: ', SIZE(obsTime)
          ! print *, 'size of raw obstime: ', SIZE(RawMWHS%obstime)
          obsTime(nobs) = RawMWHS%obstime(iobs)
          olatlon(1, nobs) = RawMWHS%latitude(iobs)
          olatlon(2, nobs) = RawMWHS%longitude(iobs)

          DO ichan = 5, this%numVars
            obsData(nobs, ichan) = RawMWHS%obsvalue(iobs, ichan - 4)
            obsErrs(nobs, ichan) = RawMWHS%obserr(ichan - 4)
          END DO

          cloud_flag(nobs) = RawMWHS%cloud_flag(iobs)
          obsTime(nobs) = RawMWHS%obstime(iobs)
          olatlon(1, nobs) = RawMWHS%latitude(iobs)
          olatlon(2, nobs) = RawMWHS%longitude(iobs)

          ! Set for four satellite angles
          Angles(1, nobs) = RawMWHS%sat_zenith(iobs)
          Angles(2, nobs) = RawMWHS%sat_azi(iobs)
          Angles(3, nobs) = RawMWHS%sol_zenith(iobs)
          Angles(4, nobs) = RawMWHS%sol_azi(iobs)

          iscanpos(nobs) = RawMWHS%iscanpos(iobs)

        END IF
      END DO

      PRINT *, 'check before Destroy_RawMWHS'
      DEALLOCATE (org_olatlon, org_obsHght)
      DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided)
      CALL RawMWHS%Destroy_RawMWHS()
    END DO

    ! Set for angles
    DO ichan = 1, 4
      obsData(:, ichan) = Angles(ichan, :)
    END DO

    ! Remove missing values
    this%numObs = nobs

    ALLOCATE (this%obsHght(this%numObs))
    this%obsHght = -1.0E8
    this%obsTime = obsTime(1:this%numObs)
    this%obsData = obsData(1:this%numObs, :)
    this%obsData_org = this%obsData
    this%obsErrs = obsErrs(1:this%numObs, :)
    this%Angles = Angles(:, 1:this%numObs)
    this%olatlon = olatlon(:, 1:this%numObs)
    this%cloud_flag = cloud_flag(1:this%numObs)
    this%iscanpos = iscanpos(1:this%numObs)

    this%varNames(1) = 'sat_zenith'
    this%varNames(2) = 'sat_azi'
    this%varNames(3) = 'sol_zenith'
    this%varNames(4) = 'sol_azi'
    this%varNames(5:this%numVars) = 'tbb'

    DEALLOCATE (obsTime, obsData, obsErrs, Angles, olatlon, cloud_flag, iscanpos)
    ! END OF reading

    ! Read coef files
    fileunit = 101
    this%npred = 8
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", filename_bias_coefs)
    filename_bias_coefs = TRIM(filename_bias_coefs)//'/Satellite/'//'VARBC.'//TRIM(this%platform_name)//'-'//TRIM(this%inst_name)
    PRINT *, 'Read_bias_coefs_WRFDA: ', TRIM(filename_bias_coefs)

    IF (TRIM(this%bias_scheme) .EQ. 'ScanAirs') THEN
      CALL Read_bias_coefs_WRFDA(fileunit, filename_bias_coefs, this%npred, pred_use, this%BCOEF)
      PRINT *, 'shape of pred_use: ', SHAPE(pred_use)
      PRINT *, 'Check coefs values: ', SHAPE(this%bcoef), MAXVAL(this%bcoef), MINVAL(this%bcoef)

      DO i = 1, SIZE(pred_use, 1)
        DO j = 1, SIZE(pred_use, 2)
          IF (pred_use(i, j) .EQ. 0) this%BCOEF(i, j) = ZERO
        END DO
      END DO
      IF (ALLOCATED(pred_use)) DEALLOCATE (pred_use)
    END IF

    ! Initialize auxiliary variables
    ! ALLOCATE(this%tb_bakAtmodel(X%sg%num_cell, X%sg%tSlots, this%numVars-4))
    ! this%tb_bakAtmodel = ZERO

    PRINT *, 'mwhs_read is successfully run'

  END SUBROUTINE mwhs_read

  SUBROUTINE destructor(this)
    CLASS(ReadWriteMWHS_t) :: this

    CALL this%destroy()
    IF (ALLOCATED(this%radFiles)) DEALLOCATE (this%radFiles)

  END SUBROUTINE destructor

END MODULE ReadWriteMWHS_m
