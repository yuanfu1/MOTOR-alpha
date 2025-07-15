!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ReadWriteHIRAS
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE ReadWriteHIRAS_m
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
  USE RawHIRAS_m, ONLY: RawHIRAS_t
  USE Dp_utils_m
  USE Satellite_utils_m
  USE RawSatellite_BC_m
  USE RadBase_m, ONLY: RadBase_t
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t

  IMPLICIT NONE

  TYPE, EXTENDS(RadBase_t) :: ReadWriteHIRAS_t
    
    CHARACTER(LEN=1024) :: inst_name='hiras', platform_name='fy3_4'
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: radFiles ! radiance(tbb) observation files

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => hiras_setup
    PROCEDURE, PUBLIC  :: ObsIngest => hiras_read
    PROCEDURE, PUBLIC  :: ObsQC => hiras_qc

    PROCEDURE, PUBLIC :: ObsPrepareForSg 
    PROCEDURE, PUBLIC :: OutputForThinning
    PROCEDURE, PUBLIC :: destructor 
  END TYPE ReadWriteHIRAS_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    
    CLASS(ReadWriteHIRAS_t) :: this
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

    CLASS(ReadWriteHIRAS_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    INTEGER(i_kind) :: istatus

    CHARACTER(len=200) :: filename
    CHARACTER(len=1024) :: outputDir
        
    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    filename = TRIM(outputDir) //'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    CALL this%OutputForThinningBase(state, mpObs, filename, this%inst_name, this%platform_name)

  END SUBROUTINE OutputForThinning

  SUBROUTINE hiras_qc(this)
    CLASS(ReadWriteHIRAS_t) :: this

  END SUBROUTINE hiras_qc

  SUBROUTINE hiras_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ReadWriteHIRAS_t) :: this
    INTEGER(i_kind), PARAMETER :: nfile_max = 1000
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    CHARACTER(LEN=1024), ALLOCATABLE :: satpaths(:)
    INTEGER(i_kind) :: numVars
    INTEGER(i_kind) :: i, j, k
    CHARACTER(LEN=1024) :: satinfo_file
    CHARACTER(LEN=2) :: chan_str
    INTEGER(i_kind) :: ichan
    INTEGER(i_kind) :: nfiles_rad, ifiles_rad
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: rad_file, rad_filenames
    INTEGER :: istat
    INTEGER(i_kind) :: file_id = 18
    INTEGER(i_kind) :: istatus
    INTEGER :: nchans

    ! (1) Read the namelist file and get necessary parameters
    this%configFile = configFile

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans)
    nchans = this%nchans

    ALLOCATE (this%chan_lists(nchans), this%rttov_chan_lists(nchans), this%ifuse(nchans), this%oversea(nchans))

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans, chan_lists=this%chan_lists, &
    rttov_chan_lists=this%rttov_chan_lists, ifuse=this%ifuse, only_over_sea=this%oversea)

    this%numVars = this%nchans + 4

    istatus = yaml_get_var(TRIM(configFile), 'FY3-HIRAS',       'satpath',         satpaths)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-HIRAS',     'rad_files',         rad_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-HIRAS',   'bias_scheme', this%bias_scheme)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-HIRAS', 'bias_predictors', this%bias_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)

    ! Read in an interpolation option:
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! (2) Prepare those obs files for later reading
    ALLOCATE(rad_filenames(nfile_max))

    ifiles_rad = 1
    DO WHILE(.TRUE.)
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(rad_file(1)), iostat = istat)
      read(file_id, *, end=200) rad_filenames(ifiles_rad)
      ! PRINT *,'radFiles: ', TRIM(rad_filenames(ifiles_rad))
      ifiles_rad = ifiles_rad + 1
      IF (istat /= 0) THEN
        PRINT *, 'ifiles_rad not exit'
        CYCLE
      END IF
    END DO 
    200 CONTINUE
    CLOSE(file_id)

    nfiles_rad = ifiles_rad - 1

    this%numfiles = nfiles_rad
    ! this%numfiles = 1 ! Test only
    IF (.NOT. ALLOCATED(this%radFiles)) &
    ALLOCATE(this%radFiles(this%numFiles))
    DO i = 1, this%numfiles
      this%radFiles(i) = TRIM(satpaths(1))//'/'//rad_filenames(i)
      ! PRINT *, this%radFiles(i)
    END DO

    IF (ALLOCATED(rad_file)) DEALLOCATE(rad_file)
    IF (ALLOCATED(rad_filenames)) DEALLOCATE(rad_filenames)
    DEALLOCATE(satpaths)    
    ! PRINT *, "hiras_setup successfully run"

  END SUBROUTINE hiras_setup

  SUBROUTINE hiras_read(this, X)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    USE domainCheck_m, ONLY: domainCheck_t
    IMPLICIT NONE
    CLASS(ReadWriteHIRAS_t) :: this
    TYPE(State_t) :: X
    TYPE(RawHIRAS_t) :: RawHIRAS
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    INTEGER(i_kind), ALLOCATABLE :: pred_use(:,:)

    INTEGER(i_kind) :: num_obs_total, nobs, ichan, i, j, iobs
    CHARACTER(LEN=2)  :: chan_name
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: nRawObs_total, nRawObs, ivar
    INTEGER(i_kind) :: nv,nc,numgrd(2),fileunit
    CHARACTER(len=1024) :: filename_bias_coefs
    REAL(r_kind) :: domain(2,2)
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :)         ! observation errors
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), iscanpos(:)
    REAL(r_kind), ALLOCATABLE :: Angles(:, :)
    REAL(r_kind), ALLOCATABLE :: org_olatlon(:,:), org_obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: Hidxtgt(:, :), Tidxtgt(:, :), Vidxtgt(:, :, :)
    REAL(r_kind), ALLOCATABLE :: Hcoetgt(:, :), Tcoetgt(:, :), Vcoetgt(:, :, :)
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: maskValided(:)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domainCheck

    ! Get the latlon grid information:
    CALL getDomain_latlon(this%configFile, nv,nc,numgrd,domain)
    domain = domain * degree2radian
    
    ! domain(1,1) = MINVAL(X%sg%cell_cntr(1,:))
    ! domain(2,1) = MAXVAL(X%sg%cell_cntr(1,:))
    ! domain(1,2) = MINVAL(X%sg%cell_cntr(2,:))
    ! domain(2,2) = MAXVAL(X%sg%cell_cntr(2,:))
    ! PRINT *, 'get domain ', domain(1,1), domain(2,1), domain(1,2), domain(2,2)

    nRawObs_total = 0
    IF (this%numFiles > 0) THEN
      CALL RawHIRAS%RawHIRAS_setup(this%configFile)
      PRINT *, 'RawHIRAS_read: ', TRIM(this%radFiles(1))
      CALL RawHIRAS%RawHIRAS_read(TRIM(this%radFiles(1)))

      PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawHIRAS%obsvalue,RawHIRAS%nchans,RawHIRAS%ifuse)
      CALL Select_InDomain(RawHIRAS%obsvalue,RawHIRAS%nRawObs,RawHIRAS%latitude,RawHIRAS%longitude,    &
                           domain(1,1), domain(2,1), domain(1,2), domain(2,2))
      ! CALL Select_InTimeWin(RawHIRAS%obsvalue,RawHIRAS%nRawObs,RawHIRAS%obstime,X%sg%tt(1), X%sg%tt(X%sg%tSlots))
      PRINT *, 'END Select OBS'

      ! TODO: Need to set num_obs_total for each channel
      DO iobs = 1, RawHIRAS%nRawObs
        IF ( ANY(RawHIRAS%obsvalue(iobs, :) .LT. missing ) .AND. &
        RawHIRAS%latitude(iobs) .LT. missing .AND. &
        RawHIRAS%longitude(iobs) .LT. missing ) THEN 
          nRawObs_total = nRawObs_total + 1
        END IF
      END DO

      CALL RawHIRAS%Destroy_RawHIRAS()
      nRawObs_total = nRawObs_total * this%numFiles
    END IF
      
    IF(X%sg%isBaseProc()) PRINT *, 'nRawObs_total = ', nRawObs_total
    this%numObs = nRawObs_total

    ! For satellite observations, number of variables is : nchans + four angles
    this%numVars = RawHIRAS%nchans + 4 
    ! ObsType: 'fy3_4_hiras'
    ! VarType: 'tbb'

    this%obsType = TRIM(RawHIRAS%platform_name)//'_'//TRIM(RawHIRAS%inst_name)
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
      CALL RawHIRAS%RawHIRAS_setup(this%configFile)
      ! PRINT *, 'RawHIRAS_read: ', TRIM(this%radFiles(i))
      CALL RawHIRAS%RawHIRAS_read(TRIM(this%radFiles(i)))

      PRINT *, i, ' file has been read in '

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawHIRAS%obsvalue,RawHIRAS%nchans,RawHIRAS%ifuse)
      nRawObs = SIZE(RawHIRAS%latitude,1)
      ALLOCATE(org_olatlon(2, nRawObs))
      ALLOCATE(org_obsHght(nRawObs))
      ALLOCATE (maskValided(nRawObs))
      org_olatlon(1,:) = RawHIRAS%latitude(:)
      org_olatlon(2,:) = RawHIRAS%longitude(:)
      org_obsHght(:) = -1.0E8
      CALL domainCheck%validation(X%sg, nRawObs, org_olatlon, org_obsHght, &
                          RawHIRAS%obstime, numValided, maskValided, &
                          Hidxtgt, Vidxtgt, Tidxtgt, &
                          Hcoetgt, Vcoetgt, Tcoetgt, &
                          1, nhstencil, 'Satellite', &
                          domain(1,1), domain(2,1), domain(1,2), domain(2,2))

      DO iobs = 1, nRawObs

        IF (maskValided(iobs) .GT. 0) THEN

          nobs = nobs + 1
          obsTime(nobs) = RawHIRAS%obstime(iobs)
          olatlon(1,nobs) = RawHIRAS%latitude(iobs)
          olatlon(2,nobs) = RawHIRAS%longitude(iobs)

          DO ichan = 5, this%numVars
            obsData(nobs,ichan) = RawHIRAS%obsvalue(iobs, ichan-4)
            obsErrs(nobs,ichan) = RawHIRAS%obserr(ichan-4)
          END DO

          cloud_flag(nobs) = RawHIRAS%cloud_flag(iobs)
          obsTime(nobs)    = RawHIRAS%obstime(iobs)
          olatlon(1,nobs) = RawHIRAS%latitude(iobs)
          olatlon(2,nobs) = RawHIRAS%longitude(iobs)

          ! Set for four satellite angles
          Angles(1, nobs) = RawHIRAS%sat_zenith(iobs)
          Angles(2, nobs) = RawHIRAS%sat_azi(iobs)
          Angles(3, nobs) = RawHIRAS%sol_zenith(iobs)
          Angles(4, nobs) = RawHIRAS%sol_azi(iobs)

          iscanpos(nobs) = RawHIRAS%iscanpos(iobs)

        END IF
      END DO

      print *, 'check before Destroy_RawHIRAS'
      DEALLOCATE(org_olatlon, org_obsHght)
      DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided)
      CALL RawHIRAS%Destroy_RawHIRAS()
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

    ! Read coef files
    fileunit = 101
    this%npred = 8
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", filename_bias_coefs)
    filename_bias_coefs=TRIM(filename_bias_coefs)//'/Satellite/'//'VARBC.'//TRIM(this%platform_name)//'-'//TRIM(this%inst_name)
    PRINT *, 'Read_bias_coefs_WRFDA: ', TRIM(filename_bias_coefs)

    IF ( TRIM(this%bias_scheme) .EQ. 'ScanAirs') THEN
      CALL Read_bias_coefs_WRFDA(fileunit, filename_bias_coefs, this%npred, pred_use, this%BCOEF)
      PRINT *, 'shape of pred_use: ', shape(pred_use)
      PRINT *, 'Check coefs values: ', shape(this%bcoef), maxval(this%bcoef), minval(this%bcoef)
    
      DO i = 1, SIZE(pred_use,1)
        DO j = 1, SIZE(pred_use,2)
          IF (pred_use (i, j) .EQ. 0) this%BCOEF(i,j) = ZERO
        END DO
      END DO
      IF (ALLOCATED(pred_use)) DEALLOCATE(pred_use)
    END IF

    ! Initialize auxiliary variables 
    ! ALLOCATE(this%tb_bakAtmodel(X%sg%num_cell, X%sg%tSlots, this%numVars-4))
    ! this%tb_bakAtmodel = ZERO

    PRINT *, 'hiras_read is successfully run'

  END SUBROUTINE hiras_read

  SUBROUTINE destructor(this)
    CLASS(ReadWriteHIRAS_t) :: this

    CALL this%destroy()
    IF (ALLOCATED(this%radFiles)) DEALLOCATE (this%radFiles)

  END SUBROUTINE destructor

END MODULE ReadWriteHIRAS_m
