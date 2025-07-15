!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ReadWriteMWTS
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
MODULE ReadWriteMWTS_m
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
  USE RawMWTS_m, ONLY: RawMWTS_t
  USE Dp_utils_m
  USE Satellite_utils_m
  USE RawSatellite_BC_m
  USE RadBase_m, ONLY: RadBase_t
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t

  IMPLICIT NONE

  TYPE, EXTENDS(RadBase_t) :: ReadWriteMWTS_t

    CHARACTER(LEN=1024) :: inst_name = 'mwts2', platform_name = 'fy3_4'
    REAL(r_kind), ALLOCATABLE :: land(:)
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: radFiles ! radiance(tbb) observation files

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => mwts_setup
    PROCEDURE, PUBLIC  :: ObsIngest => mwts_read
    PROCEDURE, PUBLIC  :: ObsQC => mwts_qc

    PROCEDURE, PUBLIC :: ObsPrepareForSg
    PROCEDURE, PUBLIC :: OutputForThinning
    PROCEDURE, PUBLIC :: destructor
  END TYPE ReadWriteMWTS_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)

    CLASS(ReadWriteMWTS_t) :: this
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

    CLASS(ReadWriteMWTS_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs

    CHARACTER(len=200) :: filename
    CHARACTER(len=1024) :: outputDir
    INTEGER :: istatus

    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    filename = TRIM(outputDir)//'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    PRINT *, 'OUTPUT TO = ', TRIM(filename)
    CALL this%OutputForThinningBase(state, mpObs, filename, this%inst_name, this%platform_name)

  END SUBROUTINE OutputForThinning

  SUBROUTINE mwts_qc(this)
    CLASS(ReadWriteMWTS_t) :: this

  END SUBROUTINE mwts_qc

  SUBROUTINE mwts_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ReadWriteMWTS_t) :: this
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

    IF (yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'satellite', platform_name) .NE. 0) STOP
    IF (yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'instrument', inst_name) .NE. 0) STOP
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

    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'satpath', satpaths)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'rad_files', rad_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'bias_scheme', this%bias_scheme)
    istatus = yaml_get_var(TRIM(configFile), 'FY3-MWTS', 'bias_predictors', this%bias_predictors)
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
    ! PRINT *, "mwts_setup successfully run"

  END SUBROUTINE mwts_setup

  SUBROUTINE mwts_read(this, X)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    USE domainCheck_m, ONLY: domainCheck_t
    IMPLICIT NONE
    CLASS(ReadWriteMWTS_t) :: this
    TYPE(State_t) :: X
    TYPE(RawMWTS_t) :: RawMWTS
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    INTEGER(i_kind), ALLOCATABLE :: pred_use(:, :)

    INTEGER(i_kind) :: num_obs_total, nobs, ichan, i, j, iobs
    CHARACTER(LEN=2)  :: chan_name
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: nRawObs_total, nRawObs, ivar
    INTEGER(i_kind) :: nv, nc, numgrd(2), fileunit
    CHARACTER(len=1024) :: filename_bias_coefs
    REAL(r_kind) :: domain(2, 2), domainL(2, 2), dlat, dlon
    REAL(r_kind), ALLOCATABLE :: bakvalue(:, :)
    REAL(r_kind), ALLOCATABLE :: satzenith(:)
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :)         ! observation errors
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), iscanpos(:), land(:)
    REAL(r_kind), ALLOCATABLE :: Angles(:, :)
    REAL(r_kind), ALLOCATABLE :: org_olatlon(:, :), org_obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: Hidxtgt(:, :), Tidxtgt(:, :), Vidxtgt(:, :, :)
    REAL(r_kind), ALLOCATABLE :: Hcoetgt(:, :), Tcoetgt(:, :), Vcoetgt(:, :, :)
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: maskValided(:)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domainCheck

    ! Get the latlon grid information:
    ! CALL getDomain_latlon(this%configFile, nv,nc,numgrd,domain)
    ! domain = domain * degree2radian

    domain(1, 1) = MINVAL(X%sg%cell_cntr(1, :))
    domain(2, 1) = MAXVAL(X%sg%cell_cntr(1, :))
    domain(1, 2) = MINVAL(X%sg%cell_cntr(2, :))
    domain(2, 2) = MAXVAL(X%sg%cell_cntr(2, :))

    dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)
    ! PRINT *, 'dlat/dlon: ', dlat/degree2radian, dlon/degree2radian
    ! We do not believe those obs fall on bdys
    domain(1, 1) = domain(1, 1) + dlat * 0.5D0
    domain(2, 1) = domain(2, 1) - dlat * 0.5D0
    domain(1, 2) = domain(1, 2) + dlon * 0.5D0
    domain(2, 2) = domain(2, 2) - dlon * 0.5D0
    ! PRINT *, 'get domain ', domain(1,1), domain(2,1), domain(1,2), domain(2,2)

    nRawObs_total = 0
    IF (this%numFiles > 0) THEN
      DO i = 1, this%numFiles
        CALL RawMWTS%RawMWTS_setup(this%configFile, this%idxFile)
        PRINT *, 'RawMWTS_read: ', TRIM(this%radFiles(i))
        CALL RawMWTS%RawMWTS_read(TRIM(this%radFiles(i)))

        PRINT *, 'START Select OBS'
        CALL Select_InDomain(RawMWTS%obsvalue, RawMWTS%nRawObs, RawMWTS%latitude, RawMWTS%longitude, &
                             domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))
        ! CALL Select_InTimeWin(RawMWTS%obsvalue,RawMWTS%nRawObs,RawMWTS%obstime,X%sg%tt(1), X%sg%tt(X%sg%tSlots))
        PRINT *, 'END Select OBS'

        ! TODO: Need to set num_obs_total for each channel
        DO iobs = 1, RawMWTS%nRawObs
          IF (ANY(RawMWTS%obsvalue(iobs, :) .LT. missing) .AND. &
              RawMWTS%latitude(iobs) .LT. missing .AND. &
              RawMWTS%longitude(iobs) .LT. missing) THEN
            nRawObs_total = nRawObs_total + 1
          END IF
        END DO

        CALL RawMWTS%Destroy_RawMWTS()
      END DO
    END IF

    IF (X%sg%isBaseProc()) PRINT *, 'nRawObs_total = ', nRawObs_total
    this%numObs = nRawObs_total

    ! For satellite observations, number of variables is : nchans + four angles
    this%numVars = RawMWTS%nchans + 4
    ! ObsType: 'fy3_5_mwts3'
    ! VarType: 'tbb'

    this%obsType = TRIM(RawMWTS%platform_name)//'_'//TRIM(RawMWTS%inst_name)
    this%numObs = nRawObs_total

    ALLOCATE (this%varNames(this%numVars))
    ALLOCATE (obsData(nRawObs_total, this%numVars))
    ALLOCATE (obsErrs(nRawObs_total, this%numVars))
    ALLOCATE (Angles(4, nRawObs_total))
    ALLOCATE (olatlon(2, nRawObs_total))
    ALLOCATE (obsTime(nRawObs_total))
    ALLOCATE (cloud_flag(nRawObs_total))
    ALLOCATE (iscanpos(nRawObs_total))
    ALLOCATE (land(nRawObs_total))

    obsTime = missing
    obsData = missing
    obsErrs = missing
    Angles = missing
    olatlon = missing
    cloud_flag = missing
    iscanpos = missing
    land = missing

    nobs = 0
    DO i = 1, this%numFiles

      !!! Read twice to reduce memory requirement !!!
      ! e.g. 10 obs files, each file has 14 2D arrays, each 2D array has a size of 2074*2074
      ! At the same time ,there are 11 obsbase types
      CALL RawMWTS%RawMWTS_setup(this%configFile, this%idxFile)
      ! PRINT *, 'RawMWTS_read: ', TRIM(this%radFiles(i))
      CALL RawMWTS%RawMWTS_read(TRIM(this%radFiles(i)))

      PRINT *, i, ' file has been read in '

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawMWTS%obsvalue,RawMWTS%nchans,RawMWTS%ifuse)
      nRawObs = SIZE(RawMWTS%latitude, 1)
      ALLOCATE (org_olatlon(2, nRawObs))
      ALLOCATE (org_obsHght(nRawObs))
      ALLOCATE (maskValided(nRawObs))
      org_olatlon(1, :) = RawMWTS%latitude(:)
      org_olatlon(2, :) = RawMWTS%longitude(:)
      org_obsHght(:) = -1.0E8
      CALL domainCheck%validation(X%sg, nRawObs, org_olatlon, org_obsHght, &
                                  RawMWTS%obstime, numValided, maskValided, &
                                  Hidxtgt, Vidxtgt, Tidxtgt, &
                                  Hcoetgt, Vcoetgt, Tcoetgt, &
                                  1, nhstencil, 'Satellite', &
                                  domain(1, 1), domain(2, 1), domain(1, 2), domain(2, 2))

      DO iobs = 1, nRawObs

        IF (maskValided(iobs) .GT. 0) THEN

          nobs = nobs + 1
          obsTime(nobs) = RawMWTS%obstime(iobs)
          olatlon(1, nobs) = RawMWTS%latitude(iobs)
          olatlon(2, nobs) = RawMWTS%longitude(iobs)

          DO ichan = 5, this%numVars
            obsData(nobs, ichan) = RawMWTS%obsvalue(iobs, ichan - 4)
            obsErrs(nobs, ichan) = RawMWTS%obserr(ichan - 4)
          END DO

          cloud_flag(nobs) = RawMWTS%cloud_flag(iobs)
          land(nobs) = RawMWTS%landseamask(iobs)
          obsTime(nobs) = RawMWTS%obstime(iobs)
          olatlon(1, nobs) = RawMWTS%latitude(iobs)
          olatlon(2, nobs) = RawMWTS%longitude(iobs)

          ! Set for four satellite angles
          Angles(1, nobs) = RawMWTS%sat_zenith(iobs)
          Angles(2, nobs) = RawMWTS%sat_azi(iobs)
          Angles(3, nobs) = RawMWTS%sol_zenith(iobs)
          Angles(4, nobs) = RawMWTS%sol_azi(iobs)

          iscanpos(nobs) = RawMWTS%iscanpos(iobs)

        END IF
      END DO

      PRINT *, 'check before Destroy_RawMWTS'
      DEALLOCATE (org_olatlon, org_obsHght)
      DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided)
      CALL RawMWTS%Destroy_RawMWTS()
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
    this%land = land(1:this%numObs)
    this%iscanpos = iscanpos(1:this%numObs)

    this%varNames(1) = 'sat_zenith'
    this%varNames(2) = 'sat_azi'
    this%varNames(3) = 'sol_zenith'
    this%varNames(4) = 'sol_azi'
    this%varNames(5:this%numVars) = 'tbb'

    DEALLOCATE (obsTime, obsData, obsErrs, Angles, olatlon, cloud_flag, iscanpos, land)
    ! END OF reading

    ! Read coef files
    fileunit = 101
    this%npred = 8
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", filename_bias_coefs)
    filename_bias_coefs = TRIM(filename_bias_coefs)//'/Satellite/'//'VARBC.'//TRIM(this%platform_name)//'-'//TRIM(this%inst_name)
    ! filename_bias_coefs=TRIM(filename_bias_coefs)//'/Satellite/'//'VARBC.fy3_5-mwts3'
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

    ! ! Initialize auxiliary variables
    ! ! ALLOCATE(this%tb_bakAtmodel(X%sg%num_cell, X%sg%tSlots, this%numVars-4))
    ! ! this%tb_bakAtmodel = ZERO
    ! ALLOCATE(bakvalue(this%numobs, this%numVars))
    ! BLOCK
    !   USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t
    !   TYPE(rttov_nl_sp_obsloc_t) :: rttov_nl_sp_obsloc
    !   REAL(r_kind), ALLOCATABLE :: pres(:,:), t(:,:), q(:,:), u(:,:), v(:,:), tskin(:), elevation(:)
    !   REAL(r_kind), ALLOCATABLE :: snowc(:), soiltype(:), landmask(:)
    !   LOGICAL, ALLOCATABLE :: PassDomainCheck(:)

    !   CALL rttov_nl_sp_obsloc%initialize(this%configFile, X, 'mwts3', 'fy3_5')

    !   ALLOCATE(pres(this%numobs, X%sg%vLevel))
    !   ALLOCATE(t(this%numobs, X%sg%vLevel))
    !   ALLOCATE(q(this%numobs, X%sg%vLevel))
    !   ALLOCATE(u(this%numobs, X%sg%vLevel))
    !   ALLOCATE(v(this%numobs, X%sg%vLevel))
    !   ALLOCATE(tskin(this%numobs))
    !   ALLOCATE(elevation(this%numobs))
    !   ALLOCATE(snowc(this%numobs))
    !   ALLOCATE(soiltype(this%numobs))
    !   ALLOCATE(landmask(this%numobs))
    !   ALLOCATE(PassDomainCheck(this%numobs))

    !   CALL Prof_bkg2obs(X, this%numobs, this%olatlon, this%obsHght, this%obsTime, &
    !                 pres=pres, t=t, q=q, u=u, v=v, tskin=tskin, elevation=elevation, &
    !                 snowc=snowc, soiltype=soiltype, landmask=landmask, PassDomainCheck=PassDomainCheck)

    !   DO ichan = 5, this%numVars

    !     CALL TB_bkgAtobs(X, this%rttov_nl_sp_obsloc, pres, t, q, u, v, tskin, elevation, &
    !                       snowc, soiltype, landmask, PassDomainCheck, i)

    !     WHERE (this%tb_inv .GE. ABS(missing) - 1.0 .AND. this%obsData(:,1) .GE. ABS(missing) - 1.0)
    !       bakvalue(:,ichan) = missing
    !     ELSEWHERE
    !       bakvalue(:,ichan) = this%tb_inv + this%obsData(:,1)
    !     END WHERE

    !   END DO
    !   DEALLOCATE(pres, t, q, u, v, tskin, elevation, snowc, soiltype, landmask, PassDomainCheck)
    ! END BLOCK
    ! ! Call for the calculation of cloud_flag (CLW) for MWTS、MWHS
    ! ! Note that the channel index should be the real channel index + 4 (angles)
    ! ! Note that fortran uses radian in cos/sin calculation
    ! ! satzenith = this%RadBase(1)%Angles(1, :) * radian2degree
    ! CALL Calc_MWCLWCheck(this%obsData(:,5), bakvalue(:,5), &
    ! this%obsData(:,6), bakvalue(:,6), &
    ! this%obsData(:,9), bakvalue(:,9), &
    ! this%Angles(1, :), this%land, this%cloud_flag)

    ! DO ichan = 1, this%numVars
    !   this%cloud_flag = this%cloud_flag
    ! END DO
    ! DEALLOCATE(bakvalue)
    ! DEALLOCATE(satzenith)
    ! MWTS does NOT have the 89GHz channel, so this SI is not available for now
    ! CALL Calc_MWScatteringIndex(this%obsvalue(:,1), this%bakvalue(:,1), &
    ! this%obsvalue(:,2), this%bakvalue(:,2), &
    ! this%obsvalue(:,9), this%bakvalue(:,9), &
    ! this%obsvalue(:,15), this%bakvalue(:,15), this%landseamask, this%SI)

    PRINT *, 'mwts_read is successfully run'

  END SUBROUTINE mwts_read

  SUBROUTINE destructor(this)
    CLASS(ReadWriteMWTS_t) :: this

    CALL this%destroy()
    IF (ALLOCATED(this%radFiles)) DEALLOCATE (this%radFiles)
    IF (ALLOCATED(this%land)) DEALLOCATE (this%land)

  END SUBROUTINE destructor

END MODULE ReadWriteMWTS_m
