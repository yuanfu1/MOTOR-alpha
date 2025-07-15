!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DP.ReadWriteAGRI
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu, Ting Shu
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2023/7/5, @GBA-MWF, Shenzhen
!   Added wavelet-based satellite data decomposition by Ting Shu (shuting@gbamwf.com), 2023/08/09, @GBA-MWF, Shenzhen 
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE ReadWriteAGRI_m
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
  USE RawAGRI_m, ONLY: RawAGRI_t
  USE Dp_utils_m
  USE RawSatellite_BC_m
  USE RadBase_m, ONLY: RadBase_t
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t

  IMPLICIT NONE

  TYPE, EXTENDS(RadBase_t) :: ReadWriteAGRI_t
    
    CHARACTER(LEN=1024) :: inst_name='agri', platform_name='fy4_1'
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: radFiles ! radiance(tbb) observation files
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: angleFiles ! satellite angle observation files
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: cmFiles   ! L2 cloud mask files

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => agri_setup
    PROCEDURE, PUBLIC  :: ObsIngest => agri_read
    PROCEDURE, PUBLIC  :: ObsQC => agri_qc

    PROCEDURE, PUBLIC :: ObsPrepareForSg 
    PROCEDURE, PUBLIC :: OutputForThinning
    PROCEDURE, PUBLIC :: destructor 
  END TYPE ReadWriteAGRI_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    
    CLASS(ReadWriteAGRI_t) :: this
    TYPE(State_t) :: X
    
    INTEGER(i_kind) :: i, j, istatus

    CALL this%initialize(X, this%bias_scheme, this%bias_predictors, this%CloudyMode, &
                         this%use_cloudy_predictors, this%cloudy_predictors, &
                         this%camin, this%camax, this%obsErrClr, this%obsErrCld, this%BTlim)
    CALL this%ObsPrepareForSgBase(X, this%configFile, this%inst_name, this%platform_name)

  END SUBROUTINE ObsPrepareForSg

  SUBROUTINE OutputForThinning(this, state, mpObs)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    USE mo_netcdf, only: NcDataset, NcDimension
    IMPLICIT NONE

    CLASS(ReadWriteAGRI_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    INTEGER(i_kind) :: istatus

    CHARACTER(len=200) :: filename
    CHARACTER(len=1024) :: outputDir
        
    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    filename = TRIM(outputDir) //'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    PRINT *, 'OUTPUT TO: ', TRIM(filename)
    CALL this%OutputForThinningBase(state, mpObs, filename, this%inst_name, this%platform_name)

  END SUBROUTINE OutputForThinning

  SUBROUTINE agri_qc(this)
    CLASS(ReadWriteAGRI_t) :: this

  END SUBROUTINE agri_qc

  SUBROUTINE agri_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ReadWriteAGRI_t) :: this
    INTEGER(i_kind), PARAMETER :: nfile_max = 1000
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    INTEGER(i_kind) :: numVars
    INTEGER(i_kind) :: i, j, k
    CHARACTER(LEN=1024) :: satinfo_file
    CHARACTER(LEN=2) :: chan_str
    INTEGER(i_kind) :: ichan
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: fdi_file, geo_file, clm_file
    INTEGER(i_kind) :: nfiles_fdi, nfiles_geo, nfiles_clm, ifiles_fdi, ifiles_geo, ifiles_clm
    CHARACTER(LEN=1024), ALLOCATABLE :: satpaths(:)
    CHARACTER(LEN=1024), ALLOCATABLE, DIMENSION(:) :: fdi_filenames, geo_filenames, clm_filenames
    INTEGER :: istat
    INTEGER(i_kind) :: file_id = 18
    INTEGER(i_kind) :: istatus
    REAL(r_kind) :: t1, t2
    INTEGER :: nchans
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: ifuse

    ! (1) Read the namelist file and get necessary parameters
    this%configFile = configFile
    ! 2023-08-30, Ting Shu
    istatus = yaml_get_var(TRIM(configFile), 'Wavelet', 'useWave', this%useWaveFlag)
    CALL CPU_TIME(t1)
    IF(this%useWaveFlag) CALL this%sateDecom%initialize(configFile)
    CALL CPU_TIME(t2)
    ! 2023-10-18, TS, calculate running time
    ! PRINT *, "Initialize running time of wavelet-based satellite thinning:", t2-t1

    istatus = yaml_get_var(TRIM(this%configFile), 'FY4-AGRI', 'satellite', this%platform_name)
    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans)
    nchans = this%nchans

    ALLOCATE (this%chan_lists(nchans), this%rttov_chan_lists(nchans), this%ifuse(nchans), this%oversea(nchans))
    ALLOCATE (this%Err_in(nchans), this%bias_in(nchans), this%camin(nchans), this%camax(nchans))
    ALLOCATE (this%obsErrClr(nchans), this%obsErrCld(nchans), this%BTlim(nchans))

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans, this%chan_lists, this%rttov_chan_lists, &
          this%ifuse, this%Err_in, this%bias_in, this%camin, this%camax, this%obsErrClr, this%obsErrCld, this%BTlim, this%oversea)
    
    IF (yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'ChanList', ifuse) .NE. 0) THEN
      PRINT *, '====== Determine chan usage from the satinfo file ======'
    ELSE
      PRINT *, '====== Determine chan usage from the YAML file ======'
      PRINT *, ifuse
      this%ifuse = ifuse
    END IF
          
    IF(ALLOCATED(ifuse)) DEALLOCATE(ifuse)
    
    this%numVars = this%nchans + 4

    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',       'satpath',         satpaths)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',     'fdi_files',         fdi_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',     'geo_files',         geo_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',     'clm_files',         clm_file)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',   'bias_scheme', this%bias_scheme)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'bias_predictors', this%bias_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'use_cloudy_predictors', this%use_cloudy_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'cloudy_predictors', this%cloudy_predictors)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)
    istatus = yaml_get_var(TRIM(configFile), 'RTTOV', 'rttov_clouds', this%CloudyMode)
    istatus = yaml_get_var(TRIM(this%configFile), 'FY4-AGRI', 'normalize_pred', this%norm_pred)
    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI',   'Source_of_bkg', this%Source_of_bkg)
    istatus = yaml_get_var(TRIM(this%configFile), 'FY4-AGRI', 'Thin_w_bkg', this%Thin_w_bkg)
    istatus = yaml_get_var(TRIM(this%configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! Read in an interpolation option:
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! (2) Prepare those obs files for later reading
    ALLOCATE(fdi_filenames(nfile_max), geo_filenames(nfile_max), clm_filenames(nfile_max))

    ifiles_fdi = 1
    DO WHILE(.TRUE.)
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(fdi_file(1)), iostat = istat)
      read(file_id, *, end=200) fdi_filenames(ifiles_fdi)
      ! PRINT *,'radFiles: ', TRIM(fdi_filenames(ifiles_fdi))
      ifiles_fdi = ifiles_fdi + 1
      IF (istat /= 0) THEN
        PRINT *, 'fdi_file not exit'
        CYCLE
      END IF
    END DO 
    200 CONTINUE
    CLOSE(file_id)
    
    ifiles_geo = 1
    DO WHILE(.TRUE.)
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(geo_file(1)), iostat = istat)
      read(file_id, *, end=201) geo_filenames(ifiles_geo)
      ! PRINT *,'angleFiles: ', TRIM(geo_filenames(ifiles_geo))
      ifiles_geo = ifiles_geo + 1
      IF (istat /= 0) THEN
        PRINT *, 'geo_file not exit'
        CYCLE
      END IF
    END DO 
    201 CONTINUE
    CLOSE(file_id)
    
    ifiles_clm = 1
    DO WHILE(.TRUE.)
      OPEN (unit=file_id, file=TRIM(satpaths(1))//'/'//TRIM(clm_file(1)), iostat = istat)
      read(file_id, *, end=202) clm_filenames(ifiles_clm)
      ! PRINT *,'cmFiles: ', TRIM(clm_filenames(ifiles_clm))
      ifiles_clm = ifiles_clm + 1
      IF (istat /= 0) THEN
        PRINT *, 'clm_file not exit'
        CYCLE
      END IF
    END DO 
    202 CONTINUE
    CLOSE(file_id)

    nfiles_fdi = ifiles_fdi - 1
    nfiles_geo = ifiles_geo - 1
    nfiles_clm = ifiles_clm - 1
    ! TODO: remove those unmatched files in SHELL
    IF ( nfiles_fdi .NE. nfiles_geo .OR. nfiles_fdi .NE. nfiles_clm) THEN
      PRINT *, 'Numbers of FDI/GEO/CLM files should be identical'
      STOP 
    END IF 
    this%numfiles = nfiles_fdi
    ! this%numfiles = 1 ! Test only
    IF (.NOT. ALLOCATED(this%radFiles)) &
    ALLOCATE(this%radFiles(this%numFiles), this%angleFiles(this%numFiles), this%cmFiles(this%numFiles))
    DO i = 1, this%numfiles
      this%radFiles(i) = TRIM(satpaths(1))//'/'//fdi_filenames(i)
      this%angleFiles(i) = TRIM(satpaths(1))//'/'//geo_filenames(i)
      this%cmFiles(i) = TRIM(satpaths(1))//'/'//clm_filenames(i)
      ! PRINT *, this%radFiles(i), this%angleFiles(i), this%cmFiles(i)
    END DO

    IF (ALLOCATED(fdi_file)) DEALLOCATE(fdi_file)
    IF (ALLOCATED(geo_file)) DEALLOCATE(geo_file)
    IF (ALLOCATED(clm_file)) DEALLOCATE(clm_file)
    IF (ALLOCATED(fdi_filenames)) DEALLOCATE(fdi_filenames)
    IF (ALLOCATED(geo_filenames)) DEALLOCATE(geo_filenames)
    IF (ALLOCATED(clm_filenames)) DEALLOCATE(clm_filenames)  
    DEALLOCATE(satpaths)  
    PRINT *, "agri_setup successfully run"

  END SUBROUTINE agri_setup

  SUBROUTINE agri_read(this, X)
    USE AdvanceTime_m, ONLY: Time_GMT_to_Unix
    ! 2023-08-31, Ting Shu, get larger latlon domain
    USE satelliteW_m, ONLY: get_latlon
    USE domainCheck_m, ONLY: domainCheck_t
    IMPLICIT NONE
    CLASS(ReadWriteAGRI_t) :: this
    TYPE(State_t) :: X
    ! CLASS(RawAGRI_t), ALLOCATABLE :: RawAGRI
    TYPE(RawAGRI_t) :: RawAGRI
    REAL(r_single), ALLOCATABLE :: sat_zenith(:), sat_azi(:)
    INTEGER(i_kind), ALLOCATABLE :: pred_use(:,:)

    INTEGER(i_kind) :: nobs, ichan, i, j, iobs
    CHARACTER(LEN=2)  :: chan_name
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: nRawObs_total, nRawObs, ivar
    INTEGER(i_kind) :: nv,nc,numgrd(2),fileunit
    CHARACTER(len=1024) :: filename_bias_coefs
    REAL(r_kind) :: domain(2,2), domainL(2,2), dlat, dlon
    REAL(r_kind) :: tt1, tt2
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
    REAL(r_kind), ALLOCATABLE :: bcoef(:,:)
    INTEGER(i_kind) :: npred_total

    ! 2023-08-31, Ting Shu, get larger latlon domain for wavelet
    !IF (this%useWaveFlag) THEN
    !  BLOCK
    !    INTEGER(i_kind) :: status
    !    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: latlon_times
    !    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: latlon_domain, fineLat, fineLon

        ! Get latitude and longitude time and domain from yaml
    !    status = yaml_get_var(TRIM(this%configFile), "Wavelet", "latlon_times", latlon_times)
    !    status = yaml_get_var(TRIM(this%configFile), "latlon_grid", "domain_latlon", latlon_domain)

        ! Get finest grid latitude and longitude
    !    CALL get_latlon(this%mgStart, this%mgEnd, latlon_times, latlon_domain, fineLat, fineLon)
    !    domain = RESHAPE((/fineLat(1), fineLat(SIZE(fineLat)), fineLon(1), fineLon(SIZE(fineLon))/), (/2, 2/))
    !    domain = domain * degree2radian

    !    dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    !    dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)
        ! PRINT *, 'dlat/dlon: ', dlat/degree2radian, dlon/degree2radian
        ! We do not believe those obs fall on bdys
    !    domain(1,1) = domain(1,1) + dlat*0.5D0
    !    domain(2,1) = domain(2,1) - dlat*0.5D0
    !    domain(1,2) = domain(1,2) + dlon*0.5D0
    !    domain(2,2) = domain(2,2) - dlon*0.5D0
        
    !    IF(ALLOCATED(latlon_times)) DEALLOCATE(latlon_times, STAT=status)
    !    IF(ALLOCATED(latlon_domain)) DEALLOCATE(latlon_domain, STAT=status)
    !    IF(ALLOCATED(fineLat)) DEALLOCATE(fineLat, STAT=status)
    !    IF(ALLOCATED(fineLon)) DEALLOCATE(fineLon, STAT=status)
    !  END BLOCK
    !ELSE
      ! Get the latlon grid information:
      ! CALL getDomain_latlon(this%configFile, nv, nc, numgrd, domainL)
      ! domainL = domainL * degree2radian
    domain(1,1) = MINVAL(X%sg%cell_cntr(1, :))
    domain(2,1) = MAXVAL(X%sg%cell_cntr(1, :))
    domain(1,2) = MINVAL(X%sg%cell_cntr(2, :))
    domain(2,2) = MAXVAL(X%sg%cell_cntr(2, :))

    dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)
    ! PRINT *, 'dlat/dlon: ', domain(1,1), dlat/degree2radian, dlon/degree2radian
      ! We do not believe those obs fall on bdys
    domain(1,1) = domain(1,1) + dlat*0.5D0
    domain(2,1) = domain(2,1) - dlat*0.5D0
    domain(1,2) = domain(1,2) + dlon*0.5D0
    domain(2,2) = domain(2,2) - dlon*0.5D0
    ! PRINT *, 'get domain ', domain(1,1), domain(2,1), domain(1,2), domain(2,2)
    ! END IF
    
    ! IF (this%platform_name == 'fy4_1') THEN
    !   ALLOCATE (RawAGRI_4a_t:: RawAGRI)
    !   PRINT *, '===== FY4A IS USED ====='
    ! ELSEIF (this%platform_name == 'fy4_2') THEN
    !   ALLOCATE (RawAGRI_4b_t:: RawAGRI)
    !   PRINT *, '===== FY4B IS USED ====='
    ! ELSE
    !   PRINT*, '===== Platform of FY4 is not defined ====='
    !   STOP
    ! END IF

    nRawObs_total = 0
    CALL cpu_time(tt1)
    IF (this%numFiles > 0) THEN
      CALL RawAGRI%RawAGRI_setup(this%configFile)
      ! PRINT *, 'RawAGRI_read: ', TRIM(this%radFiles(i))
      CALL RawAGRI%RawAGRI_read(TRIM(this%radFiles(1)), TRIM(this%angleFiles(1)), TRIM(this%cmFiles(1)))

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawAGRI%obsvalue,RawAGRI%nchans,RawAGRI%ifuse)
      CALL Select_InDomain(RawAGRI%obsvalue,RawAGRI%nRawObs,RawAGRI%latitude,RawAGRI%longitude,    &
                           domain(1,1), domain(2,1), domain(1,2), domain(2,2))
      ! CALL Select_InTimeWin(RawAGRI%obsvalue,RawAGRI%nRawObs,RawAGRI%obstime,X%sg%tt(1), X%sg%tt(X%sg%tSlots))
      ! PRINT *, 'END Select OBS'

      ! TODO: Need to set num_obs_total for each channel
      DO iobs = 1, RawAGRI%nRawObs
        IF ( ANY(RawAGRI%obsvalue(iobs, :) .LT. missing ) .AND. &
        RawAGRI%latitude(iobs) .LT. missing .AND. &
        RawAGRI%longitude(iobs) .LT. missing ) THEN 
          nRawObs_total = nRawObs_total + 1
        END IF
      END DO

      CALL RawAGRI%Destroy_RawAGRI()
      nRawObs_total = nRawObs_total * this%numFiles
    END IF
    CALL cpu_time(tt2)
    PRINT *, 'first read time is ', tt2-tt1
      
    ! IF(X%sg%isBaseProc()) PRINT *, 'nRawObs_total = ', nRawObs_total

    ! ! For satellite observations, number of variables is : nchans + four angles
    ! ! this%numVars = RawAGRI%nchans + 4 
    ! ! ObsType: 'fy4_1_agri'
    ! ! VarType: 'tbb'

    this%obsType = TRIM(RawAGRI%platform_name)//'_'//TRIM(RawAGRI%inst_name)
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
    CALL cpu_time(tt1)
    DO i = 1, this%numFiles

      CALL RawAGRI%RawAGRI_setup(this%configFile)
      ! PRINT *, 'RawAGRI_read: ', TRIM(this%radFiles(i))
      CALL RawAGRI%RawAGRI_read(TRIM(this%radFiles(i)), TRIM(this%angleFiles(i)), TRIM(this%cmFiles(i)))
      PRINT *, i, ' file has been read in '

      ! PRINT *, 'START Select OBS'
      ! CALL Select_active_channels(RawAGRI%obsvalue,RawAGRI%nchans,RawAGRI%ifuse)
      nRawObs = SIZE(RawAGRI%latitude,1)
      ALLOCATE(org_olatlon(2, nRawObs))
      ALLOCATE(org_obsHght(nRawObs))
      ALLOCATE (maskValided(nRawObs))
      org_olatlon(1,:) = RawAGRI%latitude(:)
      org_olatlon(2,:) = RawAGRI%longitude(:)
      org_obsHght(:) = -1.0E8
      CALL domainCheck%validation(X%sg, nRawObs, org_olatlon, org_obsHght, &
                          RawAGRI%obstime, numValided, maskValided, &
                          Hidxtgt, Vidxtgt, Tidxtgt, &
                          Hcoetgt, Vcoetgt, Tcoetgt, &
                          1, nhstencil, 'Satellite', &
                          domain(1,1), domain(2,1), domain(1,2), domain(2,2))

      DO iobs = 1, nRawObs

        IF (maskValided(iobs) .GT. 0) THEN
          nobs = nobs + 1
          obsTime(nobs) = RawAGRI%obstime(iobs)
          olatlon(1,nobs) = RawAGRI%latitude(iobs)
          olatlon(2,nobs) = RawAGRI%longitude(iobs)

          DO ichan = 5, this%numVars
            obsData(nobs,ichan) = RawAGRI%obsvalue(iobs, ichan-4)
            obsErrs(nobs,ichan) = RawAGRI%obserr(ichan-4)
          END DO

          cloud_flag(nobs) = RawAGRI%cloud_flag(iobs)
          obsTime(nobs)    = RawAGRI%obstime(iobs)
          olatlon(1,nobs) = RawAGRI%latitude(iobs)
          olatlon(2,nobs) = RawAGRI%longitude(iobs)

          ! Set for four satellite angles
          Angles(1, nobs) = RawAGRI%sat_zenith(iobs)
          Angles(2, nobs) = RawAGRI%sat_azi(iobs)
          Angles(3, nobs) = RawAGRI%sol_zenith(iobs)
          Angles(4, nobs) = RawAGRI%sol_azi(iobs)

          iscanpos(nobs) = RawAGRI%iscanpos(iobs)
        END IF

      END DO

      DEALLOCATE(org_olatlon, org_obsHght)
      DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided)
      CALL RawAGRI%Destroy_RawAGRI()
    END DO
    CALL cpu_time(tt2)
    PRINT *, 'second read time is ', tt2-tt1

    ! Set for angles
    DO ichan = 1, 4
      obsData(:,ichan) = Angles(ichan, :)
    END DO

    ! Remove missing values
    this%numObs = nobs
    print *, 'numObs = ', this%numObs

    ALLOCATE(this%obsHght(this%numObs))
    ALLOCATE(this%ca_mean(this%numObs, this%numVars))
    this%ca_mean = missing
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
      IF ( this%use_cloudy_predictors ) THEN
        npred_total = this%npred+SIZE(this%cloudy_predictors)
      ELSE
        npred_total = this%npred
      END IF
      PRINT *, 'npred_total = ', npred_total
      CALL Read_bias_coefs_WRFDA(fileunit, filename_bias_coefs, npred_total, pred_use, bcoef)
      PRINT *, 'shape of pred_use: ', shape(pred_use)
      PRINT *, 'Check coefs values: ', shape(bcoef), maxval(bcoef), minval(bcoef)
    
      DO i = 1, SIZE(pred_use,1)
        DO j = 1, SIZE(pred_use,2)
          IF (pred_use (i, j) .EQ. 0) bcoef(i,j) = ZERO
        END DO
      END DO

      ALLOCATE(this%BCOEF(SIZE(bcoef,1), npred_total))
      this%BCOEF(:, 1:npred_total) = bcoef

      PRINT *, 'Clear predictors: ', this%npred
      IF ( this%use_cloudy_predictors ) this%npred = npred_total
      PRINT *, 'Clear + Cloudy predictors: ', this%npred

      IF (ALLOCATED(pred_use)) DEALLOCATE(pred_use)
      IF (ALLOCATED(bcoef)) DEALLOCATE(bcoef)
    END IF

    ! PRINT *, 'agri_read is successfully run'

  END SUBROUTINE agri_read

  SUBROUTINE destructor(this)
    CLASS(ReadWriteAGRI_t) :: this

    CALL this%destroy()
    IF (ALLOCATED(this%cmFiles)) DEALLOCATE (this%cmFiles)
    IF (ALLOCATED(this%angleFiles)) DEALLOCATE (this%angleFiles)
    IF (ALLOCATED(this%radFiles)) DEALLOCATE (this%radFiles)
    ! PRINT *, 'Everything is released'

  END SUBROUTINE destructor

END MODULE ReadWriteAGRI_m
