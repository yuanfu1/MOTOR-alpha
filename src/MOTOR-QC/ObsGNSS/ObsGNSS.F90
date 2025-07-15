!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.OBSGNSS
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yongjian Huang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yongjian Huang, 2024/9/9, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite GNSS-RO data structure.
MODULE ObsGNSS_m
  USE ObsBase_m, ONLY: ObsBase_t
  USE kinds_m, ONLY: r_kind, i_kind
  USE parameters_m, ONLY: degree2radian, missing, invalid, machineEps
  USE State_m, ONLY: State_t
  USE YAMLRead_m, ONLY: yaml_get_var
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE domainCheck_m, ONLY: domainCheck_t
  USE AdvanceTime_m, ONLY: Get_Time_GMT_to_Unix
  USE ObsField_m, ONLY: ObsField_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE mpObs_m, ONLY: mpObs_t
  USE slint, ONLY: slint_init, tgt_grid

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsGNSS_t
    INTEGER(i_kind) :: mgStart, mgEnd
    INTEGER(i_kind) :: ts_current
    LOGICAL :: set_all_gnss_to_fcst_time

  CONTAINS

    PROCEDURE, PUBLIC  :: ObsInitial => gnss_setup
    PROCEDURE, PUBLIC  :: ObsIngest => gnss_read
    PROCEDURE, PUBLIC  :: ObsQC => gnss_qc
    PROCEDURE :: ObsForward => gnssForward
    PROCEDURE :: ObsTangent => gnssTangent
    PROCEDURE :: ObsAdjoint => gnssAdjoint

    PROCEDURE, PUBLIC :: ObsPrepareForSg
    ! PROCEDURE, PUBLIC :: ObsThinning

    PROCEDURE, PRIVATE :: filelist_domain_check, count_obs_in_filelist, &
      read_data_from_filelist, gen_final_obs, set_to_obsbase, gnss_refrac_qc

  END TYPE ObsGNSS_t

CONTAINS

  SUBROUTINE gnss_setup(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsGNSS_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind), OPTIONAL :: idxFile

    INTEGER(i_kind) :: i, istatus, vLevel, nObs

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileList(:)
    INTEGER(i_kind), ALLOCATABLE :: timeArray(:)

    CALL this%obsCommon%initialize(configFile)

    ! data root dir
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'input_dir_GNSSRO', obsFileDir)
    PRINT *, 'IO.input_dir_GNSSRO:', TRIM(obsFileDir)

    ! num of files
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_GNSSRO', obsFileList)
    this%numFiles = UBOUND(obsFileList, 1)
    WRITE (*, 2) this%numFiles
2   FORMAT("HYJ++: Number of gnss ro obs files: ", I3)

    ! filenames
    ALLOCATE (this%fileNames(this%numFiles))
    DO i = 1, this%numFiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileList(i))
      PRINT *, "HYJ++: name of gnss ro obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

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

    istatus = yaml_get_var(TRIM(configFile), 'analysis_para', 'end_time', timeArray)
    WRITE (*, 23) timeArray(1: UBOUND(timeArray, 1))
23   FORMAT('END TIME from yaml: ', 10I6)

    this%ts_current = Get_Time_GMT_to_Unix(timeArray)

    istatus = yaml_get_var(TRIM(configFile), 'IO', 'set_all_gnss_to_fcst_time', this%set_all_gnss_to_fcst_time)
    
    IF (istatus .EQ. 0) THEN
      WRITE (*, 24) this%set_all_gnss_to_fcst_time
  24   FORMAT('SET ALL GNSS TO FCRECAST TIME: ', L2  )
    ELSE
      this%set_all_gnss_to_fcst_time = .TRUE.
      WRITE (*, 24) this%set_all_gnss_to_fcst_time
  25  FORMAT('SET ALL GNSS TO FCRECAST TIME: ', L2  )
    END IF

    DEALLOCATE(timeArray)

  END SUBROUTINE gnss_setup

  SUBROUTINE gnss_read(this, X)
    ! inputs
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: X

    ! for step 1
    LOGICAL, ALLOCATABLE :: in_domain_file_mask(:)

    ! for step 2
    INTEGER(i_kind) :: nobs

    ! for step 3
    REAL(r_kind), ALLOCATABLE :: time_tp(:), lat_tp(:), lon_tp(:), alt_tp(:), refrac(:)
    LOGICAL, ALLOCATABLE :: qc_flag(:)

    ! for step 4
    REAL(r_kind), ALLOCATABLE :: time_tp_fnl(:), lat_tp_fnl(:), lon_tp_fnl(:), alt_tp_fnl(:), refrac_fnl(:)

    INTEGER(i_kind) :: nobs_valid

    !! 1. domain check filterd file list
    CALL this%filelist_domain_check(X, in_domain_file_mask)

    !! 2. count num obs of filtered file
    CALL this%count_obs_in_filelist(in_domain_file_mask, nobs)

    !! 3. read data
    CALL this%read_data_from_filelist(in_domain_file_mask, nobs, time_tp, lat_tp, lon_tp, alt_tp, refrac, qc_flag)

    !! 4. final obs for obsbase
    CALL this%gen_final_obs(X, nobs, time_tp, lat_tp, lon_tp, alt_tp, refrac, qc_flag, time_tp_fnl, lat_tp_fnl, lon_tp_fnl, alt_tp_fnl, refrac_fnl, nobs_valid)

    !! 5. set to obsbase
    CALL this%set_to_obsbase(nobs_valid, time_tp_fnl, lat_tp_fnl, lon_tp_fnl, alt_tp_fnl, refrac_fnl)

    IF (ALLOCATED(time_tp)) DEALLOCATE (time_tp)
    IF (ALLOCATED(lat_tp)) DEALLOCATE (lat_tp)
    IF (ALLOCATED(lon_tp)) DEALLOCATE (lon_tp)
    IF (ALLOCATED(alt_tp)) DEALLOCATE (alt_tp)
    IF (ALLOCATED(refrac)) DEALLOCATE (refrac)
    IF (ALLOCATED(qc_flag)) DEALLOCATE (qc_flag)
    IF (ALLOCATED(time_tp_fnl)) DEALLOCATE (time_tp_fnl)
    IF (ALLOCATED(lat_tp_fnl)) DEALLOCATE (lat_tp_fnl)
    IF (ALLOCATED(lon_tp_fnl)) DEALLOCATE (lon_tp_fnl)
    IF (ALLOCATED(alt_tp_fnl)) DEALLOCATE (alt_tp_fnl)
    IF (ALLOCATED(refrac_fnl)) DEALLOCATE (refrac_fnl)

  END SUBROUTINE gnss_read

  SUBROUTINE filelist_domain_check(this, X, in_domain_file_mask)
    ! INPUTS
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: in_domain_file_mask(:)

    ! for nc operations
    TYPE(NcDataset)   :: nc
    TYPE(NcVariable)   :: nc_var

    ! for domain checks
    REAL(r_kind) :: minLat, minLon, maxLat, maxLon
    REAL(r_kind) :: obs_lat, obs_lon, obs_lat_rad, obs_lon_rad
    INTEGER(i_kind) :: i

    ALLOCATE (in_domain_file_mask(this%numFiles))

    in_domain_file_mask = .FALSE.

    DO i = 1, this%numFiles
      nc = NcDataset(this%fileNames(i), "r")

      !! read lat lon from nc file and prepare for domain check
      IF ((INDEX(this%fileNames(i), 'metop') .GT. 0) .OR. (INDEX(this%fileNames(i), 'spire') .GT. 0)) THEN

        nc_var = nc%getVariable('lat'); CALL nc_var%getData(obs_lat)
        nc_var = nc%getVariable('lon'); CALL nc_var%getData(obs_lon)

      ELSE IF ((INDEX(this%fileNames(i), 'cosmic-2') .GT. 0)) THEN

        CALL nc%getAttribute('lat', obs_lat)
        CALL nc%getAttribute('lon', obs_lon)
      ELSE IF ((INDEX(this%fileNames(i), 'tianmu') .GT. 0)) THEN
        CALL nc%getAttribute('lat', obs_lat)
        CALL nc%getAttribute('lon', obs_lon)
      END IF

      ! convert to rad
      obs_lat_rad = obs_lat * degree2radian
      obs_lon_rad = obs_lon * degree2radian

      !! read domain boundary
      minLat = MINVAL(X%sg%cell_cntr(1, :))
      maxLat = MAXVAL(X%sg%cell_cntr(1, :))
      minLon = MINVAL(X%sg%cell_cntr(2, :))
      maxLon = MAXVAL(X%sg%cell_cntr(2, :))

      !! preliminary in domain check
      !! points of profile may still out of the boundary
      IF (obs_lat_rad .GE. minLat .AND. obs_lat_rad .LE. maxLat .AND. &
          obs_lon_rad .GE. minLon .AND. obs_lon_rad .LE. maxLon) THEN
        in_domain_file_mask(i) = .TRUE.
      END IF

      CALL nc%CLOSE()
    END DO

  END SUBROUTINE filelist_domain_check

  !*********************************************************
  !! 
  !!
  !!
  !!
  !!
  !*********************************************************
  SUBROUTINE count_obs_in_filelist(this, in_domain_file_mask, nobs)
    CLASS(ObsGNSS_t) :: this
    LOGICAL, DIMENSION(:), INTENT(IN) :: in_domain_file_mask

    INTEGER(i_kind), INTENT(OUT) :: nobs

    TYPE(NcDataset)   :: nc
    TYPE(NcDimension)  :: dim_vert
    INTEGER(i_kind) :: nobs_pfl
    INTEGER(i_kind) :: i

    nobs = 0
    DO i = 1, this%numFiles
      nobs_pfl = 0
      IF (in_domain_file_mask(i)) THEN
        nc = NcDataset(this%fileNames(i), "r")

        IF ((INDEX(this%fileNames(i), 'metop') .GT. 0) .OR. (INDEX(this%fileNames(i), 'spire') .GT. 0)) THEN

          ! get num pts
          dim_vert = nc%getDimension('dim_lev1b')
          nobs_pfl = dim_vert%getLength()

        ELSE IF ((INDEX(this%fileNames(i), 'cosmic-2') .GT. 0)) THEN

          dim_vert = nc%getDimension('MSL_alt')
          nobs_pfl = dim_vert%getLength()
        ELSE IF ((INDEX(this%fileNames(i), 'tianmu') .GT. 0)) THEN
          ! get num pts
          dim_vert = nc%getDimension('dim_lev2a')
          nobs_pfl = dim_vert%getLength()
        END IF
        nobs = nobs + nobs_pfl

        CALL nc%CLOSE()
      END IF

    END DO
  END SUBROUTINE count_obs_in_filelist

  SUBROUTINE read_data_from_filelist(this, in_domain_file_mask, nobs, time_tp, lat_tp, lon_tp, alt_tp, refrac, qc_flag)
    ! INPUTS OUTPUTS
    CLASS(ObsGNSS_t) :: this
    LOGICAL, DIMENSION(:), INTENT(IN) :: in_domain_file_mask
    INTEGER(i_kind), INTENT(IN) :: nobs

    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: time_tp(:), lat_tp(:), lon_tp(:), alt_tp(:), refrac(:)
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: qc_flag(:)

    ! for netcdf operations
    TYPE(NcDataset)   :: nc
    TYPE(NcVariable)   :: nc_time_tp, nc_lat_tp, nc_lon_tp, nc_alt_tp, nc_refrac
    TYPE(NcVariable) :: nc_var_tmp
    TYPE(NcDimension)  :: dim_vert


    ! tmp data
    REAL(r_kind), ALLOCATABLE :: tmp_1d(:), tmp_alt(:), tmp_refrac(:)
    LOGICAL, ALLOCATABLE :: tmp_mask(:)
    INTEGER(i_kind) :: year, month, day, hour, minute, second, ts_unix
    REAL(r_kind) :: lat_fix, lon_fix
    INTEGER(i_kind) :: i
    INTEGER(i_kind) :: obs_idx, obs_len
    LOGICAL, ALLOCATABLE :: prod_qc(:)
    INTEGER(i_kind) :: blacklist_flag



    ALLOCATE (time_tp(nobs), lat_tp(nobs), lon_tp(nobs), alt_tp(nobs), refrac(nobs), qc_flag(nobs))
    obs_idx = 1
    obs_len = 0
    DO i = 1, this%numFiles

      IF (in_domain_file_mask(i)) THEN
        nc = NcDataset(this%fileNames(i), "r")
        PRINT *, 'read gnss ro file: ', TRIM(this%fileNames(i))
        IF ((INDEX(this%fileNames(i), 'metop') .GT. 0) .OR. (INDEX(this%fileNames(i), 'spire') .GT. 0)) THEN

          ! use the same time
          nc_var_tmp = nc%getVariable('year'); CALL nc_var_tmp%getData(year)
          nc_var_tmp = nc%getVariable('month'); CALL nc_var_tmp%getData(month)
          nc_var_tmp = nc%getVariable('day'); CALL nc_var_tmp%getData(day)
          nc_var_tmp = nc%getVariable('hour'); CALL nc_var_tmp%getData(hour)
          nc_var_tmp = nc%getVariable('minute'); CALL nc_var_tmp%getData(minute)
          nc_var_tmp = nc%getVariable('second'); CALL nc_var_tmp%getData(second)

          IF (this%set_all_gnss_to_fcst_time .EQV. .TRUE. ) THEN
            ts_unix = this%ts_current
          ELSE
            ts_unix = Get_Time_GMT_to_Unix((/year, month, day, hour, minute, second/))
          ENDIF
          nc_lat_tp = nc%getVariable('lat_tp')
          nc_lon_tp = nc%getVariable('lon_tp')
          nc_alt_tp = nc%getVariable('alt_refrac')
          nc_refrac = nc%getVariable('refrac')

          dim_vert = nc%getDimension('dim_lev1b')
          obs_len = dim_vert%getLength()

          time_tp(obs_idx:obs_idx + obs_len - 1) = ts_unix

          ALLOCATE (tmp_1d(obs_len))

          CALL nc_lat_tp%getData(tmp_1d)
          lat_tp(obs_idx:obs_idx + obs_len - 1) = tmp_1d
          CALL nc_lon_tp%getData(tmp_1d)
          lon_tp(obs_idx:obs_idx + obs_len - 1) = tmp_1d
          CALL nc_alt_tp%getData(tmp_alt)
          alt_tp(obs_idx:obs_idx + obs_len - 1) = tmp_alt
          CALL nc_refrac%getData(tmp_refrac)
          refrac(obs_idx:obs_idx + obs_len - 1) = tmp_refrac

          !! no blacklist
          blacklist_flag = 0

        ELSE IF ((INDEX(this%fileNames(i), 'cosmic-2') .GT. 0)) THEN

          ! use the same time
          CALL nc%getAttribute('year', year)
          CALL nc%getAttribute('month', month)
          CALL nc%getAttribute('day', day)
          CALL nc%getAttribute('hour', hour)
          CALL nc%getAttribute('minute', minute)
          CALL nc%getAttribute('second', second)

          ts_unix = Get_Time_GMT_to_Unix((/year, month, day, hour, minute, second/))

          nc_lat_tp = nc%getVariable('Lat')
          nc_lon_tp = nc%getVariable('Lon')
          nc_alt_tp = nc%getVariable('MSL_alt')
          nc_refrac = nc%getVariable('Ref')

          dim_vert = nc%getDimension('MSL_alt')
          obs_len = dim_vert%getLength()

          time_tp(obs_idx:obs_idx + obs_len - 1) = ts_unix

          ALLOCATE (tmp_1d(obs_len))

          CALL nc_lat_tp%getData(tmp_1d)
          lat_tp(obs_idx:obs_idx + obs_len - 1) = tmp_1d
          CALL nc_lon_tp%getData(tmp_1d)
          lon_tp(obs_idx:obs_idx + obs_len - 1) = tmp_1d
          CALL nc_alt_tp%getData(tmp_alt)
          tmp_alt  = tmp_alt * 1000D0
          alt_tp(obs_idx:obs_idx + obs_len - 1) = tmp_alt
          CALL nc_refrac%getData(tmp_refrac)
          refrac(obs_idx:obs_idx + obs_len - 1) = tmp_refrac

          CALL nc%getAttribute('irs', blacklist_flag)

          IF ( blacklist_flag == -1) THEN
            blacklist_flag = -1
          END IF 

        ELSE IF ((INDEX(this%fileNames(i), 'tianmu') .GT. 0)) THEN
          ! use the same time
          CALL nc%getAttribute('year', year)
          CALL nc%getAttribute('month', month)
          CALL nc%getAttribute('day', day)
          CALL nc%getAttribute('hour', hour)
          CALL nc%getAttribute('minute', minute)
          CALL nc%getAttribute('second', second)

          ts_unix = Get_Time_GMT_to_Unix((/year, month, day, hour, minute, second/))

          CALL nc%getAttribute('lat', lat_fix)
          CALL nc%getAttribute('lon', lon_fix)
          nc_alt_tp = nc%getVariable('MSL_alt_Ref')
          nc_refrac = nc%getVariable('Ref')

          dim_vert = nc%getDimension('dim_lev2a')
          obs_len = dim_vert%getLength()

          time_tp(obs_idx:obs_idx + obs_len - 1) = ts_unix

          ALLOCATE (tmp_1d(obs_len), tmp_alt(obs_len), tmp_refrac(obs_len))


          lat_tp(obs_idx:obs_idx + obs_len - 1) = lat_fix
   
          lon_tp(obs_idx:obs_idx + obs_len - 1) = lon_fix
          CALL nc_alt_tp%getData(tmp_alt)
          tmp_alt  = tmp_alt * 1000D0
          alt_tp(obs_idx:obs_idx + obs_len - 1) = tmp_alt
          CALL nc_refrac%getData(tmp_refrac)
          refrac(obs_idx:obs_idx + obs_len - 1) = tmp_refrac

          ! tmp store the setting attr
          CALL nc%getAttribute('setting', blacklist_flag)
          IF (blacklist_flag .EQ. 1) THEN
            blacklist_flag = -2
          END IF 
        END IF

        ALLOCATE(tmp_mask(obs_len))
 
        CALL this%gnss_refrac_qc(tmp_alt, tmp_refrac, prod_qc, obs_len, blacklist_flag, tmp_mask)
        qc_flag(obs_idx:obs_idx + obs_len - 1) = tmp_mask

        obs_idx = obs_idx + obs_len

        IF (ALLOCATED(tmp_1d)) DEALLOCATE (tmp_1d)
        IF (ALLOCATED(tmp_mask)) DEALLOCATE (tmp_mask)
        IF (ALLOCATED(tmp_alt)) DEALLOCATE (tmp_alt)
        IF (ALLOCATED(tmp_refrac)) DEALLOCATE (tmp_refrac)
        CALL nc%CLOSE()
      END IF

    END DO

    PRINT *, 'HYJ+++: total obs: ', obs_idx

  END SUBROUTINE read_data_from_filelist

  SUBROUTINE gen_final_obs(this, X, nobs, time_tp, lat_tp, lon_tp, alt_tp, refrac, qc_flag,  time_tp_fnl, lat_tp_fnl, lon_tp_fnl, alt_tp_fnl, refrac_fnl, nobs_valid)
    ! INPUT OUTPUS
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: X

    INTEGER(i_kind), INTENT(IN) :: nobs
    REAL(r_kind), INTENT(IN) :: time_tp(:), lat_tp(:), lon_tp(:), alt_tp(:), refrac(:)
    LOGICAL, INTENT(IN) :: qc_flag(:)

    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: time_tp_fnl(:), lat_tp_fnl(:), lon_tp_fnl(:), alt_tp_fnl(:), refrac_fnl(:)
    INTEGER(i_kind), INTENT(OUT) :: nobs_valid


    LOGICAL, ALLOCATABLE :: final_obs_mask(:)
    INTEGER(i_kind) :: i, j

    REAL(r_kind) :: minLat, minLon, maxLat, maxLon
    REAL(r_kind) :: obs_lat, obs_lon, obs_lat_rad, obs_lon_rad

    ALLOCATE (final_obs_mask(nobs))

    !! read domain boundary
    minLat = MINVAL(X%sg%cell_cntr(1, :))
    maxLat = MAXVAL(X%sg%cell_cntr(1, :))
    minLon = MINVAL(X%sg%cell_cntr(2, :))
    maxLon = MAXVAL(X%sg%cell_cntr(2, :))

    final_obs_mask(:) = .FALSE.

    nobs_valid = 0
    DO i = 1, nobs

      obs_lat_rad = lat_tp(i) * degree2radian
      obs_lon_rad = lon_tp(i) * degree2radian

      ! domain check
      IF (obs_lat_rad .GE. minLat .AND. obs_lat_rad .LE. maxLat .AND. &
          obs_lon_rad .GE. minLon .AND. obs_lon_rad .LE. maxLon) THEN
        final_obs_mask(i) = .TRUE.
      END IF

      ! value check
      IF ((refrac(i) .LT. 0 .OR. refrac(i) .GT. 500)) THEN
        final_obs_mask(i) = .FALSE.
      END IF

      ! qc check
      IF ( .NOT. qc_flag(i) ) THEN
        final_obs_mask(i) = .FALSE.
      END IF 

      IF (final_obs_mask(i)) THEN
        nobs_valid = nobs_valid + 1
      END IF
    END DO

    ! PRINT *, 'HYJ+++, nobs_valid, nobs', nobs_valid, nobs
    ALLOCATE (time_tp_fnl(nobs_valid), lat_tp_fnl(nobs_valid), lon_tp_fnl(nobs_valid), alt_tp_fnl(nobs_valid), refrac_fnl(nobs_valid))

    j = 1
    DO i = 1, nobs

      IF (final_obs_mask(i)) THEN

        !! TODO: STUB : BEST DATA INPUT FOR TEST
        ! IF ( alt_tp(i) .GT. 10000 .AND. alt_tp(i) .LT. 20000 ) THEN

        time_tp_fnl(j) = time_tp(i)
        lat_tp_fnl(j) = lat_tp(i) * degree2radian
        lon_tp_fnl(j) = lon_tp(i) * degree2radian
        alt_tp_fnl(j) = alt_tp(i)
        refrac_fnl(j) = refrac(i)

        j = j + 1
        ! END IF
      END IF

    END DO

    ! PRINT *, 'HYJ+++ CHECKS', 'j: ', j, 'nobs_valid: ', nobs_valid

    DEALLOCATE (final_obs_mask)

  END SUBROUTINE gen_final_obs



  !> 
  SUBROUTINE set_to_obsbase(this, nobs_valid, time_tp_fnl, lat_tp_fnl, lon_tp_fnl, alt_tp_fnl, refrac_fnl)
    ! INPUT OUTPUS
    CLASS(ObsGNSS_t) :: this
    REAL(r_kind), INTENT(IN) :: time_tp_fnl(:), lat_tp_fnl(:), lon_tp_fnl(:), alt_tp_fnl(:), refrac_fnl(:)
    INTEGER(i_kind), INTENT(IN) :: nobs_valid

    INTEGER(i_kind) :: i
    INTEGER(i_kind) :: numObs(this%numVars)

    this%numVars = 1

    ALLOCATE (this%varNames(this%numVars))
    this%varNames(1) = "refractivity"

    this%obsType = "GNSSRO"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius = 5.0D0
    this%radius(3, :) = 5.0D2
    this%radius(4, :) = 5.0D1
    this%sizeInc = 10.0D0

    this%correlation_threshold = 0.25D0
    this%qcThreshold(:) = 3.0D0

    this%numObs = nobs_valid

    ALLOCATE (this%obsData(this%numObs, this%numVars), &
              this%obsErrs(this%numObs, this%numVars), &
              this%olatlon(2, this%numObs), &
              this%obsTime(this%numObs), this%obsHght(this%numObs), &
              this%land(this%numObs), this%StNames(this%numObs))

    this%land(:) = 1.0D0
    this%StNames(:) = 'gnssro'

    this%obsHght(:) = alt_tp_fnl
    this%obsTime(:) = time_tp_fnl
    this%olatlon(1, :) = lat_tp_fnl
    this%olatlon(2, :) = lon_tp_fnl
    this%obsData(:, 1) = refrac_fnl
    

    DO i = 1, this%numObs
      IF ( alt_tp_fnl(i) .LT. 10000 ) THEN
        this%obsErrs(i, 1) = refrac_fnl(i) * (0.01D0 + (1e4_r_kind - alt_tp_fnl(i)) * 0.04e4_r_kind)
      ELSE
        this%obsErrs(i, 1) = refrac_fnl(i) * 0.01D0
      END IF

      ! PRINT*, 'HYJ+++, this%obsErrs(i, 1) = ', this%obsErrs(i, 1), & 
      !     this%olatlon(1, i), this%olatlon(2, i), this%obsHght(i), &
      !     this%obsTime(i)
    END DO


  END SUBROUTINE set_to_obsbase

  !!!!!!!!!!!!!!!!!

  SUBROUTINE gnss_qc(this)
    CLASS(ObsGNSS_t) :: this
  END SUBROUTINE gnss_qc

  SUBROUTINE ObsPrepareForSg(this, X)
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: X
  END SUBROUTINE ObsPrepareForSg

  !> @brief QC for refractivity
  !! 1. use the product quality flag
  !! 2. discard the rising ro obs below 5km, for low quality reason
  !! 3. discard the obs above 30km
  !! 4. discard the obs between super refraction layer and super refraction layer + 250m
  !! 5. outlier: it may be done in the obsthinning procedure
  !! 6. gradient checks:
  !!    6.1 dN/dz > -0.05 N m^-1
  !!    6.2 dN/dz < 1e-6 N m^-1
  !!    6.3 d^2N/dz^2 < 1e-4 N^-2 m^-2
  !!    6.4 dN/dZ != 0
  !! 7. max min checks, refractivity should in range [0, 500]
  !! 8. discard the profile if 30% obs are discarded
  SUBROUTINE gnss_refrac_qc(this, alt_pfl, refrac_pfl, qc_flag_pfl, nobs_pfl, blacklist_flag, valid_mask)
    CLASS(ObsGNSS_t) :: this
    REAL(r_kind), INTENT(IN) :: alt_pfl(:), refrac_pfl(:)
    LOGICAL, ALLOCATABLE, INTENT(IN) :: qc_flag_pfl(:)
    INTEGER(i_kind), INTENT(IN) :: nobs_pfl
    INTEGER(i_kind), INTENT(IN) :: blacklist_flag !> -1 for lt 5000, -2 for 10000

    LOGICAL,  ALLOCATABLE, INTENT(OUT):: valid_mask(:)

    ! locals
    INTEGER(i_kind) :: k
    REAL(r_kind) :: dNdz1, dNdz2, dNdz3, d2Ndz2, dNdz
    INTEGER(i_kind) :: nobs_valid, hgt_drop_cnt

    ! allocate the mask array and set to ture for initilization
    ALLOCATE(valid_mask(nobs_pfl))
    valid_mask(:) = .TRUE.
    valid_mask(1:2) = .FALSE.
    valid_mask(nobs_pfl-1:nobs_pfl) = .FALSE.

    !! [for cond 1] if there is qc flag for observation
    IF (ALLOCATED(qc_flag_pfl)) THEN
      valid_mask(:) = qc_flag_pfl
    END IF 

    hgt_drop_cnt = 0
    !! [for cond 2] 
    IF (blacklist_flag .EQ. -1) THEN
      DO k = 1, nobs_pfl
        IF (alt_pfl(k) .LT. 5000) THEN
          valid_mask(k) = .FALSE.
          hgt_drop_cnt = hgt_drop_cnt+1
        END IF 
      END DO
    END IF

    IF (blacklist_flag .EQ. -2) THEN
      DO k = 1, nobs_pfl
        IF (alt_pfl(k) .LT. 10000) THEN
          valid_mask(k) = .FALSE.
          hgt_drop_cnt = hgt_drop_cnt+1
        END IF 
      END DO
    END IF

    !! [for cond 3]
    DO k = 1, nobs_pfl
      IF (alt_pfl(k) .GT. 30000 ) THEN
        valid_mask(k) = .FALSE.
        hgt_drop_cnt = hgt_drop_cnt+1
      END IF 
    END DO

    !! [for cond 4]
    !! TODO: DON'T KNOW HOW TO FIND THE SUPER REFRACTION LAYER

    !! [for cond 5]
    !! in thinning
    
    !! [for cond 6]
    DO k = 2, nobs_pfl - 1
      ! use 2 point central gradient
      IF ((alt_pfl(k + 1) - alt_pfl(k - 1)) .NE. 0) THEN
        dNdz = (refrac_pfl(k + 1) - refrac_pfl(k - 1)) / (alt_pfl(k+1) - alt_pfl(k-1))
      ELSE
        ! PRINT*, 'HYJ+++: DNDZ ZERO'
        valid_mask(k) = .FALSE.
      ENDIF
      ! use 5 point gradient
      ! dNdz = (refrac_pfl(k - 2) - 8 * refrac_pfl(k - 1) +  8 * refrac_pfl(k + 1) - refrac_pfl(k + 2))/ (3 * (alt_pfl(k+2)-alt_pfl(k-2)))

      ! PRINT*, 'HYJ+++ dNdz', dNdz

      IF ( valid_mask(k) .AND. & ! check only if the mask is true
         (  dNdz .LE. -0.05 .OR. &  ! cond 6.1
            dNdz .GE. -1E-6 .OR. &  ! cond 6.2 
            dNdz .EQ. 0 )) THEN                                         ! cond 6.4
            ! PRINT*, 'HYJ+++: dndz : ', dNdz
            valid_mask(k) = .FALSE.
      END IF 
    END DO

   
    DO k = 2, nobs_pfl - 1
      ! use 3 point central 2-order gradient
      IF ((alt_pfl(k + 1) - alt_pfl(k - 1)) .NE. 0) THEN
        d2Ndz2 = (refrac_pfl(k -1) - 2.0D0 * refrac_pfl(k) + refrac_pfl(k + 1)) / (((alt_pfl(k + 1) - alt_pfl(k - 1))**2)/4.0D0)
        IF ( valid_mask(k) .AND. (d2Ndz2 .GE. 1e-4)) THEN
          ! PRINT*, 'HYJ+++: d2Ndz2 : ', d2Ndz2
          valid_mask(k) = .FALSE.
        END IF 
      ELSE
        ! PRINT*, 'HYJ+++: d2Ndz2 :  ZERO'
        valid_mask(k) = .FALSE.
      END IF
      ! PRINT*, 'HYJ+++ d2Ndz2', d2Ndz2
      ! USE 5 POINT Second-Order gradient
      ! d2Ndz2 = ( -refrac_pfl(k + 2) + 16 * refrac_pfl(k + 1) - 30*refrac_pfl(k) + 16 * refrac_pfl(k - 1) - refrac_pfl(k - 2))/ (12 * (((alt_pfl(k+2)-alt_pfl(k-2)/4)**2)))
    

    END DO

    !! [for cond 7]
    DO k = 1, nobs_pfl
      IF (valid_mask(k) .AND. (refrac_pfl(k) .GT. 500.0D0 .OR. refrac_pfl(k) .LE. 0) ) THEN
        valid_mask(k) = .FALSE.
      END IF
    END DO

    nobs_valid = 0
    DO k = 1, nobs_pfl
      IF (valid_mask(k)) nobs_valid = nobs_valid + 1
    END DO

    PRINT*, 'valid obs:', nobs_valid,  'obs filtered by altitude condition: ',  hgt_drop_cnt, 'total obs of the pfl', nobs_pfl

    !! [for cond 7] set all valid_mask to false to discard the profile.
    ! 
    IF (1.0D0 - DBLE(nobs_valid)/ DBLE(nobs_pfl-hgt_drop_cnt) >= 0.3) THEN
      valid_mask(:) = .FALSE.
    END IF

  END SUBROUTINE gnss_refrac_qc


  FUNCTION gnssForward(this, X, O) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'GNSS: This forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION gnssForward

  FUNCTION gnssTangent(this, dX, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'GNSS: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION gnssTangent

  FUNCTION gnssAdjoint(this, dY, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(ObsGNSS_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'GNSS: This adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION gnssAdjoint
END MODULE ObsGNSS_m
