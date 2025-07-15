!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2021-11-24, created by Jiongming Pang for Sound Obsercations.
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang.
!                     2022-01-25, adjusted the parameters of raduis fro SuperObs by Jiongming Pang.
!
!   Created by Jiongming Pang for Sound obs (pang.j.m@hotmail.com), 2021/11/24, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/25, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/09/16, @GBA-MWF, Shenzhen
!     For modified soundIngest by removing hard coded obs variables and read them in from yaml
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a sound data structure.
MODULE ObsSound_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m
  USE slint, ONLY: slint_init, src_grid, tgt_grid
  USE mpObs_m, ONLY: mpObs_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE ncReadVar_m
  USE WriteVar_m
  USE AdvanceTime_m
  USE conversions_m
  USE parameters_m, ONLY: degree2radian

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsSound_t

  CONTAINS
    PROCEDURE :: ObsInitial => soundInitial
    PROCEDURE :: ObsIngest => soundIngest
    PROCEDURE :: ObsQC => soundQC
  END TYPE ObsSound_t

CONTAINS

  SUBROUTINE soundInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsSound_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    INTEGER :: i, istatus

    ! read the obs sound files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Sound', obsFileDir)

    !CALL namelist_read(TRIM(configFile), "obsFileList_Sound", obsFileName)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Sound', obsFileName)
    this%numfiles = UBOUND(obsFileName, 1)
    PRINT *, "amount of sound obs files: ", this%numfiles

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numfiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of sound obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! Read in an interpolation option:
    !CALL namelist_read(TRIM(configFile), "interpolation",this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! read the var names of obs files in configFile
    this%numVars = 5   !1!5
    ALLOCATE (this%varNames(this%numVars))

    this%varNames(4) = "uwnd"
    this%varNames(5) = "vwnd"
    this%varNames(1) = "temp"
    this%varNames(2) = "temp_dew"
    this%varNames(3) = "pres"
!  /////////
!     ! Read in analysis variable names:
!     istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsVarList_Surface', this%varNames)
!     IF (istatus .NE. 0 .OR. UBOUND(this%varNames,1) .EQ. 0) THEN
!       ! To read model state varList intead:
!       istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', this%varNames)
!     ENDIF
!     IF (UBOUND(this%varNames,1) .EQ. 0) THEN
!       WRITE(*,16)
! 16    FORMAT('There is neither obsVarList_Surface no modelState variables, Check your yaml and rerun!')
!       STOP
!     END IF
!     WRITE(*,11) UBOUND(this%varNames,1),this%varNames(1:UBOUND(this%varNames,1))(1:5)
! 11  FORMAT('SoundInitial - Obs variable lists: NumVar: ',I2,1X,20A5,' proc: ',I2)

!     ! Set the num of variables:
!     this%numVars = UBOUND(this%varNames,1)
!  /////////

    this%obsType = "SOUND"
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
    this%qcThreshold = 500.0D0 ! Temporary holders

  END SUBROUTINE soundInitial

  SUBROUTINE soundIngest(this, X)
    IMPLICIT NONE

    CLASS(ObsSound_t) :: this
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
    LOGICAL :: fileExist(this%numFiles) ! Yuanfu Xie added a safeguard for missing file 2022-05-31
    INTEGER, ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER :: num_obs_total_unique, j
    CHARACTER(LEN=1024) :: fntxt_Idtime_a, fntxt_Idtime_u

    CHARACTER(LEN=5)  :: obsTypeh, siteId
    CHARACTER(LEN=20) :: nhead
    INTEGER :: ctime(4), sitenum, ios, ilevels, imm, isite, iobs, tmp, idew
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
    REAL(r_kind) :: missingValue = 999999.0D0

    ! Get the total amount of observations in multi-files. -Jiongming Pang
    num_dataFile = this%numFiles
    ALLOCATE (num_sitein(num_dataFile))
    ALLOCATE (num_obsin(num_dataFile))
    num_obsin(:) = 0
    DO ifile = 1, num_dataFile
      ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      INQUIRE (FILE=TRIM(this%fileNames(ifile)), EXIST=fileExist(ifile))
      IF (fileExist(ifile)) THEN
        WRITE (*, 1) TRIM(this%fileNames(ifile))
1       FORMAT('ObsSound: sounding data file: ', A, ' is found and reading...')
      ELSE
        WRITE (*, 3) TRIM(this%fileNames(ifile))
3       FORMAT('ObsSound: This sound data file does not exist but the run continues, CHECK! ', A)
        CYCLE
      END IF
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
        READ (101, '(3x,a5,5f10.2,i6)', iostat=ios) siteId, rlat, rlon, relv, rday, rhm, ilevels
        IF (ios /= 0) EXIT
        num_obsin(ifile) = num_obsin(ifile) + ilevels
        DO j = 1, ilevels
          READ (101, *)
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file ', ifile, ':', num_obsin(ifile)
    END DO
    num_obs_total = SUM(num_obsin)
    WRITE (*, 2) num_obs_total, X%mpddSub%myrank
2   FORMAT('Total amount of sounding Obs in all files: ', I8, ' pc', I2)

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
      IF (.NOT. fileExist(ifile)) CYCLE  ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a5,5i8)', iostat=ios) obsTypeh, ctime(1:4), sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        CLOSE (101)
        CYCLE
      END IF
      DO isite = 1, sitenum
        READ (101, '(3x,a5,5f10.2,i6)', iostat=ios) siteId, rlat, rlon, relv, rday, rhm, ilevels
        IF (ios /= 0) EXIT
        DO j = 1, ilevels
          iobs = iobs + 1
          ! READ(101, '(10f10.2)') p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
          READ (101, *) p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d

          ! Wind conversion:
          IF (this%getVarIdx('uwnd') .GT. 0 .AND. this%getVarIdx('vwnd') .GT. 0) THEN
            IF (ABS(dd - missingValue) .LT. 0.01D0 .OR. ABS(ff - missingValue) .LT. 0.01D0) THEN
              obsData(iobs, this%getVarIdx('uwnd')) = this%missing   ! uwnd
              obsData(iobs, this%getVarIdx('vwnd')) = this%missing   ! vwnd
            ELSE
              CALL wind_to_uv(ff, dd, obsData(iobs, this%getVarIdx('uwnd')), &
                              obsData(iobs, this%getVarIdx('vwnd')))
            END IF
          END IF

          ! Temperature check:
          IF (this%getVarIdx('temp') .GT. 0) THEN
            IF (ABS(t - missingValue) .GT. 0.01D0 .AND. t .LT. 200.0D0) &
              obsData(iobs, this%getVarIdx('temp')) = t + 273.15
            IF (ABS(t - missingValue) .LT. 0.01D0 .OR. t .GE. 200.0D0) &
              obsData(iobs, this%getVarIdx('temp')) = this%missing
          END IF

          ! Dew point:
          idew = 0
          IF (this%getVarIdx('dewp') .GT. 0) THEN
            idew = this%getVarIdx('dewp')
          ELSE IF (this%getVarIdx('qvapor') .GT. 0) THEN
            idew = this%getVarIdx('qvapor')
          END IF
          IF (idew .GT. 0) THEN
            IF (ABS(td - missingValue) .GT. 0.01D0 .AND. td .LT. 200.0D0) &
              obsData(iobs, idew) = td + 273.15
            IF (ABS(td - missingValue) .LT. 0.01D0 .OR. td .GE. 200.0D0) &
              obsData(iobs, idew) = this%missing
          END IF

          ! pressure:
          IF (this%getVarIdx('pres') .GT. 0) THEN
            IF (ABS(p - missingValue) .GT. 0.01D0 .AND. p .LT. 1500.0D0 * 100.0D0) obsData(iobs, this%getVarIdx('pres')) = p * 100.0D0
            IF (ABS(p - missingValue) .LT. 0.01D0 .OR. p .GE. 1500.0D0 * 100.0D0) obsData(iobs, this%getVarIdx('pres')) = this%missing
          END IF

          obspres(iobs) = obsData(iobs, this%getVarIdx('pres'))

          IF (ABS(hgh - missingValue) .LT. 0.01D0) THEN
            obsHght(iobs) = this%missing
          ELSE
            obsHght(iobs) = hgh
          END IF

          imm = INT(rhm - ctime(4) * 100)
          CALL Time_GMT_to_Unix((/ctime(1:4), imm, 0/), tmp)
          ! IF (ABS(ttd - missingValue) .LT. 0.01D0) THEN
          obsTime(iobs) = tmp
          ! ELSE
          ! obsTime(iobs) = tmp + ttd
          ! END IF

          ! IF (ABS(lat_d - missingValue) .LT. 0.01D0 .OR. ABS(lon_d - missingValue) .LT. 0.01D0) THEN
          olatlon(1, iobs) = rlat   ! lat
          olatlon(2, iobs) = rlon   ! lon
          ! ELSE
          !   olatlon(1, iobs) = rlat + lat_d   ! lat
          !   olatlon(2, iobs) = rlon + lon_d   ! lon
          ! END IF

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
    ! WHERE (ABS(obsHght - this%missing) .LT. 0.01D0) index_unique_tmp = 0
    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    PRINT *, "the amount of all sounding input obs and the unique: ", num_obs_total, num_obs_total_unique

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

    ! Modified by Zilong to convert the pressure observations to density
    BLOCK
      USE parameters_m, ONLY: dry_air_gas_const
      USE conversions_m, ONLY: Td_to_qvapor

      IF (this%numObs .GT. 0 .AND. X%getVarIdx('qvapor') /= 0) THEN
        ! this%varNames(2) = "qvapor" ! Remove this hardcoded varNames setting
        DO j = 1, this%numObs
          IF (this%obsData(j, idew) == this%missing .OR. &
              this%obsData(j, this%getVarIdx('pres')) == this%missing) THEN
            this%obsData(j, idew) = this%missing
          ELSE
            this%obsData(j, idew) = Td_to_qvapor(this%obsData(j, idew), &
                                                 this%obsData(j, this%getVarIdx('pres')))
          END IF
        END DO
      END IF

      IF (this%numObs .GT. 0 .AND. X%getVarIdx('rho') /= 0 .OR. X%getVarIdx('rho_ctl') /= 0) THEN
        ! this%varNames(3) = "rho"
        DO j = 1, this%numObs
          IF (this%obsData(j, this%getVarIdx('pres')) == this%missing .OR. &
              this%obsData(j, this%getVarIdx('temp')) == this%missing) THEN
            this%obsData(j, this%getVarIdx('pres')) = this%missing
          ELSE
            this%obsData(j, this%getVarIdx('pres')) = &
              this%obsData(j, this%getVarIdx('pres')) / &
              dry_air_gas_const / this%obsData(j, this%getVarIdx('temp'))
          END IF
        END DO
      END IF

    END BLOCK

    this%presIdx = 0
    this%tempIdx = 0
    this%qvaporIdx = 0

    DO i = 1, this%numVars
      IF (this%varNames(i) == 'pres') this%presIdx = i
      IF (this%varNames(i) == 'temp') this%tempIdx = i
      IF (this%varNames(i) == 'qvapor') this%qvaporIdx = i
    END DO

    ! Added by Zilong for temporary scaling
    IF (this%tempIdx /= 0) this%qcThreshold(this%tempIdx) = 20.0D0
    IF (this%presIdx /= 0) this%qcThreshold(this%presIdx) = 1000.0D0
    ! IF (this%qvaporIdx /= 0) this%qcThreshold(this%qvaporIdx) = 0.006D0
    ! this%qcThreshold(4) = 10.0D0
    ! this%qcThreshold(5) = 10.0D0

    ! this%qcThreshold(4) = 10.0D0
    ! this%qcThreshold(5) = 10.0D0
    PRINT *, "DONE soundIngest ", this%varNames(1:this%numVars)

  END SUBROUTINE soundIngest

  SUBROUTINE soundQC(this)
    CLASS(ObsSound_t) :: this
  END SUBROUTINE soundQC

END MODULE ObsSound_m
