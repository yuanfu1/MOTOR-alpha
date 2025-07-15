!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2021-11-24, created by Jiongming Pang for Vwpw Obsercations.
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang.
!                     2022-01-25, adjusted the parameters of raduis fro SuperObs by Jiongming Pang.
!                     2022-05-27, modified vwpwIngest to deal with the new format by Jiongming Pang.
!
!   Created by Jiongming Pang for Vwpw obs (pang.j.m@hotmail.com), 2021/11/24, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/25, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/5/27, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a vwpw data structure.
MODULE ObsVwpw_m
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

  TYPE, EXTENDS(ObsBase_t) :: ObsVwpw_t

  CONTAINS
    PROCEDURE :: ObsInitial => vwpwInitial
    PROCEDURE :: ObsIngest => vwpwIngest
    PROCEDURE :: ObsQC => vwpwQC
  END TYPE ObsVwpw_t

CONTAINS

  SUBROUTINE vwpwInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsVwpw_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    INTEGER :: i, istatus

    ! read the obs vwpw files names
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'input_dir_Vwpw', obsFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Vwpw', obsFileName)
    this%numfiles = UBOUND(obsFileName, 1)
    PRINT *, "amount of vwpw obs files: ", this%numfiles

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numfiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of vwpw obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! Read in an interpolation option:
    !CALL namelist_read(TRIM(configFile), "interpolation", this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! Read in analysis variable names:
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsVarList_Surface', this%varNames)
    IF (istatus .NE. 0 .OR. UBOUND(this%varNames, 1) .EQ. 0) THEN
      ! To read model state varList intead:
      istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', this%varNames)
    END IF
    IF (UBOUND(this%varNames, 1) .EQ. 0) THEN
      WRITE (*, 16)
16    FORMAT('There is neither obsVarList_Surface no modelState variables, Check your yaml and rerun!')
      STOP
    END IF
    WRITE (*, 11) UBOUND(this%varNames, 1), this%varNames(1:UBOUND(this%varNames, 1)) (1:5)
11  FORMAT('VWPWIngest - Obs variable lists: NumVar: ', I2, 1X, 20A5, ' proc: ', I2)

    this%numVars = UBOUND(this%varNames, 1)
    !this%numVars = 2
    !ALLOCATE(this%varNames(this%numVars))
    ! read the var names of obs files in configFile
    !this%varNames(1) = "uwnd"
    !this%varNames(2) = "vwnd"

    this%obsType = "VWPW"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ! this%missing = -888888.000
    this%missing = 999999999.0D0
    this%invalid = 3.4D38

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius = 5.0D2        ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(3, :) = 5.0D2   ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.25D0 ! No influence under this threshold value

    ! Threshold values for QC:
    IF (this%getVarIdx('uwnd') .GT. 0) this%qcThreshold(this%getVarIdx('uwnd')) = 500.0D0
    IF (this%getVarIdx('vwnd') .GT. 0) this%qcThreshold(this%getVarIdx('vwnd')) = 500.0D0

    this%qcThreshold(1) = 6.0D0
    this%qcThreshold(2) = 6.0D0
  END SUBROUTINE vwpwInitial

  SUBROUTINE vwpwIngest(this, X)
    IMPLICIT NONE
    CLASS(ObsVwpw_t) :: this
    TYPE(State_t) :: X

    ! local variables
    INTEGER :: nf_op, nc_id, var_id, i, num_obs, num_dims, nf90_max_va_dims
    INTEGER :: ifile, num_obs_total, num_dataFile, ind_obsin_beg, ind_obsin_end
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_obsin

    ! local variables, for get the unique obs in multi-files
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), obsErrs(:, :), olatlon(:, :)
    REAL(r_kind), ALLOCATABLE :: obsHght(:), land(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:), num_sitein(:), num_obsin(:)
    CHARACTER(LEN=20), ALLOCATABLE :: stNames(:)
    CHARACTER(LEN=10), ALLOCATABLE :: timestr(:)
    CHARACTER(LEN=30), ALLOCATABLE :: Id_time(:), Id_time_unique(:)
    LOGICAL :: fileExist(this%numFiles)  ! Yuanfu Xie added a safeguard for missing file 2022-05-31
    LOGICAL, ALLOCATABLE :: index_unique_judge(:)
    INTEGER, ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER :: num_obs_total_unique, j
    CHARACTER(LEN=1024) :: fntxt_Idtime_a, fntxt_Idtime_u

    CHARACTER(LEN=5)  :: obsTypeh, siteId
    CHARACTER(LEN=20) :: chead1, chead2
    INTEGER :: ctime(5), sitenum, ios, ilevels, isite, iobs, tmp, istn, iyear, imonth, iday, ihh, imm
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, rpres, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
    REAL(r_kind) :: p, u, v, rh, q, vwu, hd, vd
    REAL(r_kind) :: missingValue = -888888.0D0

    ! Profiles have wind data only:
    IF (this%getVarIdx('uwnd') .EQ. 0 .OR. this%getVarIdx('vwnd') .EQ. 0) THEN
      this%numObs = 0
      RETURN
    END IF

    ! Get the total amount of observations in multi-files. -Jiongming Pang
    num_dataFile = this%numFiles
    ALLOCATE (num_sitein(num_dataFile))
    ALLOCATE (num_obsin(num_dataFile))
    num_obsin(:) = 0
    DO ifile = 1, num_dataFile
      ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      INQUIRE (FILE=TRIM(this%fileNames(ifile)), EXIST=fileExist(ifile))
      IF (fileExist(ifile)) THEN
        WRITE (*, 4) TRIM(this%fileNames(ifile))
4       FORMAT('ObsProfiler: profiler data file: ', A, ' is found and reading...')
      ELSE
        WRITE (*, 3) TRIM(this%fileNames(ifile))
3       FORMAT('ObsProfiler: This profiler data file does not exist but the run continues, CHECK! ', A)
        CYCLE
      END IF
      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a20,4x,5i4,2x,a7,i7)', iostat=ios) chead1, ctime(1:5), chead2, sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        num_sitein(ifile) = 0
        CLOSE (101)
        CYCLE
      END IF
      num_sitein(ifile) = sitenum
      READ (101, *)
      DO isite = 1, num_sitein(ifile)
        READ (101, '(i5,4x,a5,2f9.3,6i5,2f14.3)') &
          istn, siteId, rlat, rlon, iyear, imonth, iday, ihh, imm, ilevels, relv, rpres
        num_obsin(ifile) = num_obsin(ifile) + ilevels
        DO j = 1, ilevels
          READ (101, *)
        END DO
      END DO
      CLOSE (101)
      WRITE (*, 1) ifile, num_obsin(ifile)
1     FORMAT('Number of profiler Obs in file ', I2 ' is: ', I4)
    END DO
    num_obs_total = SUM(num_obsin)
    WRITE (*, 2) num_obs_total
2   FORMAT('Total amount of profiler Obs in all files: ', I6)

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

    ! obs file loop
    iobs = 0
    DO ifile = 1, num_dataFile
      IF (.NOT. fileExist(ifile)) CYCLE ! Yuanfu Xie added a safeguard for missing file 2022-05-31

      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a20,4x,5i4,2x,a7,i7)', iostat=ios) chead1, ctime(1:5), chead2, sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        CLOSE (101)
        CYCLE
      END IF
      READ (101, *)
      DO isite = 1, sitenum
        ! READ(101, '(i5,4x,a5,2f9.3,3i5,2f14.3)') &
        !     istn, siteId, rlat, rlon, ihh, imm, ilevels, relv, rpres
        READ (101, '(i5,4x,a5,2f9.3,6i5,2f14.3)') &
          istn, siteId, rlat, rlon, iyear, imonth, iday, ihh, imm, ilevels, relv, rpres
        DO j = 1, ilevels
          iobs = iobs + 1
          READ (101, '(11f14.3)') p, u, v, vwu, hgh, t, rh, td, q, hd, vd
          obsData(iobs, this%getVarIdx('uwnd')) = u
          obsData(iobs, this%getVarIdx('vwnd')) = v
          obsHght(iobs) = hgh

          CALL Time_GMT_to_Unix((/iyear, imonth, iday, ihh, imm, 0/), obsTime(iobs))

          olatlon(1, iobs) = rlat   ! lat
          olatlon(2, iobs) = rlon   ! lon

          land(iobs) = 1.0D0
          obsErrs(iobs, :) = this%missing * (-1.0D0)   !-999999999.0D0
          stNames(iobs) = siteId
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file ', ifile, ':', num_obsin(ifile)
    END DO

    WHERE (ABS(obsData - missingValue) .LT. 0.01D0) obsData = this%missing
    IF (this%getVarIdx('uwnd') .GT. 0) &
      WHERE (ABS(obsData(:, this%getVarIdx('uwnd')) - this%missing) .LT. 0.01D0) &
      obsData(:, this%getVarIdx('vwnd')) = this%missing
    IF (this%getVarIdx('vwnd') .GT. 0) &
      WHERE (ABS(obsData(:, this%getVarIdx('vwnd')) - this%missing) .LT. 0.01D0) &
      obsData(:, this%getVarIdx('uwnd')) = this%missing
    WHERE (ABS(obsHght - missingValue) .LT. 0.01D0) obsHght = this%missing

    ! get the index of the used elements
    ALLOCATE (index_unique_tmp(num_obs_total))
    DO i = 1, num_obs_total
      index_unique_tmp(i) = i
    END DO
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
    PRINT *, "varNames: ", this%varNames

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

    PRINT *, "DONE vwpwIngest"

  END SUBROUTINE vwpwIngest

  SUBROUTINE vwpwQC(this)
    CLASS(ObsVwpw_t) :: this
  END SUBROUTINE vwpwQC

END MODULE ObsVwpw_m
