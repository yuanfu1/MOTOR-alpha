!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-07-15, created by Jiongming Pang for GWST Obsercations.
!
!   Created by Jiongming Pang for GWST obs (pang.j.m@hotmail.com), 2021/07/15, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, @GBA-MWF, Shenzhen for adding obs
!     forward, tangent and adjoint operators.
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a GWST data structure.
MODULE ObsGWST_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m
  USE slint, ONLY: slint_init, src_grid, tgt_grid
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE ncReadVar_m
  USE WriteVar_m
  USE AdvanceTime_m
  USE conversions_m
  USE parameters_m, ONLY: degree2radian, missing

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsGWST_t

  CONTAINS
    PROCEDURE :: ObsInitial => GWSTInitial
    PROCEDURE :: ObsIngest => GWSTIngest
    PROCEDURE :: ObsForward => GWSTForward
    PROCEDURE :: ObsTangent => GWSTTangent
    PROCEDURE :: ObsAdjoint => GWSTAdjoint
    PROCEDURE :: ObsQC => GWSTQC
  END TYPE ObsGWST_t

CONTAINS

  SUBROUTINE GWSTInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsGWST_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    INTEGER :: i, istatus

    ! read the obs GWST files names
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'input_dir_GWST', obsFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_GWST', obsFileName)
    this%numfiles = UBOUND(obsFileName, 1)
    PRINT *, "amount of GWST obs files: ", this%numfiles

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numfiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of GWST obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! Read in an interpolation option:
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ! read the var names of obs files in configFile
    this%numVars = 8
    ALLOCATE (this%varNames(this%numVars))
    this%varNames(1) = "temp"
    this%varNames(2) = "temp_dew"
    this%varNames(3) = "pres"
    this%varNames(4) = "uwnd"
    this%varNames(5) = "vwnd"
    this%varNames(6) = "pcpa01hr"
    this%varNames(7) = "pcpa06hr"
    this%varNames(8) = "pcpa24hr"

    this%obsType = "GWST"
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

    this%correlation_threshold = 0.25D0 ! No influence under this threshold value

    ! Threshold values for QC:
    this%qcThreshold(:) = 500.0D0
    this%qcThreshold(3) = 1500.0D0 * 100.0D0  ! for psfc (stationPressure)
    this%qcThreshold(6) = 300.0D0
    this%qcThreshold(7) = 600.0D0
    this%qcThreshold(8) = 1000.0D0

  END SUBROUTINE GWSTInitial

  SUBROUTINE GWSTIngest(this, X)
    IMPLICIT NONE
    CLASS(ObsGWST_t) :: this
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

    CHARACTER(LEN=5)  :: obsTypeh
    CHARACTER(LEN=7)  :: siteID
    CHARACTER(LEN=10) :: siteIDcall
    CHARACTER(LEN=20) :: chead1, chead2
    INTEGER :: ctime(4), sitenum, ios, ilevels, isite, iobs, tmp, istn, iyear, imonth, iday, ihh, imm
    REAL(r_kind) :: rlat, rlon, relv, rday, rhhmm, rpres, rid, hgh, t, td, winddir, windspd, ttd, lat_d, lon_d
    REAL(r_kind) :: p, u, v, rh, q, vwu, hd, vd, pcpa01hr, pcpa06hr, pcpa24hr, dtmp(10)
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
        WRITE (*, 4) TRIM(this%fileNames(ifile))
4       FORMAT('ObsGWST: GWST data file: ', A, ' is found and reading...')
      ELSE
        WRITE (*, 3) TRIM(this%fileNames(ifile))
3       FORMAT('ObsGWST: This GWST data file does not exist but the run continues, CHECK! ', A)
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
      CLOSE (101)
      WRITE (*, 1) ifile, num_sitein(ifile)
1     FORMAT('Number of GWST Obs in file ', I2 ' is: ', I4)
    END DO
    num_obs_total = SUM(num_sitein)
    WRITE (*, 2) num_obs_total
2   FORMAT('Total amount of GWST Obs in all files: ', I6)

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

    obsData = missing

    ! obs file loop
    iobs = 0
    DO ifile = 1, num_dataFile
      IF (.NOT. fileExist(ifile)) CYCLE ! Yuanfu Xie added a safeguard for missing file 2022-05-31

      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a5,5i8)', iostat=ios) obsTypeh, ctime(1:4), sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        CLOSE (101)
        CYCLE
      END IF

      DO isite = 1, sitenum
        READ (101, '(a7, 1x, a10, 5f10.2)') &
          siteID, siteIDcall, rlat, rlon, relv, rday, rhhmm
        READ (101, '(8f10.2)') &
          p, t, td, winddir, windspd, pcpa01hr, pcpa06hr, pcpa24hr

        iobs = iobs + 1

        IF (ABS(t - missingValue) .GT. 0.01D0 .AND. t .LT. this%qcThreshold(1)) &
          obsData(iobs, 1) = t + 273.15D0
        IF (ABS(td - missingValue) .GT. 0.01D0 .AND. td .LT. this%qcThreshold(2)) &
          obsData(iobs, 2) = td + 273.15D0
        IF (ABS(p - missingValue) .GT. 0.01D0 .AND. p .LT. this%qcThreshold(3)) &
          obsData(iobs, 3) = p * 100.0D0
        IF (ABS(pcpa01hr - missingValue) .GT. 0.01D0 .AND. pcpa01hr .LT. this%qcThreshold(6)) &
          obsData(iobs, 6) = pcpa01hr
        IF (ABS(pcpa06hr - missingValue) .GT. 0.01D0 .AND. pcpa06hr .LT. this%qcThreshold(7)) &
          obsData(iobs, 7) = pcpa06hr
        IF (ABS(pcpa24hr - missingValue) .GT. 0.01D0 .AND. pcpa24hr .LT. this%qcThreshold(8)) &
          obsData(iobs, 8) = pcpa24hr

        IF (ABS(windspd - missingValue) .LT. 0.01D0 .OR. &
            ABS(winddir - missingValue) .LT. 0.01D0 .OR. &
            windspd .LT. 0.0D0 .OR. winddir .LT. 0.0D0 .OR. &
            windspd .GT. 100.0D0 .OR. winddir .GT. 360.0D0) THEN

          obsData(iobs, 4) = missing   ! uwnd
          obsData(iobs, 5) = missing   ! vwnd

        ELSE

          CALL wind_to_uv(windspd, winddir, obsData(iobs, 4), obsData(iobs, 5))

        END IF

        IF (ABS(relv - missingValue) .GT. 0.01D0 .AND. relv .LT. 50.0D0) THEN
          obsHght(iobs) = relv
        ELSE
          obsHght(iobs) = 10.0D0
        END IF

        iyear = ctime(1)
        imonth = ctime(2)
        iday = rday * 1
        ihh = rhhmm / 100
        imm = rhhmm - ihh * 100
        CALL Time_GMT_to_Unix((/iyear, imonth, iday, ihh, imm, 0/), obsTime(iobs))

        olatlon(1, iobs) = rlat   ! lat
        olatlon(2, iobs) = rlon   ! lon

        land(iobs) = 1.0D0
        obsErrs(iobs, :) = missing * (-1.0D0)   !-999999999.0D0
        stNames(iobs) = siteID
      END DO

      CLOSE (101)
      PRINT *, 'Number of Obs in file ', ifile, ':', num_obsin(ifile)

    END DO

    !WHERE ( ABS(obsData-missingValue) .LT. 0.01D0 ) obsData = missing
    !WHERE ( ABS(obsData(:,1)-missing) .LT. 0.01D0 ) obsData(:,2) = missing
    !WHERE ( ABS(obsData(:,2)-missing) .LT. 0.01D0 ) obsData(:,1) = missing
    !WHERE ( ABS(obsHght-missingValue) .LT. 0.01D0 ) obsHght = missing

    ! get the index of the used elements
    ALLOCATE (index_unique_tmp(num_obs_total))
    DO i = 1, num_obs_total
      index_unique_tmp(i) = i
    END DO
    WHERE (ABS(obsHght - missing) .LT. 0.01D0) index_unique_tmp = 0
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
        MAXVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. missing)), &
        MINVAL(this%obsData(:, i), MASK=(this%obsData(:, i) .NE. missing))
    END DO

    !PRINT *, 'obsTime:'
    !PRINT *, this%obsTime

    !PRINT *, "varNames: ", this%varNames
    !DO i=1, this%numVars
    !  PRINT *, 'temp, temp_dew, pres, uwnd, vwnd, pcpa01hr, pcpa06hr, pcpa24hr:'
    !  PRINT *, this%obsData(i,:)
    !END DO

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

    PRINT *, "DONE GWSTIngest"

    this%qcThreshold(1) = 8.0D0
    this%qcThreshold(3) = 1000.0D0
    this%qcThreshold(4) = 8.0D0
    this%qcThreshold(5) = 8.0D0

  END SUBROUTINE GWSTIngest

  FUNCTION GWSTForward(this, X, O) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsGWST_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'GWST: This forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION GWSTForward

  FUNCTION GWSTTangent(this, dX, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsGWST_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'GWST: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION GWSTTangent

  FUNCTION GWSTAdjoint(this, dY, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(ObsGWST_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'GWST: This adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION GWSTAdjoint

  SUBROUTINE GWSTQC(this)
    CLASS(ObsGWST_t) :: this
  END SUBROUTINE GWSTQC

END MODULE ObsGWST_m
