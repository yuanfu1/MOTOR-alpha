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
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/1/14, @GBA-MWF, Shenzhen
!     For adding a routine to fill in missing obs height information.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, @GBA-MWF, Shenzhen for adding obs
!     forward, tangent and adjoint operators.
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
  USE parameters_m, ONLY: degree2radian, missing

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsSound_t

    LOGICAL :: soundFill  ! Hydrostatic fill missing heights as well wind and temp and dew

  CONTAINS
    PROCEDURE :: ObsInitial => soundInitial
    PROCEDURE :: ObsIngest => soundIngest
    PROCEDURE :: ObsForward => soundForward
    PROCEDURE :: ObsTangent => soundTangent
    PROCEDURE :: ObsAdjoint => soundAdjoint
    PROCEDURE :: ObsQC => soundQC

    PROCEDURE, PUBLIC :: fillinMissing    ! Yuanfu Xie added on 2023/01/14
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
    PRINT *, 'Sound directory read status: ', istatus
    istatus = yaml_get_var(configFile, 'IO', 'SoundFill', this%soundFill)
    PRINT *, 'Sound Filling read status: ', istatus, this%soundFill
    IF (istatus .NE. 0) this%soundFill = .FALSE.  ! If users do not specify fill option, it does not fill

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

    this%obsType = "SOUND"
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
    this%qcThreshold(5) = 500.0D0
    this%qcThreshold(1) = 500.0D0
    this%qcThreshold(2) = 500.0D0
    this%qcThreshold(3) = 1500.0D0 * 100.0D0

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
    INTEGER :: ctime(4), sitenum, ios, ilevels, imm, isite, iobs, tmp
    REAL(r_kind) :: rlat, rlon, relv, rday, rhm, p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
    REAL(r_kind) :: missingValue = 999999.0D0

    ! Debugging:
!#define TRACE_OBS
    INTEGER(i_kind) :: jj

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

        ! ! Yuanfu Xie added a check to skip those sites not inside analysis domain: 2023-01-11
        ! IF (rlat*degree2radian .GE. MINVAL(X%sg%cell_cntr(1,:)) .AND. &
        !     rlat*degree2radian .LE. MAXVAL(X%sg%cell_cntr(1,:)) .AND. &
        !     rlon*degree2radian .GE. MINVAL(X%sg%cell_cntr(2,:)) .AND. &
        !     rlon*degree2radian .LE. MAXVAL(X%sg%cell_cntr(2,:)) ) THEN
        num_obsin(ifile) = num_obsin(ifile) + ilevels
        ! END IF
        DO j = 1, ilevels
          READ (101, *)
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file in preprocess ', ifile, ':', num_obsin(ifile), X%mpddGlob%myrank
    END DO
    num_obs_total = SUM(num_obsin)
    WRITE (*, 2) num_obs_total, X%sg%mpddInfo_sg%myrank
2   FORMAT('Total amount of sounding Obs in all files: ', I8, ' pc:', I4)

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
    obspres = missing
    DO ifile = 1, num_dataFile
      IF (.NOT. fileExist(ifile)) CYCLE  ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      OPEN (101, FILE=TRIM(this%fileNames(ifile)), STATUS='old', ACTION='READ')
      READ (101, '(a5,5i8)', iostat=ios) obsTypeh, ctime(1:4), sitenum
      IF (ios /= 0) THEN
        PRINT *, "ERROR in reading head of file:", TRIM(this%fileNames(ifile))
        CLOSE (101)
        CYCLE
      END IF

      ! Reading data from all sites:
      DO isite = 1, sitenum
        READ (101, '(3x,a5,5f10.2,i6)', iostat=ios) siteId, rlat, rlon, relv, rday, rhm, ilevels
        IF (ios /= 0) EXIT

        ! ! Yuanfu Xie added a check to skip those sites not inside analysis domain: 2023-01-11
        ! IF (rlat*degree2radian .GE. MINVAL(X%sg%cell_cntr(1,:)) .AND. &
        !     rlat*degree2radian .LE. MAXVAL(X%sg%cell_cntr(1,:)) .AND. &
        !     rlon*degree2radian .GE. MINVAL(X%sg%cell_cntr(2,:)) .AND. &
        !     rlon*degree2radian .LE. MAXVAL(X%sg%cell_cntr(2,:)) ) THEN

        DO j = 1, ilevels
          iobs = iobs + 1
          ! READ(101, '(10f10.2)') p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d
          READ (101, *) p, rid, hgh, t, td, dd, ff, ttd, lat_d, lon_d

          IF (ABS(dd - missingValue) .LT. 0.01D0 .OR. ABS(ff - missingValue) .LT. 0.01D0) THEN
            obsData(iobs, 4) = missing   ! uwnd
            obsData(iobs, 5) = missing   ! vwnd
          ELSE
            CALL wind_to_uv(ff, dd, obsData(iobs, 4), obsData(iobs, 5))
          END IF
          IF (ABS(t - missingValue) .GT. 0.01D0 .AND. t .LT. 200.0D0) obsData(iobs, 1) = t + 273.15
          IF (ABS(t - missingValue) .LT. 0.01D0 .OR. t .GE. 200.0D0) obsData(iobs, 1) = missing

          IF (ABS(td - missingValue) .GT. 0.01D0 .AND. td .LT. 200.0D0) obsData(iobs, 2) = td + 273.15
          IF (ABS(td - missingValue) .LT. 0.01D0 .OR. td .GE. 200.0D0) obsData(iobs, 2) = missing

          IF (ABS(p - missingValue) .GT. 0.01D0 .AND. p .LT. 1500.0D0 * 100.0D0) obsData(iobs, 3) = p * 100.0D0
          IF (ABS(p - missingValue) .LT. 0.01D0 .OR. p .GE. 1500.0D0 * 100.0D0) obsData(iobs, 3) = missing

          obspres(iobs) = obsData(iobs, 3)

          IF (hgh .LT. missingValue) THEN
            obsHght(iobs) = hgh
          ELSE
            obsHght(iobs) = missing
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
          obsErrs(iobs, :) = missing * (-1.0D0)   !-999999999.0D0
          stNames(iobs) = siteId
        END DO

        ! Filling in the missing values at this site:
        IF (this%soundFill) THEN
#ifdef TRACE_OBS
          ! Debugging: using the station name found in the ObsBase.F90 from obsThinning
          IF (X%mpddGlob%myrank .EQ. 2 .AND. TRIM(siteId) .EQ. "98328") THEN
            DO jj = 1, 50 ! ilevels
              WRITE (*, 3333) jj, iobs, ilevels, iobs - ilevels + 1, olatlon(:, iobs - ilevels + jj), &
                TRIM(siteId), &
                obsData(iobs - ilevels + jj, 1:2), obsData(iobs - ilevels + jj, 4:5), &
                obsData(iobs - ilevels + jj, 3), obsHght(iobs - ilevels + jj), obsTime(iobs - ilevels + jj) - X%sg%tt(1)
            END DO
          END IF
3333      FORMAT('Fillin before J: ', I4, ' iobs: ', I6, I4, I6, ' oll: ', 2D12.5, &
                 ' Stn: ', A, ' obs: ', 4D12.4, ' p-H: ', 2D12.4, ' oTm: ', D14.6)
#endif

          CALL this%fillinMissing(X, TRIM(siteId), X%mpddGlob%myrank, &
                                  missing, ilevels, obsData(iobs - ilevels + 1:iobs, 1:5), num_obs, &
                                  obsHght(iobs - ilevels + 1:iobs), obsTime(iobs - ilevels + 1:iobs), &
                                  olatlon(1:2, iobs - ilevels + 1:iobs), obsErrs(iobs - ilevels + 1:iobs, 1:5), &
                                  land(iobs - ilevels + 1:iobs))

          ! Remove those unused data:
          obspres(iobs - ilevels + 1:iobs) = missing
          iobs = iobs - ilevels + num_obs
          obspres(iobs - num_obs + 1:iobs) = obsData(iobs - num_obs + 1:iobs, 3)
          stNames(iobs - num_obs + 1:iobs) = siteId

#ifdef TRACE_OBS
          ! Debugging: using the station name found in the ObsBase.F90 from obsThinning
          IF (X%mpddGlob%myrank .EQ. 2 .AND. TRIM(siteId) .EQ. "98328") THEN
            DO jj = 1, 10 !num_obs
              WRITE (*, 333) jj, iobs, ilevels, iobs - num_obs + jj, num_obs, &
                olatlon(:, iobs - num_obs + jj), TRIM(siteId), &
                obsData(iobs - num_obs + jj, 1:2), obsData(iobs - num_obs + jj, 4:5), &
                obsData(iobs - num_obs + jj, 3), obsHght(iobs - num_obs + jj), obsTime(iobs - ilevels + jj) - X%sg%tt(1)
333           FORMAT('Fillin after J: ', I4, ' iobs: ', I6, I4, I6, ' nobs:', I6, ' oll: ', 2D12.5, &
                     ' Stn: ', A, ' obs: ', 4D12.4, ' p-H: ', 2D12.4, ' Otm: ', D14.6)
            END DO
          END IF
#endif

        END IF ! If filled requested

!         ELSE  ! for sites not inside domain:

! #ifdef DEBUG
!           ! Debugging: Yuanfu Xie: 2023-01-11:
!           WRITE(*,228) isite,X%mpddGlob%myrank
! 228       FORMAT('Sounding ingest skip site: ',I4,' pc',I3)
! #endif

!           DO j = 1, ilevels
!             READ (101, *)
!           END DO
!         END IF ! For sites inside the domain
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in file ingesting: ', ifile, ':', num_obsin(ifile), X%mpddGlob%myrank
    END DO

    ! get the index of the used elements
    ALLOCATE (index_unique_tmp(num_obs_total))
    DO i = 1, num_obs_total
      index_unique_tmp(i) = i
    END DO
    WHERE (ABS(obspres - missing) .LT. 0.01D0) index_unique_tmp = 0
    ! WHERE (ABS(obsHght - missing) .LT. 0.01D0) index_unique_tmp = 0
    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    PRINT *, "the amount of all sounding input obs and the unique: ", &
      num_obs_total, iobs, num_obs_total_unique, x%mpddGlob%myrank

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

    PRINT *, "DONE soundIngest"

    ! Modified by Zilong to convert the pressure observations to density
    BLOCK
      USE parameters_m, ONLY: dry_air_gas_const
      USE conversions_m, ONLY: Td_to_qvapor

      IF (X%getVarIdx('qvapor') /= 0) THEN
        this%varNames(2) = "qvapor"

        ! IF (TRIM(this%StNames(j)) == '57687') PRINT *, 'SOUND of 57687:'
        DO j = 1, this%numObs
          IF (this%obsData(j, 2) == missing .OR. this%obsData(j, 3) == missing .OR. this%obsHght(j) > 18000.0) THEN
            this%obsData(j, 2) = missing
          ELSE
            this%obsData(j, 2) = Td_to_qvapor(this%obsData(j, 2), this%obsData(j, 3))
          END IF

        END DO
      END IF

      IF (X%getVarIdx('rho') /= 0 .OR. X%getVarIdx('rho_ctl') /= 0) THEN
        this%varNames(3) = "rho"
        DO j = 1, this%numObs
          IF (this%obsData(j, 3) == missing .OR. this%obsData(j, 1) == missing) THEN
            this%obsData(j, 3) = missing
          ELSE
            this%obsData(j, 3) = this%obsData(j, 3) / dry_air_gas_const / this%obsData(j, 1)
          END IF
        END DO
      END IF

    END BLOCK

    this%presIdx = 0
    this%tempIdx = 0
    this%qvaporIdx = 0

    DO i = 1, this%numVars
      IF (this%varNames(i) == "pres") this%presIdx = i
      IF (this%varNames(i) == "temp") this%tempIdx = i
      IF (this%varNames(i) == "qvapor") this%qvaporIdx = i
    END DO

    ! Added by Zilong for temporary scaling
    ! IF (this%qvaporIdx /= 0) this%obsData(:, this%qvaporIdx) = this%obsData(:, this%qvaporIdx)*1000.0D0
    ! IF (this%presIdx /= 0) this%obsData(:, this%presIdx) = this%obsData(:, this%presIdx)/100.0D0
    IF (this%tempIdx /= 0) this%qcThreshold(this%tempIdx) = 8.0D0
    IF (this%presIdx /= 0) this%qcThreshold(this%presIdx) = 1000.0D0
    IF (this%qvaporIdx /= 0) this%qcThreshold(this%qvaporIdx) = 0.003D0
    this%qcThreshold(4) = 8.0D0
    this%qcThreshold(5) = 8.0D0

    ! this%qcThreshold(4) = 10.0D0
    ! this%qcThreshold(5) = 10.0D0

  END SUBROUTINE soundIngest

  FUNCTION soundForward(this, X, O) RESULT(Y)
    IMPLICIT NONE

    CLASS(ObsSound_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j, k

    Y = O%zeroCopy()
    DO i = 1, UBOUND(Y%ObsFields, 1)
      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
            Y%ObsFields(i)%values(k) = X%fields(j)%Get_Value(Y%ObsFields(i)%idx(k))
          END DO
        END IF
      END DO
    END DO
  END FUNCTION soundForward

  FUNCTION soundTangent(this, dX, X) RESULT(Y)
    IMPLICIT NONE

    CLASS(ObsSound_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j, k

    PRINT *, 'Sound: Testing tangent operator comparing to ctl2State code, the code stops here'
    ! DO i = 1, UBOUND(Y%ObsFields, 1)
    !   DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !     IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
    !       DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
    !         Y%ObsFields(i)%values(k) = X%fields(j)%Get_Value(Y%ObsFields(i)%idx(k))
    !       END DO
    !     END IF
    !   END DO
    ! END DO
    STOP
  END FUNCTION soundTangent

  FUNCTION soundAdjoint(this, dY, X) RESULT(dX)
    IMPLICIT NONE

    CLASS(ObsSound_t) :: this
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
  END FUNCTION soundAdjoint

  ! There are many levels in sounding data missing height information, except those mandate levels
  ! This code intends to fill the missing heights along with the missing temperature and dewpoint
  ! This is designed and developed by Yuanfu Xie

  ! Assume:
  !   r_obs array hold
  !     temperature at 1st dim;
  !     dewpoint at 2nd;
  !     pressure at 3rd.
  SUBROUTINE fillinMissing(this, x, site, rank, missing, n_obs, r_obs, n_use, h_obs, t_obs, llobs, error, lands)
    IMPLICIT NONE
    CLASS(ObsSound_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: n_obs, rank
    INTEGER(i_kind), INTENT(OUT) :: n_use
    INTEGER(i_kind), INTENT(INOUT) :: t_obs(n_obs)
    REAL(r_kind), INTENT(IN) :: missing
    REAL(r_kind), INTENT(INOUT) :: r_obs(n_obs, 5), &
                                   h_obs(n_obs), llobs(1:2, n_obs), error(n_obs, 1:5), lands(n_obs)

    CHARACTER(LEN=*) site

    ! Local variables:
    INTEGER(i_kind), PARAMETER :: nPolynomial = 25, nx = nPolynomial + 1
    INTEGER(i_kind) :: i, j, k, idx, numFit, nfill, lowHeight, uppHeight, lowBound, mask(n_obs)
    REAL(r_kind) :: coef(nx), fit_range(2), xxx, temp(5), qvapor(n_obs)

    REAL(r_kind), PARAMETER :: dryHydro = 50000.0D0   ! This is threshold value for no qvapor above this level

    ! QR decomposition:
    INTEGER :: job, info, jpvt(nx)
    DOUBLE PRECISION, ALLOCATABLE :: cmatrix(:, :)
    DOUBLE PRECISION :: fitting(n_obs, 2), xmin, xmax
    DOUBLE PRECISION :: aux(nx), work(nx), qy(n_obs), qty(n_obs), sol(nx), xb(n_obs), rsd(n_obs)

!#define DEBUG_SOUNDING
    ! Debugging:
    CHARACTER(LEN=8), PARAMETER :: site_debug = '59981' !'59644' ! Note need to check if this station fall in the processor, myrank!
    INTEGER(i_kind), PARAMETER :: myrank_debug = 2, varid_debug = 5 ! 1 temp; 2 dew

    ! Merge obs with identical pressure levels:
    DO j = 1, n_obs
      DO i = j + 1, n_obs
        IF (ABS(r_obs(j, 3) - r_obs(i, 3)) .LT. 1.0D-8 .AND. r_obs(j, 3) .LT. missing) THEN
          ! For identical pressure levels:
          DO k = 1, 5  ! Merge the obs fields
            IF (k .NE. 3) THEN ! Pressure values do not change
              IF (r_obs(i, k) .LT. missing .AND. r_obs(j, k) .LT. missing) THEN
                IF (h_obs(j) .LT. missing) THEN ! Trust the obs with height information
                  r_obs(i, k) = r_obs(j, k)
                ELSE
                  r_obs(i, k) = 0.5D0 * (r_obs(i, k) + r_obs(j, k))
                END IF
              ELSE IF (r_obs(j, k) .LT. missing) THEN
                r_obs(i, k) = r_obs(j, k)
              END IF
            END IF
            IF (error(i, k) .LT. missing .AND. error(j, k) .LT. missing) THEN
              error(i, k) = 0.5D0 * (error(i, k) + error(j, k))
            ELSE IF (error(j, k) .LT. missing) THEN
              error(i, k) = error(j, k)
            END IF
          END DO
          h_obs(i) = MIN(h_obs(j), h_obs(i))
          t_obs(i) = MIN(t_obs(j), t_obs(i))
          llobs(1, i) = MIN(llobs(1, j), llobs(1, i))
          llobs(2, i) = MIN(llobs(2, j), llobs(2, i))
          lands(i) = MIN(lands(j), lands(i))

          r_obs(j, 1:5) = missing ! void the duplicated obs
          h_obs(j) = missing ! Void the height obs as well
        END IF
      END DO
    END DO

    ! Remove unused levels:
    n_use = 0
    DO j = 1, n_obs
      ! Valid pressure:
      IF (r_obs(j, 3) .LT. missing) THEN
        n_use = n_use + 1
        r_obs(n_use, 1:5) = r_obs(j, 1:5)
        h_obs(n_use) = h_obs(j)
        t_obs(n_use) = t_obs(j)
        llobs(1:2, n_use) = llobs(1:2, j)
        error(n_use, 1:5) = error(j, 1:5)
        lands(n_use) = lands(j)
      ELSE
#ifdef DEBUG_SOUNDING
        IF (h_obs(j) .LT. missing) &
          WRITE (*, 11) j, n_obs, r_obs(j, 3), h_obs(j), rank, TRIM(site)
11      FORMAT('obs height - missing pressure: ', I3, I4, ' HghtPres: ', 2D12.4, ' pc:', I2, ' stn: ', A)
#endif
      END IF
    END DO

    ! Possible sorting according to the pressure levels: from ground to the top
    DO j = 1, n_use
      DO i = 1, n_use - j
        IF (r_obs(i, 3) .LT. r_obs(i + 1, 3)) THEN
          temp(1:5) = r_obs(i, 1:5)
          r_obs(i, 1:5) = r_obs(i + 1, 1:5)
          r_obs(i + 1, 1:5) = temp(1:5)

          temp(1) = h_obs(i)
          h_obs(i) = h_obs(i + 1)
          h_obs(i + 1) = temp(1)

          temp(1) = t_obs(i)
          t_obs(i) = t_obs(i + 1)
          t_obs(i + 1) = temp(1)

          temp(1:2) = llobs(1:2, i)
          llobs(1:2, i) = llobs(1:2, i + 1)
          llobs(1:2, i + 1) = temp(1:2)

          temp(1:5) = error(i, 1:5)
          error(i, 1:5) = error(i + 1, 1:5)
          error(i + 1, 1:5) = temp(1:5)

          temp(1) = lands(i)
          lands(i) = lands(i + 1)
          lands(i + 1) = temp(1)
        END IF
      END DO
    END DO

    ! Fill in temperature and dew point:
    ! DO idx=1,2
    DO idx = 1, 5
      IF (idx .EQ. 3) CYCLE ! idx=3 is the pressure field, no fillin needed

      numFit = 0
      nfill = 0
      fitting = 0.0D0
      fit_range(1) = 1.0D10 ! Range of the fitting in vertical
      fit_range(2) = 0.0D0
      DO j = 1, n_use
        IF (r_obs(j, idx) .LT. missing) THEN
          ! Count fitting points:
          numFit = numFit + 1
          fitting(numFit, 1) = LOG(r_obs(j, 3))
          fitting(numFit, 2) = r_obs(j, idx)
          IF (r_obs(j, 3) .GT. fit_range(2)) fit_range(2) = r_obs(j, 3)
          IF (r_obs(j, 3) .LT. fit_range(1)) fit_range(1) = r_obs(j, 3)
        ELSE
          nfill = nfill + 1
        END IF
      END DO
      ! Cannot fit with too little obs:
      IF (numFit .LT. nx) THEN
#ifdef TRACE_OBS
        WRITE (*, 1) TRIM(site), idx, numFit, nPolynomial, n_use, n_obs, r_obs(1, idx), rank
1       FORMAT('FillingMissing: ', A, ' var: ', I1, ' - too little data: ', 2I2, &
               ' in use: ', I3, ' total: ', I3, ' values: ', D10.2, ' pc', I2)
#endif
        CYCLE
      END IF

      ! Fitting the profile:
      ALLOCATE (cmatrix(numFit, nx))
      xmin = LOG(fit_range(1)) ! MINVAL(fitting(1:numFit,1))
      xmax = LOG(fit_range(2)) ! MAXVAL(fitting(1:numFit,1))
      DO j = 1, numFit
        DO i = 0, nPolynomial
          cmatrix(j, i + 1) = ((fitting(j, 1) - xmin) / (xmax - xmin))**i
        END DO
      END DO
      ! Solve the coefficients:
      job = 1 ! Pivoting
      jpvt = 0
      CALL dqrdc(cmatrix, numFit, numFit, nx, aux, jpvt, work, job)
      job = 110 ! Just compute the xbb only
      rsd = 0.0D0
      CALL dqrsl(cmatrix, numFit, numFit, nx, aux, fitting(:, 2), qy, qty, sol, rsd, xb, job, info)

      ! Pivot back:
      DO i = 1, nx
        coef(jpvt(i)) = sol(i)
      END DO

#ifdef DEBUG_SOUNDING
      WRITE (*, 2) nPolynomial, info, sol(1:nx)
      WRITE (*, 3) nPolynomial, info, coef(1:nx)
2     FORMAT('Least square solution: ', I2, ' info: ', I3, ' sol: ', 30D10.2)
3     FORMAT('Least square pivoted : ', I2, ' info: ', I3, ' sol: ', 30D10.2)
      WRITE (*, 4) idx, MAXVAL(rsd), MAXVAL(fitting(:, 2)), x%mpddGlob%myrank
4     FORMAT('Max residual of ', I1, ' is ', D12.4, ' Max values: ', D12.4, ' pc: ', I2)

      ! Save the profile to plot:
      IF (x%mpddGlob%myrank .EQ. myrank_debug) THEN
        IF (idx .EQ. varid_debug .AND. TRIM(site) .EQ. TRIM(site_debug)) THEN
          DO i = 1, numFit
            xb(i) = 0.0D0
            DO j = 1, nx
              xb(i) = xb(i) + coef(j) * ((fitting(i, 1) - xmin) / (xmax - xmin))**(j - 1)
            END DO

            !WRITE(11,5) i,(fitting(i,1)-xmin)/(xmax-xmin), &
            !  fitting(i,2)/MAXVAL(fitting(:,2)),xb(i)/MAXVAL(fitting(:,2))
            WRITE (11, 5) i, fitting(i, 1), LOG(fitting(i, 1)), &
              xb(i), fitting(i, 2)
5           FORMAT(I4, 4F12.4)
          END DO

          nfill = 0
          DO i = 1, n_use
            IF (r_obs(i, idx) .GT. missing - 1.0D0) THEN
              nfill = nfill + 1
              PRINT *, 'Missing temDew: ', i, n_use, idx, r_obs(i, 3), r_obs(i, idx)
              r_obs(i, idx) = 190.0D0
            END IF
          END DO
          PRINT *, 'Total of filling: ', nfill, ' number of used: ', n_use
          !xmin = MINVAL(LOG(r_obs(1:n_use,3)))
          !xmax = MAXVAL(LOG(r_obs(1:n_use,3)))
          !print*,'MMX2: ',xmin,xmax,rank,idx,MAXVAL(r_obs(1:n_use,3)),MINVAL(r_obs(1:n_use,3))
          DO i = 1, n_use
            xb(i) = 0.0D0
            DO j = 1, nx
              xxx = LOG(r_obs(i, 3))
              xb(i) = xb(i) + coef(j) * ((xxx - xmin) / (xmax - xmin))**(j - 1)
            END DO

            ! No extrapolation of the polynomial:
            IF (r_obs(i, 3) .LT. fit_range(1) .OR. r_obs(i, 3) .GT. fit_range(2)) CYCLE

            IF (r_obs(i, idx) .GT. missing - 1.0D0) THEN
              WRITE (12, 6) i, r_obs(i, 3), LOG(r_obs(i, 3)), xb(i)
            ELSE
              WRITE (12, 6) i, r_obs(i, 3), LOG(r_obs(i, 3)), xb(i), r_obs(i, idx)
            END IF
6           FORMAT(I4, 4F12.4)
          END DO

        END IF
      END IF
#endif

      ! Filling in temperature and dewpoint:
      DO i = 1, n_use
        ! No extrapolation of the polynomial:
        IF (r_obs(i, 3) .LT. fit_range(1) .OR. r_obs(i, 3) .GT. fit_range(2)) CYCLE

        ! No replacement for dewpoint if the obs is valid:
        IF (idx .EQ. 2 .AND. r_obs(i, idx) .LT. missing) CYCLE

        xb(i) = 0.0D0
        DO j = 1, nx
          xxx = LOG(r_obs(i, 3))
          xb(i) = xb(i) + coef(j) * ((xxx - xmin) / (xmax - xmin))**(j - 1)
        END DO
        ! Replace all temperature or dewpoint obs with the polynomial fitting:
        r_obs(i, idx) = xb(i)
      END DO

      DEALLOCATE (cmatrix)
    END DO

    ! Heights fill in:

    ! Masking the vertical levels with no fill: 0; need fill: 1; and with height: 2
    mask = 0
    qvapor = 0.0D0
    DO i = 1, n_use
      IF (r_obs(i, 3) .LT. dryHydro) THEN
        ! No qvapor is needed for hydrostatic height:
        IF (r_obs(i, 1) .LT. missing) THEN
          qvapor(i) = 0.0D0
          IF (h_obs(i) .GT. missing - 1.0D0) THEN
            mask(i) = 1
          ELSE
            mask(i) = 2
          END IF
        END IF
      ELSE
        ! Qvapor is needed for a hydrostatic height:
        IF (r_obs(i, 1) .LT. missing .AND. r_obs(i, 2) .LT. missing) THEN
          qvapor(i) = Td_to_qvapor(r_obs(i, 2), r_obs(i, 3))
          IF (h_obs(i) .GT. missing - 1.0D0) THEN
            mask(i) = 1
          ELSE
            mask(i) = 2
          END IF
        END IF
      END IF
    END DO

    ! Filling in heights: assuming the pressure levels have been sorted above
    lowHeight = 0
    uppHeight = -1
    lowBound = 0
    nfill = 0
    DO i = 1, n_use
      ! Find the consecutive fillin interval:
      IF (mask(i) .EQ. 1) THEN
        IF (uppHeight .LT. lowHeight) THEN
          lowHeight = i; uppHeight = i  ! First encounter of fillin section
        END IF
        ! Hitting the top:
        IF (i .EQ. n_use) THEN
          IF (lowBound .GT. 0) THEN
            CALL hydroHeightUp(i - lowBound + 1, r_obs(lowBound:i, 3), &
                               r_obs(lowBound:i, 1), qvapor(lowBound:i), h_obs(lowBound:i))
            nfill = nfill + n_use - lowBound
          END IF
        ELSE
          uppHeight = i
        END IF
      ELSE IF (mask(i) .EQ. 0) THEN ! No fill level
        IF (lowBound .GT. 0) THEN
          CALL hydroHeightUp(i - lowBound, r_obs(lowBound:i - 1, 3), &
                             r_obs(lowBound:i - 1, 1), qvapor(lowBound:i - 1), h_obs(lowBound:i - 1)) ! no fill in current level
          nfill = nfill + uppHeight - lowHeight + 1
        END IF
        ! Reset:
        lowBound = 0
        lowHeight = 0
        uppHeight = -1
      ELSE ! mask = 2
        ! Filling the section:
        IF (uppHeight .GE. lowHeight) THEN
          CALL hydroHeightDown(i - lowHeight + 1, r_obs(lowHeight:i, 3), &
                               r_obs(lowHeight:i, 1), qvapor(lowHeight:i), h_obs(lowHeight:i))
          nfill = nfill + uppHeight - lowHeight + 1
        END IF
        lowHeight = i + 1; uppHeight = i; lowBound = i
      END IF
    END DO
#ifdef TRACE_OBS
    WRITE (*, 7) nfill, TRIM(site), llobs(1:2, 1), x%mpddGlob%myrank, x%sg%gLevel
7   FORMAT('fillinSoundHeights - total fillin: ', I4, ' SiteName: ', A, ' obsLL: ', 2D12.4, ' pc: ', I2, ' Glvl: ', I2)
#endif

  END SUBROUTINE fillinMissing

  SUBROUTINE soundQC(this)
    CLASS(ObsSound_t) :: this
  END SUBROUTINE soundQC

END MODULE ObsSound_m
