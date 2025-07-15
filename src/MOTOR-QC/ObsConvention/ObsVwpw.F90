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
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, @GBA-MWF, Shenzhen for adding obs
!     forward, tangent and adjoint operators.
!   Modified by Yongjian Huang (huang.yongjian0709@gmail.com), 2025/2/5, @GBA-MWf, Shenzhen for adding QC
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
  USE parameters_m, ONLY: degree2radian, missing
  USE FastMCDDetector_m
  USE ObsRaw_m
  USE domainCheck_m

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsVwpw_t

  CONTAINS
    PROCEDURE :: ObsInitial => vwpwInitial
    PROCEDURE :: ObsIngest => vwpwIngest
    PROCEDURE :: ObsForward => vwpwForward
    PROCEDURE :: ObsTangent => vwpwTangent
    PROCEDURE :: ObsAdjoint => vwpwAdjoint
    PROCEDURE :: vwpwQC

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

    ! read the var names of obs files in configFile
    this%numVars = 2
    ALLOCATE (this%varNames(this%numVars))
    this%varNames(1) = "uwnd"
    this%varNames(2) = "vwnd"

    this%obsType = "VWPW"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius = 5.0D2        ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(3, :) = 5.0D2  ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1  ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.25D0 ! No influence under this threshold value

    ! Threshold values for QC:
    !! set big enough values ignore the obs thinning qc
    this%qcThreshold(1) = 100.0D0
    this%qcThreshold(2) = 100.0D0

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
    TYPE(ObsRaw_t) :: rawObs, bcgObs
    TYPE(FastMCDDetector_t) :: outlier_detector
    INTEGER(i_kind), ALLOCATABLE :: qc_mask(:)

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
          obsData(iobs, 1) = u
          obsData(iobs, 2) = v
          obsHght(iobs) = hgh

          CALL Time_GMT_to_Unix((/iyear, imonth, iday, ihh, imm, 0/), obsTime(iobs))

          olatlon(1, iobs) = rlat   ! lat
          olatlon(2, iobs) = rlon   ! lon

          land(iobs) = 1.0D0
          obsErrs(iobs, :) = missing * (-1.0D0)   !-999999999.0D0
          stNames(iobs) = siteId
        END DO
      END DO
      CLOSE (101)
      PRINT *, 'Number of Obs in VWPW file ', ifile, ':', num_obsin(ifile)
    END DO

    WHERE (ABS(obsData - missingValue) .LT. 0.01D0) obsData = missing
    WHERE (ABS(obsData(:, 1) - missing) .LT. 0.01D0) obsData(:, 2) = missing
    WHERE (ABS(obsData(:, 2) - missing) .LT. 0.01D0) obsData(:, 1) = missing
    WHERE (ABS(obsHght - missingValue) .LT. 0.01D0) obsHght = missing

    ! get the index of the used elements
    ALLOCATE (index_unique_tmp(num_obs_total))
    DO i = 1, num_obs_total
      index_unique_tmp(i) = i
    END DO
    WHERE (ABS(obsHght - missing) .LT. 0.01D0) index_unique_tmp = 0
    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    PRINT *, "before VWPW QC: ", num_obs_total, num_obs_total_unique

    ALLOCATE (qc_mask(num_obs_total_unique))
    CALL this%vwpwQC(X, num_obs_total_unique, obsData(index_unique, :), olatlon(:, index_unique) * degree2radian, &
                     obsTime(index_unique), obsHght(index_unique), qc_mask)

    !! re-mask USING QC
    j = 1
    DO i = 1, num_obs_total
      IF (index_unique_tmp(i) .NE. 0) THEN
        IF (qc_mask(j) .EQ. 0) THEN
          index_unique_tmp(i) = 0
        END IF
        j = j + 1
      END IF
    END DO

    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    DEALLOCATE (index_unique)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    PRINT *, "After VWPW QC: ", num_obs_total, num_obs_total_unique

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

    this%qcThreshold(1) = 10.0D0 !25.0D0
    this%qcThreshold(2) = 10.0D0 !25.0D0
    PRINT *, "DONE vwpwIngest QC Threshold values for wind: ", this%qcThreshold(1:2)

  END SUBROUTINE vwpwIngest


  FUNCTION vwpwForward(this, X, O) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsVwpw_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'VWPW: This forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION vwpwForward

  FUNCTION vwpwTangent(this, dX, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(ObsVwpw_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'VWPW: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION vwpwTangent

  FUNCTION vwpwAdjoint(this, dY, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(ObsVwpw_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'VWPW: This adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION vwpwAdjoint

  SUBROUTINE vwpwQC(this, state, num_unique, obsData, olatlon, obsTime, obsHght, maskValided)

    CLASS(ObsVwpw_t) :: this

    TYPE(State_t), INTENT(IN)  :: state
    REAL(r_kind), INTENT(IN) :: obsData(:, :), olatlon(:, :), obsHght(:)
    INTEGER(i_kind), INTENT(IN) :: obsTime(:)
    INTEGER(i_kind), INTENT(IN) :: num_unique

    INTEGER(i_kind), INTENT(INOUT) :: maskValided(num_unique)

    ! Thinning parameters:
    REAL(r_kind), ALLOCATABLE :: weights(:, :, :, :)       ! Observation weight
    REAL(r_kind), ALLOCATABLE :: wghtobs(:, :, :, :)       ! Weighted obs: sum w*obs
    REAL(r_kind), ALLOCATABLE :: wghterr(:, :, :, :)       ! Weighted err: sum w*err
    INTEGER, ALLOCATABLE :: maskQC(:)

    ! Local variables:
    CHARACTER(LEN=20) :: varname
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domain
    TYPE(FastMCDDetector_t) :: outlier_detector

    ! Loop variables:
    INTEGER(i_kind) :: io, iu, iv, i, j, k, istatus
    REAL(r_kind) :: bkgdAtObs, gaussian, t1, t2, t11, t22, t111, t222, bkgd

    ! Consider to save these variables in ObsBase:
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: &
      idxGrdValided(:, :), idxHgtValided(:, :, :), idxTimValided(:, :)
    REAL(r_kind), ALLOCATABLE :: &
      coeGrdValided(:, :), coeHgtValided(:, :, :), coeTimValided(:, :), &
      forward(:, :, :) ! vertical, horizontal and time
    REAL(r_kind), ALLOCATABLE :: valided_obsData(:, :), valided_bkgdata(:, :), &
                                 valided_olatlon(:, :), valided_height(:), valided_tempor(:)

    ! Validating the obs and saving their interpolation coefficients:
    CALL domain%validation(state%sg, num_unique, olatlon, obsHght, &
                           obsTime, numValided, maskValided, &
                           idxGrdValided, idxHgtValided, idxTimValided, &
                           coeGrdValided, coeHgtValided, coeTimValided, &
                           this%interpolation, nhstencil, this%obsType)

    IF (numValided .NE. 0) THEN
      ! Allocate valided parameter arrays:
      ALLOCATE (valided_obsData(numValided, this%numVars), &
                valided_bkgdata(numValided, this%numVars), &
                valided_olatlon(2, numValided), &
                valided_height(numValided), &
                valided_tempor(numValided))

      ALLOCATE (forward(2, nhstencil, 2))  ! vertical, horizontal, temporal

      CALL CPU_TIME(t1)

      ! Pass the original arrays to the valided arrays:
      i = 0
      DO io = 1, num_unique
        IF (maskValided(io) .GT. 0) THEN
          i = i + 1
          valided_obsData(i, :) = obsData(io, 1:this%numVars)
          valided_olatlon(:, i) = olatlon(:, io)
          valided_height(i) = obsHght(io)
          valided_tempor(i) = obsTime(io)

        END IF
      END DO

      DO io = 1, numValided
        DO iv = 1, this%numVars
          varname = TRIM(this%varNames(iv))
          bkgdAtObs = 0D0
          DO i = 1, 2 ! 2 time frames of interpolation
            DO j = 1, nhstencil ! 3 horizontal interpolation points
              DO k = 1, 2 ! 2 vertical levels
                forward(k, j, i) = this%GetForwardValue(state, varname, &
                                                        idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv)
                ! PRINT*, 'HYJ+++, fwd:', forward(k, j, i)r
                ! Interpolate the background at observation:
                bkgdAtObs = bkgdAtObs + forward(k, j, i) * &
                            coeHgtValided(k, j, io) * coeGrdValided(j, io) * coeTimValided(i, io)
              END DO
            END DO
          END DO
          valided_bkgdata(io, iv) = bkgdAtObs
          ! PRINT*, 'HYJ+++, bkg:', valided_bkgdata(io, iv), 'obs', valided_obsData(io, iv)
        END DO
      END DO

      ! maskQC greater than zero IS VALID
      ALLOCATE (maskQC(numValided))
      CALL outlier_detector%detect('FAST-MCD', valided_obsData, valided_bkgdata, 0D0, 0D0, maskQC)

      j = 1
      DO i = 1, num_unique
        IF (maskValided(i) .GT. 0) THEN
          maskValided(i) = maskQC(j)
          j = j + 1
        END IF
      END DO

      PRINT*, 'HYJ+++', j, numValided

      DEALLOCATE (valided_obsData, &
                  valided_bkgdata, &
                  valided_olatlon, &
                  valided_height, &
                  valided_tempor, maskQC)
    ELSE
      ! ALLOCATE(maskQC(num_unique))
      maskValided(:) = 1
    END IF

  END SUBROUTINE vwpwQC

END MODULE ObsVwpw_m
