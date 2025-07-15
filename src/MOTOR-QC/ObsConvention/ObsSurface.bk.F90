!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2021-09-01, created by Yuanfu Xie for SYNOP Observations.
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang.
!                     2022-01-25, adjusted the parameters of raduis fro SuperObs by Jiongming Pang.
!                     2022-01-25, expanded to read more variables by Jiongming Pang.
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2021/9/01, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/25, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/06/17 For limiting obs only on own processor
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/07/25 For holding off a given set of obs
!     specified in a soloSite.yaml file.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/08/01 For changing obs varList to
!   this%varNames so that it can read in obs related to these analysis variables.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a surface data structure.
MODULE ObsSurface_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m
  USE slint, ONLY: slint_init, tgt_grid
  USE mpObs_m, ONLY: mpObs_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE ncReadVar_m
  USE WriteVar_m
  USE conversions_m
  USE parameters_m, ONLY: degree2radian

  ! Yuanfu Xie added Geobox for checking obs:
  USE geoTools_m, ONLY: GeoBox_t

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsSurface_t
    CHARACTER(LEN=20), ALLOCATABLE :: obsVars(:)    ! Obs variable names

  CONTAINS
    PROCEDURE :: ObsInitial => sfcInitial
    PROCEDURE :: ObsIngest => sfcIngest
    PROCEDURE :: ObsQC => sfcQC

    PROCEDURE, PUBLIC :: GetForwardValue
  END TYPE ObsSurface_t

CONTAINS

  SUBROUTINE sfcInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsSurface_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: obsFileDir!, static_dir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    CHARACTER(LEN=20), ALLOCATABLE :: obsVarList(:)
    INTEGER :: i, istatus

    ! SoloSite variables:
    CHARACTER(LEN=1024) :: soloFile
    LOGICAL :: soloSite

    ! read the obs surface files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Surface', obsFileDir)
!     istatus = yaml_get_var(configFile, 'IO', 'static_dir', static_dir)
!     WRITE(*,1) TRIM(obsFileDir),TRIM(static_dir)
! 1   FORMAT('SFC_obsFileDir: ',A,' STATIC_DIR: ',A)

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
11  FORMAT('SfcIngest - Obs variable lists: NumVar: ', I2, 1X, 20A5, ' proc: ', I2)

    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Surface', obsFileName)
    this%numFiles = UBOUND(obsFileName, 1)
    WRITE (*, 2) this%numFiles
2   FORMAT("Number of surface obs files: ", I3)

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numFiles
      this%fileNames(i) = TRIM(obsFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of surface obs file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! Read in an interpolation option:
    !CALL namelist_read(TRIM(configFile), "interpolation",this%interpolation)
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

    ! Read the soloSite.yaml information for holding off the obs in that file: by Yuanfu Xie 2022-07-25
    ! soloFile = TRIM(static_dir)//"/"//'soloSite.yaml'
    istatus = yaml_get_var(TRIM(configFile), 'Verify', 'SoloFile', soloFile) ! Changed by Zilong

    INQUIRE (FILE=TRIM(soloFile), EXIST=soloSite)
    IF (soloSite) THEN
      WRITE (*, 8) TRIM(soloFile)
8     FORMAT('sfcInitial - soloSite file: ', A)
      istatus = yaml_get_var(TRIM(soloFile), 'sites', 'name', this%soloNames)
      IF (istatus .EQ. 0) THEN
        istatus = yaml_get_var(TRIM(soloFile), 'sites', 'latitude', this%soloLat)
        IF (istatus .EQ. 0) THEN
          istatus = yaml_get_var(TRIM(soloFile), 'sites', 'longitude', this%soloLon)
          IF (istatus .NE. 0) THEN
            WRITE (*, 3) TRIM(soloFile)
3           FORMAT('Station longitudes of SoloSite file cannot be read: ', A)
            STOP
          END IF
          WRITE (*, 7)
7         FORMAT('sfcInitial - found soloSite.yaml file and successfully read it in...')
        ELSE
          WRITE (*, 4) TRIM(soloFile)
4         FORMAT('Station latitudes of SoloSite file cannot be read: ', A)
          STOP
        END IF
      ELSE
        WRITE (*, 5) TRIM(soloFile)
5       FORMAT('Station names of SoloSite file cannot be read: ', A)
        STOP
      END IF
    ELSE
      WRITE (*, 6)
6     FORMAT('srcInitial - there is no soloSite file under static to hold off obs!')
    END IF
    ! Debugging
!     DO i=1,UBOUND(this%soloNames,1)
!       WRITE(*,10) i,this%soloNames(i),this%soloLat(i),this%soloLon(i)
! 10    FORMAT('Solo station: ',I3,1X,A,2D12.4)
!     END DO

    ! set the var names of obs files into this%obsVars(:)
    ! Modified by Yuanfu Xie 2022/08/01 to use analysis variables:
    this%numVars = UBOUND(this%varNames, 1)
    ALLOCATE (this%obsVars(this%numVars))

    ! Set obsVars to read by Yuanfu Xie 2022-08-01:
    ! Temporarily hard code these var names defined in the NC files,
    ! Consider these in a NC table file in the future:
    DO i = 1, this%numVars
      SELECT CASE (this%varNames(i))
      CASE ('temp')
        this%obsVars(i) = "temperature"
      CASE ('qvapor')
        this%obsVars(i) = "relHumidity"
      CASE ('wspd', 'uwnd')
        this%obsVars(i) = "windSpeed"
      CASE ('cdir', 'sdir', 'vwnd')
        this%obsVars(i) = "windDir"
      CASE ('pres')
        this%obsVars(i) = "stationPressure"
      CASE ('pcpa')
        this%obsVars(i) = "precipAccum"
      CASE DEFAULT
        WRITE (*, 12) TRIM(this%varNames(i))
12      FORMAT('SfcIngest - Obs variable is not defined: ', A, ' Check and rerun!')
        STOP
      END SELECT
    END DO

    this%obsType = "SYNOP"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ! this%missing = -99999.9
    ! this%missing = -9999.0D0
    this%missing = 999999999.0D0  !-9999.0D0
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
    this%qcThreshold = 500.0D0 ! Temporary holders
    ! For more general thresholds setting, the following is temporarily on hold: Yuanfu Xie 2022-08-01
    !this%qcThreshold(4) = 1500.0D0*100.0D0  ! for psfc (stationPressure)
    !this%qcThreshold(2) =  200.0D0
  END SUBROUTINE sfcInitial

  SUBROUTINE sfcIngest(this, X)
    USE NETCDF
    IMPLICIT NONE

    CLASS(ObsSurface_t) :: this
    TYPE(State_t) :: X

    ! local variables
    INTEGER :: nf_op, nc_id, var_id, i, num_obs, num_dims, nf90_max_va_dims, iuv, ii
    INTEGER :: ifile, num_obs_total, num_dataFile, ind_obsin_beg, ind_obsin_end
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_obsin
    INTEGER, DIMENSION(nf90_max_var_dims) :: var_dim_ids

    ! local variables, for get the unique obs in multi-files
    TYPE(ncStr_t), ALLOCATABLE  ::  varstr

    REAL(r_kind), ALLOCATABLE :: obsData(:, :), obsErrs(:, :), olatlon(:, :), obsDataWsWd(:, :)
    REAL(r_kind), ALLOCATABLE :: obsHght(:), land(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    CHARACTER(LEN=6), ALLOCATABLE :: stNames(:)
    LOGICAL :: fileExist(this%numFiles) ! Yuanfu Xie added a safeguard for missing file 2022-05-31
    LOGICAL, ALLOCATABLE :: index_unique_judge(:)
    INTEGER, ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER :: num_obs_total_unique, j, soloLen
    CHARACTER(LEN=1024) :: fntxt_Idtime_a, fntxt_Idtime_u
    REAL(r_kind) :: missingValue = -9999.0D0

    ! Yuanfu Xie added geobox for checking obs inside processor's domain:
    TYPE(GeoBox_t)  :: GeoBox

    ! Yuanfu Xie added uv indices for options selecting wind obs as uv:
    INTEGER(i_kind) :: uIdx, vIdx, cIdx, sIdx, wIdx
    REAL(r_kind) :: wspd, wdir

    ! Get the total amount of observations in multi-files. added for test by Jiongming Pang
    num_dataFile = this%numFiles
    ALLOCATE (ind_obsin(num_dataFile))
    ! Yuanfu Xie added initialization to ind_obsin 2022-02-19:
    ind_obsin = 0
    DO ifile = 1, num_dataFile
      ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      INQUIRE (FILE=TRIM(this%fileNames(ifile)), EXIST=fileExist(ifile))
      IF (fileExist(ifile)) THEN
        WRITE (*, 4) TRIM(this%fileNames(ifile))
4       FORMAT('ObsSurface: surface data file: ', A, ' is found and reading...')
      ELSE
        WRITE (*, 3) TRIM(this%fileNames(ifile))
3       FORMAT('ObsSurface: This surface data file does not exist but the run continues, CHECK! ', A)
        CYCLE
      END IF
      nf_op = nf90_open(TRIM(this%fileNames(ifile)), nf90_nowrite, nc_id)
      ! Yuanfu Xie added safeguard 2022-02-19:
      IF (nf_op .EQ. NF90_NOERR) THEN
        nf_op = nf90_inq_varid(nc_id, TRIM(this%obsVars(1)), var_id)
        nf_op = nf90_inquire_variable(nc_id, var_id, ndims=num_dims)
        nf_op = nf90_inquire_variable(nc_id, var_id, dimids=var_dim_ids(:num_dims))
        ! observation variable should only have one dimension
        nf_op = nf90_inquire_dimension(nc_id, var_dim_ids(1), len=num_obs)
        nf_op = nf90_close(nc_id)
        ind_obsin(ifile) = num_obs
        WRITE (*, 1) num_obs, TRIM(this%fileNames(ifile))
1       FORMAT('Number obs read in: ', I8, /, ' from file: ', A)
      END IF
    END DO
    num_obs_total = SUM(ind_obsin)
    ! this%numObs = num_obs_total
    ! Yuanfu Xie modified the print for specifically stated as surface 2022-02-19:
    PRINT *, 'Total amount of surface Obs in all files: ', num_obs_total
    ! Yuanfu Xie added a case of no surface obs to return to caller 2022-02-19:
    IF (num_obs_total .EQ. 0) THEN
      this%numObs = 0
      RETURN
    END IF

    ! Allocated memory for arrays for all data from multi-nc-files.
    ! added by Jiongming Pang
    ALLOCATE (obsData(num_obs_total, this%numVars))
    ALLOCATE (obsErrs(num_obs_total, this%numVars))
    ALLOCATE (olatlon(2, num_obs_total))
    ALLOCATE (obsTime(num_obs_total))
    ALLOCATE (obsHght(num_obs_total))
    ALLOCATE (land(num_obs_total))
    ALLOCATE (stNames(num_obs_total))
    ALLOCATE (obsDataWsWd(num_obs_total, 2))

    obsData = this%missing

    ind_obsin_end = 0
    ! obs file loop
    DO ifile = 1, num_dataFile

      IF (.NOT. fileExist(ifile)) CYCLE  ! Yuanfu Xie added a safeguard for missing file 2022-05-31

      ind_obsin_beg = ind_obsin_end + 1
      ind_obsin_end = ind_obsin_beg + ind_obsin(ifile) - 1
      PRINT *, "ind_obsin_beg and ind_obsin_end:", ind_obsin_beg, ind_obsin_end
      PRINT *, "amount in calculation and file:", ind_obsin_end - ind_obsin_beg + 1, ind_obsin(ifile)

      ! Open the nc file
      nf_op = nf90_open(TRIM(this%fileNames(ifile)), nf90_nowrite, nc_id)
      ! Yuanfu Xie added safeguard 2022-02-19:
      IF (nf_op .EQ. NF90_NOERR) THEN
        ! get variables one by one
        DO i = 1, this%numVars
          ! get the id of the current variable
#ifdef DEBUG
          WRITE (*, 22) i, TRIM(this%varNames(i)), TRIM(this%obsVars(i)), ifile
22        FORMAT('SfcIngest - Var name of: ', I2, 1X, A, ' obsVar: ', A, ' from file: ', I2)
#endif
          nf_op = nf90_inq_varid(nc_id, TRIM(this%obsVars(i)), var_id)

          ! Get variable
          nf_op = nf90_get_var(nc_id, var_id, obsData(ind_obsin_beg:ind_obsin_end, i))

          ! Get variable error
          nf_op = nf90_inq_varid(nc_id, TRIM(this%obsVars(i))//'QCD', var_id)

          IF (nf_op .EQ. nf90_noerr) THEN
            nf_op = nf90_get_var(nc_id, var_id, obsErrs(ind_obsin_beg:ind_obsin_end, i))
          ELSE
            ! this variable does not have qcd4
            obsErrs(ind_obsin_beg:ind_obsin_end, i) = this%missing * (-1.0D0)   !-9999.0D0
          END IF

        END DO

        ! Convert "windSpeed" and "windDir" to "u10m" and "v10m" if uv are requested;
        ! otherwise, wspd and wdir, modified by Yuanfu Xie 2022-08-01:
        ! Search uv if requested:
        uIdx = 0
        vIdx = 0
        cIdx = 0
        sIdx = 0
        wIdx = 0
        DO i = 1, this%numVars
          IF (TRIM(this%varNames(i)) .EQ. 'uwnd') uIdx = i
          IF (TRIM(this%varNames(i)) .EQ. 'vwnd') vIdx = i
          IF (TRIM(this%varNames(i)) .EQ. 'cdir') cIdx = i
          IF (TRIM(this%varNames(i)) .EQ. 'sdir') sIdx = i
          IF (TRIM(this%varNames(i)) .EQ. 'wspd') wIdx = i
        END DO
        IF ((uIdx .EQ. 0 .AND. vIdx .NE. 0) .OR. (uIdx .NE. 0 .AND. vIdx .EQ. 0)) THEN
          WRITE (*, 15) uIdx, vIdx
15        FORMAT('SfcIngest - User did not provided both uwnd and vwnd in obsVarList_Surface', 2I2, /, &
                 'Edit your yaml and rerun!')
          STOP
        END IF

        ! For UV options:
        IF (uIdx .NE. 0 .AND. vIdx .NE. 0) THEN
          ! Convert wind to uv:
          DO iuv = ind_obsin_beg, ind_obsin_end
            IF (ABS(obsData(iuv, uIdx)) .GT. ABS(this%missing) .OR. &
                ABS(obsData(iuv, vIdx)) .GT. ABS(this%missing) .OR. &
                obsData(iuv, vIdx) .LT. 0.0D0 .OR. &
                obsData(iuv, vIdx) .GT. 500.0D0) THEN
              obsData(iuv, uIdx) = this%invalid
              obsData(iuv, vIdx) = this%invalid
            ELSE
              wspd = obsData(iuv, uIdx)  ! Note that when uwnd is selected, the NC reads wind speed in uwnd position
              wdir = obsData(iuv, vIdx)  ! Note that when vwnd is selected, the NC reads wind direc in vwnd position
              CALL wind_to_uv(wspd, wdir, obsData(iuv, uIdx), obsData(iuv, vIdx))
              wspd = obsErrs(iuv, uIdx)
              wdir = obsErrs(iuv, vIdx)
              CALL wind_to_uv(wspd, wdir, obsErrs(iuv, uIdx), obsErrs(iuv, vIdx))
            END IF
          END DO
        ELSE
          ! For wind dir and spd analysis option:
          IF (cIdx .NE. 0 .AND. sIdx .NE. 0 .AND. wIdx .NE. 0) THEN

            ! Convert wind to uv:
            DO iuv = ind_obsin_beg, ind_obsin_end
              IF (ABS(obsData(iuv, cIdx)) .GT. ABS(this%missing) .OR. &
                  ABS(obsData(iuv, wIdx)) .GT. ABS(this%missing) .OR. &
                  obsData(iuv, cIdx) .LT. 0.0D0 .OR. &
                  obsData(iuv, wIdx) .GT. 500.0D0) THEN
                obsData(iuv, cIdx) = this%invalid
                obsData(iuv, sIdx) = this%invalid
                obsData(iuv, wIdx) = this%invalid
              ELSE
                obsData(iuv, cIdx) = DCOS(obsData(iuv, cIdx))
                obsData(iuv, sIdx) = DSIN(obsData(iuv, sIdx))
                obsData(iuv, wIdx) = obsData(iuv, wIdx) * obsData(iuv, wIdx) ! Speed square
                obsErrs(iuv, cIdx) = DCOS(obsErrs(iuv, cIdx))
                obsErrs(iuv, sIdx) = DSIN(obsErrs(iuv, sIdx))
                obsErrs(iuv, wIdx) = obsErrs(iuv, wIdx) * obsErrs(iuv, wIdx)
              END IF
            END DO

          ELSE IF ((uIdx .EQ. 0 .AND. vIdx .NE. 0) .OR. &
                   (uIdx .NE. 0 .AND. vIdx .EQ. 0) .OR. &
                   (cIdx .EQ. 0 .AND. wIdx .NE. 0) .OR. &
                   (cIdx .NE. 0 .AND. wIdx .EQ. 0)) THEN
            WRITE (*, 17) uIdx, vIdx, cIdx, sIdx, wIdx
17          FORMAT('SfcIngest - incomplete wind vars: u: ', I2, ' v: ', I2, ' cdir: ', I2, ' sdir: ', I2, ' wspd: ', I2, /, &
                   ' Check your yaml file and rerun! ')
            STOP
          END IF
        END IF

        ! get observation latitude and longtitude
        nf_op = nf90_inq_varid(nc_id, 'latitude', var_id)
        nf_op = nf90_get_var(nc_id, var_id, olatlon(1, ind_obsin_beg:ind_obsin_end))
        nf_op = nf90_inq_varid(nc_id, 'longitude', var_id)
        nf_op = nf90_get_var(nc_id, var_id, olatlon(2, ind_obsin_beg:ind_obsin_end))

        ! get observation time
        nf_op = nf90_inq_varid(nc_id, 'observationTime', var_id)
        nf_op = nf90_get_var(nc_id, var_id, obsTime(ind_obsin_beg:ind_obsin_end))

        ! get topography
        nf_op = nf90_inq_varid(nc_id, 'elevation', var_id)
        nf_op = nf90_get_var(nc_id, var_id, obsHght(ind_obsin_beg:ind_obsin_end))

        ! observation data file does not have land mask information
        land(ind_obsin_beg:ind_obsin_end) = 1.0D0 ! Assume all land for testing

        ! close nc file
        nf_op = nf90_close(nc_id)
      END IF
    END DO ! obs file loop

    WHERE (ABS(obsData - missingValue) .LT. 0.01D0) obsData = this%missing

    ! get the stationId
    ind_obsin_end = 0
    DO ifile = 1, num_dataFile
      PRINT *, "obs file", ifile
      ind_obsin_beg = ind_obsin_end + 1
      ind_obsin_end = ind_obsin_beg + ind_obsin(ifile) - 1
      PRINT *, "ind_obsin_beg and ind_obsin_end:", ind_obsin_beg, ind_obsin_end
      ALLOCATE (varstr)
      CALL varstr%ncGetStr(TRIM(this%fileNames(ifile)), "stationId")
      stNames(ind_obsin_beg:ind_obsin_end) = varstr%ncVar(:)
      DEALLOCATE (varstr)
    END DO

    olatlon = olatlon * degree2radian   ! convert degree to radian. - Zilong
    GeoBox = GeoBox_t(X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob)

    ! PRINT*, 'X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob: ', &
    !           X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob
    ! get the index of unique elements
    ! Yuanfu Xie modified this uniqueness check to count data on own processor: 2022-06-17
    ALLOCATE (index_unique_judge(num_obs_total))
    ALLOCATE (index_unique_tmp(num_obs_total))
    index_unique_judge = .TRUE.
    index_unique_tmp = 0

    ! DO i = 1, num_obs_total
    !   ! Yuanfu Xie checks if the obs is inside the processor:
    !   IF (.NOT. GeoBox%inBox(olatlon(:,i))) THEN
    !     index_unique_judge(i) = .FALSE.
    !     index_unique_tmp(i) = 0
    !   END IF
    ! END DO

    ! Only use national stations for 3Dvar, Zilong
    IF (X%sg%vLevel == 1) THEN
      FORALL (i=1:num_obs_total, &
              (.NOT. GeoBox%inBoxElem(olatlon(1, i), olatlon(2, i), X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob)))
        index_unique_judge(i) = .FALSE.
        index_unique_tmp(i) = 0
      END FORALL
    ELSE
      FORALL (i=1:num_obs_total, &
              (.NOT. GeoBox%inBoxElem(olatlon(1, i), olatlon(2, i), X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob)) .OR. &
              ! stNames(i) (1:1) /= '5' .OR. &
              obsHght(i) > 1E4 &
              )
        index_unique_judge(i) = .FALSE.
        index_unique_tmp(i) = 0
      END FORALL
    END IF

    DO i = 1, num_obs_total
      ! Checks duplicated obs is in the soloSite:
      IF (index_unique_judge(i)) THEN
        DO j = i + 1, num_obs_total ! Yuanfu Xie changed this loop to i+1 instead of 1 saving some
          IF (.NOT. index_unique_judge(j)) CYCLE

          IF (obsTime(i) == obsTime(j)) THEN
            IF (stNames(i) == stNames(j)) THEN
              index_unique_judge(j) = .FALSE.
              index_unique_tmp(j) = 0
            END IF
          END IF

        END DO
        index_unique_tmp(i) = i
      END IF
    END DO

    ! Add soloSite blacklist to hold off some obs: Yuanfu Xie 2022-07-26
    IF (UBOUND(this%soloNames, 1) .GT. 0) THEN
      soloLen = LEN(TRIM(this%soloNames(1)))
      DO i = 1, num_obs_total
        IF (.NOT. index_unique_judge(i)) CYCLE

        ! Checks if the obs is in the soloSite:
        DO ii = 1, UBOUND(this%soloNames, 1)
          IF (TRIM(stNames(i) (1:soloLen)) .EQ. TRIM(this%soloNames(ii))) THEN
            index_unique_judge(i) = .FALSE.
            index_unique_tmp(i) = 0
#ifdef DEBUG
            ! WRITE(*,13) i,ii,TRIM(stNames(i))
#endif
13          FORMAT('ObsSurface obs station name: ', 2I6, 2X, A)
            EXIT
          END IF
        END DO
      END DO
    END IF

    num_obs_total_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_obs_total_unique))
    index_unique = PACK(index_unique_tmp, index_unique_tmp /= 0)

    WRITE (*, 12) num_obs_total, num_obs_total_unique, X%mpddGlob%myrank
12  FORMAT('SfcIngest - the amount of all input obs and accepted obs: ', 2I8, ' on proc: ', I3)

    ! ! get the unique elements Commented by Zilong, the ALLOCATABLE variable do not need to be initialized
    ! Not sure if the modern fortran has default initialization; otherwise, one has to initialize them

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
    this%olatlon = olatlon(:, index_unique)  ! convert degree to radian. - JMPang
    this%obsTime = obsTime(index_unique)
    this%obsHght = obsHght(index_unique)
    this%land = land(index_unique)
    this%StNames = stNames(index_unique)

    DO i = 1, this%numObs
      DO j = 1, this%numVars
        IF (this%obsData(i, j) .GE. this%qcThreshold(j)) this%obsData(i, j) = this%missing
      END DO
    END DO

    ! deallocate the local arrays
    IF (ALLOCATED(obsData)) DEALLOCATE (obsData)
    IF (ALLOCATED(obsErrs)) DEALLOCATE (obsErrs)
    IF (ALLOCATED(olatlon)) DEALLOCATE (olatlon)
    IF (ALLOCATED(obsTime)) DEALLOCATE (obsTime)
    IF (ALLOCATED(obsHght)) DEALLOCATE (obsHght)
    IF (ALLOCATED(land)) DEALLOCATE (land)
    IF (ALLOCATED(stNames)) DEALLOCATE (stNames)
    IF (ALLOCATED(index_unique_judge)) DEALLOCATE (index_unique_judge)
    IF (ALLOCATED(index_unique_tmp)) DEALLOCATE (index_unique_tmp)
    IF (ALLOCATED(index_unique)) DEALLOCATE (index_unique)
    IF (ALLOCATED(obsDataWsWd)) DEALLOCATE (obsDataWsWd)

    PRINT *, "DONE sfcIngest"

    ! Modified by Zilong to convert the pressure observations to density
    ! 1-temp, 2-RH, 4-pres
    BLOCK
      USE parameters_m, ONLY: dry_air_gas_const
      USE conversions_m, ONLY: RH_to_qvapor, uv_to_wind
      REAL(r_kind) :: u, v

      INTEGER(i_kind) :: iq, ip, it

      ! RH -> qvapor
      IF (X%getVarIdx('qvapor') /= 0) THEN
        iq = this%getVarIdx('qvapor')
        ip = this%getVarIdx('pres')
        it = this%getVarIdx('temp')
        IF (iq .NE. 0 .AND. ip .NE. 0 .AND. it .NE. 0) THEN
          DO j = 1, this%numObs
            IF (this%obsData(j, ip) == this%missing .OR. this%obsData(j, it) == this%missing &
                .OR. this%obsData(j, iq) == this%missing) THEN
              this%obsData(j, iq) = this%missing
            ELSE
              this%obsData(j, iq) = RH_to_qvapor(this%obsData(j, iq), this%obsData(j, it), this%obsData(j, ip)) * 1000.0D0
            END IF
          END DO
        ELSE
          WRITE (*, 20) ip, iq, it
20        FORMAT('sfcIngest - As the analysis requests qvapor, sfc analysis variables must have pres, qvapor and temp present in your yaml', /, &
                 'Please check and rerun! ', 3I2)
          STOP
        END IF
      END IF

      ! pres -> rho
      IF (X%getVarIdx('rho') /= 0 .OR. X%getVarIdx('rho_ctl') /= 0) THEN
        iq = this%getVarIdx('rho')
        ip = this%getVarIdx('pres')
        it = this%getVarIdx('temp')
        DO j = 1, this%numObs
          IF (this%obsData(j, ip) == this%missing .OR. this%obsData(j, it) == this%missing) THEN
            this%obsData(j, iq) = this%missing
          ELSE
            this%obsData(j, iq) = this%obsData(j, ip) / dry_air_gas_const / this%obsData(j, it)
          END IF
        END DO
      END IF

      ! Pressure in hPa:
      DO j = 1, this%numObs
        IF (this%obsData(j, ip) .NE. this%missing) this%obsData(i, ip) = this%obsData(i, ip) / 100.0D0
      END DO

    END BLOCK

    DO i = 1, this%numVars
      WRITE (*, 31) i, this%numVars, TRIM(this%varNames(i)), this%missing, this%invalid, X%mpddGlob%myrank
31    FORMAT('SfcIngest valid obs:', I2, ' of total', I2, ' var ', A6, ' missing/invalid vals:', 2D14.6, ' pc', I2)
      WRITE (*, 32) i, TRIM(this%varNames(i)), &
        MAXVAL(this%obsData(:, i), &
               !MASK=(this%obsData(:, i) .ne. this%missing .and. this%obsData(:, i) .le. this%invalid)), &
               MASK=(this%obsData(:, i) .NE. this%missing)), &
        MINVAL(this%obsData(:, i), &
               MASK=(this%obsData(:, i) .NE. this%missing .AND. this%obsData(:, i) .LE. this%invalid)), &
        COUNT(ABS(this%obsData(:, i)) .LT. this%missing), TRIM(this%obsType), X%mpddGlob%myrank
32    FORMAT('SfcIngest obs range/counts: ', I2, ' ', A6, 2D14.6, ' count:', I6, ' obsType ', A, ' pc', I2)
    END DO

    this%presIdx = 0
    this%tempIdx = 0
    this%qvaporIdx = 0

    DO i = 1, this%numVars
      IF (this%varNames(i) == 'pres') this%presIdx = i
      IF (this%varNames(i) == 'temp') this%tempIdx = i
      IF (this%varNames(i) == 'qvapor') this%qvaporIdx = i
    END DO

    ! Added by Zilong for temporary scaling

    ! IF(this%qvaporIdx /=0 ) this%obsData(:, this%qvaporIdx) = this%obsData(:, this%qvaporIdx)*1000.0D0
    ! IF(this%presIdx /=0 ) this%obsData(:, this%presIdx) = this%obsData(:, this%presIdx)/100.0D0
    IF (this%tempIdx /= 0) this%qcThreshold(this%tempIdx) = 20.0D0
    IF (this%presIdx /= 0) this%qcThreshold(this%presIdx) = 1000.0D0
    ! IF (this%tempIdx /= 0) this%obsData(:, this%tempIdx) = this%obsData(:, this%tempIdx)-1.8D0
    ! IF (this%qvaporIdx /= 0) this%qcThreshold(this%qvaporIdx) = 0.008D0
    ! this%qcThreshold(5) = 10.0D0
    ! this%qcThreshold(6) = 10.0D0

    ! BLOCK
    !   USE RectangleGridUtility_m, ONLY: RectangleGridUtility_t

    !   INTEGER(i_kind), ALLOCATABLE :: idxtgt(:, :)
    !   REAL(r_kind), ALLOCATABLE :: coetgt(:, :)
    !   REAL(r_kind), ALLOCATABLE :: topo(:)
    !   REAL(r_kind), ALLOCATABLE :: cell_cntr(:, :)
    !   REAL(r_kind), ALLOCATABLE :: ObsTopoAt0Level(:)
    !   TYPE(RectangleGridUtility_t) :: regular

    !   IF (X%sg%isBaseProc()) THEN
    !     ALLOCATE (topo(X%sg%num_icell_global))
    !     ALLOCATE (ObsTopoAt0Level(this%numObs))
    !     ALLOCATE (cell_cntr(2, X%sg%num_icell_global))
    !   END IF
    !   CALL X%sg%aggrGridReal1D(X%sg%topo, topo, X%sg%num_icell_global)
    !   CALL X%sg%aggrGridReal2D(X%sg%cell_cntr, cell_cntr, [2, X%sg%num_icell_global])

    !   IF (X%sg%isBaseProc()) THEN
    !     ALLOCATE (idxtgt(4, this%numObs), coetgt(4, this%numObs))
    !     CALL regular%RectangleHorizontalIntp(TRANSPOSE(cell_cntr), X%sg%num_icell_global, &
    !                                          TRANSPOSE(this%olatlon), this%numObs, &
    !                                          idxtgt, coetgt, X%sg%mpddInfo_sg%myrank)
    !     DO i = 1, this%numObs
    !       DO j = 1, 4
    !         IF (idxtgt(j, i) /= 0)
    !         ObsTopoAt0Level(i) = topo(idxtgt(1, i))*coetgt(1, i) + topo(idxtgt(2, i))*coetgt(2, i) + &
    !                              topo(idxtgt(3, i))*coetgt(3, i) + topo(idxtgt(4, i))*coetgt(4, i)
    !       END DO
    !     END DO

    !     DEALLOCATE (idxtgt, coetgt, cell_cntr, ObsTopoAt0Level)
    !   END IF

    !   IF (X%sg%isBaseProc()) THEN
    !     DO j = 1, this%numObs
    !       PRINT *, 'Height: ', this%StNames(j), this%obsHght(j), ObsTopoAt0Level(j)
    !     END DO
    !   END IF
    ! END BLOCK

    ! Add by Zilong to set the Surface obs always on the ground.
    ! this%obsHght = -2.0E7
  END SUBROUTINE sfcIngest

  FUNCTION GetForwardValue(this, X, varnametmp, vIdx, hIdx, tIdx, iv) RESULT(val)
    CLASS(ObsSurface_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: vIdx, hIdx, tIdx
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iv
    REAL(r_kind) :: val
    CHARACTER(LEN=*) :: varnametmp

    ! Local variables:
    INTEGER(i_kind) :: iu, iw
    REAL(r_kind) :: u, v, r

    IF (TRIM(varnametmp) .NE. 'cdir' .AND. &
        TRIM(varnametmp) .NE. 'sdir' .AND. &
        TRIM(varnametmp) .NE. 'wspd') THEN
      val = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(vIdx, hIdx, tIdx)
    ELSE
      ! For wind dir and spd analysis:
      iu = X%getVarIdx('uwnd')
      iw = X%getVarIdx('vwnd')
      IF (iu .GT. 0 .AND. iw .GT. 0) THEN
        u = X%Fields(iu)%DATA(vIdx, hIdx, tIdx)
        v = X%Fields(iw)%DATA(vIdx, hIdx, tIdx)
        r = DSQRT(u * u + v * v)
        SELECT CASE (TRIM(varnametmp))
        CASE ('cdir')
          val = -v / r
        CASE ('sdir')
          val = -u / r
        CASE ('wspd')
          val = r * r
        END SELECT
      ELSE
        WRITE (*, 3)
3       FORMAT('Surface - GetForwardValue: Model does not have uwnd/vwnd, STOP')
        STOP
      END IF
    END IF
  END FUNCTION GetForwardValue

  SUBROUTINE sfcQC(this)
    CLASS(ObsSurface_t) :: this
  END SUBROUTINE sfcQC

END MODULE ObsSurface_m
