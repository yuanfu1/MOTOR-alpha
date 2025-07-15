!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.satellite_utils
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/06/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE Satellite_QC_utils_m
  USE kinds_m
  USE parameters_m
  
CONTAINS

  SUBROUTINE Calc_average(ObsData, ObsData_avg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: ObsData(:)
    REAL(r_kind), INTENT(OUT) :: ObsData_avg
    REAL(r_kind), ALLOCATABLE :: tmp(:)
    INTEGER :: numVars, ichan

    ! Calculate the average of obsData for thinning without bkg
    tmp = PACK(ObsData, ObsData<400.0)
    IF (SIZE(tmp,1) > 0) THEN
      obsData_avg = SUM(tmp)/REAL(SIZE(tmp,1))
    ELSE 
      obsData_avg = missing
    END IF
    DEALLOCATE(tmp)
  END SUBROUTINE Calc_average

  SUBROUTINE Open_NcFile_from_DP(filename, numVars, numObs, nc, olatlon, ObsTime)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(i_kind), INTENT(OUT) :: numVars, numObs
    REAL(r_kind), ALLOCATABLE, INTENT(OUT), OPTIONAL :: olatlon(:,:)
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ObsTime(:)
    TYPE(NcDataset), INTENT(OUT)   :: nc
    type(NcDimension)  :: dim
    TYPE(NcVariable)  :: var
    INTEGER :: Nlat, Nlon, ichan, ivar
    REAL(r_kind), ALLOCATABLE :: tmp(:)

    nc = NcDataset(TRIM(filename), "r")
    dim = nc%getDimension("Nlongitude"); Nlon = dim%getLength()
    dim = nc%getDimension("Nlatitude"); Nlat = dim%getLength()
    dim = nc%getDimension("NOBS"); numObs = dim%getLength()
    dim = nc%getDimension("Nchannels"); numVars = dim%getLength()

    numVars = numVars + 4

    ALLOCATE(tmp(numObs))
    IF (PRESENT(olatlon)) THEN
      ALLOCATE(olatlon(2, numObs))
      var = nc%getVariable("lat"); CALL var%getData(tmp)
      olatlon(1,:) = tmp
      var = nc%getVariable("lon"); CALL var%getData(tmp)
      olatlon(2,:) = tmp
      olatlon = olatlon * degree2radian
    END IF

    IF (PRESENT(obsTime)) THEN
      ALLOCATE(obsTime(numObs))
      var = nc%getVariable("ObsTime"); CALL var%getData(tmp)
      ObsTime = tmp
    END IF
    DEALLOCATE(tmp)

  END SUBROUTINE Open_NcFile_from_DP

  SUBROUTINE Close_NcFile_from_DP(nc)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    TYPE(NcDataset), INTENT(INOUT)   :: nc

    CALL nc%close()

  END SUBROUTINE Close_NcFile_from_DP

  SUBROUTINE Read_OBS_from_DP_EachSlot(nc, numObs, nchans, chan_lists, ObsData, ObsErrs, InvData)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable, hasVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: numObs, nchans
    INTEGER(i_kind), INTENT(IN) :: chan_lists(:)
    REAL(r_kind), INTENT(INOUT) :: ObsData(:,:), ObsErrs(:,:), InvData(:,:)
    TYPE(NcDataset), INTENT(INOUT)   :: nc
    TYPE(NcVariable)  :: var
    INTEGER :: ichan
    CHARACTER(10) :: chan_str, chan_name, Err_name, Inv_name
    REAL(r_kind), ALLOCATABLE :: tmp(:)

    ALLOCATE(tmp(numObs))

    var = nc%getVariable("sat_zenith"); CALL var%getData(tmp)
    PRINT *, 'Check shape ', SHAPE(obsdata), SHAPE(tmp)
    ObsData(:, 1) = tmp
    var = nc%getVariable("sat_azi"); CALL var%getData(tmp)
    ObsData(:, 2) = tmp
    var = nc%getVariable("sol_zenith"); CALL var%getData(tmp)
    ObsData(:, 3) = tmp
    var = nc%getVariable("sol_azi"); CALL var%getData(tmp)
    ObsData(:, 4) = tmp
    ObsData(:, 1:4) = ObsData(:, 1:4) * degree2radian

    ObsErrs(:, 1:4) = missing
    InvData(:, 1:4) = 0.0D0

    DO ichan = 5, nchans ! 1~4 are satellite angles

      IF (chan_lists(ichan-4) < 10) THEN
        WRITE (chan_str, '(i1)') chan_lists(ichan-4)
        chan_name = "cor_TBB00"//TRIM(chan_str)
        Err_name = "ERROR00"//TRIM(chan_str)
        Inv_name = "org_Inv00"//TRIM(chan_str)
      ELSE
        IF (chan_lists(ichan-4) < 100) THEN
          WRITE (chan_str, '(i2)') chan_lists(ichan-4)
          chan_name = "cor_TBB0"//TRIM(chan_str)
          Err_name = "ERROR0"//TRIM(chan_str)
          Inv_name = "org_Inv0"//TRIM(chan_str)
        ELSE
          IF (chan_lists(ichan-4) < 1000) THEN
            WRITE (chan_str, '(i3)') chan_lists(ichan-4)
            chan_name = "cor_TBB"//TRIM(chan_str)
            Err_name = "ERROR"//TRIM(chan_str)
            Inv_name = "org_Inv"//TRIM(chan_str)
          END IF
        END IF
      END IF

      PRINT *, ' Read in ', TRIM(chan_name)
      IF ( hasVariable(nc, TRIM(chan_name)) ) THEN
        var = nc%getVariable(TRIM(chan_name)); CALL var%getData(tmp)
        ObsData(:, ichan) = tmp

        var = nc%getVariable(TRIM(Err_name)); CALL var%getData(tmp)
        ObsErrs(:, ichan) = tmp

        var = nc%getVariable(TRIM(Inv_name)); CALL var%getData(tmp)
        InvData(:, ichan) = tmp
      ELSE
        ObsData(:, ichan) = missing
        ObsErrs(:, ichan) = missing
        InvData(:, ichan) = missing
      END IF
    END DO
    DEALLOCATE(tmp)

  END SUBROUTINE Read_OBS_from_DP_EachSlot

  SUBROUTINE Read_OBS_from_DP(filename, ntimes, chan_lists, nchans, numObs, olatlon, ObsTime, obsData, obsErrs, tb_inv)
    USE mo_netcdf, only: NcDataset
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(i_kind), INTENT(IN) :: ntimes, chan_lists(:)
    INTEGER(i_kind), INTENT(OUT) :: nchans, numObs
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: olatlon(:,:)
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT) :: ObsTime(:)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: ObsData(:,:), ObsErrs(:,:), tb_inv(:,:)
    ! LOCAL VARS
    TYPE(NcDataset)    :: nc
    INTEGER :: numObsEachSlot
    REAL(r_kind), ALLOCATABLE    :: ObsDataEachSlot(:,:)
    INTEGER(i_kind), ALLOCATABLE :: ObsTimeEachSlot(:) 
    REAL(r_kind), ALLOCATABLE    :: olatlonEachSlot(:,:)
    REAL(r_kind), ALLOCATABLE    :: ObsErrsEachSlot(:,:)
    REAL(r_kind), ALLOCATABLE    :: INVDataEachSlot(:,:)
    INTEGER :: it, tt1, tt2, IdxStart, IdxEnd
    REAL :: t1, t2
    CHARACTER(LEN=3) :: tSlotsStr
    CHARACTER(LEN=200) :: filenameEachSlot
    LOGICAL :: file_exists

    numObs = 0
    DO it = 0, ntimes-1 
      WRITE(tSlotsSTR, "(I1)") it+1
      filenameEachSlot = TRIM(filename) // '_'// TRIM(tSlotsStr) // '.nc'
      
      INQUIRE(file=TRIM(filenameEachSlot), exist=file_exists)
      IF (.NOT. file_exists) THEN
        numObsEachSlot = 0
      ELSE
        CALL Open_NcFile_from_DP(TRIM(filenameEachSlot), nchans, numObsEachSlot, nc)
        CALL Close_NcFile_from_DP(nc)
      END IF
      numObs = numObs + numObsEachSlot
    END DO

    IF (numObs < 1) nchans = 4  ! 4 is the number of satellite angles, nchans is actually 0

    ALLOCATE(olatlon(2, numObs))
    ALLOCATE(obsTime(numObs))
    ALLOCATE(obsData(numObs, nchans))
    ALLOCATE(obsErrs(numObs, nchans))
    ALLOCATE(tb_inv(numObs, nchans))
    IdxStart = 1
    IdxEnd = 0

    IF (numObs < 1) RETURN

    DO it = 0, ntimes-1 
      WRITE(tSlotsSTR, "(I1)") it+1
      filenameEachSlot = TRIM(filename) // '_'// TRIM(tSlotsStr) // '.nc'
      PRINT *, 'filenameEachSlot = ', TRIM(filenameEachSlot)
      CALL Open_NcFile_from_DP(TRIM(filenameEachSlot), nchans, numObsEachSlot, nc, olatlonEachSlot, ObsTimeEachSlot)
      PRINT *, 'Check file infos from DP Each slot: ', it, nchans, numObsEachSlot
      IF (numObsEachSlot .EQ. 0) CYCLE
      ALLOCATE(ObsDataEachSlot(numObsEachSlot, nchans))
      ALLOCATE(ObsErrsEachSlot(numObsEachSlot, nchans))
      ALLOCATE(INVDataEachSlot(numObsEachSlot, nchans))
      CALL Read_OBS_from_DP_EachSlot(nc, numObsEachSlot, nchans, chan_lists, ObsDataEachSlot, ObsErrsEachSlot, INVDataEachSlot)
      CALL Close_NcFile_from_DP(nc)
      IdxEnd = IdxEnd + numObsEachSlot
      olatlon(:, IdxStart:IdxEnd) = olatlonEachSlot
      obsTime(IdxStart:IdxEnd) = ObsTimeEachSlot
      obsData(IdxStart:IdxEnd, :) = ObsDataEachSlot
      obsErrs(IdxStart:IdxEnd, :) = ObsErrsEachSlot
      tb_inv(IdxStart:IdxEnd, :)  = INVDataEachSlot
      IdxStart = IdxStart + numObsEachSlot
      ! PRINT *, 'Read_OBS_from_DP_EachSlot: ', maxval(ObsDataEachSlot(:,6)), minval(ObsDataEachSlot(:,6))
      DEALLOCATE(olatlonEachSlot, ObsTimeEachSlot, ObsDataEachSlot, ObsErrsEachSlot, INVDataEachSlot)
    END DO
    ! PRINT *, 'Read_OBS_from_DP: ', SHAPE(obsData), MAXVAL(tb_inv(:,6)), MINVAL(tb_inv(:,6))

  END SUBROUTINE Read_OBS_from_DP

END MODULE Satellite_QC_utils_m
