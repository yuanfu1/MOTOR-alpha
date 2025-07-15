!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.DP_utils
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
MODULE DP_utils_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE mpObs_m, ONLY: mpObs_t
  
CONTAINS

  SUBROUTINE MPI_GATHER_Interface(sg, array_in, array_out)
  
    IMPLICIT NONE 
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: array_in(:)
    REAL(r_kind), INTENT(INOUT) :: array_out(:)
    INTEGER(i_kind) :: nobs
    INTEGER(i_kind) :: i
    ! For mpi
    INTEGER(i_kind), ALLOCATABLE :: ncount_group(:), disp_group(:)
    REAL(r_kind), ALLOCATABLE :: idxArray(:)
    INTEGER(i_kind) :: icount

    INCLUDE "mpif.h"
    nobs = UBOUND(array_in, 1)
    ! PRINT *, 'MPI_GATHER_Interface: nobs = ', nobs

    ALLOCATE (ncount_group(sg%mpddInfo_sg%nproc))
    ALLOCATE (disp_group(sg%mpddInfo_sg%nproc))

    ncount_group = 0
    disp_group = 0

    ! PRINT *, 'nproc ', sg%mpddInfo_sg%nProc
    ! 将不同进程上的观测个数汇总到主线程 0
    CALL MPI_GATHER(nobs, 1, MPI_INTEGER4, &
                    ncount_group, 1, MPI_INTEGER4, 0, &
                    sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    ! PRINT *, 'check myrank MPI_GATHER_Interface: ', sg%mpddInfo_sg%myrank,sg%isBaseProc(), SHAPE(array_in)
    IF (sg%isBaseProc()) THEN
      ! PRINT *, "ncount_group is ", ncount_group
      ALLOCATE (idxArray(SUM(ncount_group)))
      idxArray = ZERO
      disp_group = 0
      FORALL (i=2:sg%mpddInfo_sg%nproc) disp_group(i) = SUM(ncount_group(1:i - 1)) !Gatherv 时，每个进程所需的位移
      ! PRINT *, 'check out MPI_GATHER_Interface:', sg%mpddInfo_sg%myrank, 'disp_group=', disp_group, 'ncount_group=', ncount_group
    END IF

    CALL sg%mpddInfo_sg%bcast(ncount_group)
    CALL sg%mpddInfo_sg%bcast(disp_group)

    ! 将所有进程的innovation值汇总，存在idxArray中
    CALL MPI_GATHERV(array_in, SIZE(array_in), MPI_DOUBLE_PRECISION, &
                    idxArray, ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                    sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    IF (sg%isActiveProc() .AND. (.NOT. sg%isBaseProc())) &
      ALLOCATE (idxArray(SUM(ncount_group)))
    CALL sg%mpddInfo_sg%bcast(idxArray)
    array_out = idxArray

    IF (ALLOCATED (idxArray)) DEALLOCATE(idxArray)

    DEALLOCATE(ncount_group, disp_group)
    ! PRINT *, 'MPI_GATHER_Interface is successfully called'

  END SUBROUTINE MPI_GATHER_Interface

  SUBROUTINE MPI_GATHER_Interface_Int(sg, array_in, array_out)
  
    IMPLICIT NONE 
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(IN) :: array_in(:)
    INTEGER(i_kind), INTENT(INOUT) :: array_out(:)
    INTEGER(i_kind) :: nobs
    INTEGER(i_kind) :: i
    ! For mpi
    INTEGER(i_kind), ALLOCATABLE :: ncount_group(:), disp_group(:)
    INTEGER(i_kind), ALLOCATABLE :: idxArray(:)
    INTEGER(i_kind) :: icount

    INCLUDE "mpif.h"
    nobs = UBOUND(array_in, 1)
    ! PRINT *, 'MPI_GATHER_Interface_Int: nobs = ', nobs, 'isBaseProc ', sg%isBaseProc()

    ALLOCATE (ncount_group(sg%mpddInfo_sg%nproc))
    ALLOCATE (disp_group(sg%mpddInfo_sg%nproc))

    ncount_group = 0
    disp_group = 0

    ! PRINT *, 'nproc ', sg%mpddInfo_sg%nProc
    ! 将不同进程上的观测个数汇总到主线程 0
    ! print *, 'before MPI_GATHER ', sg%isBaseProc()
    CALL MPI_GATHER(nobs, 1, MPI_INTEGER4, &
                    ncount_group, 1, MPI_INTEGER4, 0, &
                    sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    ! PRINT* , 'before isBaseProc ',sg%isBaseProc()
    IF (sg%isBaseProc()) THEN
      ! PRINT *, "ncount_group is ", ncount_group
      ALLOCATE (idxArray(SUM(ncount_group)))
      idxArray = ZERO
      disp_group = 0
      FORALL (i=2:sg%mpddInfo_sg%nproc) disp_group(i) = SUM(ncount_group(1:i - 1)) !Gatherv 时，每个进程所需的位移
      ! PRINT *, 'check out MPI_GATHER_Interface_Int:', sg%mpddInfo_sg%myrank, 'disp_group=', disp_group, 'ncount_group=', ncount_group
    END IF

    CALL sg%mpddInfo_sg%bcast(ncount_group)
    CALL sg%mpddInfo_sg%bcast(disp_group)

    ! 将所有进程的innovation值汇总，存在idxArray中
    CALL MPI_GATHERV(array_in, SIZE(array_in), MPI_INTEGER4, &
                    idxArray, ncount_group, disp_group, MPI_INTEGER4, 0, &
                    sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    IF (sg%isActiveProc() .AND. (.NOT. sg%isBaseProc())) &
      ALLOCATE (idxArray(SUM(ncount_group)))
    CALL sg%mpddInfo_sg%bcast(idxArray)
    array_out = idxArray

    IF (ALLOCATED (idxArray)) DEALLOCATE(idxArray)

    DEALLOCATE(ncount_group, disp_group)
    ! PRINT *, 'MPI_GATHER_Interface_Int is successfully called'

  END SUBROUTINE MPI_GATHER_Interface_Int

  SUBROUTINE Write_OBS_BC(numObs, olatlon, ObsData, ObsData_nobc, BkgAtRawobs, filename)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable, NcGroup
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    ! TYPE(SingleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(IN) :: numObs
    REAL(r_kind), INTENT(IN) :: olatlon(:,:)
    CHARACTER(*), INTENT(IN) :: filename
    REAL(r_kind), INTENT(IN) :: ObsData(:,:), ObsData_nobc(:,:), BkgAtRawobs(:,:)
    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1, dim2, dim3
    TYPE(NcVariable)  :: var
    TYPE(NcGroup)     :: grp
    INTEGER :: Nlat, Nlon, Nchans
    REAL(r_kind), ALLOCATABLE    :: out_array(:), tmp(:)

    ! IF (.not. sg%isBaseProc()) RETURN

    nc = NcDataset(TRIM(filename), "w")
    Nlat = numObs !this%olatlon
    Nlon = numObs

    dim1 = nc%setDimension("Nlongitude", Nlon)
    dim2 = nc%setDimension("Nlatitude", Nlat)
    dim3 = nc%setDimension("NOBS", numObs)

    ALLOCATE(out_array(numObs))

    ! Reverse latitude for plotting
    var = nc%setVariable("lat", "f32", (/dim2/))
    call var%setFillValue(9999.0)
    WHERE (olatlon(1,:) .NE. missing)  
      out_array = olatlon(1,:) * radian2degree
    ELSEWHERE 
      out_array = 9999.0
    END WHERE
    call var%setData(out_array(numObs:1:-1))

    var = nc%setVariable("lon", "f32", (/dim1/))
    call var%setFillValue(9999.0)
    WHERE ( olatlon(2,:) .NE. missing )  
      out_array = olatlon(2,:) * radian2degree
    ELSEWHERE
      out_array = 9999.0
    END WHERE
    call var%setData(out_array(numObs:1:-1))
  
    var = nc%setVariable("tbb01_bak", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(BkgAtRawobs(numObs:1:-1, 1))

    var = nc%setVariable("tbb02_bak", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(BkgAtRawobs(numObs:1:-1, 2))

    var = nc%setVariable("tbb01_nobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(ObsData_nobc(numObs:1:-1, 1))

    var = nc%setVariable("tbb02_nobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(ObsData_nobc(numObs:1:-1, 2))

    var = nc%setVariable("tbb01_bc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(obsData(numObs:1:-1, 1))

    var = nc%setVariable("tbb02_bc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    call var%setData(obsData(numObs:1:-1, 2))

    var = nc%setVariable("tbb01_omb_bc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = obsData(numObs:1:-1, 1)-BkgAtRawobs(numObs:1:-1, 1)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    var = nc%setVariable("tbb02_omb_bc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = obsData(numObs:1:-1, 2)-BkgAtRawobs(numObs:1:-1, 2)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    var = nc%setVariable("tbb01_omb_nobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = ObsData_nobc(numObs:1:-1, 1)-BkgAtRawobs(numObs:1:-1, 1)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    var = nc%setVariable("tbb02_omb_nobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = ObsData_nobc(numObs:1:-1, 2)-BkgAtRawobs(numObs:1:-1, 2)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    var = nc%setVariable("tbb01_obs_bcmnobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = obsData(numObs:1:-1, 1)-ObsData_nobc(numObs:1:-1, 1)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    var = nc%setVariable("tbb02_obs_bcmnobc", "f32", (/dim3/))
    call var%setFillValue(9999.0)
    ALLOCATE(tmp(numObs))
    tmp = obsData(numObs:1:-1, 2)-ObsData_nobc(numObs:1:-1, 2)
    call var%setData(tmp)
    DEALLOCATE(tmp)

    ! close the dataset
    call nc%close()
    DEALLOCATE(out_array)

  END SUBROUTINE Write_OBS_BC

  SUBROUTINE Open_NcFile_for_thinning(numVars, numObs, olatlon, Angles, ObsTime, filename, nc, dim3)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: numVars, numObs
    CHARACTER(*), INTENT(IN) :: filename
    REAL(r_kind), INTENT(IN) :: olatlon(:,:)
    REAL(r_kind), INTENT(IN) :: Angles(:,:)
    INTEGER(i_kind), INTENT(IN) :: ObsTime(:)
    TYPE(NcDataset), INTENT(OUT)   :: nc
    TYPE(NcVariable)  :: var
    TYPE(NcDimension) :: dim1, dim2, dim4
    TYPE(NcDimension), INTENT(OUT) :: dim3
    INTEGER :: Nlat, Nlon, ichan, ivar
    CHARACTER(10) :: chan_str, chan_name
    REAL(r_kind), ALLOCATABLE :: out_array(:)

    nc = NcDataset(TRIM(filename), "w")
    Nlat = numObs !this%olatlon
    Nlon = numObs

    dim1 = nc%setDimension("Nlongitude", Nlon)
    dim2 = nc%setDimension("Nlatitude", Nlat)
    dim3 = nc%setDimension("NOBS", numObs)
    dim4 = nc%setDimension("Nchannels", numVars-4)

    ALLOCATE(out_array(numObs))

    ! Reverse latitude for plotting
    var = nc%setVariable("lat", "f32", (/dim2/))
    ! CALL var%setFillValue(9999.0)
    WHERE (olatlon(1,:) .LE. missing)  
      out_array = olatlon(1,:) * radian2degree
    ELSEWHERE 
      out_array = missing
    END WHERE
    CALL var%setData(out_array)

    var = nc%setVariable("lon", "f32", (/dim1/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( olatlon(2,:) .LE. missing )  
      out_array = olatlon(2,:) * radian2degree
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array)    

    var = nc%setVariable("ObsTime", "i32", (/dim3/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( ObsTime .NE. missing )  
      out_array = ObsTime
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array)   
    ! print *, 'check obstime = ', maxval(out_array), minval(out_array) 

    var = nc%setVariable("sat_zenith", "f32", (/dim3/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( Angles(:,1) .LE. missing )  
      out_array = Angles(:,1) * radian2degree
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array)   

    var = nc%setVariable("sat_azi", "f32", (/dim3/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( Angles(:,2) .LE. missing )  
      out_array = Angles(:,2) * radian2degree
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array) 

    var = nc%setVariable("sol_zenith", "f32", (/dim3/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( Angles(:,3) .LE. missing )  
      out_array = Angles(:,3) * radian2degree
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array) 

    var = nc%setVariable("sol_azi", "f32", (/dim3/))
    ! CALL var%setFillValue(9999.0)
    WHERE ( Angles(:,4) .LE. missing )  
      out_array = Angles(:,4) * radian2degree
    ELSEWHERE
      out_array = missing
    END WHERE
    CALL var%setData(out_array) 

    DEALLOCATE(out_array)

  END SUBROUTINE Open_NcFile_for_thinning

  SUBROUTINE Close_NcFile_for_thinning(nc)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    TYPE(NcDataset), INTENT(INOUT)   :: nc

    CALL nc%close()

  END SUBROUTINE Close_NcFile_for_thinning

  SUBROUTINE Write_OBS_for_thinning(nc, dim3, numObs, nchans, chan_lists, ObsData, &
                                    ObsData_raw, ObsErrs, InvData, ClrData, ca_mean)
    USE mo_netcdf, only: NcDataset, NcDimension, NcVariable
    ! Refer to https://github.com/schaefed/mo_netcdf for usage
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: numObs, nchans
    INTEGER(i_kind), INTENT(IN) :: chan_lists(:)
    REAL(r_kind), INTENT(IN) :: ObsData(:,:), ObsData_raw(:,:), ObsErrs(:,:), InvData(:,:), ClrData(:,:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: ca_mean(:,:)
    TYPE(NcDataset), INTENT(INOUT)   :: nc
    TYPE(NcDimension), INTENT(IN) :: dim3
    TYPE(NcVariable)  :: var
    INTEGER :: ichan
    CHARACTER(10) :: chan_str, chan_name, raw_name, Err_name, ca_name, Inv_name, Clr_name

    ! PRINT *, 'Write_OBS_for_thinning: ', nchans, ' channels: ', chan_lists
    DO ichan = 1, nchans

      IF (chan_lists(ichan) < 10) THEN
        WRITE (chan_str, '(i1)') chan_lists(ichan)
        chan_name = "cor_TBB00"//TRIM(chan_str)
        raw_name = "org_TBB00"//TRIM(chan_str)
        Err_name = "ERROR00"//TRIM(chan_str)
        ca_name = "ca_mean00"//TRIM(chan_str)
        Inv_name = "org_Inv00"//TRIM(chan_str)
        Clr_name = "Clr00"//TRIM(chan_str)
      ELSE
        IF (chan_lists(ichan) < 100) THEN
          WRITE (chan_str, '(i2)') chan_lists(ichan)
          chan_name = "cor_TBB0"//TRIM(chan_str)
          raw_name = "org_TBB0"//TRIM(chan_str)
          Err_name = "ERROR0"//TRIM(chan_str)
          ca_name = "ca_mean0"//TRIM(chan_str)
          Inv_name = "org_Inv0"//TRIM(chan_str)
          Clr_name = "Clr0"//TRIM(chan_str)
        ELSE
          IF (chan_lists(ichan) < 1000) THEN
            WRITE (chan_str, '(i3)') chan_lists(ichan)
            chan_name = "cor_TBB"//TRIM(chan_str)
            raw_name = "org_TBB"//TRIM(chan_str)
            Err_name = "ERROR"//TRIM(chan_str)
            ca_name = "ca_mean"//TRIM(chan_str)
            Inv_name = "org_Inv"//TRIM(chan_str)
            Clr_name = "Clr"//TRIM(chan_str)
          END IF
        END IF
      END IF

      IF ( .NOT. ALL(ABS(ObsData(:, ichan)) .GT. 500)) THEN
        PRINT *, ' Write out ', TRIM(chan_name)
        var = nc%setVariable(TRIM(chan_name), "f32", (/dim3/))
        ! CALL var%setFillValue(9999.0)
        CALL var%setData(ObsData(:, ichan))

        var = nc%setVariable(TRIM(raw_name), "f32", (/dim3/))
        ! CALL var%setFillValue(9999.0)
        CALL var%setData(ObsData_raw(:, ichan))

        var = nc%setVariable(TRIM(Err_name), "f32", (/dim3/))
        ! CALL var%setFillValue(9999.0)
        CALL var%setData(ObsErrs(:, ichan))

        IF (PRESENT(ca_mean)) THEN
          var = nc%setVariable(TRIM(ca_name), "f32", (/dim3/))
          ! CALL var%setFillValue(9999.0)
          CALL var%setData(ca_mean(:, ichan))
        END IF

        var = nc%setVariable(TRIM(Inv_name), "f32", (/dim3/))
        ! CALL var%setFillValue(9999.0)
        CALL var%setData(InvData(:, ichan))

        var = nc%setVariable(TRIM(Clr_name), "f32", (/dim3/))
        ! CALL var%setFillValue(9999.0)
        CALL var%setData(ClrData(:, ichan))
      END IF
    END DO

  END SUBROUTINE Write_OBS_for_thinning

    ! SUBROUTINE Select_active_channels(obsvalue,nchans,ifuse)
    !   ! This subroutine should be CALLed at the very beginning
    !   ! avoid unnecessary QC procedures
    !   IMPLICIT NONE
    !   REAL(r_kind),    INTENT(INOUT)  :: obsvalue(:, :)
    !   INTEGER(i_kind), INTENT(IN)     :: nchans
    !   INTEGER(i_kind), INTENT(IN)     :: ifuse(:)
    !   !local arguments:
    !   INTEGER(i_kind)  :: ichan
    !   INTEGER(i_kind)  :: ichan_sel
  
    !   IF (SIZE(obsvalue) .LE. 0) RETURN

    !   ichan_sel=0
    !   DO ichan = 1, nchans
    !     IF (ifuse(ichan) .LT. 1) THEN
    !       obsvalue(:, ichan) = missing
    !       ichan_sel=ichan_sel+1
    !     ENDIF
    !   END DO
    !   ! PRINT *, 'check blacklist channels: ', ichan_sel
  
    ! END SUBROUTINE Select_active_channels

    SUBROUTINE Select_Clear(inst_name,obsData,cloud_flag, ichn)
      ! This subroutine should be CALLed at the very beginning
      ! avoid unnecessary QC procedures
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN)     :: inst_name
      REAL(r_kind),     INTENT(INOUT)  :: obsData(:)
      REAL(r_kind),     INTENT(IN)     :: cloud_flag(:)
      INTEGER(i_kind), INTENT(IN)      :: ichn
      !local arguments:
      REAL(r_kind)    :: r
      INTEGER(i_kind) :: iobs,icld
  
      icld=0
      DO iobs = 1, SIZE(obsData,1)
        SELECT CASE(TRIM(inst_name))
        CASE('agri')
          IF (cloud_flag(iobs) .NE. 3.0) THEN
            obsData(iobs) = missing
            icld=icld+1
          ENDIF
        CASE('giirs')
          IF (cloud_flag(iobs) .GT. 0.0) THEN
            obsData(iobs) = missing
            icld=icld+1
          ENDIF
        END SELECT
      END DO
      r=real(icld)/real(SIZE(obsData,1))
      print *,'Select_Clear: cld ratio is ',r,' ', icld,' out of ', SIZE(obsData,1), ' are clouds for channel ', ichn
  
  END SUBROUTINE Select_Clear
  
  SUBROUTINE QC_BoundaryCheck(obsvalue,num_obs,iscanpos)
    IMPLICIT NONE
    REAL(r_kind),    INTENT(INOUT)  :: obsvalue(:, :)
    INTEGER(i_kind), INTENT(IN)     :: num_obs
    INTEGER(i_kind), INTENT(IN)     :: iscanpos(:)
    !local arguments:
    INTEGER(i_kind) :: iobs, iqc
    
    iqc = 0
    DO iobs = 1, num_obs
        IF (( MOD(iscanpos(iobs),32) == 0 ) .OR. ( MOD(iscanpos(iobs),32) <= 5 ) .OR. ( MOD(iscanpos(iobs),32) >= 28 )) THEN !<detector>
          obsvalue(iobs, :) = missing
          iqc = iqc + 1
          ! print *,iobs,'iscanpos(iobs):',iscanpos(iobs)
        ENDIF  
        
    END DO 
    PRINT *, iqc, 'out of ', num_obs, ' observations are removed by QC_BoundaryCheck'

  END SUBROUTINE QC_BoundaryCheck

  SUBROUTINE Select_InDomain(obsvalue,num_obs,lat,lon,minLat, maxLat, minLon, maxLon)
    ! This subroutine should be CALLed at the very beginning
    ! avoid unnecessary QC procedures
    IMPLICIT NONE
    REAL(r_kind),    INTENT(INOUT)  :: obsvalue(:, :)
    INTEGER(i_kind), INTENT(IN)     :: num_obs
    REAL(r_kind),    INTENT(IN)     :: lat(:),lon(:)
    REAL(r_kind),    INTENT(IN)     :: minLat, maxLat, minLon, maxLon
    !local arguments:
    INTEGER(i_kind) :: iobs,iqc

    IF (SIZE(obsvalue) .LE. 0) RETURN

    iqc=0
    DO iobs = 1, num_obs
      
      IF ((lat(iobs) .LT. minLat) .OR. (lat(iobs) .GT. maxLat) .OR. &      
            (lon(iobs) .LT. minLon) .OR. (lon(iobs) .GT. maxLon)) THEN      
              obsvalue(iobs, :) = missing
              iqc=iqc+1
      ENDIF
      
    END DO
  
    PRINT *,iqc, 'out of ', num_obs,' are removed by Select_InDomain lat:',minLat/degree2radian,maxlat/degree2radian,  &
    ' lon:',minLon/degree2radian,maxLon/degree2radian

  END SUBROUTINE Select_InDomain

  SUBROUTINE Select_InTimeWin(obsvalue,num_obs,obstime,minTime, maxTime)
    ! This subroutine should be CALLed at the very beginning
    ! avoid unnecessary QC procedures
    IMPLICIT NONE
    REAL(r_kind),    INTENT(INOUT)  :: obsvalue(:, :)
    INTEGER(i_kind), INTENT(IN)     :: num_obs
    INTEGER(i_kind), INTENT(IN)     :: obstime(:)
    REAL(r_kind),    INTENT(IN)     :: minTime, maxTime
    !local arguments:
    INTEGER(i_kind) :: iobs,iqc

    IF (SIZE(obsvalue) .LE. 0) RETURN
    
    iqc=0
    DO iobs = 1, num_obs
      
      IF((obstime(iobs) .LT. minTime) .OR. (obstime(iobs) .GT. maxTime)) THEN  
        obsvalue(iobs, :) = missing
        iqc=iqc+1
      ENDIF
      
    ENDDO
    PRINT *,iqc, 'out of ', num_obs, ' are removed by Select_InTimeWin ',obstime(1),minTime,maxTime
    
  END SUBROUTINE Select_InTimeWin

  FUNCTION COUNT_MISS_REAL_1D(data) RESULT(num_miss)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: data(:)
    INTEGER :: num_miss

    num_miss = COUNT(ABS(data - missing) .LT. 1.0)

  END FUNCTION COUNT_MISS_REAL_1D

  FUNCTION ISMISS_REAL(data) RESULT(ifmiss)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: data
    LOGICAL :: ifmiss

    IF (ABS(data - missing) .LT. 1.0) THEN
      ifmiss = .TRUE.
    ELSE
      ifmiss = .FALSE.
    END IF

  END FUNCTION ISMISS_REAL

  SUBROUTINE Select_active_channels(obsData, ifuse, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    INTEGER(i_kind), INTENT(IN) :: ifuse
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn

    IF (ifuse .LT. 1.0) THEN
      qcflag = 'Select_active_channels'
      obsData = missing
      PRINT *, SIZE(obsData,1), 'out of ', SIZE(obsData,1), ' observations are removed by Select_active_channels for Channnel:', ichn
    ELSE 
      PRINT *, ichn, ' channel is active after Select_active_channels'
    END IF

  END SUBROUTINE Select_active_channels

  SUBROUTINE Select_sea_channels(obsData, oversea, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    INTEGER(i_kind), INTENT(IN) :: oversea
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn

    IF ( oversea .LT. 1.0) THEN
      qcflag = 'Select_sea_channels'
      obsData = missing
      PRINT *, SIZE(obsData,1), 'out of ', SIZE(obsData,1), ' observations are removed by Select_sea_channels for Channnel:', ichn
    ELSE 
      PRINT *, ichn, ' channel is active after Select_sea_channels'
    END IF

  END SUBROUTINE Select_sea_channels

  SUBROUTINE QC_InitialCheck(obsData, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs

    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF (ABS(obsData(iobs) - missing) .LE. 1.0) THEN
        qcflag(iobs) = 'QC_InitialCheck'
        qcnum = qcnum + 1
      END IF
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_InitialCheck for Channnel:', ichn
  END SUBROUTINE QC_InitialCheck

  SUBROUTINE QC_AbsoluteCheck(obsData, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF (obsData(iobs) .GT. 350.0 .OR. obsData(iobs) .LT. 180.0) THEN
        qcflag(iobs) = 'QC_AbsoluteCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_AbsoluteCheck for Channnel:', ichn
  END SUBROUTINE QC_AbsoluteCheck

  SUBROUTINE QC_GeoCheck(obsData, lat, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), lat(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF (lat(iobs) .LT. 0.0) THEN
        qcflag(iobs) = 'QC_GeoCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_GeoCheck for Channnel:', ichn
  END SUBROUTINE QC_GeoCheck

  SUBROUTINE QC_MWCLWCheck(obsData, lat, lon, CldIndex, qcflag, ichn)
    !  Any radiances corresponding to C>0.5 are considered to be 
    !  affected by thick clouds and are therefore excluded.
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), lat(:), lon(:)
    REAL(r_kind), INTENT(IN) :: CldIndex(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind)    :: iobs
    INTEGER :: file_unit, ios
    CHARACTER(len=100) :: filename
    
  !     !!!!!!!!!!!!!!!!!!!!!!!!!!!! Set the filename
  !   filename = 'MWCLWCheck.txt'
  ! ! Open the file
  !   file_unit = 10
  !   OPEN(UNIT=file_unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
  !   IF (ios /= 0) THEN
  !     PRINT *, 'Error opening file: ', ios
  !     RETURN
  !   END IF
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
          !!!!!!!!!!!!!!!!!!!!!!!Write clwp to file
      ! WRITE(file_unit, '(I5, F10.5, F10.5, F10.5)') iobs, lat(iobs) * radian2degree, lon(iobs) * radian2degree, CldIndex(iobs)
      IF ( CldIndex(iobs) > 0.5 ) THEN
        qcflag(iobs) = 'QC_MWCLWCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 
        ! !!!!!!!!!!!!!!!!!!!!!!!!!Close the file
    CLOSE(file_unit)
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWCLWCheck for Channnel:', ichn

  END SUBROUTINE QC_MWCLWCheck

!!!new
  SUBROUTINE QC_MWTS3LWPestimate(obsData, BTo_238, BTo_314,lat,lon, zenith,landseamask, qcflag, ichn)
    !  Any radiances corresponding to clwp>0.15 are considered to be 
    ! Refer to DAVID et al.(2021)
    ! RawMWTS3.
    ! landseamask: 1 land, 2 continentalwater, 0 sea, 5 boundary

    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    REAL(r_kind), INTENT(IN) :: zenith(:)
    REAL(r_kind), INTENT(IN) :: BTo_238(:), BTo_314(:), landseamask(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    REAL(r_kind), INTENT(IN) :: lat(:), lon(:)
    REAL(r_kind) :: clwp
    REAL(r_kind) :: SI
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs

    INTEGER :: file_unit, ios
    CHARACTER(len=100) :: filename
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! Set the filename
  !   filename = 'landmask.txt'
  ! ! Open the file
  !   file_unit = 10
  !   OPEN(UNIT=file_unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
  !   IF (ios /= 0) THEN
  !     PRINT *, 'Error opening file: ', ios
  !     RETURN
  !   END IF
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      clwp = 8.24 - cos(zenith(iobs)) * (2.539 - 1.744 * cos(zenith(iobs)) ) + 0.754 * log((285 - BTo_238(iobs))) - 2.265 * log((285 - BTo_314(iobs)))
      SI =  BTo_238(iobs) -  BTo_314(iobs)
      !!!!!!!!!!!!!!!!!!!!!!!!Write clwp to file
      ! WRITE(file_unit, '(I5, F10.5)') iobs, landseamask(iobs)  !!!landamsk
      ! WRITE(file_unit, '(I5, F10.5, F10.5, F10.5, F10.5)') iobs, lat(iobs)*radian2degree, lon(iobs)*radian2degree, clwp,SI
      ! ! PRINT *, "lat=",lat*radian2degree
      IF (  (ABS(landseamask(iobs) -1) .LT. 0.1 .AND. SI .GE. 3.0)  .OR. &
      (ABS(landseamask(iobs) -0) .LT. 0.1 .AND. clwp .LT. 0.35) ) THEN 

        qcflag(iobs) = 'QC_MWTS3LWPestimate'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 

    PRINT *,  'obsdata = ', maxval(ObsData), minval(ObsData)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!Close the file
    ! CLOSE(file_unit)
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWTS3LWPestimate for Channnel:', ichn

  END SUBROUTINE QC_MWTS3LWPestimate

  SUBROUTINE QC_MWCLWPguessCheck(obsData, pres, qcloud, qcflag, ichn)
    !  Any radiances corresponding to clwp>0.2 are considered to be 
    !  affected by thick clouds and are therefore excluded.
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    REAL(r_kind), INTENT(IN) :: pres(:,:), qcloud(:,:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    REAL(r_kind) :: clwp
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      CALL Calc_MWclwp_guess(pres(iobs, :), qcloud(iobs, :), clwp)
      IF ( clwp > 0.2 ) THEN
        qcflag(iobs) = 'QC_MWCLWPguessCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *,  'obsdata = ', maxval(obsData), minval(obsData)
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWCLWPguessCheck for Channnel:', ichn

  END SUBROUTINE QC_MWCLWPguessCheck

  SUBROUTINE QC_MWScatteringIndex(obsData, SI, qcflag, ichn)
    ! Any radiances with SI>=1 are excluded
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    REAL(r_kind), INTENT(IN) :: SI(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( SI(iobs) .GE. 4.0 ) THEN
        qcflag(iobs) = 'QC_MWScatteringIndex'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *,  'obsData =  ', maxval(obsData), minval(obsData)
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWScatteringIndex for Channnel:', ichn

  END SUBROUTINE QC_MWScatteringIndex

!  SUBROUTINE QC_MWHS2_ScatteringIndex(obsData, BTo_89, BTo_166, qcflag, ichn)
!     ! Any radiances with SI>=1 are excluded
!     IMPLICIT NONE
!     REAL(r_kind), INTENT(INOUT) :: obsData(:)
!     REAL(r_kind), INTENT(IN) :: BTo_89(:), BTo_166(:)
!     REAL(r_kind), INTENT(IN) :: zenith(:)
!     REAL(r_kind) :: SI, IWP , omega
!     CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
!     INTEGER(i_kind), INTENT(IN)     :: ichn
!     INTEGER(i_kind)    :: qcnum
!     INTEGER(i_kind) :: iobs
    
!     qcnum = 0
!     DO iobs = 1, SIZE(obsData,1)
!       SI = BTo_89(iobs) - BTo_166(iobs)
!       omega = BTo_89(iobs) / BTo_166(iobs)
!       IWP = cos(zenith(iobs)) * omega * SI

!       IF ( IWP .GE. 0.02 ) THEN   !差异值 IWP 大于或等于 0.02
!         qcflag(iobs) = 'QC_MWHS2_ScatteringIndex_IWP' !质控标志 qcflag(iobs) 设置为 'QC_MWHS2_ScatteringIndex'，表示这个观测点的数据受到了散射的影响。
!         qcnum = qcnum + 1
!         obsData(iobs) = missing
!       END IF  
!     END DO 
!     PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWHS2_ScatteringIndex_IWP for Channnel:', ichn

!   END SUBROUTINE QC_MWHS2_ScatteringIndex

  SUBROUTINE Calc_MWCLWCheck(BTo_238, BTm_238, BTo_314, BTm_314, BTo_528, BTm_528, zenith, sfctype, CldIndex)
    ! Refer to Zhu et al.(2016)
    ! This subroutine should be called before obsThinning.
    ! RawMWTS2 has zenith and landmask data.
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: BTo_238(:), BTm_238(:), BTo_314(:), BTm_314(:), BTo_528(:), BTm_528(:)
    REAL(r_kind), INTENT(IN) :: zenith(:), sfctype(:)
    REAL(r_kind), INTENT(INOUT) :: CldIndex(:)
    REAL(r_kind) :: c0, c1 = 0.754, c2 = -2.265 
    REAL(r_kind) :: w1 = 1.0D0/0.3D0, w2 = 1.0D0/1.8D0
    REAL(r_kind) :: CLWo, CLWm, deltaCLW
    INTEGER(i_kind) :: iobs
    LOGICAL :: isflg

    ! cloud-free FOVs and optically thin clouds are assimilated.
    ! CldIndex is to exclude those optically thick clouds.
    ! CldIndex uses CLW and 52.8GHz channel.
    ! Typically, for MW channels, only those FOVs over ocean surface are assimilated.
    ! PRINT *, 'Calc_MWCLWCheck'
    ! PRINT *, SIZE(BTo_528,1)
    DO iobs = 1, SIZE(BTo_528,1)   
      isflg = (BTo_238(iobs) < 999) .AND. (BTm_238(iobs) < 999) .AND. (BTo_314(iobs) < 999) .AND.(BTm_314(iobs) < 999) &
      .AND.  (BTO_528(iobs) < 999) .AND. (BTm_528(iobs) < 999) .AND. (zenith(iobs) < 999) .AND. (sfctype(iobs) < 999)
      CldIndex(iobs) = 0.0D0
      IF (isflg) THEN
        c0 = 8.240 - (2.622 - 1.846 * cos(zenith(iobs))) * cos(zenith(iobs))  ! Note that cos function uses a radian value
        CLWo = cos(zenith(iobs)) * (c0 + c1 * log(285 - BTo_238(iobs)) + c2 * log(285 - BTo_314(iobs)))
        CLWm = cos(zenith(iobs)) * (c0 + c1 * log(285 - BTm_238(iobs)) + c2 * log(285 - BTm_314(iobs)))
        deltaCLW = CLWo - CLWm
        ! PRINT *, 'Check deltaCLW: ', deltaCLW

        IF ( ABS(sfctype(iobs) -3) <0.1 ) THEN
          ! Over ocean
          ! PRINT *, 'Check deltaCLW: ', deltaCLW
          CldIndex(iobs) = (w1 * deltaCLW)**2 + (w2 * (BTo_528(iobs) - BTm_528(iobs)))**2
        ELSE
          ! Over land, etc
          CldIndex(iobs) = (w1 * 0.6)**2 + (w2 * (BTo_528(iobs) - BTm_528(iobs)))**2
          ! PRINT *, 'Check values: ', CldIndex(iobs), (w1 * 0.6)**2, (w2 * (BTo_528(iobs) - BTm_528(iobs)))**2
        END IF
      END IF
    END DO
    ! PRINT *, 'Check CldIndex: ', MAXVAL(CldIndex), MINVAL(CldIndex)

  END SUBROUTINE Calc_MWCLWCheck

  SUBROUTINE Calc_MWScatteringIndex(BTo_238, BTm_238, BTo_314, BTm_314, BTo_544, BTm_544, BTo_89, BTm_89, sfctype, SI)
    ! Refer to Grody et al. (1999)
    ! This subroutine should be called after obsThinning.
    ! RawMWTS2 has landmask.
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: BTo_238(:), BTm_238(:)
    REAL(r_kind), INTENT(IN) :: BTo_314(:), BTm_314(:)
    REAL(r_kind), INTENT(IN) :: BTo_544(:), BTm_544(:)
    REAL(r_kind), INTENT(IN) :: BTo_89(:), BTm_89(:), sfctype(:)
    REAL(r_kind), INTENT(OUT) :: SI(:)
    REAL(r_kind) :: w3 = 0.1D0, w4 = 1.0D0/0.8D0, deltas
    INTEGER(i_kind) :: iobs

    ! MWTS2 does NOT has 89GHz channel; MWHS2 has. 
    ! MWTS2 and MWHS2 has different nadir resolutions.
    ! Interpolation is required. Another option is to use MERSI image cloud mask product.
    ! Scattering Index is to exclude precipitating clouds.
    DO iobs = 1, SIZE(BTo_238,1)
      deltas = ( 2.41 - 0.0098 * BTo_238(iobs) ) * (BTo_238(iobs) - BTm_238(iobs)) &
                + 0.454 * (BTo_314(iobs) - BTm_314(iobs)) - (BTo_89(iobs) - BTm_89(iobs))
      IF ( ABS(sfctype(iobs) -3) <0.1 ) THEN
        SI(iobs) = ( w3 * deltas )**2 + ( w4 * (BTo_544(iobs) - BTm_544(iobs)) )**2
      ELSE
        SI(iobs) = 0.8**2 + ( w4 * (BTo_544(iobs) - BTm_544(iobs)) )**2
      END IF
    END DO

  END SUBROUTINE Calc_MWScatteringIndex

  SUBROUTINE QC_MWHS2_ScatteringIndex(obsData, lat, lon, BTo_89, BTo_166, qcflag, ichn)
    ! Any radiances with SI>=1 are excluded
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), lat(:), lon(:)
    REAL(r_kind), INTENT(IN) :: BTo_89(:), BTo_166(:)
    REAL(r_kind) :: SI
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs

    INTEGER :: file_unit, ios
    CHARACTER(len=100) :: filename
    
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!! Set the filename
  !   filename = 'QC_MWHS2_ScatteringIndex.txt'
  ! ! Open the file
  !   file_unit = 10
  !   OPEN(UNIT=file_unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
  !   IF (ios /= 0) THEN
  !     PRINT *, 'Error opening file: ', ios
  !     RETURN
  !   END IF
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      SI = BTo_89(iobs) - BTo_166(iobs)
        !!!!!!!!!!!!!!!!!!!!!!!!Write clwp to file
      ! WRITE(file_unit, '(I5, F10.5)') iobs, landseamask(iobs)  !!!landamsk
      ! WRITE(file_unit, '(I5, F10.5, F10.5, F10.5)') iobs, lat(iobs)*radian2degree, lon(iobs)*radian2degree, SI
      
      IF ( SI .GE. 4 ) THEN
        qcflag(iobs) = 'QC_MWHS2_ScatteringIndex'
        qcnum = qcnum +  1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWHS2_ScatteringIndex for Channnel:', ichn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!close the file 
    ! CLOSE(file_unit)
  END SUBROUTINE QC_MWHS2_ScatteringIndex

  SUBROUTINE QC_MWHS_ScatteringIndex(obsData, lat, lon, BTo_89, BTo_166, zenith, landseamask, qcflag, ichn)
    ! landseamask: 1 land, 2 continentalwater, 0 sea, 5 boundary
    ! Bennartz,et al,1999

    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), lat(:), lon(:)
    REAL(r_kind), INTENT(IN) :: BTo_89(:), BTo_166(:)
    REAL(r_kind), INTENT(IN) ::zenith(:), landseamask(:)

    REAL(r_kind) :: IS
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs

    INTEGER :: file_unit, ios
    CHARACTER(len=100) :: filename
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!! Set the filename
  !   filename = 'QC_MWHS_ScatteringIndex.txt'
  ! ! Open the file
  !   file_unit = 10
  !   OPEN(UNIT=file_unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
  !   IF (ios /= 0) THEN
  !     PRINT *, 'Error opening file: ', ios
  !     RETURN
  !   END IF
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IS = (BTo_89(iobs) - BTo_166(iobs)) - (-39.2010 + 0.1140*cos(zenith(iobs)))
        !!!!!!!!!!!!!!!!!!!!!!!Write clwp to file
      ! WRITE(file_unit, '(I5, F10.5, F10.5, F10.5)') iobs, lat(iobs)*radian2degree, lon(iobs)*radian2degree, IS
      
      IF (  (ABS(landseamask(iobs) -1) .LT. 0.1 .AND. IS .GE. 37.5)  .OR. &
      (ABS(landseamask(iobs) -0) .LT. 0.1 .AND. IS .GE. 5) ) THEN 

        qcflag(iobs) = 'QC_MWHS_ScatteringIndex'
        qcnum = qcnum +  1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWHS_ScatteringIndex for Channnel:', ichn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!close the file 
    ! CLOSE(file_unit)
  END SUBROUTINE QC_MWHS_ScatteringIndex

  SUBROUTINE QC_MWTS2_ScatteringIndex(obsData, BTo_503, BTm_503, landseamask, qcflag, ichn)
    ! landseamask: 1 land, 2 continentalwater, 0 sea, 5 boundary
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    REAL(r_kind), INTENT(IN) :: BTo_503(:), BTm_503(:), landseamask(:)
    REAL(r_kind) :: SI
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      SI = BTo_503(iobs) - BTm_503(iobs)
      IF (  (ABS(landseamask(iobs) -0) .LT. 0.1 .AND. SI .GE. 3.0)  .OR. &
      (ABS(landseamask(iobs) -1) .LT. 0.1 .AND. SI .GE. 1) ) THEN
        qcflag(iobs) = 'QC_MWTS2_ScatteringIndex'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_MWTS2_ScatteringIndex for Channnel:', ichn

  END SUBROUTINE QC_MWTS2_ScatteringIndex

  SUBROUTINE Calc_MWclwp_guess(pres, qcloud, clwp)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: pres(:), qcloud(:)
    REAL(r_kind), INTENT(INOUT) :: clwp
    REAL(r_kind), ALLOCATABLE :: dp(:), clw(:)
    INTEGER :: k

    clwp = ZERO
    ALLOCATE(dp(SIZE(pres)), clw(SIZE(pres)))
    DO k = 1, SIZE(pres)-1
      dp(k) = pres(k) - pres(k+1)
      clw(k) = qcloud(k) * dp(k) / g
      if (pres(k) < 10000.0) clw(k) = ZERO
      clwp = clwp + clw(k)
    END DO

  END SUBROUTINE Calc_MWclwp_guess

  SUBROUTINE ObsErrModel(obsData, innovation, BTclr, camin, camax, ErrClr, ErrCld, BTlim, obsErrs, ca_mean_out)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), innovation(:), BTclr(:), obsErrs(:), ca_mean_out(:)
    REAL(r_kind), INTENT(IN) :: camin, camax, ErrClr, ErrCld, BTlim
    INTEGER(i_kind) :: iobs
    REAL(r_kind) :: BTbak, ca_mean, cx, cy
    
    DO iobs = 1, SIZE(obsData,1)

      IF (ABS(obsData(iobs) - 200.0) < 150.0) THEN

          ! Calculate cloud amount, which is the average of observed cloud amount and background cloud amount
          ! Okamoto et al. (2014)
          ! BTbak = innovation(iobs) + obsData(iobs)
          ! cx = ABS( BTclr(iobs) - BTbak )
          ! cy = ABS( BTclr(iobs) - obsData(iobs))
          ! ca_mean = 0.5 * ( cx +  cy )
          ! print *, 'check ca_mean ', BTclr(iobs), obsData(iobs), innovation(iobs), BTbak, ca_mean

          ! Harnisch et al. (2016)
          BTbak = innovation(iobs) + obsData(iobs)
          cx = MAX( 0.0, BTlim - BTbak )
          cy = MAX( 0.0, BTlim - obsData(iobs))
          ca_mean = 0.5 * ( cx +  cy )
          ca_mean_out(iobs) = ca_mean

          ! Recalculate obs error, which varies with cloud amount
          IF ( ca_mean .LT. camin .AND. ca_mean .GT. missing ) THEN
            obsErrs(iobs) = ErrClr
          ELSE IF ( ca_mean .LT. camax ) then
            obsErrs(iobs) = ErrClr + ( ErrCld - ErrClr ) * &
                ( ca_mean - camin ) / ( camax - camin)
          ELSE
            obsErrs(iobs) = ErrCld
          END IF
      ELSE
          obsErrs(iobs) = missing
      END IF
    END DO

  END SUBROUTINE ObsErrModel

  SUBROUTINE QC_FirstguessCheck(obsData, innovation, obsErrs, inv_threshold, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), innovation(:), obsErrs(:)
    REAL(r_kind), INTENT(IN) :: inv_threshold
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    print *,ichn,'innovation: ',MINVAL(innovation),MAXVAL(innovation), inv_threshold
    print *, 'obsErrs = ', MAXVAL(obsErrs), MINVAL(obsErrs)
    print *, 'obsData = ', MAXVAL(obsData), MINVAL(obsData)
    DO iobs = 1, SIZE(obsData,1)
      
      IF (ABS( innovation(iobs) - missing ) < 1.0 .OR. ABS(innovation(iobs)) > 10.0 * obsErrs(iobs) .OR. &
        ABS(innovation(iobs)) > inv_threshold ) THEN
        qcflag(iobs) = 'QC_FirstguessCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
        ! print *, 'qced: ', ABS( innovation(iobs) - missing ), ABS(innovation(iobs)), 10*obsErrs(iobs), inv_threshold
      END IF
    END DO 

    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_FirstguessCheck for Channnel:', ichn
  END SUBROUTINE QC_FirstguessCheck

  SUBROUTINE QC_SfctypeCheck(obsData, sfctype, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    CHARACTER(len=*), INTENT(IN) :: sfctype(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( TRIM(sfctype(iobs)) .EQ. 'snow' ) THEN
        qcflag(iobs) = 'QC_SfctypeCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF  
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_SfctypeCheck for Channnel:', ichn
  END SUBROUTINE QC_SfctypeCheck

  SUBROUTINE QC_LandseaCheck(obsData, sfctype, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    CHARACTER(len=*), INTENT(IN) :: sfctype(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( TRIM(sfctype(iobs)) .EQ. 'land' ) THEN
        qcflag(iobs) = 'QC_LandseaCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_LandseaCheck for Channnel:', ichn
  END SUBROUTINE QC_LandseaCheck

  !!!!first way
  SUBROUTINE Select_only_over_sea_channels(obsData, sea, sfctype,  qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:)
    INTEGER(i_kind), INTENT(IN) :: sea
    CHARACTER(len=*), INTENT(IN) :: sfctype(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn

    IF (sea .LT. 1.0)  &
      CALL QC_LandseaCheck(obsData, sfctype, qcflag, ichn)
    !   PRINT *, SIZE(obsData,1), 'out of ', SIZE(obsData,1), ' observations are removed by Select_only_over_sea_channels for Channnel:', ichn
    ! ELSE 
    !   PRINT *, ichn, ' channel is active'

  END SUBROUTINE Select_only_over_sea_channels

  !!!!second way
  ! SUBROUTINE Select_only_over_sea_channels(obsData, sea, qcflag, ichn)
  !   IMPLICIT NONE
  !   REAL(r_kind), INTENT(INOUT) :: obsData(:)
  !   CHARACTER(len=*), INTENT(IN) :: sfctype(:)
  !   INTEGER(i_kind), INTENT(IN) :: sea
  !   CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
  !   INTEGER(i_kind), INTENT(IN)     :: ichn
  !   LOGICAL, ALLOCATABLE :: only_over_sea(:)


  !   IF (sea .LT. 1.0) THEN
  !     only_over_sea(this%nchans) = .TRUE.
  !      IF (only_over_sea(this%nchans)) &
  !      CALL QC_LandseaCheck(this%obsData(:, ichan), this%sfctype, qcflag, ichn)
  !   END IF
  
  ! DEALLOCATE(only_over_sea)
  ! END SUBROUTINE Select_only_over_sea_channels

  SUBROUTINE QC_ZenithCheck(obsData, zenith, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), zenith(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( zenith(iobs) .GE. 60.0 .OR. zenith(iobs) .LE. -60.0 ) THEN
        qcflag(iobs) = 'QC_ZenithCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF 
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_ZenithCheck for Channnel:', ichn
  END SUBROUTINE QC_ZenithCheck

  SUBROUTINE QC_ScanEdgeCheck(ObsData, scanpos, Edge1, Edge2, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), scanpos(:)
    REAL(r_kind) :: Edge1, Edge2
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs

    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( scanpos(iobs) .LT. Edge1 .OR. scanpos(iobs) .GT. Edge2 ) THEN
        qcflag(iobs) = 'QC_ScanEdgeCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF 
    END DO 
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_ScanEdgeCheck for Channnel:', ichn

  END SUBROUTINE QC_ScanEdgeCheck

  SUBROUTINE QC_GlintCheck(obsData, satzen, solzen, rel_azi, qcflag, ichn)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: obsData(:), satzen(:), solzen(:), rel_azi(:)
    REAL(r_kind), ALLOCATABLE :: glint_angle(:)
    CHARACTER(len=*), INTENT(INOUT) :: qcflag(:)
    INTEGER(i_kind), INTENT(IN)     :: ichn
    INTEGER(i_kind)    :: qcnum
    INTEGER(i_kind) :: iobs
    
    ! Calculate glint angle
    ALLOCATE(glint_angle(SIZE(satzen)))

    glint_angle = COS(solzen) * COS(satzen) + &
                SIN(solzen) * SIN(satzen) * COS(rel_azi)
    glint_angle = MAX(-1.0 , MIN( glint_angle ,1.0 ))
    glint_angle = ACOS(glint_angle) * radian2degree

    qcnum = 0
    DO iobs = 1, SIZE(obsData,1)
      IF ( glint_angle(iobs) .GE. 40.0 ) THEN
        qcflag(iobs) = 'QC_GlintCheck'
        qcnum = qcnum + 1
        obsData(iobs) = missing
      END IF 
    END DO 
    DEALLOCATE(glint_angle)
    PRINT *, qcnum, 'out of ', SIZE(obsData,1), ' observations are removed by QC_GlintCheck for Channnel:', ichn
  END SUBROUTINE QC_GlintCheck

  FUNCTION Relative_Azimuth ( solazi, satazi )
    IMPLICIT NONE
    REAL(r_kind) :: solazi
    REAL(r_kind) :: satazi
    REAL(r_kind) :: relative_azimuth
  
    relative_azimuth = ABS(solazi - satazi)
    if (relative_azimuth > 180.0) then
         relative_azimuth = 360.0 - relative_azimuth
    endif
    relative_azimuth = (180.0 - relative_azimuth) * degree2radian
  
  end function Relative_Azimuth

  SUBROUTINE Read_bias_coefs_WRFDA(fileunit, filename, npred_max, pred_use, coeffs)
    IMPLICIT NONE
    ! INTEGER, PARAMETER :: npredmax = 8
    INTEGER(i_kind), INTENT(IN)  :: fileunit  ! number of predictors
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER(i_kind), INTENT(IN) :: npred_max
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT) :: pred_use(:,:)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: coeffs(:,:)
    INTEGER(i_kind) :: nchan_in, npred_in, ichan, tmp
    REAL(r_kind), ALLOCATABLE :: mean_in(:),std_in(:)
    INTEGER(i_kind), ALLOCATABLE :: bias_error_guess(:)
    CHARACTER(len=120) :: fparam

    OPEN(fileunit,file=filename,form='formatted')

    READ(fileunit, *) nchan_in, npred_in
    READ(fileunit, *)
    ALLOCATE(mean_in(npred_in), std_in(npred_in), bias_error_guess(npred_in))
    ALLOCATE(pred_use(nchan_in, npred_max), coeffs(nchan_in, npred_max))
  
    READ(fileunit,'(30F10.1)') mean_in(1:npred_in)
    READ(fileunit,'(30F10.1)') std_in(1:npred_in)
    READ(fileunit,'(30I10)')   bias_error_guess(1:npred_in)

    READ(fileunit, *) 

    DO ichan = 1, nchan_in
      IF (npred_in <= 0) cycle               
      WRITE(fparam,*) '(I4,I6,',npred_max,'I3,',npred_in,'F8.3)'
      READ(fileunit, fmt=trim(adjustl(fparam))) tmp, tmp, pred_use(ichan,:), coeffs(ichan,:)   
      ! PRINT *, 'Check varbc.out: ', ichan, ichan, pred_use(ichan,:), coeffs(ichan,:)
    ENDDO

    CLOSE(fileunit)
        
    DEALLOCATE(mean_in, std_in, bias_error_guess)
    ! DEALLOCATE(pred_use, coeffs)

  END SUBROUTINE Read_bias_coefs_WRFDA

  SUBROUTINE Read_bias_coefs_WRFDA_CMAGFS(fileunit, filename, npred, BCOEF, SCOEF)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)  :: fileunit, npred  ! number of predictors
    CHARACTER(len=*), INTENT(IN) :: filename
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: BCOEF(:,:)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT), OPTIONAL :: SCOEF(:,:,:)
    INTEGER :: satid, sensor, nchn, yy, mm, dd, nx, ic, i, j
    CHARACTER(len=17) :: st, sensorname
    INTEGER(i_kind) :: nband,nscan

    OPEN(fileunit,file=filename,form='formatted')
        
    READ(fileunit,*)
    ! READ(fileunit,'(7i6)') satid, sensor, nchn, yy, mm, dd, nx
    READ (fileunit, *) sensorname, nchn
    print *, 'Read headers: ', sensorname, nchn
    ! read air-mass dependent bias coeficient
    ALLOCATE(BCOEF(nchn,npred+1))
    !ALLOCATE(BCOR(nchanl))
    READ(fileunit,*)
    DO i=1,nchn
      READ(fileunit,*) ic, BCOEF(i,1:npred+1)
      ! print *,ic, i, BCOEF(i,1:npred+1),nchn
    ENDDO
  
    ! IF (TRIM(sensorname) .EQ. 'fy4_1-giirs') THEN
    IF ( PRESENT(SCOEF)) THEN
    ! 2.0 read scan dependent bias correction     
      READ(fileunit,'(1x,a17,2i4)') st, nband, nscan
      ! print *,st, nband, nscan
      ALLOCATE(SCOEF(nchn,nscan,nband))
      DO i=1,nchn
        DO j=1,nband
          READ(fileunit,*) ic, SCOEF(i,1:nscan,j)
        ENDDO
      ENDDO
    END IF
    i=nscan*nchn*nband
    CLOSE(fileunit)

    ! DEALLOCATE(BCOEF, SCOEF)

  END SUBROUTINE Read_bias_coefs_WRFDA_CMAGFS

  SUBROUTINE Prof_bkg2obs_fullDomain(X, numObs, olatlon, obsHght, obsTime, Angles, TotalObs, &
                                pres, t, q, tskin, angles_BC)
    USE domainCheck_m, ONLY: domainCheck_t
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: numObs
    REAL(r_kind), INTENT(IN)  :: olatlon(:, :), obsHght(:), Angles(:,:)
    INTEGER(i_kind), INTENT(IN) :: obsTime(:)
    INTEGER(i_kind), INTENT(OUT) :: TotalObs
    TYPE(domainCheck_t) :: domain
    INTEGER(i_kind) :: iobs, ivar, ichan, igrid, itime, ichantmp, io, i, j, k, ivalid
    INTEGER(i_kind), ALLOCATABLE :: Hidxtgt(:, :), Tidxtgt(:, :), Vidxtgt(:, :, :)
    REAL(r_kind), ALLOCATABLE :: Hcoetgt(:, :), Tcoetgt(:, :), Vcoetgt(:, :, :)
    ! Dimensions of 3D bak pres at obs loctions are: (nobs,X%sg%vLevel)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: pres(:, :), t(:, :), q(:, :), tskin(:), angles_BC(:,:)
    ! for local calculation on each processor
    REAL(r_kind), ALLOCATABLE  :: valid_pres(:, :), valid_t(:, :), valid_q(:, :), valid_tskin(:), valid_angles(:,:)
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: maskValided(:)
    INTEGER(i_kind) :: nhstencil, vIdx, hIdx, tIdx, idx(1)
    REAL(r_kind) :: domainLatlon(4), dlat, dlon
    REAL(r_kind), ALLOCATABLE :: forward(:)
    CHARACTER(LEN=20) :: varnametmp

    INCLUDE "mpif.h"

    ALLOCATE (maskValided(numObs))
    domainLatlon(1) = MINVAL(X%sg%cell_cntr(1, :))
    domainLatlon(2) = MAXVAL(X%sg%cell_cntr(1, :))
    domainLatlon(3) = MINVAL(X%sg%cell_cntr(2, :))
    domainLatlon(4) = MAXVAL(X%sg%cell_cntr(2, :))
    ! PRINT *, 'get domain 1: ', domainLatlon(1), domainLatlon(2), domainLatlon(3), domainLatlon(4)

    dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)

    domainLatlon(1) = domainLatlon(1) + dlat*0.5D0
    domainLatlon(2) = domainLatlon(2) - dlat*0.5D0
    domainLatlon(3) = domainLatlon(3) + dlon*0.5D0
    domainLatlon(4) = domainLatlon(4) - dlon*0.5D0
    PRINT *, 'get domain 2: ', domainLatlon(1), domainLatlon(2), domainLatlon(3), domainLatlon(4)

    PRINT *, 'start calculate Prof_bkg2obs_fullDomain'
    ALLOCATE(forward(X%sg%vLevel))
    ! Validating the obs and saving their interpolation coefficients:
    CALL domain%validation(X%sg, numObs, olatlon, obsHght, &
                          obsTime, numValided, maskValided, &
                          Hidxtgt, Vidxtgt, Tidxtgt, &
                          Hcoetgt, Vcoetgt, Tcoetgt, &
                          1, nhstencil, 'Satellite', &
                          domainLatlon(1), domainLatlon(2), domainLatlon(3), domainLatlon(4))
    PRINT *, 'Check nums: ', numObs, numValided
    ! PRINT *, 'check coefs and idx: ', 'H - ', shape(Hcoetgt), 'V - ', shape(Vcoetgt), 'T - ', shape(Tcoetgt)
    ! PRINT *, 'check coefs: ', Hcoetgt(:, 4), Tcoetgt(:, 4)

    IF (PRESENT(pres)) THEN
      ALLOCATE (valid_pres(numValided, X%sg%vLevel)); valid_pres = ZERO
    END IF
    IF (PRESENT(t)) THEN
      ALLOCATE (valid_t(numValided, X%sg%vLevel)); valid_t = ZERO
    END IF
    IF (PRESENT(q)) THEN
      ALLOCATE (valid_q(numValided, X%sg%vLevel)); valid_q = ZERO
    END IF
    IF (PRESENT(tskin)) THEN
      ALLOCATE (valid_tskin(numValided)); valid_tskin = ZERO
    END IF
    IF (PRESENT(angles_BC)) THEN
      ALLOCATE (valid_angles(4, numValided)); valid_angles = ZERO
    END IF

    ivalid = 0
    DO io = 1, numObs

      IF (maskValided(io) .GT. 0) THEN
        ivalid = ivalid + 1

        IF (PRESENT(angles_BC)) valid_angles(:, ivalid) = Angles(:, io)

        ! PRINT *, 'Check coefs: ', Hcoetgt(:,ivalid), Tcoetgt(:,ivalid)
        DO i = 1, 2 ! 2 time frames of interpolation
          DO j = 1, nhstencil ! 3 horizontal interpolation points

            ! vIdx = Vidxtgt(k,j,ivalid)
            hIdx = Hidxtgt(j, ivalid)
            tIdx = Tidxtgt(i, ivalid)

            varnametmp = 'pres'
            IF (PRESENT(pres)) THEN
              forward = X%sg%FGPres(:, hIdx, tIdx)
              valid_pres(ivalid, :) = valid_pres(ivalid, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'temp'
            IF (PRESENT(t)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              valid_t(ivalid, :) = valid_t(ivalid, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qvapor'
            IF (PRESENT(q)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              valid_q(ivalid, :) = valid_q(ivalid, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            IF (PRESENT(tskin)) THEN
              forward(1) = X%sg%tskin(hIdx, tIdx)
              valid_tskin(ivalid) = valid_tskin(ivalid) + forward(1) * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

          END DO
        END DO

      END IF

    END DO

    ! Deallocate memory:

    DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided, forward)

    BLOCK 
      ! For mpi
      INTEGER(i_kind), ALLOCATABLE :: ncount_group(:), disp_group(:)
      REAL(r_kind), ALLOCATABLE :: idxArray2D(:,:), idxArray1D(:)
      INTEGER(i_kind) :: i, icount, ilevel

      ALLOCATE (ncount_group(X%sg%mpddInfo_sg%nproc))
      ALLOCATE (disp_group(X%sg%mpddInfo_sg%nproc))

      ncount_group = 0
      disp_group = 0

      ! PRINT *, 'nproc ', X%sg%mpddInfo_sg%nProc
      ! 将不同进程上的观测个数汇总到主线程 0
      ! print *, 'Prof_bkgAtobs valided: ', numValided, ivalid
      CALL MPI_GATHER(numValided, 1, MPI_INTEGER4, &
                      ncount_group, 1, MPI_INTEGER4, 0, &
                      X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)

      TotalObs = SUM(ncount_group)
      CALL X%sg%mpddInfo_sg%bcast(TotalObs)
      PRINT *, 'TotalObs of the full domain is: ', TotalObs
      ALLOCATE (pres(TotalObs, X%sg%vLevel))
      ALLOCATE (t(TotalObs, X%sg%vLevel))
      ALLOCATE (q(TotalObs, X%sg%vLevel))
      ALLOCATE (angles_BC(4, TotalObs))
      ALLOCATE(tskin(TotalObs))

      IF (X%sg%isBaseProc()) THEN
        ! PRINT *, "ncount_group is ", SUM(ncount_group)
        ALLOCATE (idxArray2D(TotalObs, X%sg%vLevel), idxArray1D(TotalObs))
        idxArray2D = ZERO
        idxArray1D = ZERO
        disp_group = 0
        FORALL (i=2:X%sg%mpddInfo_sg%nproc) disp_group(i) = SUM(ncount_group(1:i - 1)) !Gatherv 时，每个进程所需的位移
        ! PRINT *, 'disp_group = ', disp_group
      END IF

      CALL X%sg%mpddInfo_sg%bcast(ncount_group)
      CALL X%sg%mpddInfo_sg%bcast(disp_group)

      ! 将所有进程的innovation值汇总，存在idxArray中
      IF (PRESENT(pres)) THEN
        DO ilevel = 1, X%sg%vLevel
          CALL MPI_GATHERV(valid_pres(:,ilevel), numValided, MPI_DOUBLE_PRECISION, &
                          idxArray2D(:,ilevel), ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                          X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)
        END DO
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) &
          ALLOCATE (idxArray2D(TotalObs, X%sg%vLevel))
        CALL X%sg%mpddInfo_sg%bcast(idxArray2D)
        pres = idxArray2D
        ! PRINT *, 'p is bcast'
        DEALLOCATE(valid_pres)
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) DEALLOCATE(idxArray2D)
      END IF

      IF (PRESENT(t)) THEN
        DO ilevel = 1, X%sg%vLevel
          CALL MPI_GATHERV(valid_t(:,ilevel), numValided, MPI_DOUBLE_PRECISION, &
                          idxArray2D(:,ilevel), ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                          X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)
        END DO
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) &
          ALLOCATE (idxArray2D(TotalObs, X%sg%vLevel))
        CALL X%sg%mpddInfo_sg%bcast(idxArray2D)
        t = idxArray2D
        ! PRINT *, 't is bcast'
        DEALLOCATE(valid_t)
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) DEALLOCATE(idxArray2D)
      END IF

      IF (PRESENT(q)) THEN
        DO ilevel = 1, X%sg%vLevel
          CALL MPI_GATHERV(valid_q(:,ilevel), numValided, MPI_DOUBLE_PRECISION, &
                          idxArray2D(:,ilevel), ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                          X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)
        END DO
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) &
          ALLOCATE (idxArray2D(TotalObs, X%sg%vLevel))
        CALL X%sg%mpddInfo_sg%bcast(idxArray2D)
        q = idxArray2D
        ! PRINT *, 'q is bcast'
        DEALLOCATE(valid_q)
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) DEALLOCATE(idxArray2D)
      END IF

      IF (PRESENT(angles_BC)) THEN
        DO ilevel = 1, 4
          CALL MPI_GATHERV(valid_angles(ilevel,:), numValided, MPI_DOUBLE_PRECISION, &
                          idxArray2D(:, ilevel), ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                          X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)
        END DO
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) &
          ALLOCATE (idxArray2D(TotalObs, X%sg%vLevel))
        CALL X%sg%mpddInfo_sg%bcast(idxArray2D)
        DO ilevel = 1, 4
          angles_BC(ilevel,:) = idxArray2D(:,ilevel)
        END DO
        ! PRINT *, 'Angles is bcast: ', SHAPE(angles_BC)
        DEALLOCATE(valid_angles)
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) DEALLOCATE(idxArray2D)
      END IF

      IF (PRESENT(tskin)) THEN
        CALL MPI_GATHERV(valid_tskin, numValided, MPI_DOUBLE_PRECISION, &
                        idxArray1D, ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                        X%sg%mpddInfo_sg%comm, X%sg%mpddInfo_sg%ierr)
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) &
          ALLOCATE (idxArray1D(TotalObs))
        CALL X%sg%mpddInfo_sg%bcast(idxArray1D)
        tskin = idxArray1D
        DEALLOCATE(valid_tskin)
        ! PRINT *, 'tskin is bcast'
        IF (X%sg%isActiveProc() .AND. (.NOT. X%sg%isBaseProc())) DEALLOCATE(idxArray1D)
      END IF

      IF (X%sg%isBaseProc()) DEALLOCATE (idxArray2D, idxArray1D)
      DEALLOCATE (ncount_group, disp_group)

    END BLOCK

  END SUBROUTINE Prof_bkg2obs_fullDomain

  SUBROUTINE Prof_bkg2obs(X, numObs, olatlon, obsHght, obsTime, &
                          pres, t, q, u, v, qc, qr, qi, qs, qg, tskin, elevation, &
                          snowc, soiltype, landmask, PassDomainCheck)
    USE domainCheck_m, ONLY: domainCheck_t
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: numObs
    REAL(r_kind), INTENT(IN)  :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), INTENT(IN) :: obsTime(:)
    TYPE(domainCheck_t) :: domain
    INTEGER(i_kind) :: iobs, ivar, ichan, igrid, itime, ichantmp, io, i, j, k, ivalid
    INTEGER(i_kind), ALLOCATABLE :: Hidxtgt(:, :), Tidxtgt(:, :), Vidxtgt(:, :, :)
    REAL(r_kind), ALLOCATABLE :: Hcoetgt(:, :), Tcoetgt(:, :), Vcoetgt(:, :, :)
    ! Dimensions of 3D bak pres at obs loctions are: (nobs,X%sg%vLevel)
    REAL(r_kind), INTENT(INOUT), OPTIONAL  :: pres(:, :), t(:, :), q(:, :), u(:, :), v(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL  :: qc(:,:), qr(:,:), qi(:,:), qs(:,:), qg(:,:)
    ! Dimensions of 2D bak tskin at obs loctions are: (nobs)
    REAL(r_kind), INTENT(INOUT), OPTIONAL  :: tskin(:), elevation(:), snowc(:), soiltype(:), landmask(:)
    LOGICAL, INTENT(OUT), OPTIONAL :: PassDomainCheck(:)
    INTEGER(i_kind) :: numValided
    INTEGER(i_kind), ALLOCATABLE :: maskValided(:)
    INTEGER(i_kind) :: nhstencil, vIdx, hIdx, tIdx, idx(1)
    REAL(r_kind), ALLOCATABLE :: forward(:)
    CHARACTER(LEN=20) :: varnametmp
    REAL(r_kind) :: domainLatlon(4), dlat, dlon

    ALLOCATE (maskValided(numObs))
    IF (PRESENT(PassDomainCheck)) PassDomainCheck = .TRUE.

    IF (PRESENT(pres)) pres = ZERO
    IF (PRESENT(t)) t = ZERO
    IF (PRESENT(q)) q = ZERO
    IF (PRESENT(qc)) qc = ZERO
    IF (PRESENT(qr)) qr = ZERO
    IF (PRESENT(qi)) qi = ZERO
    IF (PRESENT(qs)) qs = ZERO
    IF (PRESENT(qg)) qg = ZERO
    IF (PRESENT(u)) u = ZERO
    IF (PRESENT(v)) v = ZERO
    IF (PRESENT(tskin)) tskin = ZERO
    IF (PRESENT(elevation)) elevation = ZERO
    IF (PRESENT(snowc)) snowc = ZERO
    IF (PRESENT(soiltype)) soiltype = ZERO
    IF (PRESENT(landmask)) landmask = ZERO

    PRINT *, 'start calculate Prof_bkg2obs'
    ALLOCATE(forward(X%sg%vLevel))
    ! Validating the obs and saving their interpolation coefficients:

    domainLatlon(1) = MINVAL(X%sg%cell_cntr(1, :))
    domainLatlon(2) = MAXVAL(X%sg%cell_cntr(1, :))
    domainLatlon(3) = MINVAL(X%sg%cell_cntr(2, :))
    domainLatlon(4) = MAXVAL(X%sg%cell_cntr(2, :))
    ! PRINT *, 'get domain 1: ', domainLatlon(1), domainLatlon(2), domainLatlon(3), domainLatlon(4)

    dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)

    domainLatlon(1) = domainLatlon(1) + dlat*0.5D0
    domainLatlon(2) = domainLatlon(2) - dlat*0.5D0
    domainLatlon(3) = domainLatlon(3) + dlon*0.5D0
    domainLatlon(4) = domainLatlon(4) - dlon*0.5D0

    CALL domain%validation(X%sg, numObs, olatlon, obsHght, &
                          obsTime, numValided, maskValided, &
                          Hidxtgt, Vidxtgt, Tidxtgt, &
                          Hcoetgt, Vcoetgt, Tcoetgt, &
                          1, nhstencil, 'Satellite', &
                          domainLatlon(1), domainLatlon(2), domainLatlon(3), domainLatlon(4))

    ! PRINT *, 'check coefs and idx: ', 'H - ', shape(Hcoetgt), 'V - ', shape(Vcoetgt), 'T - ', shape(Tcoetgt)
    ! ! PRINT *, 'Prof_bkg2obs numObs = ', numObs
    ! PRINT *, 'Prof_bkg2obs numValided = ', numValided

    ivalid = 0
    DO io = 1, numObs

      IF (maskValided(io) .GT. 0) THEN
        ivalid = ivalid + 1

        idx(1:1) = MAXLOC(Hcoetgt(:, ivalid))

        ! PRINT *, 'MAXLOC idx = ', idx
        IF (PRESENT(landmask)) landmask(io) = X%sg%landmask(Hidxtgt(idx(1), ivalid))
        IF (PRESENT(soiltype)) soiltype(io) = X%sg%soil_type(Hidxtgt(idx(1), ivalid))

        ! PRINT *, 'Check coefs: ', Hcoetgt(:,ivalid), Tcoetgt(:,ivalid)
        DO i = 1, 2 ! 2 time frames of interpolation
          DO j = 1, nhstencil ! 3 horizontal interpolation points

            ! vIdx = Vidxtgt(k,j,ivalid)
            hIdx = Hidxtgt(j, ivalid)
            tIdx = Tidxtgt(i, ivalid)

            varnametmp = 'pres'
            IF (PRESENT(pres)) THEN
              forward = X%sg%FGPres(:, hIdx, tIdx)
              pres(io, :) = pres(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'temp'
            IF (PRESENT(t)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              t(io, :) = t(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qvapor'
            IF (PRESENT(t)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              q(io, :) = q(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'uwnd'
            IF (PRESENT(t)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              u(io, :) = u(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'vwnd'
            IF (PRESENT(t)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              v(io, :) = v(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qcloud'
            IF (PRESENT(qc)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              qc(io, :) = qc(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qrain'
            IF (PRESENT(qr)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              qr(io, :) = qr(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qice'
            IF (PRESENT(qi)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              qi(io, :) = qi(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qsnow'
            IF (PRESENT(qs)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              qs(io, :) = qs(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            varnametmp = 'qgraupel'
            IF (PRESENT(qg)) THEN
              forward = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(:, hIdx, tIdx)
              qg(io, :) = qg(io, :) + forward * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            IF (PRESENT(tskin)) THEN
              forward(1) = X%sg%tskin(hIdx, tIdx)
              tskin(io) = tskin(io) + forward(1) * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

            IF (PRESENT(elevation)) THEN
              forward(1) = X%sg%topo(hIdx)
              elevation(io) = elevation(io) + forward(1) * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF
            
            IF (PRESENT(snowc)) THEN
              forward(1) = X%sg%snowc(hIdx, tIdx)
              snowc(io) = snowc(io) + forward(1) * Hcoetgt(j, ivalid) * Tcoetgt(i, ivalid)
            END IF

          END DO
        END DO

      ELSE
        IF (PRESENT(PassDomainCheck)) PassDomainCheck(io) = .FALSE.
      END IF

    END DO
    print *, 'check nums 1: ', ivalid, numObs

    IF (PRESENT(qc)) WHERE(qc <= 1.0D-12) qc = 0.0D0
    IF (PRESENT(qi)) WHERE(qi <= 1.0D-12) qi = 0.0D0

    ! Deallocate memory:

    DEALLOCATE (Hidxtgt, Hcoetgt, Tidxtgt, Tcoetgt, Vidxtgt, Vcoetgt, maskValided, forward)

END SUBROUTINE Prof_bkg2obs

END MODULE DP_utils_m
