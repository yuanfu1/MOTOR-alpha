!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Obs
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu, Ting Shu
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!   Added "wavelet-based satellite data and error decomposition" by Ting Shu (shuting@gbamwf.com), 2023/08/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE RadBase_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  ! USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE Dp_utils_m
  USE Satellite_utils_m
  USE RawSatellite_BC_m
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_nl_sp_obsloc_m, ONLY: rttov_nl_sp_obsloc_t
  ! 2023-08-09, Ting Shu
  USE satelliteDecom_m, ONLY: satelliteDecom_t

  IMPLICIT NONE

  TYPE :: RadBase_t
    TYPE(rttov_nl_sp_t) :: OprRTTOV
    TYPE(rttov_nl_sp_obsloc_t) :: rttov_nl_sp_obsloc
    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind) :: nchans, npred
    INTEGER(i_kind) :: mgStart, mgEnd
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :), &      ! observation errors
                                 ca_mean(:,:)          ! cloud amount
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    CHARACTER(LEN=40) :: obsType
    INTEGER(i_kind) :: numObs, numVars, numFiles       ! Number of obs, vars and obs files
    INTEGER(i_kind) :: interpolation 
    INTEGER :: idxFile
    CHARACTER(LEN=20), ALLOCATABLE :: varNames(:)
    INTEGER(i_kind), ALLOCATABLE :: chan_lists(:)
    INTEGER(i_kind), ALLOCATABLE :: rttov_chan_lists(:)
    INTEGER(i_kind), ALLOCATABLE :: ifuse(:)
    INTEGER(i_kind), ALLOCATABLE :: oversea(:)
    INTEGER(i_kind), ALLOCATABLE :: iscanpos(:)
    REAL(r_kind), ALLOCATABLE :: BCOEF(:,:), SCOEF(:, :,:)
    REAL(r_kind), ALLOCATABLE :: Err_in(:), bias_in(:), camin(:), camax(:), obsErrClr(:), obsErrCld(:), BTlim(:)
    CHARACTER(LEN=20) :: bias_scheme = 'pdf_mode'
    CHARACTER(LEN=50), ALLOCATABLE :: bias_predictors(:), cloudy_predictors(:)
    LOGICAL :: CloudyMode = .FALSE., use_cloudy_predictors = .FALSE.
    REAL(r_kind), ALLOCATABLE :: tb_inv(:), tb_clr(:), tb_bakAtmodel(:,:)
    CHARACTER(len=10), ALLOCATABLE :: sfctype(:)
    ! INTEGER(i_kind) :: rttov_chan = 0
    LOGICAL :: Thin_w_bkg = .TRUE.
    ! LOGICAL, ALLOCATABLE :: Nchannels(:)
    CHARACTER(LEN=10) :: Source_of_bkg = 'obsavg'
    REAL(r_kind), ALLOCATABLE :: Angles(:, :), Thinned_Angles(:,:)
    INTEGER(i_kind), ALLOCATABLE :: nbins
    REAL(r_kind) :: obsData_avg = missing
    INTEGER(i_kind), ALLOCATABLE :: Thinned_Angles_hidx(:), Thinned_Angles_tidx(:)
    REAL(r_kind), ALLOCATABLE :: cloud_flag(:), obsData_org(:,:), SI(:)
    ! CHARACTER(LEN=20) :: bias_scheme = 'pdf_mode'
    LOGICAL           :: norm_pred = .FALSE.
    ! 2023-08-09, Ting Shu
    TYPE(satelliteDecom_t) :: sateDecom
    ! 2023-08-30, Ting Shu, add flag of using wavelet or not
    LOGICAL :: useWaveFlag

  CONTAINS

    PROCEDURE, PUBLIC  :: RadBase_QC1, RadBase_QC2
    PROCEDURE, PUBLIC  :: RadBase_BC
    PROCEDURE, PUBLIC :: ObsPrepareForSgBase
    PROCEDURE, PUBLIC :: OutputForThinningBase
    PROCEDURE, PUBLIC :: Initialize
    PROCEDURE, PUBLIC :: destroy
  END TYPE RadBase_t

CONTAINS

  SUBROUTINE Initialize(this, X, bias_scheme, bias_predictors, CloudyMode, &
                        use_cloudy_predictors, cloudy_predictors, &
                        camin, camax, obsErrClr, obsErrCld, BTlim)
    
    CLASS(RadBase_t) :: this
    TYPE(State_t) :: X
    ! REAL(r_kind), INTENT(IN) :: tb_bakAtmodel(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: bias_scheme
    CHARACTER(LEN=50), INTENT(IN) :: bias_predictors(:)
    LOGICAL, INTENT(IN), OPTIONAL :: CloudyMode
    LOGICAL, INTENT(IN), OPTIONAL :: use_cloudy_predictors
    CHARACTER(LEN=50), INTENT(IN), OPTIONAL :: cloudy_predictors(:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: camin(:), camax(:), obsErrClr(:), obsErrCld(:), BTlim(:)
    
    ! This is an optional holder for data prepare polymorphic subroutine
    ! this%tb_bakAtmodel = tb_bakAtmodel
    this%bias_scheme = bias_scheme
    this%bias_predictors = bias_predictors
    IF (PRESENT(CloudyMode)) this%CloudyMode = CloudyMode
    IF (PRESENT(use_cloudy_predictors)) this%use_cloudy_predictors = use_cloudy_predictors
    IF (PRESENT(cloudy_predictors)) this%cloudy_predictors = cloudy_predictors
    IF (PRESENT(camin)) THEN
      this%camin = camin
      this%camax = camax
      this%obsErrClr = obsErrClr
      this%obsErrCld = obsErrCld
      this%BTlim = BTlim
    END IF
  END SUBROUTINE Initialize

  SUBROUTINE ObsPrepareForSgBase(this, X, configFile, inst_name, platform_name)
    
    CLASS(RadBase_t) :: this
    TYPE(State_t) :: X
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(len=20), INTENT(IN) :: inst_name, platform_name
    
    INTEGER(i_kind) :: i, j, istatus, ifile
    
    this%configFile = configFile
    this%nbins = 40
    this%numVars = this%nchans + 4 
    !FORALL (i = 1: this%numVars) this%RadBase(i)%nbins = this%nbins 
    !FORALL is not work for intel
    ! IF(X%sg%isBaseProc()) PRINT *, 'Check nbins ', this%nbins

    ! This is an optional holder for data prepare polymorphic subroutine
    ALLOCATE (this%tb_bakAtmodel(X%sg%num_cell, X%sg%tSlots))
    this%tb_bakAtmodel = ZERO

    CALL this%OprRTTOV%initialize(this%configFile, X, inst_name, platform_name)
    CALL this%rttov_nl_sp_obsloc%initialize(this%configFile, X, inst_name, platform_name)

    WHERE (X%fields(X%getVarIdx('qvapor'))%data < qv_limit) X%fields(X%getVarIdx('qvapor'))%data = qv_limit

    ! 2023-08-30, Ting Shu, get flag of using wavelet or not from yaml
    ifile = yaml_get_var(TRIM(configFile), 'Wavelet', 'useWave', this%useWaveFlag)
    ifile = yaml_get_var(this%configFile, 'RTTOV', 'rttov_clouds', this%CloudyMode)

  END SUBROUTINE ObsPrepareForSgBase

  SUBROUTINE GetForwardValue_obsloc(this, rttov_nl_sp_obsloc, vLevel, rttov_ichan, pres, t, q, u, v, tskin, landmask, &
                                    elevation, soil_type, snowc, &
                                    lat, lon, SatAngles, BT, BTclr, qc, qr, qi, qs, qg)
    TYPE(RadBase_t), INTENT(INOUT) :: this
    TYPE(rttov_nl_sp_obsloc_t), INTENT(IN) :: rttov_nl_sp_obsloc
    INTEGER(i_kind), INTENT(IN)  :: vLevel
    REAL(r_kind), INTENT(IN) :: pres(:), t(:), q(:), u(:), v(:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: qc(:), qr(:), qi(:), qs(:), qg(:)
    REAL(r_kind), INTENT(IN) :: tskin, elevation, lat, lon, landmask, snowc, soil_type
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    INTEGER(i_kind), INTENT(IN) :: rttov_ichan
    REAL(r_kind), INTENT(OUT) :: BT, BTclr       ! For each instrument, each channel

    IF (PRESENT(qc) .AND. PRESENT(qi)) THEN
      CALL rttov_nl_sp_obsloc%simobsloc(vLevel, rttov_ichan, pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc, &
          lat, lon, SatAngles, BT, BTclr, &
          qc=qc, qr=qr, qi=qi, qs=qs, qg=qg)
    ELSE
      CALL rttov_nl_sp_obsloc%simobsloc(vLevel, rttov_ichan, pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc, &
          lat, lon, SatAngles, BT, BTclr)
    END IF

  END SUBROUTINE GetForwardValue_obsloc

  SUBROUTINE OutputForThinningBase(this, state, mpObs, filename, inst_name, platform_name)
    USE State_m, ONLY: State_t
    ! USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    USE mo_netcdf, only: NcDataset, NcDimension
    IMPLICIT NONE

    CLASS(RadBase_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    CHARACTER(*), INTENT(IN) :: filename, inst_name, platform_name
    INTEGER(i_kind) :: TotalObs
    REAL(r_kind), ALLOCATABLE :: pres(:,:), t(:,:), q(:,:), u(:,:), v(:,:), tskin(:), elevation(:)
    REAL(r_kind), ALLOCATABLE :: qc(:,:), qr(:,:), qi(:,:), qs(:,:), qg(:,:)
    REAL(r_kind), ALLOCATABLE :: snowc(:), soiltype(:), landmask(:)
    REAL(r_kind), ALLOCATABLE :: pres_BC(:, :), t_BC(:, :), q_BC(:, :), tskin_BC(:), angles_BC(:,:)
    LOGICAL, ALLOCATABLE :: PassDomainCheck(:)
    INTEGER :: i, j, Nqcout
    TYPE(NcDataset)    :: nc
    TYPE(NcDimension) :: dim3
    REAL(r_kind), ALLOCATABLE :: ObsData(:, :), ObsData_raw(:,:), ObsErrs(:,:), ca_mean(:,:), INVData(:,:), ClrData(:,:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    REAL(r_kind), ALLOCATABLE :: olatlon(:,:)
    REAL(r_kind), ALLOCATABLE :: tb_inv(:), tb_bak503(:)
    INTEGER(i_kind) :: ch503
    INTEGER :: ichan
    ! 2023-08-09, Ting Shu
    INTEGER(i_kind) :: fn_len
    CHARACTER(LEN=:), ALLOCATABLE :: coeff_fileName
    
    IF(state%sg%isBaseProc()) PRINT *, "START QC and BC for satellite obs. Num of obs = ", this%numobs
    IF (this%numObs < 1) RETURN
    IF (mpObs%isActiveProc()) THEN
           
      ! Calculate background profiles at obs locations
      ALLOCATE (pres(this%numobs, state%sg%vLevel))
      ALLOCATE (t(this%numobs, state%sg%vLevel))
      ALLOCATE (q(this%numobs, state%sg%vLevel))
      ALLOCATE(u(this%numobs, state%sg%vLevel))
      ALLOCATE(v(this%numobs, state%sg%vLevel))
      IF (this%CloudyMode) THEN
        ALLOCATE(qc(this%numobs, state%sg%vLevel))
        ALLOCATE(qr(this%numobs, state%sg%vLevel))
        ALLOCATE(qi(this%numobs, state%sg%vLevel))
        ALLOCATE(qs(this%numobs, state%sg%vLevel))
        ALLOCATE(qg(this%numobs, state%sg%vLevel))
      END IF
      ALLOCATE(tskin(this%numobs))
      ALLOCATE(elevation(this%numobs))
      ALLOCATE(snowc(this%numobs))
      ALLOCATE(soiltype(this%numobs))
      ALLOCATE(landmask(this%numobs))
      ALLOCATE(PassDomainCheck(this%numobs))
      ALLOCATE(tb_inv(this%numobs))
      ALLOCATE(tb_bak503(this%numobs))
      
      PRINT *, 'this%numobs = ', this%numobs
      IF (this%CloudyMode) THEN
        CALL Prof_bkg2obs(state, this%numobs, this%olatlon, this%obsHght, this%obsTime, &
                        pres=pres, t=t, q=q, u=u, v=v, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, tskin=tskin, elevation=elevation, &
                        snowc=snowc, soiltype=soiltype, landmask=landmask, PassDomainCheck=PassDomainCheck)
      ELSE
        CALL Prof_bkg2obs(state, this%numobs, this%olatlon, this%obsHght, this%obsTime, &
                        pres=pres, t=t, q=q, u=u, v=v, tskin=tskin, elevation=elevation, &
                        snowc=snowc, soiltype=soiltype, landmask=landmask, PassDomainCheck=PassDomainCheck)
      END IF
      CALL Prof_bkg2obs_fullDomain(state, this%numobs, this%olatlon, this%obsHght, this%obsTime, this%Angles, TotalObs, &
                        pres_BC, t_BC, q_BC, tskin_BC, angles_BC)

      ALLOCATE(ObsData(TotalObs, this%numVars))
      ALLOCATE(ObsData_raw(TotalObs, this%numVars))
      ALLOCATE(ObsErrs(TotalObs, this%numVars))
      IF (ALLOCATED(this%ca_mean)) ALLOCATE(ca_mean(TotalObs, this%numVars))
      ALLOCATE(INVData(TotalObs, this%numVars))
      ALLOCATE(ClrData(TotalObs, this%numVars))
      ALLOCATE(obsTime(TotalObs))
      ALLOCATE(olatlon(2, TotalObs))

      DO i = 1, this%numVars
       
        ! ObsData_raw(:,i) = this%ObsData(:,i)
        CALL MPI_GATHER_Interface(state%sg, this%ObsData(:,i), ObsData_raw(:, i))
        CALL MPI_GATHER_Interface_Int(state%sg, this%ObsTime, obsTime)
        CALL MPI_GATHER_Interface(state%sg, this%olatlon(1,:), olatlon(1,:))
        CALL MPI_GATHER_Interface(state%sg, this%olatlon(2,:), olatlon(2,:)) 

        IF (SIZE(this%obsData,1) > 0) THEN
          this%obsData_avg = SUM(ObsData_raw(:,i))/REAL(SIZE(ObsData_raw,1))
        ELSE 
          this%obsData_avg = missing
        END IF

        IF (i<5) THEN

          ! ObsData(:, i) = this%ObsData(:,i)
          ! ObsErrs(:, i) = this%obsErrs(:,i)
          InvData(:, i) = 0.0D0
          IF (ALLOCATED(this%ca_mean)) ca_mean(:,i) = 0.0D0
          ! ClrData(:,i) = this%ObsData(:,i)
          CALL MPI_GATHER_Interface(state%sg, this%ObsData(:,i), ObsData(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%obsErrs(:,i), ObsErrs(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%ObsData(:,i), ClrData(:, i))

        ELSE

          this%obsData(:,i) = this%obsData_org(:,i)
          ! CALL this%ObsPrepareForSgBase(state)

          IF (ALLOCATED(qc)) THEN
            CALL this%RadBase_QC1(i, TRIM(inst_name), this%rttov_chan_lists(i-4), Nqcout, pres=pres, qcloud=qc)
          ELSE
            CALL this%RadBase_QC1(i, TRIM(inst_name), this%rttov_chan_lists(i-4), Nqcout)
          END IF

          IF (this%CloudyMode) THEN
            CALL TB_bkgAtobs(this, state, this%rttov_nl_sp_obsloc, pres, t, q, u, v, tskin, elevation, &
                        snowc, soiltype, landmask, PassDomainCheck, i, &
                        qc, qr, qi, qs, qg)
          ELSE
            CALL TB_bkgAtobs(this, state, this%rttov_nl_sp_obsloc, pres, t, q, u, v, tskin, elevation, &
                        snowc, soiltype, landmask, PassDomainCheck, i) 
          END IF

          IF (TRIM(inst_name) .EQ. 'mwts2') THEN
            ch503 = 5
          END IF
          IF (TRIM(inst_name) .EQ. 'mwts3') THEN
            ch503 = 7
          END IF
          IF (TRIM(inst_name) .EQ. 'mwts2' .OR. TRIM(inst_name) .EQ. 'mwts3') THEN
            IF (i == ch503) tb_inv = this%tb_inv
            tb_bak503 = tb_inv + this%obsData_org(:,i)
          ELSE
            tb_bak503 = missing
          END IF

          ! CALL this%initialize(state, this%bias_scheme, this%bias_predictors, this%CloudyMode, &
          !                                 this%camin, this%camax, this%obsErrClr, this%obsErrCld)

          CALL this%RadBase_QC2(i, TRIM(inst_name), this%rttov_chan_lists(i-4), Nqcout, landmask, tb_bak503, ch503)
          ! print *, 'check sizes of ca_mean 1: ', SHAPE(this%ca_mean)
          IF (ALLOCATED(this%ca_mean)) &
          CALL MPI_GATHER_Interface(state%sg, this%ca_mean(:,i), ca_mean(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%ObsData(:,i), ObsData(:, i))
          ! print *, 'check sizes of ca_mean 2: ', SHAPE(ca_mean)
          IF (state%sg%gLevel == this%mgEnd) THEN
            ! PRINT *, '000: ', state%sg%mpddInfo_sg%myrank,TRIM(this%bias_scheme),Nqcout,SIZE(this%obsData,1), TotalObs
            IF ( TRIM(this%bias_scheme) .EQ. 'pdf_mode' .AND. Nqcout .LT. SIZE(this%obsData,1) .AND. TotalObs .GT. 100) &
              CALL this%RadBase_BC(i, state, TotalObs)
            IF ( TRIM(this%bias_scheme) .EQ. 'ScanAirs' .AND. Nqcout .LT. SIZE(this%obsData,1) .AND. TotalObs .GT. 100) THEN
              IF ( this%use_cloudy_predictors .AND. ALLOCATED(this%ca_mean) ) THEN
                CALL this%RadBase_BC(i, state, TotalObs, PassDomainCheck=PassDomainCheck, pres=pres, t=t, q=q, tskin=tskin, &
                                              pres_BC=pres_BC, t_BC=t_BC, q_BC=q_BC, tskin_BC=tskin_BC, angles_BC=angles_BC, &
                                              ca_mean = this%ca_mean(:,i), tb_obs = this%ObsData(:,i), &
                                              ca_mean_BC = ca_mean(:,i), tb_obs_BC = ObsData(:,i), &
                                              npred=this%npred, BCOEF=this%BCOEF(this%rttov_chan_lists(i-4), :), &
                                              iscanpos=this%iscanpos)
              ELSE
                CALL this%RadBase_BC(i, state, TotalObs, PassDomainCheck=PassDomainCheck, pres=pres, t=t, q=q, tskin=tskin, &
                                              pres_BC=pres_BC, t_BC=t_BC, q_BC=q_BC, tskin_BC=tskin_BC, angles_BC=angles_BC, &
                                              npred=this%npred, BCOEF=this%BCOEF(this%rttov_chan_lists(i-4), :), &
                                              iscanpos=this%iscanpos)
              END IF
            END IF
          END IF
          CALL MPI_GATHER_Interface(state%sg, this%ObsData(:,i), ObsData(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%obsErrs(:,i), ObsErrs(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%tb_inv(:), InvData(:, i))
          CALL MPI_GATHER_Interface(state%sg, this%tb_clr(:), ClrData(:, i))
          ! ObsData(:, i) = this%ObsData(:,i)
          ! ObsErrs(:, i) = this%obsErrs(:,i)
          ! InvData(:, i) = this%tb_inv(:)
          ! ClrData(:, i) = this%tb_clr(:)

        END IF

        ! CALL destroy(this)

      END DO
      DEALLOCATE(tb_bak503, tb_inv)
      
      ! For satellite data, data fall between tt1 and tt2 are all processed at tt2.
      ! That is, if tt1 < t <= tt2, then t(ObsData) = tt2.
      ! This is for wavelet convenience.
      BLOCK 
        INTEGER :: numObsEachSlot
        REAL(r_kind), ALLOCATABLE    :: ObsDataEachSlot(:,:)
        INTEGER(i_kind), ALLOCATABLE :: ObsTimeEachSlot(:) 
        REAL(r_kind), ALLOCATABLE    :: olatlonEachSlot(:,:)
        REAL(r_kind), ALLOCATABLE    :: ObsData_rawEachSlot(:,:)
        REAL(r_kind), ALLOCATABLE    :: ObsErrsEachSlot(:,:)
        REAL(r_kind), ALLOCATABLE    :: CaMeanEachSlot(:,:)
        REAL(r_kind), ALLOCATABLE    :: INVDataEachSlot(:,:)
        REAL(r_kind), ALLOCATABLE    :: ClrDataEachSlot(:,:)
        INTEGER(i_kind), ALLOCATABLE :: IdxEachSlot(:), IdxOut(:)
        INTEGER :: it, tt1, tt2, iobs
        REAL(r_kind) :: t1, t2, t3, t4
        CHARACTER(LEN=3) :: tSlotsStr
        CHARACTER(LEN=200) :: filenameEachSlot

        CALL cpu_time(t1)
        ! ALLOCATE(IdxEachSlot(this%numObs))
        ALLOCATE(IdxEachSlot(TotalObs))

        DO it = 0, state%sg%tSlots -1 
          IdxEachSlot = 0

          tt2 = state%sg%tt(it+1)
          IF (it .EQ. 0) THEN
            tt1 = tt2 - 1
          ELSE
            tt1 = state%sg%tt(it)
          END IF

          ! DO iobs = 1, this%numObs
          DO iobs = 1, TotalObs
            ! IF (this%ObsTime(iobs) .GT. tt1 .AND. this%ObsTime(iobs) .LE. tt2) IdxEachSlot(iobs) = iobs
            IF (ObsTime(iobs) .GT. tt1 .AND. ObsTime(iobs) .LE. tt2) IdxEachSlot(iobs) = iobs
          END DO

          numObsEachSlot = COUNT(IdxEachSlot .NE. 0)
          ALLOCATE(IdxOut(numObsEachSlot))

          IdxOut = PACK(IdxEachSlot, IdxEachSlot .NE. 0)
          PRINT *, 'Num of OBS for each time slot is: ', it, numObsEachSlot, TotalObs, SIZE(IdxOut), MAXVAL(IdxOut), MINVAL(IdxOut)

          ObsDataEachSlot = ObsData(IdxOut, :)
          ! ObsTimeEachSlot = this%ObsTime(IdxOut)
          ObsTimeEachSlot = ObsTime(IdxOut)
          ObsTimeEachSlot = tt2
          ! olatlonEachSlot = this%olatlon(:, IdxOut)
          olatlonEachSlot = olatlon(:, IdxOut)
          ObsData_rawEachSlot = ObsData_raw(IdxOut, :)
          ObsErrsEachSlot = ObsErrs(IdxOut, :)
          IF (ALLOCATED(this%ca_mean)) CaMeanEachSlot = ca_mean(IdxOut, :)
          INVDataEachSlot = -INVData(IdxOut, :) ! transfer B-O to O-B
          ClrDataEachSlot = ClrData(IdxOut, :)
          PRINT *, 'ObsDataEachSlot: ', MAXVAL(obsTime), MINVAL(obsTime)

          IF (mpObs%isBaseProc() .AND. state%sg%gLevel == this%mgEnd) THEN
            ! filename = "/Users/yaliwu/Desktop/MOTOR/MOTOR/output/fy4_1_agri_after_QC_BC.nc"
            WRITE(tSlotsSTR, "(I1)") it+1
            filenameEachSlot = TRIM(filename) // '_'// TRIM(tSlotsStr) // '.nc'
            ! PRINT *, 'filenameEachSlot is: ', tSlotsSTR, filenameEachSlot
            CALL Open_NcFile_for_thinning(this%numVars, numObsEachSlot , olatlonEachSlot, ObsDataEachSlot(:, 1:4), ObsTimeEachSlot, filenameEachSlot, nc, dim3)
            ! print *, 'obstime = ', maxval(this%obstime), minval(this%obstime), maxval(ObsData), minval(ObsData)
            ! PRINT *, 'INVData = ', it, maxval(INVData), minval(INVData)
            ! PRINT *, 'Err = ', it, maxval(ObsErrsEachSlot), minval(ObsErrsEachSlot)
            ! PRINT *, 'num of obs out = ', SHAPE(ObsDataEachSlot)
            IF (ALLOCATED(CaMeanEachSlot)) THEN
              CALL Write_OBS_for_thinning(nc, dim3, numObsEachSlot, this%nchans, this%chan_lists, ObsDataEachSlot(:, 5:this%numVars), &
              ObsData_rawEachSlot(:, 5:this%numVars), ObsErrsEachSlot(:, 5:this%numVars), &
              INVDataEachSlot(:, 5:this%numVars), ClrDataEachSlot(:, 5:this%numVars), CaMeanEachSlot(:, 5:this%numVars))
            ELSE
              CALL Write_OBS_for_thinning(nc, dim3, numObsEachSlot, this%nchans, this%chan_lists, ObsDataEachSlot(:, 5:this%numVars), &
              ObsData_rawEachSlot(:, 5:this%numVars), ObsErrsEachSlot(:, 5:this%numVars), &
              INVDataEachSlot(:, 5:this%numVars), ClrDataEachSlot(:, 5:this%numVars))
            END IF
            CALL Close_NcFile_for_thinning(nc)

            ! 2023-08-09, Ting Shu, wavelet-based satellite data and error decomposition
            IF (this%useWaveFlag) THEN
              IF (UBOUND(ObsDataEachSlot, 1) .EQ. 0) THEN
                PRINT *, "There is no data for this satellite slot."
              ELSE
                PRINT *, "Start wavelet-based satellite decomposition."
                fn_len = LEN(TRIM(filenameEachSlot)) - 3
                coeff_fileName = filenameEachSlot(1 : fn_len)//"_wcoeff.nc"
                CALL CPU_TIME(t3)
                ! CALL this%sateDecom%decom(this%olatlon, this%ObsTime(1)*1.0D0, ObsData(:, 5:this%numVars), ObsErrs(:, 5:this%numVars), coeff_fileName)
                CALL this%sateDecom%decom(olatlonEachSlot, ObsTimeEachSlot(1)*1.0D0, ObsDataEachSlot(:, 5:this%numVars), ObsErrsEachSlot(:, 5:this%numVars), coeff_fileName, t3)
                CALL CPU_TIME(t4)
                ! 2023-10-18, TS, calculate running time
                PRINT *, "Decomposition running time of wavelet-based satellite thinning: ", tSlotsSTR, t4-t3
                PRINT *, "Wavelet-based satellite decomposition Done."
                PRINT *, "The wavelet coefficient nc file has been written to "//TRIM(coeff_fileName)
              END IF
            END IF

          END IF

          DEALLOCATE(IdxOut, ObsDataEachSlot, ObsTimeEachSlot, olatlonEachSlot)
          DEALLOCATE(ObsData_rawEachSlot, ObsErrsEachSlot, INVDataEachSlot, ClrDataEachSlot)
          IF (ALLOCATED(coeff_fileName)) DEALLOCATE(coeff_fileName)
          IF (ALLOCATED(CaMeanEachSlot)) DEALLOCATE(CaMeanEachSlot)
        END DO

        DEALLOCATE(IdxEachSlot)
        CALL cpu_time(t2)
        PRINT *, 'TIME COSTS OF RESETTING AT ANALYSIS SLOTS: ', t2-t1

      END BLOCK

      DEALLOCATE(ObsData, ObsData_raw, ObsErrs, INVData, clrData)
      IF (ALLOCATED(this%ca_mean)) DEALLOCATE(ca_mean)
      DEALLOCATE (pres, pres_BC, t, t_BC, q, q_BC, angles_BC, u, v, tskin, tskin_BC, &
       elevation, snowc, soiltype, landmask, PassDomainCheck)
      IF (this%CloudyMode) DEALLOCATE(qc, qr, qi, qs, qg)
    END IF

  END SUBROUTINE OutputForThinningBase

  SUBROUTINE RadBase_QC1(this, ichan, inst_name, ichn, qcnum, pres, qcloud)
    CLASS(RadBase_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: ichan
    CHARACTER(len=*), INTENT(IN) :: inst_name
    INTEGER(i_kind), INTENT(IN)  :: ichn
    REAL(r_kind), INTENT(IN), OPTIONAL :: pres(:,:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: qcloud(:,:)
    INTEGER(i_kind) :: iobs
    INTEGER(i_kind), INTENT(OUT) :: qcnum
    CHARACTER(len=20), ALLOCATABLE :: qcflag(:)
    REAL(r_kind), ALLOCATABLE :: satzenith(:), solzenith(:), satazi(:), solazi(:), relaz(:)
    INTEGER(i_kind) :: i
    REAL(r_kind) :: inv_threshold

    ALLOCATE(qcflag(SIZE(this%obsData,1))) 
    qcflag = 'good'

    ! 0. Initial check of missing values

    PRINT *, '=== QC_InitialCheck ==='
    CALL QC_InitialCheck(this%obsData(:, ichan), qcflag, ichn)

    ! inv_threshold = 10.0 ! default
    IF (.NOT. this%CloudyMode) THEN
      PRINT *, '=== Select_Clear ==='
      CALL Select_Clear(TRIM(inst_name),this%obsData(:, ichan),this%cloud_flag, ichn)
    !   inv_threshold = 10.0
    ! ELSE
    !   PRINT *, '=== ObsErrModel ==='
    !   CALL ObsErrModel(this%obsData(:,ichan), this%tb_inv(:), this%tb_clr(:), &
    !                   this%camin(ichan-4), this%camax(ichan-4), this%obsErrClr(ichan-4), &
    !                   this%obsErrCld(ichan-4), this%BTlim(ichan-4), this%obsErrs(:,ichan))
    !   inv_threshold = 100.0
    !   ! PRINT *, 'Check obserr parameters: ', this%camin, this%camax, this%obsErrClr, this%obsErrCld, &
    !   ! maxval(this%obsErrs(:,ichan)), minval(this%obsErrs(:,ichan))
    END IF

    ! 1. Absolute value check
    PRINT *, '=== QC_AbsoluteCheck ==='
    CALL QC_AbsoluteCheck(this%obsData(:, ichan), qcflag, ichn)

    ! 2. Exclude data in Southern hemisphere
    PRINT *, '=== QC_AbsoluteCheck ==='
    CALL QC_GeoCheck(this%obsData(:, ichan), this%olatlon(1, :), qcflag, ichn)

    ! ! 3. First guess (innovation) check
    ! PRINT *, '=== QC_FirstguessCheck ==='
    ! CALL QC_FirstguessCheck(this%obsData(:, ichan), this%tb_inv(:), this%obsErrs(:, ichan), inv_threshold, qcflag, ichn)

    ! ! 4. Surface type check
    ! PRINT *, '=== QC_SfctypeCheck ==='
    ! CALL QC_SfctypeCheck(this%obsData(:, ichan), this%sfctype, qcflag, ichn)

    ! 6. Zenith angle check
    PRINT *, '=== QC_ZenithCheck ==='
    satzenith = this%Angles(1, :) * radian2degree
    CALL QC_ZenithCheck(this%obsData(:, ichan), satzenith, qcflag, ichn)
    ! PRINT *, '=== QC_GlintCheck ==='
    ! Note to pass unit radian for cos/sin calculation
    ! satzenith = this%Angles(1, :)
    ! solzenith = this%Angles(3, :)
    ! satazi = this%Angles(2, :) * radian2degree
    ! solazi = this%Angles(4, :) * radian2degree
    ! relaz = Relative_Azimuth(solazi, satazi)
    ! CALL QC_GlintCheck(this%obsData(:, ichan), satzenith, solzenith, relaz, qcflag, ichn)

    ! 7. CLW check (MW)
    PRINT *, '=== QC_MWCLWCheck ==='
    IF (TRIM(inst_name) .EQ. 'mwts2' .OR. TRIM(inst_name) .EQ. 'mwhs2' .OR.  &
        TRIM(inst_name) .EQ. 'mwts3' .OR. TRIM(inst_name) .EQ. 'mwhs3') THEN
        ! TRIM(inst_name) .EQ. 'mwhs3') THEN
      IF (PRESENT(qcloud)) CALL QC_MWCLWPguessCheck(this%obsData(:, ichan), pres, qcloud, qcflag, ichn)
    END IF
    
    ! CALL QC_MWCLWCheck(this%obsData(:, ichan), this%olatlon(1, :), this%olatlon(2, :), this%cloud_flag, qcflag, ichn)
    ! PRINT *, '=== QC_MWCLWCheck ==='

    ! ! 7. LWP check (MW)
    ! PRINT *, '=== QC_MWTS3CLWCheck ==='
    ! satzenith = this%Angles(1, :)
    ! IF (TRIM(inst_name) .EQ. 'mwts3' ) &
    ! CALL QC_MWTS3LWPestimate(this%obsData(:, ichan),this%obsData_org(:, 5), this%obsData_org(:, 6),this%obsData_org(:, 15),satzenith, landmask,  qcflag, ichn)
    ! PRINT *, qcnum, ' observations out of ', SIZE(this%obsData,1), ' are removed by QC procedures'

    ! 8. Scattering Index check (MW)
    PRINT *, '=== QC_MWScatteringIndex ==='
    ! IF (TRIM(inst_name) .EQ. 'mwts2') &
    ! CALL QC_MWTS2_ScatteringIndex(this%obsData(:, ichan), this%obsData_org(:, 5), tb_bak503, landmask, qcflag, ichn)  
    IF (TRIM(inst_name) .EQ. 'mwhs2' .OR. TRIM(inst_name) .EQ. 'mwhs3') &
    CALL QC_MWHS2_ScatteringIndex(this%obsData(:, ichan),this%olatlon(1, :), this%olatlon(2, :), this%obsData_org(:, 5), this%obsData_org(:, 14), qcflag, ichn)
    ! CALL QC_MWScatteringIndex(this%obsData(:, ichan), this%SI, qcflag, ichn)

    ! ! 8.5 check: select_only_over_sea_channels
    ! PRINT *, '=== Select_only_over_sea_channels ==='
    ! CALL Select_only_over_sea_channels(this%obsData(:, ichan), this%sea(ichan-4), this%sfctype,  qcflag, ichn)

    ! 9. Last check: select_active_channels
    ! IF (ALLOCATED(this%oversea)) DEALLOCATE (this%oversea)
    ! PRINT *, '=== Select_sea_channels ==='
    ! ALLOCATE(this%oversea(this%nchans))
    ! ! ALLOCATE(pred_std(UBOUND(this%bias_predictors, 1)))
    ! this%oversea = this%oversea(ichn)
    ! PRINT *, this%oversea(ichan-4)
    ! PRINT *, 'Lower bound of oversea:', LBOUND(this%oversea)
    ! PRINT *, 'Upper bound of oversea:', UBOUND(this%oversea)
    ! PRINT *, 'this%oversea size:', SIZE(this%oversea)
    ! ! CALL Select_sea_channels(this%obsData(:, ichan), this%oversea(ichan-4), qcflag, ichn)

    ! 9. Last check: select_active_channels
    CALL Select_active_channels(this%obsData(:, ichan), this%ifuse(ichan-4), qcflag, ichn)

    qcnum = COUNT_MISS_REAL_1D(this%obsData(:,ichan))
    PRINT *, qcnum, ' observations out of ', SIZE(this%obsData,1), ' are removed by QC procedures'
    DEALLOCATE(qcflag, satzenith)
    IF (ALLOCATED(relaz)) DEALLOCATE(solzenith, satazi, solazi, relaz)

  END SUBROUTINE RadBase_QC1

  SUBROUTINE RadBase_QC2(this, ichan, inst_name, ichn, qcnum, landmask, tb_bak503, ch503)
    CLASS(RadBase_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: ichan
    CHARACTER(len=*), INTENT(IN) :: inst_name
    INTEGER(i_kind), INTENT(IN)  :: ichn, ch503
    REAL(r_kind), INTENT(IN) :: landmask(:), tb_bak503(:)
    INTEGER(i_kind) :: iobs
    INTEGER(i_kind), INTENT(OUT) :: qcnum
    CHARACTER(len=20), ALLOCATABLE :: qcflag(:)
    REAL(r_kind), ALLOCATABLE :: satzenith(:), solzenith(:), satazi(:), solazi(:), relaz(:)
    INTEGER(i_kind) :: i, istatus
    REAL(r_kind) :: inv_threshold

    ! PRINT *, 'check before QC: max/min obsdata = ', MAXVAL(this%obsData), MINVAL(this%obsData)

    ALLOCATE(qcflag(SIZE(this%obsData,1)))
    qcflag = 'good'

    inv_threshold = 10.0 ! default
    IF (.NOT. this%CloudyMode) THEN
      ! PRINT *, '=== Select_Clear ==='
      ! CALL Select_Clear(TRIM(inst_name),this%obsData(:, ichan),this%cloud_flag, ichn)
      inv_threshold = 10.0
    ELSE
      PRINT *, '=== ObsErrModel ==='
      CALL ObsErrModel(this%obsData(:,ichan), this%tb_inv(:), this%tb_clr(:), &
                      this%camin(ichan-4), this%camax(ichan-4), this%obsErrClr(ichan-4), &
                      this%obsErrCld(ichan-4), this%BTlim(ichan-4), this%obsErrs(:,ichan), this%ca_mean(:,ichan))
      inv_threshold = 100.0
      ! PRINT *, 'Check obserr parameters: ', this%camin, this%camax, this%obsErrClr, this%obsErrCld, &
      ! maxval(this%obsErrs(:,ichan)), minval(this%obsErrs(:,ichan))
    END IF

    ! 4. Surface type check
    PRINT *, '=== QC_SfctypeCheck ==='
    CALL QC_SfctypeCheck(this%obsData(:, ichan), this%sfctype, qcflag, ichn)
    
    ! ! ! 5. Land/sea check
    PRINT *, '=== QC_LandseaCheck ==='

    IF (this%oversea(ichan-4)>-1) THEN  
      CALL QC_LandseaCheck(this%obsData(:, ichan), this%sfctype, qcflag, ichn)  
    END IF  

    ! IF (ichn<7) &
    !   ! PRINT *, only_over_sea(ichn)
    !   CALL QC_LandseaCheck(this%obsData(:, ichan), this%sfctype, qcflag, ichn)

    ! 7. CLW check (MW)
    PRINT *, '=== QC_MWTS3CLWCheck ==='
    satzenith = this%Angles(1, :)
    IF (TRIM(inst_name) .EQ. 'mwts3' ) &
    CALL QC_MWTS3LWPestimate(this%obsData(:, ichan),this%obsData_org(:, 5), this%obsData_org(:, 6),this%olatlon(1, :), this%olatlon(2, :),satzenith, landmask,  qcflag, ichn)


   ! 8. Scattering Index check (MW)
    PRINT *, '=== QC_MWHSScatteringIndex ==='
    satzenith = this%Angles(4, :)
    IF (TRIM(inst_name) .EQ. 'mwhs2' .OR. TRIM(inst_name) .EQ. 'mwhs3') &
    CALL QC_MWHS_ScatteringIndex(this%obsData(:, ichan),this%olatlon(1, :), this%olatlon(2, :), this%obsData_org(:, 5), this%obsData_org(:, 14), satzenith, landmask, qcflag, ichn)

    ! 8. Scattering Index check (MW)
    PRINT *, '=== QC_MWTS_ScatteringIndex ==='
    IF (TRIM(inst_name) .EQ. 'mwts2' .OR. TRIM(inst_name) .EQ. 'mwts3') &
    CALL QC_MWTS2_ScatteringIndex(this%obsData(:, ichan), this%obsData_org(:, ch503), tb_bak503, landmask, qcflag, ichn)   

    ! 3. First guess (innovation) check
    PRINT *, '=== QC_FirstguessCheck ==='
    CALL QC_FirstguessCheck(this%obsData(:, ichan), this%tb_inv(:), this%obsErrs(:, ichan), inv_threshold, qcflag, ichn)

    DEALLOCATE(qcflag)

  END SUBROUTINE RadBase_QC2

  SUBROUTINE RadBase_BC(this, ichan, X, TotalObs, PassDomainCheck, pres, t, q, tskin, &
                        pres_BC, t_BC, q_BC, tskin_BC, angles_BC, &
                        ca_mean, tb_obs, ca_mean_BC, tb_obs_BC, npred, BCOEF, SCOEF, iscanpos)
    CLASS(RadBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: ichan
    TYPE(State_t),  INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: TotalObs
    LOGICAL, INTENT(IN), OPTIONAL       :: PassDomainCheck(:)
    REAL(r_kind), INTENT(IN), OPTIONAL  :: pres(:,:), t(:,:), q(:,:), tskin(:)
    REAL(r_kind), INTENT(IN), OPTIONAL  :: pres_BC(:, :), t_BC(:, :), q_BC(:, :), tskin_BC(:), angles_BC(:,:)
    REAL(r_kind), INTENT(IN), OPTIONAL  :: ca_mean(:), tb_obs(:), ca_mean_BC(:), tb_obs_BC(:)
    INTEGER(i_kind),INTENT(IN), OPTIONAL :: npred
    INTEGER(i_kind),INTENT(IN), OPTIONAL :: iscanpos(:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: BCOEF(:), SCOEF(:,:)
    REAL(r_kind), ALLOCATABLE :: zpred(:,:)
    REAL(r_kind), ALLOCATABLE :: bias(:)
    INTEGER(i_kind) :: iobs, ivar, i, igrid, itime, nlevs, ipred, istatus
    ! REAL(r_kind), ALLOCATABLE :: ObsBias(:)
    REAL(r_kind) :: ObsBias
    REAL(r_kind), ALLOCATABLE :: pred(:,:), beta(:,:), pred_mean(:), pred_std(:), mean(:), rms(:)
    CHARACTER(len=100) :: outputDir

    ! ALLOCATE (ObsBias(this%numVars))
    ObsBias = ZERO

    PRINT *, 'USE BC scheme of ', TRIM(this%bias_scheme), this%numVars
    SELECT CASE (TRIM(this%bias_scheme))
    CASE ('pdf_mode')
      PRINT *, 'START CALL Get_pdf_mode'
      CALL Get_pdf_mode(X%sg, this%nbins, this%tb_inv, ObsBias)
      PRINT *, 'START ASSIGN BIAS VALUES'
      ! print *, 'check before bias correction ', maxval(this%obsData(:,ichan)), minval(this%obsData(:,ichan))
      ! ObsBias(ichan) = 0.0D0 ! For debug only
      PRINT *, 'Check bias values: ', ObsBias
      PRINT *, 'START TO correct the observation'
      WHERE((this%obsData(:, ichan) .GT. 100.0D0 .AND. this%obsData(:, ichan) .LT. 350.0)) &
      this%obsData(:, ichan) = this%obsData(:, ichan) + ObsBias
      ! tb_inv = B - O
      ! bias = pdf model of tb_inv
      ! The bias of O is actually related to B. That is, 
      ! if bias>0, O is colder than B, Y should be corrected as Y+|bias|.
      ! if bias<0, O is warmer than B, Y should be corrected as Y-|bias|.
      ! Thus, 
      ! Y_corrected = Y + bias
      PRINT *, 'pdf_mode is successfully run'
    CASE ('ScanAirs')
      IF (PRESENT(BCOEF) .AND. PRESENT(npred)) THEN
      
        nlevs = X%sg%vLevel

        ALLOCATE(zpred(this%numObs,npred))
        ALLOCATE(bias(this%numObs))
        zpred = 0.0D0
        bias = 0.0D0

        IF (present(pres) .AND. present(t) .AND. present(q) .AND. present(tskin)) THEN 
          IF (PRESENT(ca_mean) .AND. PRESENT(tb_obs)) THEN
            CALL Biasprep(X, PassDomainCheck, pres, t, q, tskin, pres_BC, t_BC, q_BC, tskin_BC, angles_BC,  &
              ca_mean=ca_mean, tb_obs=tb_obs, ca_mean_BC=ca_mean_BC, tb_obs_BC=tb_obs_BC, &
              num_rad=this%numObs, TotalObs=TotalObs, npred=npred, &
              njplev=X%sg%vLevel, olatlon=this%olatlon, obsHght=this%obsHght, obsTime=this%obsTime, &
              Angles=this%Angles, norm_pred=this%norm_pred, zpred=zpred)
          ELSE
            CALL Biasprep(X, PassDomainCheck, pres, t, q, tskin, pres_BC, t_BC, q_BC, tskin_BC, angles_BC,  &
              num_rad=this%numObs, TotalObs=TotalObs, npred=npred, &
              njplev=X%sg%vLevel, olatlon=this%olatlon, obsHght=this%obsHght, obsTime=this%obsTime, &
              Angles=this%Angles, norm_pred=this%norm_pred, zpred=zpred)
          END IF
        ELSE
          zpred = ZERO
        END IF

        ! Note to use a single processor or serial run. Otherwise, the out file will only be part of the whole domain
        IF ( this%use_cloudy_predictors .AND. ALLOCATED(this%ca_mean) .AND. X%sg%mpddInfo_sg%nproc == 1) THEN
        BLOCK
          ! WRITE OUT predictors for offline diagnose
            USE mo_netcdf, only: NcDataset, NcDimension, NcVariable, NcGroup
            ! Refer to https://github.com/schaefed/mo_netcdf for usage
            TYPE(NcDataset)   :: nc
            TYPE(NcDimension) :: dim1, dim2
            TYPE(NcVariable)  :: var
            TYPE(NcGroup)     :: grp
            INTEGER :: Npixels, Nchans
            REAL(r_kind), ALLOCATABLE    :: out_array(:)
            CHARACTER(len=100) :: filename = './predictors.nc'
            INTEGER :: out_chidx
            CHARACTER(len=10) :: tmpstr
        
            ! IF (.not. sg%isBaseProc()) RETURN
            out_chidx = this%chan_lists(ichan-4)
            IF (out_chidx < 10) THEN
              WRITE(tmpstr, "(I1)") out_chidx
            ELSE 
              IF (out_chidx < 100) THEN
                WRITE(tmpstr, "(I2)") out_chidx
              ELSE 
                IF (out_chidx < 1000) THEN
                  WRITE(tmpstr, "(I3)") out_chidx
                ELSE
                  IF (out_chidx < 10000) THEN
                    WRITE(tmpstr, "(I4)") out_chidx
                  END IF
                END IF
              END IF
            END IF
            filename = './predictors_'//TRIM(tmpstr)//'.nc'
            
            nc = NcDataset(TRIM(filename), "w")
            Npixels = this%numobs
            dim1 = nc%setDimension("Npixels", Npixels)
            ALLOCATE(out_array(Npixels))
        
            ! Reverse latitude for plotting
            var = nc%setVariable("lat", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            WHERE (this%olatlon(1,:) .LT. missing-1.0)  
              out_array = this%olatlon(1,:) * radian2degree
            ELSEWHERE 
              out_array = 9999.0
            END WHERE
            ! print *, 'lat = ', maxval(out_array), minval(out_array), maxval(this%latitude), minval(this%latitude)
            call var%setData(out_array)
        
            var = nc%setVariable("lon", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            WHERE ( this%olatlon(2,:) .LT. missing-1.0 )  
              out_array = this%olatlon(2,:) * radian2degree
            ELSEWHERE
              out_array = 9999.0
            END WHERE
            ! print *, 'lon = ', maxval(out_array), minval(out_array), maxval(this%longitude), minval(this%longitude)
            call var%setData(out_array)
          
            var = nc%setVariable("Norm_Thickness1000-300hPa", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            ! print *, 'tbb02 = ', maxval(this%ObsValue(:, 2)), minval(this%ObsValue(:, 2))
            call var%setData(zpred(:,2))
        
            var = nc%setVariable("Norm_Thickness200-50hPa", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,3))
        
            var = nc%setVariable("Norm_Tskin", "f32", (/dim1/)) 
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,4))
        
            var = nc%setVariable("Norm_TotalPWater", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,5))
        
            var = nc%setVariable("Norm_Satzen", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,6))

            var = nc%setVariable("Norm_Satzen2", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,7))

            var = nc%setVariable("Norm_Satzen3", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,8))

            var = nc%setVariable("Norm_Ca_mean", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,9))

            var = nc%setVariable("Norm_Tb_obs", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(zpred(:,10))

            var = nc%setVariable("Norm_INV", "f32", (/dim1/))
            call var%setFillValue(9999.0)
            call var%setData(-this%tb_inv) ! transfer B-O to O-B
        
            ! close the dataset
            call nc%close()
            DEALLOCATE(out_array)

        END BLOCK
        END IF

        IF ( PRESENT (SCOEF)) THEN
          CALL GetBias_ScanAirs(this%numObs, PassDomainCheck, npred, zpred,iscanpos,this%olatlon(1,:), BCOEF, bias, SCOEF)
        ELSE 
          CALL GetBias_ScanAirs(this%numObs, PassDomainCheck, npred, zpred,iscanpos,this%olatlon(1,:), BCOEF, bias)
        END IF

        PRINT *, 'Check bias values: ', shape(bias), maxval(bias), minval(bias)
        PRINT *, 'FINISH GetBias_ScanAirs'
      ELSE 
        bias = 0.0D0
        PRINT *, 'WARNING: ScanAirs BC is not implemented'
      END IF
      ! CASE ('static_BC')
      !   CALL static_BC(this)
      ! CASE ('var_BC')
      !   CALL var_BC(this)
      WHERE((this%obsData(:, ichan) .GT. 100.0D0 .AND. this%obsData(:, ichan) .LT. 350.0)) &
      this%obsData(:, ichan) = this%obsData(:, ichan) + bias(:)
      ! tb_inv = B-O; the output innovation file is from O-B
      ! The bias of O is actually related to B. That is, 
      ! if bias>0, O is colder than B, Y should be corrected as Y+|bias|.
      ! if bias<0, O is warmer than B, Y should be corrected as Y-|bias|.
      ! Thus, 
      ! Y_corrected = Y + bias
      ! I can NOT recall why I commented the below line. Looks like it was wrong.
      ! "This is different from the pdf_norm method, which comes from B-O statistics."
      
      IF ( ALLOCATED(zpred) ) DEALLOCATE (zpred)
      IF ( ALLOCATED(bias) ) DEALLOCATE (bias)
      PRINT *, 'FINISH ScanAirs'

    CASE ('static_BC')
      print *, 'START CALL BC_predictors of static_BC'
      ! Allocation
      nlevs = X%sg%vLevel
      ALLOCATE(pred(this%numobs, SIZE(this%bias_predictors, 1))) ! bias predictors
      ALLOCATE(beta(this%nchans, SIZE(this%bias_predictors, 1))) ! coefficients of predictors for linear regression
      ALLOCATE(pred_mean(UBOUND(this%bias_predictors, 1)))
      ALLOCATE(pred_std(UBOUND(this%bias_predictors, 1)))
      ALLOCATE(mean(SIZE(this%bias_predictors, 1)))
      ALLOCATE(rms(SIZE(this%bias_predictors, 1)))

      pred = ZERO; beta = ZERO; pred_mean = ZERO; pred_std = ZERO; mean = ZERO; rms = ZERO

      IF (present(pres) .AND. present(t) .AND. present(q) .AND. present(tskin)) THEN 
      ! Step 01: calculate predictors. Predictors are specified from YAML.
        CALL BC_predictors(this%numobs, this%configFile, X, pres, t, q, tskin, this%Angles(1, :), pred)
      ELSE
        pred = ZERO
      END IF
      DO ipred = 1, SIZE(this%bias_predictors, 1)
        ! TODO: exclude halo???
        IF (COUNT(pred(:, ipred) .GT. ZERO) .GT. ZERO) THEN
          mean(ipred) = SUM(pred(:, ipred)) / REAL(COUNT(pred(:, ipred) .GT. ZERO))
        ELSE 
          mean(ipred) = ZERO
        END IF
        ! CALL MPI_GATHER_Interface(X%sg, mean(ipred:ipred), pred_mean(ipred:ipred))
      END DO

      ! Step 02: calculate statistics of each predictor.
      DO ipred = 1, SIZE(this%bias_predictors, 1)
        ! TODO: exclude halo???
        rms(ipred) = SUM((pred(:, ipred) - pred_mean(ipred))**2)
        ! CALL MPI_GATHER_Interface(X%sg, rms(ipred:ipred), pred_std(ipred:ipred))
        IF (COUNT(pred(:, ipred) .GT. ZERO) .GT. ZERO) THEN
          pred_std(ipred) = SQRT(pred_std(ipred)) / REAL(COUNT(pred(:, ipred) .GT. ZERO))
        ELSE 
          pred_std(ipred) = ZERO
        END IF
      END DO
      
      ! Step 03: normalize each predictor according to the mean and std values.
      DO ipred = 1, SIZE(this%bias_predictors, 1)
        IF (pred_std(ipred) .GT. ZERO) THEN
          pred(:,ipred) = (pred(:,ipred) - pred_mean(ipred)) / pred_std(ipred)
        ELSE 
          pred(:,ipred) = ZERO
        END IF
      END DO

      ! Step 04: estimate coefficients of predictors for linear regression
      ! For the first predictor (constant value, 1), mode is the coefficient of the predictor (1)
      beta = ZERO
      PRINT *, 'Check beta shape ', SHAPE(beta)
      ! CALL Get_pdf_mode(X%sg, this%nbins, this%tb_inv, ObsBias) 
      ! beta(:,1) = ObsBias(:)
      ! CALL linear_least_squares(n, x, y, a, b, square_error)

      ! Step 05: estimate biases of each channel according to linear regression.
      CALL static_BC(X%sg, pred, beta, ObsBias)

      IF(ALLOCATED(pred)) DEALLOCATE(pred, pred_mean, pred_std, mean, rms)
    ! CASE ('var_BC')
      ! CALL var_BC(this)
    CASE default
      CALL Get_pdf_mode(X%sg, this%nbins, this%tb_inv, ObsBias)
      this%obsData(:, ichan) = this%obsData(:, ichan) + ObsBias
    END SELECT

    IF (X%sg%isBaseProc() ) &
    ! IF (SIZE(this%obsData(:, 1),1) > 0) &
    ! PRINT *, 'On sg', X%sg%gLevel, ': ObsBias = ', ObsBias
    
    ! PRINT *, 'Check after BC: ', MAXVAL(this%obsData), MINVAL(this%obsData)

    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)

    ! DEALLOCATE (ObsBias)
    PRINT *, 'RadBase_BC is successfully run'

  END SUBROUTINE RadBase_BC

  SUBROUTINE TB_bkgAtobs(this, X, rttov_nl_sp_obsloc, pres, t, q, u, v, tskin, elevation, &
    snowc, soiltype, landmask, PassDomainCheck, iv, &
    qc, qr, qi, qs, qg)
    CLASS(RadBase_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(rttov_nl_sp_obsloc_t), INTENT(IN) :: rttov_nl_sp_obsloc
    INTEGER(i_kind), INTENT(IN) :: iv
    REAL(r_kind), INTENT(IN) :: pres(:,:), t(:,:), q(:,:), u(:,:), v(:,:)
    REAL(r_kind), INTENT(IN) :: tskin(:), elevation(:), snowc(:), soiltype(:), landmask(:)
    LOGICAL, INTENT(IN)      :: PassDomainCheck(:)
    REAL(r_kind), INTENT(IN), OPTIONAL :: qc(:,:), qr(:,:), qi(:,:), qs(:,:), qg(:,:)
    INTEGER(i_kind) :: iobs, ivar, ichan, io, i, j, k
    REAL(r_kind) :: forward
    REAL(r_kind) :: BT, BTclr

    PRINT *, 'START TB_bkgAtobs'
    IF (this%numVars < 1) RETURN

    print *, 'check TB_bkgAtobs valid: ', COUNT(PassDomainCheck)
    IF (.NOT. ALLOCATED(this%tb_inv)) ALLOCATE (this%tb_inv(this%numObs))
    IF (.NOT. ALLOCATED(this%tb_clr)) ALLOCATE (this%tb_clr(this%numObs))
    IF (.NOT. ALLOCATED(this%sfctype)) ALLOCATE (this%sfctype(this%numObs))

    this%sfctype = 'land'
    this%tb_inv = missing 
    this%tb_clr = missing

    DO io=1,this%numObs

      IF (PassDomainCheck(io)) THEN

        IF (this%obsData(io,iv) .GE. ABS(missing) - 1.0) CYCLE
        this%tb_inv(io) = ZERO
        this%tb_clr(io) = ZERO

        IF (snowc(io) > 1.0E-5) this%sfctype(io) = 'snow'
        IF (landmask(io) < 1.0E-5) this%sfctype(io) = 'sea' ! landmask = 0, sea; landmask = 1, land

        IF (iv .LE. 4) THEN 
          this%tb_inv(io) = ZERO
          this%tb_clr(io) = ZERO
        ELSE 
          IF (PRESENT(qc) .AND. PRESENT(qi)) THEN
            CALL GetForwardValue_obsloc(this, rttov_nl_sp_obsloc, X%sg%vLevel, this%rttov_chan_lists(iv-4), pres(io,:), t(io,:), q(io,:), &
                u(io,:), v(io,:), tskin(io), landmask(io), &
                elevation(io), soiltype(io), snowc(io), &
                this%olatlon(1, io), this%olatlon(2, io), this%Angles(:,io), BT, BTclr, &
                qc(io,:), qr(io,:), qi(io,:), qs(io,:), qg(io,:))
          ELSE
            CALL GetForwardValue_obsloc(this, rttov_nl_sp_obsloc, X%sg%vLevel, this%rttov_chan_lists(iv-4), pres(io,:), t(io,:), q(io,:), &
                u(io,:), v(io,:), tskin(io), landmask(io), &
                elevation(io), soiltype(io), snowc(io), &
                this%olatlon(1, io), this%olatlon(2, io), this%Angles(:,io), BT, BTclr)
          END IF
          forward = BT

          this%tb_inv(io) = forward
          this%tb_clr(io) = BTclr

          ! IF ( forward < 150.0 ) PRINT *, 'forward = ', forward
        END IF

      ELSE 
        ! this%tb_inv(io) = this%obsData_avg
        ! this%tb_clr(io) = this%obsData_avg
        this%tb_inv(io) = missing
        this%tb_clr(io) = missing
        !    this%tb_inv(io, 1) = missing
      END IF

    END DO

    !WHERE (this%obsData .GE. ABS(missing) - 1.0 .OR. ABS(this%tb_inv) .LE. 100.0)
    WHERE (this%obsData(:,iv) .GE. ABS(missing) - 1.0 .OR. ABS(this%tb_inv) .LE. 200.0 .OR. this%tb_inv .GE. ABS(missing) - 1.0)
      this%tb_inv = missing
      this%tb_clr = missing
    ELSEWHERE
      this%tb_inv = this%tb_inv - this%obsData(:,iv) ! here, inv = B - O
    END WHERE

    print *, 'finish TB_bkgAtobs'

  END SUBROUTINE TB_bkgAtobs

  SUBROUTINE destroy(this)
    CLASS(RadBase_t) :: this
    REAL(r_kind) :: t1, t2

    IF (ALLOCATED(this%obsData)) DEALLOCATE(this%obsData)
    IF (ALLOCATED(this%obsErrs)) DEALLOCATE(this%obsErrs)
    IF (ALLOCATED(this%ca_mean)) DEALLOCATE(this%ca_mean)
    IF (ALLOCATED(this%olatlon)) DEALLOCATE(this%olatlon)
    IF (ALLOCATED(this%obsHght)) DEALLOCATE(this%obsHght)
    IF (ALLOCATED(this%obsTime)) DEALLOCATE(this%obsTime)
    IF (ALLOCATED(this%varNames)) DEALLOCATE(this%varNames)
    IF (ALLOCATED(this%Angles)) DEALLOCATE (this%Angles)
    IF (ALLOCATED(this%cloud_flag)) DEALLOCATE(this%cloud_flag)
    IF (ALLOCATED(this%obsData_org)) DEALLOCATE(this%obsData_org)

    IF (ALLOCATED(this%tb_inv)) DEALLOCATE(this%tb_inv)
    IF (ALLOCATED(this%tb_clr)) DEALLOCATE(this%tb_clr)
    IF (ALLOCATED(this%tb_bakAtmodel)) DEALLOCATE(this%tb_bakAtmodel)
    IF (ALLOCATED(this%sfctype)) DEALLOCATE(this%sfctype)
    IF (ALLOCATED(this%Thinned_Angles)) DEALLOCATE (this%Thinned_Angles)
    IF (ALLOCATED(this%Thinned_Angles_hidx)) DEALLOCATE(this%Thinned_Angles_hidx)
    IF (ALLOCATED(this%Thinned_Angles_tidx)) DEALLOCATE(this%Thinned_Angles_tidx)
    IF (ALLOCATED(this%bias_predictors)) DEALLOCATE(this%bias_predictors)
    IF (ALLOCATED(this%cloudy_predictors)) DEALLOCATE(this%cloudy_predictors)

    IF (ALLOCATED(this%chan_lists)) DEALLOCATE (this%chan_lists)
    IF (ALLOCATED(this%rttov_chan_lists)) DEALLOCATE (this%rttov_chan_lists)
    IF (ALLOCATED(this%ifuse)) DEALLOCATE (this%ifuse)
    IF (ALLOCATED(this%oversea)) DEALLOCATE (this%oversea)
    IF (ALLOCATED(this%Err_in)) DEALLOCATE (this%Err_in)
    IF (ALLOCATED(this%bias_in)) DEALLOCATE (this%bias_in)
    IF (ALLOCATED(this%camin)) DEALLOCATE (this%camin)
    IF (ALLOCATED(this%camax)) DEALLOCATE (this%camax)
    IF (ALLOCATED(this%obsErrClr)) DEALLOCATE (this%obsErrClr)
    IF (ALLOCATED(this%obsErrCld)) DEALLOCATE (this%obsErrCld)
    IF (ALLOCATED(this%BTlim)) DEALLOCATE (this%BTlim)
    IF (ALLOCATED(this%iscanpos)) DEALLOCATE (this%iscanpos)
    IF ( ALLOCATED(this%BCOEF) ) DEALLOCATE(this%BCOEF)
    IF ( ALLOCATED(this%SCOEF) ) DEALLOCATE(this%SCOEF)
    IF ( ALLOCATED(this%tb_bakAtmodel) ) DEALLOCATE(this%tb_bakAtmodel)

    ! 2023-08-30, Ting Shu
    CALL CPU_TIME(t1)
    IF(this%useWaveFlag) CALL this%sateDecom%destroy()
    CALL CPU_TIME(t2)
    PRINT *, "Destroy running time of wavelet-based satellite thinning:", t2-t1
  END SUBROUTINE

END MODULE RadBase_m
