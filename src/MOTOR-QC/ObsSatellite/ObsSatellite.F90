!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsSatellite
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu, Ting Shu
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!   Added wavelet-based satellite reconstruction by Ting Shu (shuting@gbamwf.com), 2023/08/11, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE ObsSatellite_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE satellite_utils_m
  USE Interp_utils_m
  USE Thinning_utils_m
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE EachChannel_m, ONLY: EachChannel_t
  ! 2023-08-11, Ting Shu
  USE satelliteRecon_m, ONLY: satelliteRecon_t
  USE Satellite_utils_m
  USE Satellite_QC_utils_m

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsSatellite_t
    TYPE(rttov_nl_sp_t) :: OprRTTOV
    TYPE(EachChannel_t), ALLOCATABLE :: EachChannel(:)
    INTEGER(i_kind) :: nchans
    INTEGER(i_kind) :: mgStart, mgEnd
    INTEGER(i_kind), ALLOCATABLE :: chan_lists(:)
    INTEGER(i_kind), ALLOCATABLE :: rttov_chan_lists(:)
    CHARACTER(LEN=1024) :: inst_name = 'agri', platform_name = 'fy4_1'
    ! 2023-08-11, Ting Shu
    LOGICAL :: useWaveFlag
    TYPE(satelliteRecon_t) :: satRecon

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => ObsSatellite_setup
    PROCEDURE, PUBLIC  :: ObsIngest => ObsSatellite_read
    PROCEDURE, PUBLIC  :: ObsQC => ObsSatellite_qc
    PROCEDURE, PUBLIC  :: ObsForward => ObsSatellite_fwd
    PROCEDURE, PUBLIC  :: ObsTangent => ObsSatellite_tgt
    PROCEDURE, PUBLIC  :: ObsAdjoint => ObsSatellite_adj

    PROCEDURE, PUBLIC :: ObsPrepareForSg
    PROCEDURE, PUBLIC :: ObsThinning
    FINAL :: destructor
  END TYPE ObsSatellite_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)

    CLASS(ObsSatellite_t) :: this
    TYPE(State_t) :: X

    INTEGER(i_kind) :: i, j, istatus

    this%numVars = this%nchans + 4
    CALL this%OprRTTOV%initialize(this%configFile, X, TRIM(this%inst_name), TRIM(this%platform_name))

    WHERE (X%fields(X%getVarIdx('qvapor'))%DATA < qv_limit) X%fields(X%getVarIdx('qvapor'))%DATA = qv_limit

  END SUBROUTINE ObsPrepareForSg

  SUBROUTINE ObsThinning(this, state, thinObs, mpObs, useBkg, vector)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    IMPLICIT NONE
    CLASS(ObsSatellite_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    LOGICAL, INTENT(IN) :: useBkg, vector
    TYPE(ObsSet_t) :: YAngles
    INTEGER :: i
    REAL(r_kind) :: t1, t2

    IF (state%sg%isBaseProc()) PRINT *, "START thinning satellite obs for Angles and TBB"

    IF (mpObs%isActiveProc()) THEN
      thinObs = ObsSet_t(this%configFile, mpObs)
      ALLOCATE (thinObs%ObsFields(this%numVars - 4))
      ! 2023-08-11, Ting Shu
      IF (this%useWaveFlag) THEN
        CALL CPU_TIME(t1)
        CALL this%satRecon%recon(state, thinObs, mpObs)
        CALL CPU_TIME(t2)
        PRINT *, "************Recon: running time:", t2 - t1
        PRINT *, 'Wavelet-based thinning is successfully done.'
      ELSE
        YAngles = ObsSet_t(this%configFile, mpObs)
        ALLOCATE (YAngles%ObsFields(4))

        DO i = 1, this%numVars

          WHERE (this%EachChannel(i)%obsData .GE. missing) this%EachChannel(i)%obsData = missing
          PRINT *, 'obsdata = ', MAXVAL(this%EachChannel(i)%obsData), MINVAL(this%EachChannel(i)%obsData)

          BLOCK
            TYPE(ObsSet_t) :: thinObsTemp

            CALL this%EachChannel(i)%ObsPrepareForSg(state)
            CALL this%EachChannel(i)%ObsInitial(this%configFile)

            IF (i .LE. 4) THEN ! For satellite angles's thinning
              PRINT *, this%EachChannel(i)%varnames
              CALL this%EachChannel(i)%ObsThinning(state, thinObsTemp, mpObs, .TRUE., .FALSE.)
              YAngles%ObsFields(i) = thinObsTemp%ObsFields(1)
            ELSE
              IF (ALL(ABS(this%EachChannel(i)%ObsData) > 400.0)) THEN
                ! IF all data are missing, no need to calculate tb_bakAtmodel.
                this%EachChannel(i)%tb_bakAtmodel = this%EachChannel(i)%obsData_avg
              ELSE
                ! Angles are hard to deal with, otherwise, we do NOT need the ThinningBase module.
                CALL SetObsAttrSAT(YAngles, this%EachChannel(i))
                CALL TB_bkgAtmodel(state, this%OprRTTOV, this%rttov_chan_lists(i - 4), this%EachChannel(i))
              END IF
              CALL this%EachChannel(i)%ObsThinning(state, thinObsTemp, mpObs, .TRUE., .FALSE.)
              thinObs%ObsFields(i - 4) = thinObsTemp%ObsFields(1)

              IF (SIZE(thinObsTemp%ObsFields(1)%values) > 0) THEN
                PRINT *, 'check errors 1: ', MAXVAL(thinObsTemp%ObsFields(1)%errors), MINVAL(thinObsTemp%ObsFields(1)%errors)
                thinObs%ObsFields(i - 4)%errors = thinObsTemp%ObsFields(1)%errors
              END IF
            END IF

          END BLOCK
          CALL this%EachChannel(i)%destroy_EachChannel()
        END DO
        CALL SetObsAttrSAT_Y(thinObs, YAngles, state)
        ! PRINT *, 'Check angles: ', SIZE(thinObs%ObsFields(1)%ObsAttrSat%zenangle),SIZE(thinObs%ObsFields(2)%ObsAttrSat%zenangle),SIZE(thinObs%ObsFields(3)%ObsAttrSat%zenangle), &
        ! SIZE(thinObs%ObsFields(4)%ObsAttrSat%zenangle), SIZE(thinObs%ObsFields(5)%ObsAttrSat%zenangle), SIZE(thinObs%ObsFields(6)%ObsAttrSat%zenangle)

        ! BLOCK
        !   USE Obs2State_m
        !   USE State2NC_m
        !   TYPE(State_t) :: AA
        !   CHARACTER(len=1024) :: outputDir, varName
        !   IF(yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
        !   AA = Obs2State_BaseTypeName(state%sg, YAngles)

        !   DO i = 1, SIZE(AA%Fields)
        !     AA%Fields(i)%data =  AA%Fields(i)%data / degree2radian
        !   END DO

        !   CALL Output_NC_State_AV(AA, TRIM(outputDir), "thinned_Angles", .true., .TRUE.)
        ! END BLOCK

        ! BLOCK
        !   USE Obs2State_m
        !   USE State2NC_m
        !   TYPE(State_t) :: YY
        !   CHARACTER(len=1024) :: outputDir, varname
        !   IF(yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

        !   YY = Obs2State_BaseTypeName(state%sg, thinObs_noBC)
        !   ! print *, 'WRITE OUT thinobs_noBC'
        !   CALL Output_NC_State_AV(YY, TRIM(outputDir), "obs_noBC", .true., .TRUE.)
        ! END BLOCK
      END IF
      PRINT *, 'Thinning is successfully done'
    END IF
  END SUBROUTINE ObsThinning

  SUBROUTINE ObsSatellite_qc(this)
    CLASS(ObsSatellite_t) :: this

  END SUBROUTINE ObsSatellite_qc

  SUBROUTINE ObsSatellite_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(ObsSatellite_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    LOGICAL :: istat
    INTEGER(i_kind) :: istatus
    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)

    this%configFile = configFile

    IF (yaml_get_var(TRIM(this%configFile), 'RTTOV', 'inst_name', inst_name) /= 0) STOP
    IF (yaml_get_var(TRIM(this%configFile), 'RTTOV', 'platform_name', platform_name) /= 0) STOP

    IF (PRESENT(idxFile)) THEN
      this%platform_name = TRIM(platform_name(idxFile))
      this%inst_name = TRIM(inst_name(idxFile))
    ELSE
      PRINT *, 'NO instrument is used'
      STOP
    END IF
    DEALLOCATE (platform_name, inst_name)

    CALL Get_rttov_chan_info(TRIM(this%platform_name)//"-"//TRIM(this%inst_name), this%nchans, this%chan_lists, this%rttov_chan_lists)
    this%numVars = this%nchans + 4

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)
    ! 2023-08-11, Ting Shu
    ! Initilize satellite reconstruction object
    istatus = yaml_get_var(TRIM(configFile), 'Wavelet', 'useWave', this%useWaveFlag)
    PRINT *, "This satellite data use wavelet for thinning: ", this%useWaveFlag

  END SUBROUTINE ObsSatellite_setup

  SUBROUTINE ObsSatellite_read(this, X)
    IMPLICIT NONE
    CLASS(ObsSatellite_t) :: this
    TYPE(State_t) :: X
    CHARACTER(len=200) :: filename
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), tb_inv(:, :), obsData(:, :), obsErrs(:, :)
    INTEGER(i_kind), ALLOCATABLE :: ObsTime(:)
    INTEGER(i_kind) :: ivar, istatus
    CHARACTER(len=1024) :: outputDir
    CHARACTER(LEN=3) :: tSlotsStr
    CHARACTER(LEN=200) :: filenameEachSlot
    REAL(r_kind) :: t1, t2

    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', outputDir)
    ! Read obs from MOTOR-DP
    filename = TRIM(outputDir)//'/'//TRIM(this%platform_name)//"_"//TRIM(this%inst_name)//"_after_QC_BC"
    CALL Read_OBS_from_DP(TRIM(filename), X%sg%tSlots, this%chan_lists, &
                          this%numVars, this%numObs, olatlon, ObsTime, obsData, obsErrs, tb_inv)

    ! update nchans to avoid inconsistency b.t. nchans and numVars.
    ! this may happen when obs file is absent
    this%nchans = this%numVars - 4
    PRINT *, 'nchans will be processed is ', this%nchans
    ! Conversions from MOTOR-DP to Xie's ObsThinning ObsBase type
    ALLOCATE (this%EachChannel(this%numVars))
    DO ivar = 1, this%numVars
      ALLOCATE (this%EachChannel(ivar)%ObsData(this%numObs, 1))
      ALLOCATE (this%EachChannel(ivar)%ObsErrs(this%numObs, 1))
      ALLOCATE (this%EachChannel(ivar)%obsHght(this%numObs))
      ALLOCATE (this%EachChannel(ivar)%land(this%numObs))
      ALLOCATE (this%EachChannel(ivar)%varNames(1))
      ALLOCATE (this%EachChannel(ivar)%radius(4, 1), this%EachChannel(ivar)%sizeInc(1), this%EachChannel(ivar)%qcThreshold(1))

      ! PRINT *, 'Check obsdata values: ',maxval(PACK(ObsData(:, ivar),ObsData(:, ivar) <400 )), minval(PACK(ObsData(:, ivar),ObsData(:, ivar) <400 ))
      this%EachChannel(ivar)%numObs = this%numObs
      this%EachChannel(ivar)%ObsData(:, 1) = obsData(:, ivar)
      this%EachChannel(ivar)%obsErrs(:, 1) = obsErrs(:, ivar)
      this%EachChannel(ivar)%olatlon = olatlon
      this%EachChannel(ivar)%obsTime = obsTime

      this%EachChannel(ivar)%radius = 5.0D2        ! Temporary holders
      this%EachChannel(ivar)%radius(3, :) = 5.0D2   ! Topography influence radius
      this%EachChannel(ivar)%radius(4, :) = 5.0D1   ! Temporal influence radius
      this%EachChannel(ivar)%sizeInc = 10.0D0 ! Temporary holders
      this%EachChannel(ivar)%correlation_threshold = 0.25D0   ! No influence under this threshold value
      this%EachChannel(ivar)%qcThreshold = 50.0D0  ! Temporary holders
      this%EachChannel(ivar)%land = missing
      this%EachChannel(ivar)%obsHght = -1.0E8

      this%EachChannel(ivar)%obsType = TRIM(this%platform_name)//'_'//TRIM(this%inst_name)
      IF (ivar .EQ. 1) this%EachChannel(ivar)%varNames = "sat_zenith"
      IF (ivar .EQ. 2) this%EachChannel(ivar)%varNames = "sat_azi"
      IF (ivar .EQ. 3) this%EachChannel(ivar)%varNames = "sol_zenith"
      IF (ivar .EQ. 4) this%EachChannel(ivar)%varNames = "sol_azi"
      IF (ivar .GT. 4) this%EachChannel(ivar)%varNames = 'tbb'
    END DO

    DEALLOCATE (obsData, obsErrs, tb_inv, olatlon, obsTime)

    ! 2023-08-11, Ting Shu, add satellite data loading code
    CALL CPU_TIME(t1)
    IF (this%useWaveFlag) CALL this%satRecon%ingest(TRIM(this%configFile), X%sg%tSlots, TRIM(this%platform_name), TRIM(this%inst_name))
    CALL CPU_TIME(t2)
    PRINT *, "************Ingest: running time:", t2 - t1, X%sg%gLevel

  END SUBROUTINE ObsSatellite_read

  FUNCTION ObsSatellite_fwd(this, X, O) RESULT(Y)
    CLASS(ObsSatellite_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs structure
    TYPE(ObsSet_t) :: Y
  END FUNCTION ObsSatellite_fwd

  FUNCTION ObsSatellite_tgt(this, dX, X) RESULT(Y)
    CLASS(ObsSatellite_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y
  END FUNCTION ObsSatellite_tgt

  FUNCTION ObsSatellite_adj(this, dY, X) RESULT(dX)
    CLASS(ObsSatellite_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY
  END FUNCTION ObsSatellite_adj

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(ObsSatellite_t), INTENT(INOUT) :: this
    INTEGER(i_kind) :: ivar

    IF (ALLOCATED(this%EachChannel)) THEN
    DO ivar = 1, this%numVars
      ! IF (ALLOCATED(this%EachChannel(ivar)%varNames)) DEALLOCATE(this%EachChannel(ivar)%varNames)
      ! IF (ALLOCATED(this%EachChannel(ivar)%obsData_org)) DEALLOCATE(this%EachChannel(ivar)%obsData_org)
      CALL this%EachChannel(ivar)%ObsDestroy()
    END DO
    END IF
    IF (ALLOCATED(this%EachChannel)) DEALLOCATE (this%EachChannel)
    IF (ALLOCATED(this%chan_lists)) DEALLOCATE (this%chan_lists)
    IF (ALLOCATED(this%rttov_chan_lists)) DEALLOCATE (this%rttov_chan_lists)
    ! 2023-08-11, Ting Shu
    IF (this%useWaveFlag) CALL this%satRecon%destroy()
    PRINT *, 'ObsSatellite destructor is successfully run'

  END SUBROUTINE destructor

END MODULE ObsSatellite_m
