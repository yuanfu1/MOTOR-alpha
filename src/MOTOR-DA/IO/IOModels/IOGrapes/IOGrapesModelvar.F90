!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-06-29, created by JiongmingPang for reading modelvar for BEC.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-06-29,, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE IOGrapesModelvar_m
  USE kinds_m
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE YAMLRead_m
  USE AdvanceTime_m
  USE GrapesModelvarIO_m, ONLY: GrapesModelvarIO_t
  USE InterpolateTime_m
  USE Interp1D_m
  USE Widgets_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE Interp1D_m
  USE conversions_m

  TYPE, EXTENDS(IOModel_t) :: IOGrapesModelvar_t
    CHARACTER(LEN=1024)                              :: nmlDAFilename, modvarCtlFn, hfFn
    CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE   :: modvarFns
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE     :: varNames
    INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE    :: mdate
    INTEGER(i_kind)                                  :: numFiles

  CONTAINS
    PROCEDURE, PUBLIC                :: initialize
    PROCEDURE, PRIVATE, PASS(this)   :: m_read_config_file
    PROCEDURE, PUBLIC, PASS(this)    :: m_read_bcg_into_Xm_Ens

  END TYPE IOGrapesModelvar_t

CONTAINS

  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOGrapesModelvar_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN)   :: geometry
    CHARACTER(LEN=1024), INTENT(IN)        :: configFile

    CALL this%initializeIOModel(configFile, geometry)
    this%nmlDAFilename = TRIM(configFile)
  END SUBROUTINE

  SUBROUTINE interpVerticalModelvar(sg, Source_realheight, Source_data, Target_data)
    IMPLICIT NONE
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN)    :: Source_realheight(:, :, :)   ! X x Y x Z
    REAL(r_kind), INTENT(IN)    :: Source_data(:, :, :, :)     ! X x Y x Z x T
    REAL(r_kind), INTENT(OUT)   :: Target_data(:, :, :)    ! vLevel x Horizontal x T
    INTEGER(i_kind)             :: shape_indata(4)
    INTEGER(i_kind)             :: i, j, z, t
    REAL(r_kind), ALLOCATABLE   :: Source_sigma(:)
    REAL(r_kind), ALLOCATABLE   :: Target_data_tmp(:, :, :, :)   ! X x Y x vLevel x T

    shape_indata = SHAPE(Source_data)
    ALLOCATE (Target_data_tmp(shape_indata(1), shape_indata(2), sg%vLevel, shape_indata(4)))

    PRINT *, "PJM---DEBUG, shape_Source_data:", shape_indata

    DO t = 1, shape_indata(4)
      DO i = 1, shape_indata(1)
        DO j = 1, shape_indata(2)
          CALL interp1d(Source_realheight(i, j, :), Source_data(i, j, :, t), sg%sigma, Target_data_tmp(i, j, :, t))
        END DO
      END DO
    END DO

    DO z = 1, sg%vLevel
      DO t = 1, shape_indata(4)
        Target_data(z, :, t) = RESHAPE(Target_data_tmp(:, :, z, t), (/shape_indata(1) * shape_indata(2)/))
      END DO
    END DO

    DEALLOCATE (Target_data_tmp)

  END SUBROUTINE interpVerticalModelvar

  SUBROUTINE m_read_config_file(this, ensIndx)
    CLASS(IOGrapesModelvar_t) :: this
    INTEGER(i_kind), INTENT(IN), OPTIONAL  :: ensIndx
    CHARACTER(LEN=1024), ALLOCATABLE       :: modvarFns(:)
    CHARACTER(LEN=20), ALLOCATABLE         :: inVarNames(:)
    CHARACTER(LEN=10), ALLOCATABLE         :: timestr(:)
    CHARACTER(LEN=1024)                    :: configFile
    CHARACTER(LEN=1024)                    :: inFileDir
    CHARACTER(LEN=1024)                    :: ensDirtmp
    CHARACTER(LEN=3)                       :: memid
    INTEGER(i_kind)                        :: t, inFilesNum, ifile, inTimeNum
    LOGICAL                                :: ensFlag

    ifile = yaml_get_var(this%m_configFile, 'BMat', 'ensinput_dir_model', inFileDir)
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensbakFileList', modvarFns)
    inFilesNum = UBOUND(modvarFns, 1)
    PRINT *, "inFilesNum: ", inFilesNum
    this%numFiles = inFilesNum
    ALLOCATE (this%modvarFns(inFilesNum))

    ! and get the full position of modelvar files
    IF (PRESENT(ensIndx)) THEN

      ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensFlag', ensFlag)
      IF (ensFlag) THEN
        ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensPath', ensDirtmp)
        WRITE (memid, '(A1,I2.2)') "m", ensIndx
        DO t = 1, inFilesNum
          PRINT *, 'modvarFns: ', TRIM(modvarFns(t))
          this%modvarFns(t) = TRIM(inFileDir)//"/"//TRIM(ensDirtmp)//"/"//memid//"/"//TRIM(modvarFns(t))
          PRINT *, 'modvarFns: ', TRIM(this%modvarFns(t))
        END DO
      END IF

    ELSE

      DO t = 1, inFilesNum
        PRINT *, 'modvarFns: ', TRIM(modvarFns(t))
        this%modvarFns(t) = TRIM(inFileDir)//"/"//TRIM(modvarFns(t))
        PRINT *, 'modvarFns: ', TRIM(this%modvarFns(t))
      END DO

    END IF

    ! get the full position of modelvar control files (*.ctl)
    this%modvarCtlFn = TRIM(inFileDir)//"/model.ctl"
    PRINT *, 'modvarCtlFn: ', TRIM(this%modvarCtlFn)

    ! get the full position of HeightField file
    this%hfFn = TRIM(inFileDir)//"/HeightField"
    PRINT *, 'hfFn: ', TRIM(this%hfFn)

    ! get the times of modelvar files
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensbakFileTime_modelvar', timestr)
    inTimeNum = UBOUND(modvarFns, 1)
    IF (inTimeNum .NE. inFilesNum) THEN
      PRINT *, "Error: inTimeNum is not consistent with inFilesNum, &
&      please check 'bakFileList' and 'bakFileTime_modelvar' in your yaml file!"
      STOP
    END IF

    ALLOCATE (this%mdate(6, inTimeNum))
    DO t = 1, inTimeNum
      READ (timestr(t) (1:4), "(i4)") this%mdate(1, t) !year
      READ (timestr(t) (5:6), "(i2)") this%mdate(2, t) !mon
      READ (timestr(t) (7:8), "(i2)") this%mdate(3, t) !day
      READ (timestr(t) (9:10), "(i2)") this%mdate(4, t) !hor
      this%mdate(5, t) = 0 !min
      this%mdate(6, t) = 0 !sec
      PRINT *, "GMT Time in t-step: ", t, this%mdate(:, t)
    END DO

    ! get the analysis variables
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'modelState', 'varList', this%varNames)

    IF (ALLOCATED(modvarFns)) DEALLOCATE (modvarFns)
    IF (ALLOCATED(timestr)) DEALLOCATE (timestr)

  END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm_Ens(this, Xm, sg, ensIndx)
    USE GrapesModelvarIO_m, ONLY: GrapesModelvarIO_t
    USE InterpolateTime_m
    USE ModelCoupler_m, ONLY: ModelCoupler_t

    IMPLICIT NONE

    CLASS(IOGrapesModelvar_t) :: this
    TYPE(State_t), INTENT(INOUT)            :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT)       :: sg
    TYPE(GrapesModelvarIO_t)                :: GrapesModelvarIO
    TYPE(ModelCoupler_t)                    :: ModelCoupler
    INTEGER(i_kind), INTENT(IN), OPTIONAL   :: ensIndx

    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    REAL(r_kind), ALLOCATABLE :: cellCntrDA(:, :)
    REAL(r_kind), ALLOCATABLE :: valueModel(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueModel_tmp(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueDA(:, :, :)

    INTEGER(i_kind), ALLOCATABLE :: TimeSlots_target(:)   ! for time-scale interpolation, Jiongming Pang

    INTEGER(i_kind) :: i, j, k, HorNumGridModel, nprocInterp, varIdx
    TYPE(InterpolateData_Time_t), ALLOCATABLE :: IntpTime

    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind)   :: iv, varIndx
    REAL(r_kind), ALLOCATABLE :: timeModel(:), timeDA(:)
    REAL(r_kind) :: u, v

    ASSOCIATE (sg_mp => sg)

      IF (sg_mp%isBaseProc()) THEN
        IF (PRESENT(ensIndx)) THEN
          CALL this%m_read_config_file(ensIndx)
        ELSE
          CALL this%m_read_config_file()
        END IF
      END IF

      CALL sg_mp%mpddInfo_sg%bCast(this%numFiles)
      IF (this%numFiles == 0) THEN
        PRINT *, 'No file listed in the configurations.'
        RETURN
      END IF

      IF (sg_mp%isBaseProc()) THEN
        GrapesModelvarIO = GrapesModelvarIO_t(this%modvarFns, this%modvarCtlFn, this%hfFn, this%mdate)

        ! get the time_slots_target
        ALLOCATE (TimeSlots_target(sg_mp%tSlots))
        CALL IntpLinearTimeSlots(GrapesModelvarIO%time_unix, TimeSlots_target)

        HorNumGridModel = GrapesModelvarIO%idn * GrapesModelvarIO%jdn

        ALLOCATE (cellCntrModel(2, HorNumGridModel), &
                  valueModel(sg_mp%vLevel, HorNumGridModel, sg_mp%tSlots), &
                  valueDA(sg_mp%vLevel, sg_mp%num_icell_global, sg_mp%tSlots), &
                  valueModel_tmp(sg_mp%vLevel, HorNumGridModel, GrapesModelvarIO%tdn))

        cellCntrModel(1, :) = RESHAPE(GrapesModelvarIO%lat2D, (/HorNumGridModel/))
        cellCntrModel(2, :) = RESHAPE(GrapesModelvarIO%lon2D, (/HorNumGridModel/))

      END IF

      CALL sg_mp%mpddInfo_sg%bCast(GrapesModelvarIO%kdn)

      ! Initialize the model coupler.
      CALL ModelCoupler%Initialize(this%m_configFile, this%geometry)
      CALL ModelCoupler%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)

      ! Ingest the topo and update the height in MOTOR-DA
      CALL ModelCoupler%ingest_to_topo_and_update_zHght_d(GrapesModelvarIO%topo, HorNumGridModel, sg_mp)

      IF (sg_mp%vLevel == 1) THEN
        !! Ingest the 2D background field

        IF (sg_mp%isBaseProc()) THEN
          timeModel = GrapesModelvarIO%time_unix
          timeDA = TimeSlots_target
        END IF

        ! Ingest t2m-temperature at 2m
        IF (Xm%getVarIdx('temp') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%t2m)) THEN
              valueModel_tmp = RESHAPE(GrapesModelvarIO%temp(:, :, 1, :), (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
              PRINT *, "WARNNING!!! There is no 't2m' in GrapesModelvar, use the first level temp instead! "
            ELSE
              valueModel_tmp = RESHAPE(GrapesModelvarIO%t2m, (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - temp has been loaded.', Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, 1, 1), &
            sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
        END IF

        ! Ingest uwnd at 10m
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%u10m)) THEN
              valueModel_tmp = RESHAPE(GrapesModelvarIO%uwnd(:, :, 1, :), (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
              PRINT *, "WARNNING!!! There is no 'u10m' in GrapesModelvar, use the first level uwnd instead! "
            ELSE
              valueModel_tmp = RESHAPE(GrapesModelvarIO%u10m, (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
            END IF
            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - u10m has been loaded.'
        END IF

        ! Ingest vwnd at 10m
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%v10m)) THEN
              valueModel_tmp = RESHAPE(GrapesModelvarIO%vwnd(:, :, 1, :), (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
              PRINT *, "WARNNING!!! There is no 'v10m' in GrapesModelvar, use the first level vwnd instead! "
            ELSE
              valueModel_tmp = RESHAPE(GrapesModelvarIO%v10m, (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - v10m has been loaded.'
        END IF

        ! Ingest q2m-qvapor at 2m
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%q2m)) THEN
              valueModel_tmp = RESHAPE(GrapesModelvarIO%qvapor(:, :, 1, :), (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
              PRINT *, "WARNNING!!! There is no 'q2m' in GrapesModelvar, use the first level qvapor instead! "
            ELSE
              valueModel_tmp = RESHAPE(GrapesModelvarIO%q2m, (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - q2m has been loaded.'
        END IF

        ! Ingest pressure at ground
        IF (Xm%getVarIdx('pres') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%psfc)) THEN
              valueModel_tmp = RESHAPE(GrapesModelvarIO%pres(:, :, 1, :), (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
              PRINT *, "WARNNING!!! There is no 'psfc' in GrapesModelvar, use the first level pres instead! "
            ELSE
              valueModel_tmp = RESHAPE(GrapesModelvarIO%psfc, (/1, HorNumGridModel, GrapesModelvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - psfc has been loaded.'
        END IF

        ! Ingest the u10min and v10min
        IF (Xm%getVarIdx('uwnd10min') .NE. 0 .AND. Xm%getVarIdx('uwnd') .NE. 0) THEN
          Xm%Fields(Xm%getVarIdx('uwnd10min'))%DATA = Xm%Fields(Xm%getVarIdx('uwnd'))%DATA
        END IF

        IF (Xm%getVarIdx('vwnd10min') .NE. 0 .AND. Xm%getVarIdx('vwnd') .NE. 0) THEN
          Xm%Fields(Xm%getVarIdx('vwnd10min'))%DATA = Xm%Fields(Xm%getVarIdx('vwnd'))%DATA
        END IF

        ! Ingest the u2min and v2min
        IF (Xm%getVarIdx('uwnd2min') .NE. 0 .AND. Xm%getVarIdx('uwnd') .NE. 0) THEN
          Xm%Fields(Xm%getVarIdx('uwnd2min'))%DATA = Xm%Fields(Xm%getVarIdx('uwnd'))%DATA
        END IF

        IF (Xm%getVarIdx('vwnd2min') .NE. 0 .AND. Xm%getVarIdx('vwnd') .NE. 0) THEN
          Xm%Fields(Xm%getVarIdx('vwnd2min'))%DATA = Xm%Fields(Xm%getVarIdx('vwnd'))%DATA
        END IF

      ELSE

        ! Ingest temp-temperature
        PRINT *, 'Ingestion of temp-temperature'
        IF (Xm%getVarIdx('temp') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%temp)) THEN
              PRINT *, "ERROR!!! There is no 'temp' in GrapesModelvarIO! "
              STOP
            END IF

            ! for vertical inter/extra-polation
            IF (.NOT. ALLOCATED(GrapesModelvarIO%topo) .OR. .NOT. ALLOCATED(GrapesModelvarIO%zRHght_s)) THEN
              PRINT *, "ERROR!!! There is no 'topo' or no 'zRHght_s' in GrapesModelvarIO! "
              STOP
            END IF

            CALL interpVerticalModelvar(sg_mp, GrapesModelvarIO%zRHght_s, GrapesModelvarIO%temp, valueModel_tmp(:, :, :))
            ! for time-scale interpolation
            ALLOCATE (IntpTime)
            CALL IntpTime%IntpInitialize(valueModel_tmp, GrapesModelvarIO%time_unix, TimeSlots_target)
            CALL IntpTime%IntpLinearTime(valueModel_tmp)
            valueModel = IntpTime%data_intperpolate
            DEALLOCATE (IntpTime)
          END IF

          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - temp has been loaded.'
        END IF

        ! Ingest uwnd-u wind
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%uwnd)) THEN
              PRINT *, "ERROR!!! There is no 'uwnd' in GrapesModelvarIO! "
              STOP
            END IF

            ! for vertical inter/extra-polation
            IF (.NOT. ALLOCATED(GrapesModelvarIO%topo) .OR. .NOT. ALLOCATED(GrapesModelvarIO%zRHght_s)) THEN
              PRINT *, "ERROR!!! There is no 'topo' or no 'zRHght_s' in GrapesModelvarIO! "
              STOP
            END IF
            CALL interpVerticalModelvar(sg_mp, GrapesModelvarIO%zRHght_s, GrapesModelvarIO%uwnd, valueModel_tmp(:, :, :))
            ! for time-scale interpolation
            ALLOCATE (IntpTime)
            CALL IntpTime%IntpInitialize(valueModel_tmp, GrapesModelvarIO%time_unix, TimeSlots_target)
            CALL IntpTime%IntpLinearTime(valueModel_tmp)
            valueModel = IntpTime%data_intperpolate
            DEALLOCATE (IntpTime)
          END IF

          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - uwnd has been loaded.'
        END IF

        ! Ingest vwnd-v wind
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%vwnd)) THEN
              PRINT *, "ERROR!!! There is no 'vwnd' in GrapesModelvarIO! "
              STOP
            END IF

            ! for vertical inter/extra-polation
            IF (.NOT. ALLOCATED(GrapesModelvarIO%topo) .OR. .NOT. ALLOCATED(GrapesModelvarIO%zRHght_s)) THEN
              PRINT *, "ERROR!!! There is no 'topo' or no 'zRHght_s' in GrapesModelvarIO!"
              STOP
            END IF
            CALL interpVerticalModelvar(sg_mp, GrapesModelvarIO%zRHght_s, GrapesModelvarIO%vwnd, valueModel_tmp(:, :, :))
            ! for time-scale interpolation
            ALLOCATE (IntpTime)
            CALL IntpTime%IntpInitialize(valueModel_tmp, GrapesModelvarIO%time_unix, TimeSlots_target)
            CALL IntpTime%IntpLinearTime(valueModel_tmp)
            valueModel = IntpTime%data_intperpolate
            DEALLOCATE (IntpTime)
          END IF

          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesModelvar - vwnd has been loaded.'
        END IF

        ! Ingest pres-pressure
        IF (Xm%getVarIdx('pres') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%pres)) THEN
              PRINT *, "ERROR!!! There is no 'pres' in GrapesModelvarIO! "
              STOP
            END IF

            ! for vertical inter/extra-polation
            IF (.NOT. ALLOCATED(GrapesModelvarIO%topo) .OR. .NOT. ALLOCATED(GrapesModelvarIO%zRHght_s)) THEN
              PRINT *, "ERROR!!! There is no 'topo' or no 'zRHght_s' in GrapesModelvarIO! "
              STOP
            END IF
            CALL interpVerticalModelvar(sg_mp, GrapesModelvarIO%zRHght_s, GrapesModelvarIO%pres, valueModel_tmp(:, :, :))
            ! for time-scale interpolation
            ALLOCATE (IntpTime)
            CALL IntpTime%IntpInitialize(valueModel_tmp, GrapesModelvarIO%time_unix, TimeSlots_target)
            CALL IntpTime%IntpLinearTime(valueModel_tmp)
            valueModel = IntpTime%data_intperpolate
            DEALLOCATE (IntpTime)
          END IF

          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
          Xm%Fields(Xm%getVarIdx('pres'))%DATA = Xm%Fields(Xm%getVarIdx('pres'))%DATA / 100.0D0
          PRINT *, 'IOGrapesModelvar - pres has been loaded.'
        END IF

        ! Ingest qvapor
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesModelvarIO%qvapor)) THEN
              PRINT *, "ERROR!!! There is no 'qvapor' in GrapesModelvarIO! "
              STOP
            END IF

            ! for vertical inter/extra-polation
            IF (.NOT. ALLOCATED(GrapesModelvarIO%topo) .OR. .NOT. ALLOCATED(GrapesModelvarIO%zRHght_s)) THEN
              PRINT *, "ERROR!!! There is no 'topo' or no 'zRHght_s' in GrapesModelvarIO! "
              STOP
            END IF
            CALL interpVerticalModelvar(sg_mp, GrapesModelvarIO%zRHght_s, GrapesModelvarIO%qvapor, valueModel_tmp(:, :, :))
            ! for time-scale interpolation
            ALLOCATE (IntpTime)
            CALL IntpTime%IntpInitialize(valueModel_tmp, GrapesModelvarIO%time_unix, TimeSlots_target)
            CALL IntpTime%IntpLinearTime(valueModel_tmp)
            valueModel = IntpTime%data_intperpolate
            DEALLOCATE (IntpTime)
          END IF

          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)

          varIdx = Xm%getVarIdx('qvapor')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%DATA(i, j, k) <= 1.0D-7) &
            Xm%Fields(varIdx)%DATA(i, j, k) = 1.0D-7

          Xm%Fields(Xm%getVarIdx('qvapor'))%DATA = Xm%Fields(Xm%getVarIdx('qvapor'))%DATA * 1000.0D0
          PRINT *, 'IOGrapesModelvar - qvapor has been loaded.'
        END IF

      END IF

      CALL sg_mp%mpddInfo_sg%barrier

    END ASSOCIATE

    IF (ALLOCATED(cellCntrModel)) DEALLOCATE (cellCntrModel)
    IF (ALLOCATED(valueModel)) DEALLOCATE (valueModel)
    IF (ALLOCATED(valueDA)) DEALLOCATE (valueDA)
    IF (ALLOCATED(valueModel_tmp)) DEALLOCATE (valueModel_tmp)
    IF (ALLOCATED(TimeSlots_target)) DEALLOCATE (TimeSlots_target)
    IF (ALLOCATED(this%modvarFns)) DEALLOCATE (this%modvarFns)
    IF (ALLOCATED(this%varNames)) DEALLOCATE (this%varNames)
    IF (ALLOCATED(this%mdate)) DEALLOCATE (this%mdate)


  END SUBROUTINE m_read_bcg_into_Xm_Ens

END MODULE IOGrapesModelvar_m
