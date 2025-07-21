!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2021-09-26, created by JiongmingPang for reading WRF as background.
!                     2022-01-27, modified by JiongmingPang for dealing with Multi-Var.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-09-26,, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang, 2022-01-27, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/6/11, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE IOWRF_m
  USE kinds_m
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE ncReadVar_m
  USE AdvanceTime_m
  USE WRFIO_m, ONLY: WRFIO_t
  USE InterpolateTime_m
  USE Interp1D_m
  USE Widgets_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE Interp1D_m
  USE conversions_m
  USE Filter_m, ONLY: smoothField, guidedfilter
  ! USE netcdf

  TYPE, EXTENDS(IOModel_t) :: IOWRF_t
    CHARACTER(LEN=1024) :: nmlDAFilename, DEMFileName
    CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: inFileNames
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: varNames
    INTEGER(i_kind) :: numFiles

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm_Ens

  END TYPE IOWRF_t

CONTAINS

  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOWRF_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%initializeIOModel(configFile, geometry)
    this%nmlDAFilename = TRIM(configFile)
  END SUBROUTINE

  SUBROUTINE interpVerticalWRF_unStaged(sg, WRF_topo, WRF_height, WRF_data, Target_data)
    IMPLICIT NONE
    TYPE(SingleGrid_t), INTENT(IN)  :: sg
    REAL(r_kind), INTENT(IN)  :: WRF_topo(:)        ! X x Y
    REAL(r_kind), INTENT(IN)  :: WRF_height(:, :, :, :)   ! X x Y x Z x T
    REAL(r_kind), INTENT(IN)  :: WRF_data(:, :, :, :)     ! X x Y x Z x T
    REAL(r_kind), INTENT(OUT) :: Target_data(:, :, :)    ! vLevel x Horizontal x T

    INTEGER(i_kind) :: shape_indata(4)
    INTEGER(i_kind) :: i, j, z, t
    REAL(r_kind), ALLOCATABLE :: WRF_sigma(:)
    REAL(r_kind), ALLOCATABLE :: Target_data_tmp(:, :, :, :)   ! X x Y x vLevel x T

    shape_indata = SHAPE(WRF_data)
    ALLOCATE (Target_data_tmp(shape_indata(1), shape_indata(2), sg%vLevel, shape_indata(4)))

    DO t = 1, shape_indata(4)
      DO i = 1, shape_indata(1)
        DO j = 1, shape_indata(2)
          ALLOCATE (WRF_sigma(shape_indata(3)))
          WRF_sigma = (WRF_height(j, i, :, t) - WRF_topo(i * j)) / (sg%ztop - WRF_topo(i * j)) * sg%ztop
          CALL interp1d(WRF_sigma, WRF_data(i, j, :, t), sg%sigma, Target_data_tmp(i, j, :, t))
          DEALLOCATE (WRF_sigma)
        END DO
      END DO
    END DO

    DO z = 1, sg%vLevel
      DO t = 1, shape_indata(4)
        Target_data(z, :, t) = RESHAPE(Target_data_tmp(:, :, z, t), (/shape_indata(1) * shape_indata(2)/))
      END DO
    END DO

    DEALLOCATE (Target_data_tmp)

  END SUBROUTINE interpVerticalWRF_unStaged

  SUBROUTINE m_read_config_file(this, ensIndx)
    CLASS(IOWRF_t) :: this
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: ensIndx

    CHARACTER(LEN=1024), ALLOCATABLE :: inFileNames(:)
    CHARACTER(LEN=20), ALLOCATABLE :: inVarNames(:)
    CHARACTER(LEN=1024) :: configFile
    CHARACTER(LEN=1024) :: inFileDir, inFileDirt
    CHARACTER(LEN=1024) :: ensDirtmp
    CHARACTER(LEN=10)   :: vartmp
    CHARACTER(LEN=3)    :: memid
    INTEGER(i_kind)     :: t, inFilesNum, varWRFNum, ifile
    LOGICAL             :: ensFlag

    ifile = yaml_get_var(this%m_configFile, 'IO', 'input_dir_model', inFileDir)
    ! inFileDirt = TRIM(inFileDirt)
    !CALL namelist_read(TRIM(this%nmlDAFilename), 'bakFileList', inFileNames)
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'IO', 'bakFileList', inFileNames)
    inFilesNum = UBOUND(inFileNames, 1)
    PRINT *, "inFilesNum: ", inFilesNum
    this%numFiles = inFilesNum
    ALLOCATE (this%inFileNames(inFilesNum))

    ! IF ( PRESENT(ensemble) .AND. PRESENT(ensIndx) ) THEN
    IF (PRESENT(ensIndx)) THEN
      ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensFlag', ensFlag)
      IF (ensFlag) THEN
        ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'BMat', 'ensPath', ensDirtmp)
        WRITE (memid, '(A1,I2.2)') "m", ensIndx
        inFileDir = TRIM(inFileDir)//"/"//TRIM(ensDirtmp)//"/"//memid
      END IF
      ! ELSE
      !   inFileDir = TRIM(inFileDirt)
    END IF

    DO t = 1, inFilesNum
      PRINT *, 'inFileNames: ', TRIM(inFileNames(t))
      this%inFileNames(t) = TRIM(inFileDir)//"/"//TRIM(inFileNames(t))
      PRINT *, 'inFileNames: ', TRIM(this%inFileNames(t))
    END DO

    !CALL namelist_read(TRIM(this%nmlDAFilename), 'varWRFList', this%varNames)
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'IO', 'varWRFList', this%varNames)

    this%DEMFileName = TRIM(inFileDir)//"/DL_DEMdata.nc"
    PRINT *, 'DEMFileName: ', TRIM(this%DEMFileName)

    IF (ALLOCATED(inFileNames)) DEALLOCATE (inFileNames)

  END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm_Ens(this, Xm, sg, ensIndx)
    USE WRFIO_m, ONLY: WRFIO_t
    USE InterpolateTime_m
    USE WriteVar_m
    USE ModelCoupler_m, ONLY: ModelCoupler_t

    IMPLICIT NONE

    CLASS(IOWRF_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(WRFIO_t)  :: WRFIO
    TYPE(SingleGrid_t), INTENT(INOUT)  :: sg
    TYPE(ModelCoupler_t) :: ModelCoupler
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: ensIndx

    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    REAL(r_kind), ALLOCATABLE :: cellCntrDA(:, :)
    REAL(r_kind), ALLOCATABLE :: valueModel(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueModel_tmp(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueDA(:, :, :)

    INTEGER(i_kind), ALLOCATABLE :: TimeSlots_target(:)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA(:, :)

    INTEGER(i_kind) :: i, j, k, HorNumGridModel, nprocInterp
    TYPE(InterpolateData_Time_t), ALLOCATABLE :: IntpTime

    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind)   :: iv, varIndx
    REAL(r_kind), ALLOCATABLE :: timeWRF(:), timeDA(:), timeModel(:)
    REAL(r_kind) :: u, v

    ! CHARACTER(LEN=1024) :: OUTPUT_DIR
    ! CHARACTER(LEN=1024) :: iobakoutfnnc
    ! CHARACTER(LEN=1024) :: iobakoutfntxt

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
        WRFIO = WRFIO_t(this%inFileNames, this%varNames)

        ! get the time_slots_target
        ALLOCATE (TimeSlots_target(sg_mp%tSlots))
        CALL IntpLinearTimeSlots(WRFIO%time_unix, TimeSlots_target)

        HorNumGridModel = WRFIO%nx * WRFIO%ny

        ALLOCATE (cellCntrModel(2, HorNumGridModel), &
                  valueModel(sg_mp%vLevel, HorNumGridModel, sg_mp%tSlots), &
                  valueDA(sg_mp%vLevel, sg_mp%num_icell_global, sg_mp%tSlots), &
                  valueModel_tmp(sg_mp%vLevel, HorNumGridModel, WRFIO%nt))

        cellCntrModel(1, :) = RESHAPE(WRFIO%lat2D, (/HorNumGridModel/))
        cellCntrModel(2, :) = RESHAPE(WRFIO%lon2D, (/HorNumGridModel/))

      END IF

      CALL sg%mpddInfo_sg%bCast(WRFIO%nt)
      CALL sg%mpddInfo_sg%bCast(WRFIO%nx)
      CALL sg%mpddInfo_sg%bCast(WRFIO%ny)
      CALL sg%mpddInfo_sg%bCast(WRFIO%nz)
      CALL sg%mpddInfo_sg%bCast(HorNumGridModel)
      ! CALL sg%mpddInfo_sg%bCast(this%DEMFileName)

      BLOCK
        USE mo_netcdf
        USE parameters_m
        TYPE(NcDataset) :: nc
        TYPE(NcVariable)  :: var
        TYPE(ModelCoupler_t)     :: ModelCouplerForTopo

        REAL(r_kind), ALLOCATABLE :: latTopoFromNC(:), lonTopoFromNC(:), valueTopoFromNC(:, :)
        REAL(r_kind), ALLOCATABLE :: llTopo(:, :), valueTopo(:)
        INTEGER(i_kind) :: lat_len, lon_len

        IF (sg_mp%isBaseProc()) THEN

          PRINT *, 'DEMFileName: ', TRIM(this%DEMFileName)
          nc = NcDataset(this%DEMFileName, "r")

          var = nc%getVariable("lat"); CALL var%getData(latTopoFromNC)
          var = nc%getVariable("lon"); CALL var%getData(lonTopoFromNC)
          var = nc%getVariable("height"); CALL var%getData(valueTopoFromNC)

          CALL nc%CLOSE()

          lat_len = SIZE(latTopoFromNC)
          lon_len = SIZE(lonTopoFromNC)

          PRINT *, 'lat_len: ', lat_len
          PRINT *, 'lon_len: ', lon_len

          ALLOCATE (llTopo(2, lat_len * lon_len), valueTopo(lat_len * lon_len))

          DO i = 1, lat_len
            DO j = 1, lon_len
              llTopo(1, (j - 1) * lat_len + i) = latTopoFromNC(i) * degree2radian
              llTopo(2, (j - 1) * lat_len + i) = lonTopoFromNC(j) * degree2radian
              valueTopo((j - 1) * lat_len + i) = valueTopoFromNC(j, i)
            END DO
          END DO

          BLOCK
            REAL(r_kind), ALLOCATABLE ::  tmp2D(:, :)
            IF (sg%isBaseProc()) THEN
              ALLOCATE (tmp2D(lon_len, lat_len))
              tmp2D = RESHAPE(valueTopo, [lon_len, lat_len])
              CALL guidedfilter(tmp2D, tmp2D, 5, 4.0D0, tmp2D)
              valueTopo = RESHAPE(tmp2D, [lon_len * lat_len])
            END IF

            IF (sg%isBaseProc()) DEALLOCATE (tmp2D)
          END BLOCK

          PRINT *, 'MAXVAL(llTopo)', MAXVAL(llTopo(1, :)) * radian2degree, MAXVAL(llTopo(2, :)) * radian2degree
          PRINT *, 'MINVAL(llTopo)', MINVAL(llTopo(1, :)) * radian2degree, MINVAL(llTopo(2, :)) * radian2degree
        END IF

        CALL ModelCouplerForTopo%Initialize(this%m_configFile, this%geometry)
        CALL ModelCouplerForTopo%gen_interp_coeffs(llTopo, lat_len * lon_len, sg)

        ! -----------------------------------------------------------
        ! CALL ModelCouplerForTopo%ingest_to_topo_and_update_zHght_d(valueTopo*1.0D0, lat_len*lon_len, sg_mp)
        ! -----------------------------------------------------------

        CALL ModelCouplerForTopo%interhpN%interp_singleLevel_d(valueTopo * 1.0D0, sg_mp%topo)
        CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
        CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:), tmp2D(:, :)
          INTEGER(i_kind) :: bufShape(1)
          IF (sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%num_icell_global))
          END IF

          CALL sg_mp%aggrGridReal1D(sg_mp%topo, tmp, sg_mp%num_icell_global)

          IF (sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg_mp%dimCell_global(2), sg_mp%dimCell_global(1)))

            tmp2D = RESHAPE(tmp(:), [sg_mp%dimCell_global(2), sg_mp%dimCell_global(1)])
            CALL guidedfilter(tmp2D, tmp2D, 5, 4.0D0, tmp2D)
            tmp(:) = RESHAPE(tmp2D, [sg_mp%num_icell_global])
          END IF

          CALL sg%DistGridReal1D(tmp, sg_mp%topo, [sg_mp%num_icell_global])
          IF (sg%isBaseProc()) DEALLOCATE (tmp, tmp2D)
        END BLOCK

        CALL sg_mp%update_zHght_from_topo_and_sigma

      END BLOCK

      ! Initialize the model coupler.
      CALL ModelCoupler%Initialize(this%m_configFile, this%geometry)
      CALL ModelCoupler%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)

      IF (sg_mp%isBaseProc()) THEN
        timeModel = WRFIO%time_unix
        timeDA = TimeSlots_target
      ELSE
        ALLOCATE (timeModel(WRFIO%nt), timeDA(sg%tSlots))
      END IF
      PRINT *, "ccc: ", WRFIO%nt, sg%tSlots, SHAPE(timeModel), SHAPE(timeDA)

      CALL sg%mpddInfo_sg%bCastReal1D(timeModel)
      CALL sg%mpddInfo_sg%bCastReal1D(timeDA)

      ! Modelvar only work for surface analysis
      PRINT *, "sg_mp%vLevel", sg_mp%vLevel

      ! Ingest the topo and update the height in MOTOR-DA
      ! CALL ModelCoupler%ingest_to_topo_and_update_zHght_d(WRFIO%topo, HorNumGridModel, sg_mp)

      IF (sg_mp%vLevel == 1) THEN
        !! Ingest the 2D background field

        IF (sg_mp%isBaseProc()) THEN
          timeWRF = WRFIO%time_unix
          timeDA = TimeSlots_target
        END IF

        ! Ingest t2m-temperature at 2m
        IF (Xm%getVarIdx('temp') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 't2m') < 0) THEN
              PRINT *, "ERROR!!! There is no 't2m'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%t2m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
          PRINT *, 'IOWRF - temp has been loaded.', Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, 1, 1), &
            sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
        END IF

        ! Ingest U10-temperature at 2m
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'u10m') < 0) THEN
              PRINT *, "ERROR!!! There is no 'U10'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%u10m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOWRF - u10m has been loaded.'
        END IF

        ! Ingest V10-temperature at 2m
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'v10m') < 0) THEN
              PRINT *, "ERROR!!! There is no 'V10'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%v10m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOWRF - v10m has been loaded.'
        END IF

        IF (Xm%getVarIdx('winds') .NE. 0 .AND. Xm%getVarIdx('windd') .NE. 0) THEN

          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'u10m') < 0) THEN
              PRINT *, "ERROR!!! There is no 'U10'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%u10m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('winds'))%DATA, sg_mp)
          PRINT *, 'IOWRF - u10m has been loaded.'

          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'v10m') < 0) THEN
              PRINT *, "ERROR!!! There is no 'V10'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%v10m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('windd'))%DATA, sg_mp)
          PRINT *, 'IOWRF - v10m has been loaded.'

          DO i = 1, sg_mp%num_cell
            DO j = 1, sg_mp%tSlots
              u = Xm%Fields(Xm%getVarIdx('winds'))%DATA(1, i, j)
              v = Xm%Fields(Xm%getVarIdx('windd'))%DATA(1, i, j)
              CALL uv_to_wind(u, v, Xm%Fields(Xm%getVarIdx('winds'))%DATA(1, i, j), Xm%Fields(Xm%getVarIdx('windd'))%DATA(1, i, j))
            END DO
          END DO
        END IF

        ! Ingest q2m-qvapor at 2m
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'q2m') < 0) THEN
              PRINT *, "ERROR!!! There is no 'q2m'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%q2m, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
          PRINT *, 'IOWRF - q2m has been loaded.'
        END IF

        ! Ingest pressure at ground
        IF (Xm%getVarIdx('pres') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (getStrIndex(this%varNames, 'psfc') < 0) THEN
              PRINT *, "ERROR!!! There is no 'psfc'. Please CHECK 'varWRFList' in your namelist file! "
              STOP
            END IF
            valueModel_tmp = RESHAPE(WRFIO%psfc, (/1, HorNumGridModel, WRFIO%nt/))

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeWRF, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
          PRINT *, 'IOWRF - psfc has been loaded.'
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
        PRINT *, "GrapesModelvarIO_500m_t run in vertical level other than 1  vLevel: vLevel", sg_mp%vLevel
        ALLOCATE (zRHghtModelAtDA(WRFIO%nz, sg_mp%num_cell))
        BLOCK
          REAL(r_single), ALLOCATABLE :: valueTemp(:, :)

          ALLOCATE (valueTemp(WRFIO%nz, HorNumGridModel))

          PRINT *, 'WRFIO%nt: ', WRFIO%nt

          IF (sg%isBaseProc()) THEN
            FORALL (i=1:WRFIO%nx, j=1:WRFIO%ny)
              valueTemp(:, (j - 1) * WRFIO%nx + i) = WRFIO%height(i, j, :, 1)
            END FORALL
          END IF
          CALL sg%mpddInfo_sg%barrier
          CALL ModelCoupler%ingest_to_field_data_2D(WRFIO%nz, &
                                                    valueTemp, zRHghtModelAtDA, sg_mp)          ! CALL sg%mpddInfo_sg%barrier

          PRINT *, "zRHghtModelAtDA", MAXVAL(zRHghtModelAtDA), timeModel, timeDA
          DEALLOCATE (valueTemp)
        END BLOCK

        !!! 开始进行数据插值
        BLOCK
          REAL(r_kind), ALLOCATABLE :: valueDA_mp(:, :, :)
          REAL(r_single), ALLOCATABLE :: tmpModel_u(:, :, :)
          REAL(r_single), ALLOCATABLE :: tmp_u(:, :, :)

          IF (sg_mp%isBaseProc()) THEN
            ALLOCATE (tmpModel_u(WRFIO%nz, HorNumGridModel, WRFIO%nt), &
                      tmp_u(WRFIO%nz, HorNumGridModel, sg_mp%tSlots))
          END IF

          ! Ingest temp-temperature
          IF (Xm%getVarIdx('temp') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (getStrIndex(this%varNames, 'temp') < 0) THEN
                PRINT *, "ERROR!!! There is no 'temp'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF

              ! for vertical inter/extra-polation
              IF (.NOT. ALLOCATED(sg_mp%topo) .OR. .NOT. ALLOCATED(WRFIO%height)) THEN
                PRINT *, "ERROR!!! There is no 'topo' or no 'height'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF

              CALL interpVerticalWRF_unStaged(sg_mp, sg_mp%topo, WRFIO%height, WRFIO%temp, valueModel_tmp(:, :, :))
              ! for time-scale interpolation
              ALLOCATE (IntpTime)
              CALL IntpTime%IntpInitialize(valueModel_tmp, WRFIO%time_unix, TimeSlots_target)
              CALL IntpTime%IntpLinearTime(valueModel_tmp)
              valueModel = IntpTime%data_intperpolate
              DEALLOCATE (IntpTime)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
            PRINT *, 'IOWRF - temp has been loaded.'
          END IF

          ! Ingest uwnd
          IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              FORALL (k=1:WRFIO%nt, i=1:WRFIO%nx, j=1:WRFIO%ny)
                tmpModel_u(:, (j - 1) * WRFIO%nx + i, k) = WRFIO%uwnd(i, j, :, k)
              END FORALL

              CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
              PRINT *, "End of interp1d_3D_idx3."
            END IF

            CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(WRFIO%nz, zRHghtModelAtDA, tmp_u, &
                                                               Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - uwnd has been loaded.'
          END IF

          ! Ingest vwnd
          IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              FORALL (k=1:WRFIO%nt, i=1:WRFIO%nx, j=1:WRFIO%ny)
                tmpModel_u(:, (j - 1) * WRFIO%nx + i, k) = WRFIO%vwnd(i, j, :, k)
              END FORALL

              CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
              PRINT *, "End of interp1d_3D_idx3."
            END IF

            CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(WRFIO%nz, zRHghtModelAtDA, tmp_u, &
                                                               Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - vwnd has been loaded.'
          END IF

          ! ! Ingest wwnd
          ! IF (Xm%getVarIdx('wwnd') .NE. 0) THEN
          !   IF (sg_mp%isBaseProc()) THEN
          !     FORALL (k=1:WRFIO%nt, i=1:WRFIO%nx, j=1:WRFIO%ny)
          !       tmpModel_u(:, (j - 1) * WRFIO%nx + i, k) = WRFIO%wwnd(i, j, :, k)
          !     END FORALL

          !     CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
          !     PRINT *, "End of interp1d_3D_idx3."
          !   END IF

          !   CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(WRFIO%nz, zRHghtModelAtDA, tmp_u, &
          !                                                      Xm%Fields(Xm%getVarIdx('wwnd'))%DATA, sg_mp)
          !   PRINT *, 'IOGrapes - wwnd has been loaded.'

          ! END IF

          ! ! Ingest uwnd-u wind
          ! IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          !   IF (sg_mp%isBaseProc()) THEN
          !     IF (getStrIndex(this%varNames, 'uwnd') < 0) THEN
          !       PRINT *, "ERROR!!! There is no 'uwnd'. Please CHECK 'varWRFList' in your namelist file! "
          !       STOP
          !     END IF

          !     ! for vertical inter/extra-polation
          !     IF (.NOT. ALLOCATED(WRFIO%topo) .OR. .NOT. ALLOCATED(WRFIO%height)) THEN
          !       PRINT *, "ERROR!!! There is no 'topo' or no 'height'. Please CHECK 'varWRFList' in your namelist file! "
          !       STOP
          !     END IF
          !     CALL interpVerticalWRF_unStaged(sg_mp, sg_mp%topo, WRFIO%height, WRFIO%uwnd, valueModel_tmp(:, :, :))
          !     ! for time-scale interpolation
          !     ALLOCATE (IntpTime)
          !     CALL IntpTime%IntpInitialize(valueModel_tmp, WRFIO%time_unix, TimeSlots_target)
          !     CALL IntpTime%IntpLinearTime(valueModel_tmp)
          !     valueModel = IntpTime%data_intperpolate
          !     DEALLOCATE (IntpTime)
          !   END IF

          !   CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          !   PRINT *, 'IOWRF - uwnd has been loaded.'
          ! END IF

          ! Ingest vwnd-v wind
          ! IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          !   IF (sg_mp%isBaseProc()) THEN
          !     IF (getStrIndex(this%varNames, 'vwnd') < 0) THEN
          !       PRINT *, "ERROR!!! There is no 'vwnd'. Please CHECK 'varWRFList' in your namelist file! "
          !       STOP
          !     END IF

          !     ! for vertical inter/extra-polation
          !     IF (.NOT. ALLOCATED(sg_mp%topo) .OR. .NOT. ALLOCATED(WRFIO%height)) THEN
          !       PRINT *, "ERROR!!! There is no 'topo' or no 'height'. Please CHECK 'varWRFList' in your namelist file! "
          !       STOP
          !     END IF
          !     CALL interpVerticalWRF_unStaged(sg_mp, sg_mp%topo, WRFIO%height, WRFIO%vwnd, valueModel_tmp(:, :, :))
          !     ! for time-scale interpolation
          !     ALLOCATE (IntpTime)
          !     CALL IntpTime%IntpInitialize(valueModel_tmp, WRFIO%time_unix, TimeSlots_target)
          !     CALL IntpTime%IntpLinearTime(valueModel_tmp)
          !     valueModel = IntpTime%data_intperpolate
          !     DEALLOCATE (IntpTime)
          !   END IF

          !   CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          !   PRINT *, 'IOWRF - vwnd has been loaded.'
          ! END IF

          ! Ingest pres-pressure
          IF (Xm%getVarIdx('pres') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (getStrIndex(this%varNames, 'pres') < 0) THEN
                PRINT *, "ERROR!!! There is no 'pres'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF

              ! for vertical inter/extra-polation
              IF (.NOT. ALLOCATED(sg_mp%topo) .OR. .NOT. ALLOCATED(WRFIO%height)) THEN
                PRINT *, "ERROR!!! There is no 'topo' or no 'height'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF
              CALL interpVerticalWRF_unStaged(sg_mp, sg_mp%topo, WRFIO%height, WRFIO%pres, valueModel_tmp(:, :, :))
              ! for time-scale interpolation
              ALLOCATE (IntpTime)
              CALL IntpTime%IntpInitialize(valueModel_tmp, WRFIO%time_unix, TimeSlots_target)
              CALL IntpTime%IntpLinearTime(valueModel_tmp)
              valueModel = IntpTime%data_intperpolate
              DEALLOCATE (IntpTime)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
            PRINT *, 'IOWRF - pres has been loaded.'
          END IF

          ! Ingest qvapor
          IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (getStrIndex(this%varNames, 'qvapor') < 0) THEN
                PRINT *, "ERROR!!! There is no 'pres'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF

              ! for vertical inter/extra-polation
              IF (.NOT. ALLOCATED(sg_mp%topo) .OR. .NOT. ALLOCATED(WRFIO%height)) THEN
                PRINT *, "ERROR!!! There is no 'topo' or no 'height'. Please CHECK 'varWRFList' in your namelist file! "
                STOP
              END IF
              CALL interpVerticalWRF_unStaged(sg_mp, sg_mp%topo, WRFIO%height, WRFIO%qvapor, valueModel_tmp(:, :, :))
              ! for time-scale interpolation
              ALLOCATE (IntpTime)
              CALL IntpTime%IntpInitialize(valueModel_tmp, WRFIO%time_unix, TimeSlots_target)
              CALL IntpTime%IntpLinearTime(valueModel_tmp)
              valueModel = IntpTime%data_intperpolate
              DEALLOCATE (IntpTime)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
            PRINT *, 'IOWRF - qvapor has been loaded.'
          END IF
        END BLOCK
      END IF

      CALL sg_mp%mpddInfo_sg%barrier

    END ASSOCIATE

    IF (ALLOCATED(cellCntrModel)) DEALLOCATE (cellCntrModel)
    IF (ALLOCATED(valueModel)) DEALLOCATE (valueModel)
    IF (ALLOCATED(valueDA)) DEALLOCATE (valueDA)
    IF (ALLOCATED(valueModel_tmp)) DEALLOCATE (valueModel_tmp)
    IF (ALLOCATED(TimeSlots_target)) DEALLOCATE (TimeSlots_target)
    IF (ALLOCATED(this%inFileNames)) DEALLOCATE (this%inFileNames)
    ! DEALLOCATE (this%inFileNames)

  END SUBROUTINE m_read_bcg_into_Xm_Ens

END MODULE IOWRF_m
