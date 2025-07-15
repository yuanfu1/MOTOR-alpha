MODULE IOECM_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE ECMIO_m, ONLY: ECMIO_t
  USE YAMLRead_m
  USE Interp1D_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE AdvanceTime_m
  USE parameters_m
  USE Filter_m, ONLY: smoothField, guidedfilter

  IMPLICIT NONE

  TYPE, EXTENDS(IOModel_t) :: IOECM_t
    CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: inFileNames
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: varNames
    INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE :: mdate
    INTEGER(i_kind) :: numFiles
    TYPE(ECMIO_t) :: ECMIO
    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm
    PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
  END TYPE IOECM_t

CONTAINS

  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOECM_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%initializeIOModel(configFile, geometry)
    PRINT *, "configFile: ", configFile
    this%m_configFile = TRIM(configFile)
  END SUBROUTINE

  SUBROUTINE m_read_config_file(this)
    CLASS(IOECM_t) :: this
    CHARACTER(LEN=1024), ALLOCATABLE  :: ECMFileNames(:)
    CHARACTER(LEN=20), ALLOCATABLE    :: inVarNames(:)
    CHARACTER(LEN=1024)               :: inputECMDir
    CHARACTER(LEN=1024)               :: startTime, endTime
    CHARACTER(LEN=10), ALLOCATABLE    :: timestr(:)
    INTEGER(i_kind)                   :: t, inFilesNum, ifile, inTimeNum
    PRINT *, "Reading ECM data...", this%m_configFile
    ifile = yaml_get_var(TRIM(this%m_configFile), 'IO', 'ModelFileName', ECMFileNames)
    ifile = yaml_get_var(TRIM(this%m_configFile), 'IO', 'ModelFileTime', timestr)
    ifile = yaml_get_var(TRIM(this%m_configFile), 'IO', 'input_dir_model', inputECMDir)
    ifile = yaml_get_var(TRIM(this%m_configFile), 'analysis_para', 'start_time', startTime)
    ifile = yaml_get_var(TRIM(this%m_configFile), 'analysis_para', 'end_time', endTime)
    ifile = yaml_get_var(TRIM(this%m_configFile), 'IO', 'varECMList', inVarNames)

    this%varNames = inVarNames
    inFilesNum = UBOUND(ECMFileNames, 1)
    PRINT *, "inFilesNum: ", inFilesNum
    this%numFiles = inFilesNum
    ALLOCATE (this%inFileNames(inFilesNum))

    ! get the full position of ECM nc files
    DO t = 1, inFilesNum
      PRINT *, 'ECMFileNames: ', TRIM(ECMFileNames(t))
      this%inFileNames(t) = TRIM(inputECMDir)//"/"//TRIM(ECMFileNames(t))
      PRINT *, 'ECMFileNames: ', TRIM(this%inFileNames(t))
    END DO

    ! get the times of ECM nc files
    ifile = yaml_get_var(TRIM(this%m_configFile), 'IO', 'ModelFileTime', timestr)
    inTimeNum = UBOUND(timestr, 1)
    IF (inTimeNum .NE. inFilesNum) THEN
      PRINT *, "Error: inTimeNum is not consistent with inFilesNum, please check 'ModelFileName' and 'ModelFileTime' in your yaml file!"
      STOP
    END IF

    PRINT *, "timestr: ", timestr

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

    IF (ALLOCATED(ECMFileNames)) DEALLOCATE (ECMFileNames)
    IF (ALLOCATED(timestr)) DEALLOCATE (timestr)
    IF (ALLOCATED(inVarNames)) DEALLOCATE (inVarNames)

  END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    USE InterpolateTime_m
    USE WriteVar_m
    IMPLICIT NONE
    CLASS(IOECM_t) :: this
    TYPE(ECMIO_t) :: ECMIO
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT)  :: sg
    TYPE(ModelCoupler_t) :: ModelCoupler
    TYPE(InterpolateData_Time_t), ALLOCATABLE :: IntpTime
    CHARACTER(LEN=20) :: varName
    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    REAL(r_kind), ALLOCATABLE :: cellCntrDA(:, :)
    REAL(r_kind), ALLOCATABLE :: valueModel(:, :, :)
    REAL(r_kind), ALLOCATABLE :: valueModel_tmp(:, :, :)
    REAL(r_kind), ALLOCATABLE :: valueDA(:, :, :)
    REAL(r_kind), ALLOCATABLE :: timeModel(:), timeDA(:)
    INTEGER(i_kind), ALLOCATABLE :: TimeSlots_target(:)
    INTEGER(i_kind) :: i, j, k, HorNumGridModel, nprocInterp
    INTEGER(i_kind)   :: iv, varIndx

    ASSOCIATE (sg_mp => sg)

      IF (sg_mp%isBaseProc()) THEN
        CALL this%m_read_config_file()
      END IF

      CALL sg_mp%mpddInfo_sg%bCast(this%numFiles)
      IF (this%numFiles == 0) THEN
        PRINT *, 'No file listed in the configurations.'
        RETURN
      END IF
      PRINT *, '*************************Number of files to be read: ', this%numFiles
      PRINT *, this%mdate(:, 1), this%mdate(:, 2)

      IF (sg_mp%isBaseProc()) THEN
        ECMIO = ECMIO_t(this%inFileNames, this%varNames, this%mdate)

        ALLOCATE (TimeSlots_target(sg_mp%tSlots))
        CALL IntpLinearTimeSlots(ECMIO%time_unix, TimeSlots_target)
        PRINT *, 'ECMIO%time_unix: ', ECMIO%time_unix

        HorNumGridModel = ECMIO%nx * ECMIO%ny

        ALLOCATE (cellCntrModel(2, HorNumGridModel), &
                  valueModel(1, HorNumGridModel, sg_mp%tSlots), &
                  valueDA(sg_mp%vLevel, sg_mp%num_icell_global, sg_mp%tSlots), &
                  valueModel_tmp(1, HorNumGridModel, ECMIO%nt))
        cellCntrModel(1, :) = RESHAPE(ECMIO%lat2D, (/HorNumGridModel/))
        cellCntrModel(2, :) = RESHAPE(ECMIO%lon2D, (/HorNumGridModel/))
      END IF
      CALL sg%mpddInfo_sg%bCast(HorNumGridModel)
      CALL sg%mpddInfo_sg%bCast(ECMIO%nt)
      CALL sg%mpddInfo_sg%bCast(ECMIO%nx)
      CALL sg%mpddInfo_sg%bCast(ECMIO%ny)

      BLOCK
        TYPE(ModelCoupler_t)     :: ModelCouplerForTopo
        ! REAL(r_kind), ALLOCATABLE :: llTopo(:, :), valueTopo(:)

        ! ALLOCATE (llTopo(2, ECMIO%ny * ECMIO%nx), valueTopo(ECMIO%ny * ECMIO%nx))

        ! DO i = 1, ECMIO%ny
        !   DO j = 1, ECMIO%nx
        !     llTopo(1, (j - 1) * ECMIO%ny + i) = ECMIO%lat(i) * degree2radian
        !     llTopo(2, (j - 1) * ECMIO%ny + i) = ECMIO%lon(j) * degree2radian
        !     valueTopo((j - 1) * ECMIO%ny + i) = ECMIO%lat(j, i)
        !   END DO
        ! END DO
        !TODO
        CALL ModelCouplerForTopo%Initialize(this%m_configFile, this%geometry)
        CALL ModelCouplerForTopo%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)
        ! elCoupler%ingest_to_topo_and_update_zHght_d(ECMIO%topo, HorNumGridModel, sg_mp)
        ! CALL ModelCouplerForTopo%interhpN%interp_singleLevel_d(valueTopo * 1.0D0, sg_mp%topo)
        ! CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
        ! CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)
        ! -----------------------------------------------------------

        PRINT *, 'here'
      END BLOCK

      CALL ModelCoupler%Initialize(this%m_configFile, this%geometry)
      CALL ModelCoupler%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)

      IF (sg_mp%isBaseProc()) THEN
        timeModel = ECMIO%time_unix
        timeDA = TimeSlots_target
      ELSE
        ALLOCATE (timeModel(ECMIO%nt), timeDA(sg%tSlots))
      END IF
      PRINT *, "ccc: ", ECMIO%nt, sg%tSlots, SHAPE(timeModel), SHAPE(timeDA)

      CALL sg%mpddInfo_sg%bCastReal1D(timeModel)
      CALL sg%mpddInfo_sg%bCastReal1D(timeDA)
      ! PRINT *, 'aaa: ', timeModel
      ! PRINT*, 'bbb: ', timeDA

      ! STOP

      ! postvar only work for surface analysis
      PRINT *, "sg_mp%vLevel", sg_mp%vLevel

      IF (Xm%getVarIdx('temp') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%t2m)) THEN
            PRINT *, "ERROR!!! There is no 't2m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%t2m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('temp'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('temp'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK

        PRINT *, 'ECM_NC_file - t2m has been loaded.', Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%q2m)) THEN
            PRINT *, "ERROR!!! There is no 'q2m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%q2m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK

        PRINT *, 'ECM_NC_file - q2m has been loaded.', Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('pcpa') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%precipitation)) THEN
            PRINT *, "ERROR!!! There is no 'precipitation' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%precipitation, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pcpa'))%DATA, sg_mp)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('pcpa'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('pcpa'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK

        PRINT *, 'ECM_NC_file - precipitation has been loaded.', Xm%Fields(Xm%getVarIdx('pcpa'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%u10m)) THEN
            PRINT *, "ERROR!!! There is no 'u10m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%u10m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        
        PRINT *, 'ECM_NC_file - u10m has been loaded.', Xm%Fields(Xm%getVarIdx('uwnd'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%v10m)) THEN
            PRINT *, "ERROR!!! There is no 'v10m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%v10m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK

        PRINT *, 'ECM_NC_file - v10m has been loaded.', Xm%Fields(Xm%getVarIdx('vwnd'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('u200m') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%u200m)) THEN
            PRINT *, "ERROR!!! There is no 'u200m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%u200m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('u200m'))%DATA, sg_mp)

        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('u200m'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('u200m'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        PRINT *, 'ECM_NC_file - u200m has been loaded.', Xm%Fields(Xm%getVarIdx('u200m'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('v200m') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%u200m)) THEN
            PRINT *, "ERROR!!! There is no 'u200m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%u200m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('v200m'))%DATA, sg_mp)
        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('v200m'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('v200m'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        PRINT *, 'ECM_NC_file - u200m has been loaded.', Xm%Fields(Xm%getVarIdx('v200m'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('u100m') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%u100m)) THEN
            PRINT *, "ERROR!!! There is no 'u100m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%u100m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('u100m'))%DATA, sg_mp)
         BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('u100m'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('u100m'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        PRINT *, 'ECM_NC_file - u100m has been loaded.', Xm%Fields(Xm%getVarIdx('u100m'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('v100m') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%v100m)) THEN
            PRINT *, "ERROR!!! There is no 'v100m' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%v100m, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('v100m'))%DATA, sg_mp)
        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('v100m'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('v100m'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        PRINT *, 'ECM_NC_file - v100m has been loaded.', Xm%Fields(Xm%getVarIdx('v100m'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      IF (Xm%getVarIdx('pres') .NE. 0) THEN
        IF (sg_mp%isBaseProc()) THEN
          IF (.NOT. ALLOCATED(ECMIO%pressure)) THEN
            PRINT *, "ERROR!!! There is no 'pressure' in ECM_NC_file! "
            STOP
          ELSE
            valueModel_tmp = RESHAPE(ECMIO%pressure, (/1, HorNumGridModel, ECMIO%nt/))
          END IF
          CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
        END IF
        CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
        BLOCK
          REAL(r_kind), ALLOCATABLE :: tmp(:, :, :), tmp2D(:, :)
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%vLevel, sg%num_icell_global, sg%tSlots))
          END IF
          CALL sg%aggrGridRealForFieldGrid(Xm%Fields(Xm%getVarIdx('pres'))%DATA, &
                                         tmp, [sg%vLevel, sg%num_icell_global, sg%tSlots])
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg%dimCell_global(2), sg%dimCell_global(1)))

            DO i = 1, sg%tSlots
              DO j = 1, sg%vLevel
                tmp2D = reshape(tmp(j, :, i), [sg%dimCell_global(2), sg%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 3, 0.4D0, tmp2D)
                tmp(j, :, i) = reshape(tmp2D, [sg%num_icell_global])
              END DO
            END DO  
          END IF

          CALL sg%DistGridRealForFieldGrid(tmp, Xm%Fields(Xm%getVarIdx('pres'))%DATA, [sg%vLevel, sg%num_icell_global, sg%tSlots])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK
        PRINT *, 'ECM_NC_file - pressure has been loaded.', Xm%Fields(Xm%getVarIdx('pres'))%DATA(1, 1, 1), &
          sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      END IF

      ! DO iv = 1, ECMIO%numVars
      !   varName = this%varNames(iv)
      !   varIndx = Xm%getVarIdx(varName)
      !   IF (varIndx .NE. 0) THEN
      !     IF (sg_mp%isBaseProc()) THEN
      !       IF (.NOT. ALLOCATED(ECMIO%var2D)) THEN
      !         PRINT *, "ERROR!!! There is no '"//TRIM(varName)//"' in ECM_NC_file! "
      !         STOP
      !       ELSE
      !         valueModel_tmp = RESHAPE(ECMIO%var2D, (/1, HorNumGridModel, ECMIO%nt/))
      !       END IF
      !       CALL interp1d_3D_idx3(timeModel, valueModel_tmp, timeDA, valueModel)
      !       CALL ModelCoupler%ingest_to_field_data_MT( sg_mp%vLevel, valueModel, Xm%Fields(varIndx), sg_mp)
      !       PRINT *, TRIM(varName)//' has been loaded.', Xm%Fields(varIndx)%DATA(1, 1, 1), &
      !         sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
      !         ! STOP
      !     END IF
      !   END IF
      ! END DO
    END ASSOCIATE

    IF (ALLOCATED(cellCntrModel)) DEALLOCATE (cellCntrModel)
    IF (ALLOCATED(valueModel)) DEALLOCATE (valueModel)
    IF (ALLOCATED(valueDA)) DEALLOCATE (valueDA)
    IF (ALLOCATED(valueModel_tmp)) DEALLOCATE (valueModel_tmp)
    IF (ALLOCATED(TimeSlots_target)) DEALLOCATE (TimeSlots_target)
    IF (ALLOCATED(this%inFileNames)) DEALLOCATE (this%inFileNames)
    ! Deallocate memories

  END SUBROUTINE m_read_bcg_into_Xm

END MODULE IOECM_m
