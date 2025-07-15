!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-08-09, created by JiongmingPang for reading Modelvar.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-08-09, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE IOGrapesModelvar_500m_m
    USE kinds_m
    USE IOModel_m, ONLY: IOModel_t
    USE State_m, ONLY: State_t
    USE geometry_m, ONLY: geometry_t
    USE singleGrid_m, ONLY: singleGrid_t
    USE YAMLRead_m
    USE AdvanceTime_m
    USE GrapesModelvarIO_500m_m, ONLY: GrapesModelvarIO_500m_t
    USE InterpolateTime_m
    USE Widgets_m
    USE ModelCoupler_m, ONLY: ModelCoupler_t
    USE MGOpts_m, ONLY: Calc_sg_vars, update_restrictionOfStatics
    USE Interp1D_m
    USE conversions_m
    USE parameters_m
    USE Filter_m, ONLY: smoothField, guidedfilter

    TYPE, EXTENDS(IOModel_t) :: IOGrapesModelvar_500m_t
      CHARACTER(LEN=1024) :: nmlDAFilename, ModelvarCtlFn, DEMFileName
      CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: ModelvarFns
      CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: varNames
      INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE :: mdate
      INTEGER(i_kind) :: numFiles
  
    CONTAINS
      PROCEDURE, PUBLIC :: initialize
      PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
      PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm
  
    END TYPE IOGrapesModelvar_500m_t
  
  CONTAINS
  
    SUBROUTINE initialize(this, configFile, geometry)
      IMPLICIT NONE
      CLASS(IOGrapesModelvar_500m_t) :: this
      TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
      CHARACTER(LEN=1024), INTENT(IN) :: configFile
  
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

    SUBROUTINE m_read_config_file(this)
      CLASS(IOGrapesModelvar_500m_t) :: this
  
      CHARACTER(LEN=1024), ALLOCATABLE :: ModelvarFns(:)
      CHARACTER(LEN=20), ALLOCATABLE :: inVarNames(:)
      CHARACTER(LEN=12), ALLOCATABLE :: timestr(:)
      CHARACTER(LEN=1024) :: configFile
      CHARACTER(LEN=1024) :: inFileDir
      CHARACTER(LEN=1024) :: ensDirtmp
      CHARACTER(LEN=3)    :: memid
      INTEGER(i_kind)     :: t, inFilesNum, ifile, inTimeNum
      LOGICAL             :: ensFlag
  
      ifile = yaml_get_var(this%m_configFile, 'IO', 'input_dir_model', inFileDir)
      ifile = yaml_get_var(this%m_configFile, 'IO', 'ModelFileName', ModelvarFns)
      inFilesNum = UBOUND(ModelvarFns, 1)
      PRINT *, "inFilesNum: ", inFilesNum
      this%numFiles = inFilesNum
      ALLOCATE (this%ModelvarFns(inFilesNum))
  
      ! get the full position of Modelvar files
      DO t = 1, inFilesNum
        PRINT *, 'ModelvarFns: ', TRIM(ModelvarFns(t))
        this%ModelvarFns(t) = TRIM(inFileDir)//"/"//TRIM(ModelvarFns(t))
        PRINT *, 'ModelvarFns: ', TRIM(this%ModelvarFns(t))
      END DO
  
      ! get the full position of Modelvar control files (post.ctl)
      this%ModelvarCtlFn = TRIM(inFileDir)//"/modelvar.ctl"
      PRINT *, 'ModelvarCtlFn: ', TRIM(this%ModelvarCtlFn)
  
      ! get the full position of dem files (DL_DEMdata.csv)
      this%DEMFileName = TRIM(inFileDir)//"/DL_DEMdata.nc"
      PRINT *, 'DEMFileName: ', TRIM(this%DEMFileName)

      ! get the times of Modelvar files
      ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'IO', 'ModelFileTime', timestr)
      inTimeNum = UBOUND(timestr, 1)
      IF (inTimeNum .NE. inFilesNum) THEN
        PRINT *, "Error: inTimeNum is not consistent with inFilesNum, &
  &      please check 'ModelFileName' and 'ModelFileTime' in your yaml file!"
        STOP
      END IF
  
      ALLOCATE (this%mdate(6, inTimeNum))
      DO t = 1, inTimeNum
        READ (timestr(t) (1:4), "(i4)") this%mdate(1, t) !year
        READ (timestr(t) (5:6), "(i2)") this%mdate(2, t) !mon
        READ (timestr(t) (7:8), "(i2)") this%mdate(3, t) !day
        READ (timestr(t) (9:10), "(i2)") this%mdate(4, t) !hor
        READ (timestr(t) (11:12), "(i2)") this%mdate(5, t) !min
        !this%mdate(5, t) = 0 !min
        this%mdate(6, t) = 0 !sec
        PRINT *, "GMT Time in t-step: ", t, this%mdate(:, t)
      END DO
      
      ! get the analysis variables
      ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'modelState', 'varList', this%varNames)

      IF (ALLOCATED(ModelvarFns)) DEALLOCATE (ModelvarFns)
      IF (ALLOCATED(timestr)) DEALLOCATE (timestr)
  
    END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    USE GrapesModelvarIO_500m_m, ONLY: GrapesModelvarIO_500m_t
    USE InterpolateTime_m
    USE ModelCoupler_m, ONLY: ModelCoupler_t

    IMPLICIT NONE

    CLASS(IOGrapesModelvar_500m_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(GrapesModelvarIO_500m_t) :: GrapesModelvarIO_500m
    TYPE(ModelCoupler_t)     :: ModelCoupler

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
    REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA(:, :)

    ASSOCIATE (sg_mp => sg)
      CALL this%m_read_config_file()
      CALL sg_mp%mpddInfo_sg%bCast(this%numFiles)

      IF (this%numFiles == 0) THEN
        PRINT *, 'No file listed in the configurations.'
        STOP
      END IF

      PRINT *, '*************************Number of files to be read: ', this%numFiles
      PRINT *, this%mdate(:, 1), this%mdate(:, 2)
      IF (sg_mp%isBaseProc()) THEN
        GrapesModelvarIO_500m = GrapesModelvarIO_500m_t(this%ModelvarFns, this%ModelvarCtlFn, this%mdate, sg_mp%vLevel)

        ! get the time_slots_target=-jnk
        ALLOCATE (TimeSlots_target(sg_mp%tSlots))
        CALL IntpLinearTimeSlots(GrapesModelvarIO_500m%time_unix, TimeSlots_target)
        PRINT *, 'TimeSlots_target: ', TimeSlots_target
        PRINT *, 'GrapesModelvarIO_500m%time_unix: ', GrapesModelvarIO_500m%time_unix

        HorNumGridModel = GrapesModelvarIO_500m%idn * GrapesModelvarIO_500m%jdn

        ALLOCATE (cellCntrModel(2, HorNumGridModel), &
                  valueModel(GrapesModelvarIO_500m%kdn, HorNumGridModel, sg_mp%tSlots), &
                  valueDA(sg_mp%vLevel, sg_mp%num_icell_global, sg_mp%tSlots), &
                  valueModel_tmp(GrapesModelvarIO_500m%kdn, HorNumGridModel, GrapesModelvarIO_500m%tdn))
        cellCntrModel(1, :) = RESHAPE(GrapesModelvarIO_500m%lat2D, (/HorNumGridModel/))
        cellCntrModel(2, :) = RESHAPE(GrapesModelvarIO_500m%lon2D, (/HorNumGridModel/))

      END IF
      CALL sg%mpddInfo_sg%bCast(HorNumGridModel)
      CALL sg%mpddInfo_sg%bCast(GrapesModelvarIO_500m%tdn)
      CALL sg%mpddInfo_sg%bCast(GrapesModelvarIO_500m%kdn)
      CALL sg%mpddInfo_sg%bCast(GrapesModelvarIO_500m%idn)
      CALL sg%mpddInfo_sg%bCast(GrapesModelvarIO_500m%jdn)

      ! Ingest the topo and update the height in MOTOR-DA
      ! CALL ModelCoupler%ingest_to_topo_and_update_zHght_d(GrapesModelvarIO_500m%zs * 1.0D0, HorNumGridModel, sg_mp)

      BLOCK
        USE mo_netcdf
        USE parameters_m
        TYPE(NcDataset) :: nc
        TYPE(NcVariable)  :: var
        TYPE(ModelCoupler_t)     :: ModelCouplerForTopo

        REAL(r_kind), ALLOCATABLE :: latTopoFromNC(:), lonTopoFromNC(:), valueTopoFromNC(:, :)
        REAL(r_kind), ALLOCATABLE :: llTopo(:, :), valueTopo(:)
        INTEGER(i_kind) :: lat_len, lon_len

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
          IF(sg%isBaseProc()) THEN
            ALLOCATE (tmp2D(lon_len, lat_len))
                tmp2D = reshape(valueTopo, [lon_len, lat_len])
                CALL guidedfilter(tmp2D, tmp2D, 5, 4.0D0, tmp2D)
                valueTopo = reshape(tmp2D, [lon_len*lat_len])
          END IF

          IF(sg%isBaseProc())DEALLOCATE ( tmp2D)
        END BLOCK

        PRINT *, 'MAXVAL(llTopo)', MAXVAL(llTopo(1, :)) * radian2degree, MAXVAL(llTopo(2, :)) * radian2degree
        PRINT *, 'MINVAL(llTopo)', MINVAL(llTopo(1, :)) * radian2degree, MINVAL(llTopo(2, :)) * radian2degree

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
          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp(sg%num_icell_global))
          END IF

          CALL sg_mp%aggrGridReal1D( sg_mp%topo, tmp, sg_mp%num_icell_global)

          IF(sg_mp%isBaseProc()) THEN
            ALLOCATE (tmp2D(sg_mp%dimCell_global(2), sg_mp%dimCell_global(1)))

                tmp2D = reshape(tmp(:), [sg_mp%dimCell_global(2), sg_mp%dimCell_global(1)])
                CALL guidedfilter(tmp2D, tmp2D, 5, 4.0D0, tmp2D)
                tmp(:) = reshape(tmp2D, [sg_mp%num_icell_global])
          END IF

          CALL sg%DistGridReal1D( tmp, sg_mp%topo, [sg_mp%num_icell_global])          
          IF(sg%isBaseProc())DEALLOCATE (tmp, tmp2D)
        END BLOCK

        CALL sg_mp%update_zHght_from_topo_and_sigma

      END BLOCK

      ! Initialize the model coupler.
      CALL ModelCoupler%Initialize(this%m_configFile, this%geometry)
      CALL ModelCoupler%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)

      IF (sg_mp%isBaseProc()) THEN
        timeModel = GrapesModelvarIO_500m%time_unix
        timeDA = TimeSlots_target
      ELSE
        ALLOCATE (timeModel(GrapesModelvarIO_500m%tdn), timeDA(sg%tSlots))
      END IF
      PRINT *, "ccc: ", GrapesModelvarIO_500m%tdn, sg%tSlots, SHAPE(timeModel), SHAPE(timeDA)

      CALL sg%mpddInfo_sg%bCastReal1D(timeModel)
      CALL sg%mpddInfo_sg%bCastReal1D(timeDA)

      ! Modelvar only work for surface analysis
      PRINT *, "sg_mp%vLevel", sg_mp%vLevel
      IF (sg_mp%vLevel == 1) THEN
        ! Ingest surface_vars
        BLOCK
          REAL(r_kind), ALLOCATABLE :: data_tmp(:, :, :)
          IF (sg%isBaseProc()) THEN
            ALLOCATE (data_tmp(1, HorNumGridModel, GrapesModelvarIO_500m%tdn))
            DEALLOCATE(valueModel)
                ALLOCATE (valueModel(1, HorNumGridModel, sg_mp%tSlots)) 
         END IF

          ! Ingest rainnc
          IF (Xm%getVarIdx('rainnc') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%rainnc)) THEN
                PRINT *, "ERROR!!! There is no 'rainnc' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%rainnc, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('rainnc'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - rainnc has been loaded.'
          END IF

          ! IF (Xm%getVarIdx('pcpa') .NE. 0) THEN
          !   IF (sg_mp%isBaseProc()) THEN
          !     IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%rainnc)) THEN
          !       PRINT *, "ERROR!!! There is no 'rainnc' in GrapesModelvarIO_500m! "
          !       STOP
          !     ELSE
          !       data_tmp = RESHAPE(GrapesModelvarIO_500m%rainnc, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
          !     END IF
          !     CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
          !   END IF

          !   CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pcpa'))%DATA, sg_mp)
          !   PRINT *, 'IOGrapes - rainnc has been loaded.'
          ! END IF

          ! Ingest t2m
          IF (Xm%getVarIdx('t2m') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%t2m)) THEN
                PRINT *, "ERROR!!! There is no 't2m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%t2m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('t2m'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - t2m has been loaded.'
          END IF

          ! Ingest temp
          IF (Xm%getVarIdx('temp') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%t2m)) THEN
                PRINT *, "ERROR!!! There is no 't2m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%t2m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - temp has been loaded.'
          END IF


          ! Ingest q2m
          IF (Xm%getVarIdx('q2m') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%q2m)) THEN
                PRINT *, "ERROR!!! There is no 'q2m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%q2m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('q2m'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - q2m has been loaded.'
          END IF

          ! Ingest q2m
          IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%q2m)) THEN
                PRINT *, "ERROR!!! There is no 'q2m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%q2m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - q2m has been loaded.'
          END IF

          ! Ingest u10m
          IF (Xm%getVarIdx('u10m') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%u10m)) THEN
                PRINT *, "ERROR!!! There is no 'u10m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%u10m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('u10m'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - u10m has been loaded.'
          END IF

          ! Ingest uwnd
          IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%u10m)) THEN
                PRINT *, "ERROR!!! There is no 'u10m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%u10m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - uwnd has been loaded.'
          END IF

          ! Ingest v10m
          IF (Xm%getVarIdx('v10m') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%v10m)) THEN
                PRINT *, "ERROR!!! There is no 'v10m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%v10m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('v10m'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - v10m has been loaded.'
          END IF

          ! Ingest vwnd
          IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%v10m)) THEN
                PRINT *, "ERROR!!! There is no 'v10m' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%v10m, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - vwnd has been loaded.'
          END IF

          ! Ingest vis
          IF (Xm%getVarIdx('vis') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%vis)) THEN
                PRINT *, "ERROR!!! There is no 'vis' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%vis, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vis'))%DATA, sg_mp, 1)
            PRINT *, 'IOGrapes - vis has been loaded.'
          END IF

          ! Ingest psl
          IF (Xm%getVarIdx('psl') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%psl)) THEN
                PRINT *, "ERROR!!! There is no 'psl' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%psl, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('psl'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - psl has been loaded.'
          END IF

          ! Ingest press
          IF (Xm%getVarIdx('pres') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%psl)) THEN
                PRINT *, "ERROR!!! There is no 'psl' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%psl, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
            Xm%Fields(Xm%getVarIdx('pres'))%DATA = Xm%Fields(Xm%getVarIdx('pres'))%DATA * 100.0_r_kind  ! convert to Pa
            PRINT *, 'IOGrapes - pres has been loaded.'
          END IF

          ! Ingest cldt
          IF (Xm%getVarIdx('cldt') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesModelvarIO_500m%cldt)) THEN
                PRINT *, "ERROR!!! There is no 'cldt' in GrapesModelvarIO_500m! "
                STOP
              ELSE
                data_tmp = RESHAPE(GrapesModelvarIO_500m%cldt, (/1, HorNumGridModel, GrapesModelvarIO_500m%tdn/))
              END IF
              CALL interp1d_3D_idx3(timeModel, data_tmp, timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('cldt'))%DATA, sg_mp)
            PRINT *, 'IOGrapes - cldt has been loaded.'
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

          IF (ALLOCATED(data_tmp)) DEALLOCATE (data_tmp)
        END BLOCK
      ELSE
        PRINT *, "GrapesModelvarIO_500m_t run in vertical level other than 1  vLevel: vLevel", sg_mp%vLevel
        ALLOCATE (zRHghtModelAtDA(GrapesModelvarIO_500m%kdn-1, sg_mp%num_cell))
        BLOCK
        REAL(r_single), ALLOCATABLE :: valueTemp(:, :)

        ALLOCATE (valueTemp(GrapesModelvarIO_500m%kdn-1, HorNumGridModel))

        PRINT *, 'GrapesModelvarIO_500m%tdn: ', GrapesModelvarIO_500m%tdn

        IF (sg%isBaseProc()) THEN
              FORALL (i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              valueTemp(:, (j - 1) * GrapesModelvarIO_500m%idn + i) = GrapesModelvarIO_500m%zRHght_u(i, j, :)
              END FORALL
        END IF
        CALL sg%mpddInfo_sg%barrier
        CALL ModelCoupler%ingest_to_field_data_2D(GrapesModelvarIO_500m%kdn-1, &
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
            ALLOCATE (tmpModel_u(GrapesModelvarIO_500m%kdn-1, HorNumGridModel, GrapesModelvarIO_500m%tdn), &
                    tmp_u(GrapesModelvarIO_500m%kdn-1, HorNumGridModel, sg_mp%tSlots))
          END IF

        ! Ingest uwnd
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%uwnd(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - uwnd has been loaded.'
        END IF

        ! Ingest vwnd
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%vwnd(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - vwnd has been loaded.'
        END IF

        ! Ingest temp
        IF (Xm%getVarIdx('temp') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%temp(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - temp has been loaded.'
        END IF

        ! Ingest qvapor
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%qvapor(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - qvapor has been loaded.'
        END IF

        ! Ingest pres
        IF (Xm%getVarIdx('pres') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%pres(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
          ! Xm%Fields(Xm%getVarIdx('pres'))%DATA = Xm%Fields(Xm%getVarIdx('pres'))%DATA * 100.0_r_kind  ! convert to Pa
          PRINT *, 'IOGrapes - pres has been loaded.'
        END IF

        ! Ingest wwnd
        IF (Xm%getVarIdx('wwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:GrapesModelvarIO_500m%tdn, i=1:GrapesModelvarIO_500m%idn, j=1:GrapesModelvarIO_500m%jdn)
              tmpModel_u(:, (j - 1) * GrapesModelvarIO_500m%idn + i, k) = GrapesModelvarIO_500m%wwnd(i, j, :, k)
            END FORALL

            CALL interp1d(timeModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3."
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(GrapesModelvarIO_500m%kdn-1, zRHghtModelAtDA, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('wwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - wwnd has been loaded.'
        END IF
        END BLOCK


        DEALLOCATE (zRHghtModelAtDA)
        CALL sg_mp%mpddInfo_sg%barrier

      END IF
      
      CALL Calc_sg_vars(this%m_configFile, Xm, sg_mp)
      CALL update_restrictionOfStatics(this%m_configFile, Xm, sg_mp)

    END ASSOCIATE

    IF (ALLOCATED(cellCntrModel)) DEALLOCATE (cellCntrModel)
    IF (ALLOCATED(valueModel)) DEALLOCATE (valueModel)
    IF (ALLOCATED(valueDA)) DEALLOCATE (valueDA)
    IF (ALLOCATED(valueModel_tmp)) DEALLOCATE (valueModel_tmp)
    IF (ALLOCATED(TimeSlots_target)) DEALLOCATE (TimeSlots_target)
    IF (ALLOCATED(this%ModelvarFns)) DEALLOCATE (this%ModelvarFns)
    IF (ALLOCATED(this%varNames)) DEALLOCATE (this%varNames)
    IF (ALLOCATED(this%mdate)) DEALLOCATE (this%mdate)

  END SUBROUTINE m_read_bcg_into_Xm
  
END MODULE IOGrapesModelvar_500m_m
  