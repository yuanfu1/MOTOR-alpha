!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-08-09, created by JiongmingPang for reading postvar.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-08-09, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE IOGrapesPostvar_m
  USE kinds_m
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE YAMLRead_m
  USE AdvanceTime_m
  USE GrapesPostvarIO_m, ONLY: GrapesPostvarIO_t
  USE InterpolateTime_m
  USE Widgets_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE MGOpts_m, ONLY: Calc_sg_vars, update_restrictionOfStatics
  USE Interp1D_m
  USE conversions_m
  USE parameters_m

  TYPE, EXTENDS(IOModel_t) :: IOGrapesPostvar_t
    CHARACTER(LEN=1024) :: nmlDAFilename, postvarCtlFn, DEMFileName
    CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: postvarFns
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: varNames
    INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE :: mdate
    INTEGER(i_kind) :: numFiles

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm

  END TYPE IOGrapesPostvar_t

CONTAINS

  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOGrapesPostvar_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%initializeIOModel(configFile, geometry)
    this%nmlDAFilename = TRIM(configFile)
  END SUBROUTINE

  SUBROUTINE m_read_config_file(this)
    CLASS(IOGrapesPostvar_t) :: this

    CHARACTER(LEN=1024), ALLOCATABLE :: postvarFns(:)
    CHARACTER(LEN=20), ALLOCATABLE :: inVarNames(:)
    CHARACTER(LEN=10), ALLOCATABLE :: timestr(:)
    CHARACTER(LEN=1024) :: configFile
    CHARACTER(LEN=1024) :: inFileDir
    CHARACTER(LEN=1024) :: ensDirtmp
    CHARACTER(LEN=3)    :: memid
    INTEGER(i_kind)     :: t, inFilesNum, ifile, inTimeNum
    LOGICAL             :: ensFlag

    ifile = yaml_get_var(this%m_configFile, 'IO', 'input_dir_model', inFileDir)
    ifile = yaml_get_var(this%m_configFile, 'IO', 'ModelFileName', postvarFns)
    inFilesNum = UBOUND(postvarFns, 1)
    PRINT *, "inFilesNum: ", inFilesNum
    this%numFiles = inFilesNum
    ALLOCATE (this%postvarFns(inFilesNum))

    ! get the full position of postvar files
    DO t = 1, inFilesNum
      PRINT *, 'postvarFns: ', TRIM(postvarFns(t))
      this%postvarFns(t) = TRIM(inFileDir)//"/"//TRIM(postvarFns(t))
      PRINT *, 'postvarFns: ', TRIM(this%postvarFns(t))
    END DO

    ! get the full position of postvar control files (post.ctl)
    this%postvarCtlFn = TRIM(inFileDir)//"/post.ctl"
    PRINT *, 'postvarCtlFn: ', TRIM(this%postvarCtlFn)

    ! get the full position of dem files (DL_DEMdata.csv)
    this%DEMFileName = TRIM(inFileDir)//"/DL_DEMdata.nc"
    PRINT *, 'DEMFileName: ', TRIM(this%DEMFileName)

    ! get the times of postvar files
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
      this%mdate(5, t) = 0 !min
      this%mdate(6, t) = 0 !sec
      PRINT *, "GMT Time in t-step: ", t, this%mdate(:, t)
    END DO

    ! get the analysis variables
    ifile = yaml_get_var(TRIM(this%nmlDAFilename), 'modelState', 'varList', this%varNames)

    IF (ALLOCATED(postvarFns)) DEALLOCATE (postvarFns)
    IF (ALLOCATED(timestr)) DEALLOCATE (timestr)

  END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    USE GrapesPostvarIO_m, ONLY: GrapesPostvarIO_t
    USE InterpolateTime_m
    USE ModelCoupler_m, ONLY: ModelCoupler_t

    IMPLICIT NONE

    CLASS(IOGrapesPostvar_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(GrapesPostvarIO_t) :: GrapesPostvarIO
    TYPE(ModelCoupler_t)     :: ModelCoupler

    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    REAL(r_kind), ALLOCATABLE :: cellCntrDA(:, :)
    REAL(r_kind), ALLOCATABLE :: valueModel(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueModel_tmp(:, :, :)   ! for time-scale interpolation, Jiongming Pang
    REAL(r_kind), ALLOCATABLE :: valueDA(:, :, :)

    INTEGER(i_kind), ALLOCATABLE :: TimeSlots_target(:)   ! for time-scale interpolation, Jiongming Pang

    INTEGER(i_kind) :: i, j, k, HorNumGridModel, nprocInterp
    TYPE(InterpolateData_Time_t), ALLOCATABLE :: IntpTime

    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind)   :: iv, varIndx
    REAL(r_kind), ALLOCATABLE :: timeModel(:), timeDA(:)
    REAL(r_kind) :: u, v
    REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA(:, :, :)

    ! CHARACTER(LEN=1024) :: OUTPUT_DIR
    ! CHARACTER(LEN=1024) :: iobakoutfnnc
    ! CHARACTER(LEN=1024) :: iobakoutfntxt

    ASSOCIATE (sg_mp => sg)

      ! IF (sg_mp%isBaseProc()) THEN
      CALL this%m_read_config_file()
      ! END IF

      CALL sg_mp%mpddInfo_sg%bCast(this%numFiles)
      IF (this%numFiles == 0) THEN
        PRINT *, 'No file listed in the configurations.'
        RETURN
      END IF
      PRINT *, '*************************Number of files to be read: ', this%numFiles
      PRINT *, this%mdate(:, 1), this%mdate(:, 2)

      IF (sg_mp%isBaseProc()) THEN
        GrapesPostvarIO = GrapesPostvarIO_t(this%postvarFns, this%postvarCtlFn, this%DEMFileName, this%mdate)

        ! get the time_slots_target=-jnk
        ALLOCATE (TimeSlots_target(sg_mp%tSlots))
        CALL IntpLinearTimeSlots(GrapesPostvarIO%time_unix, TimeSlots_target)
        PRINT *, 'GrapesPostvarIO%time_unix: ', GrapesPostvarIO%time_unix

        HorNumGridModel = GrapesPostvarIO%idn * GrapesPostvarIO%jdn

        ALLOCATE (cellCntrModel(2, HorNumGridModel), &
                  valueModel(GrapesPostvarIO%kdn, HorNumGridModel, sg_mp%tSlots), &
                  valueDA(sg_mp%vLevel, sg_mp%num_icell_global, sg_mp%tSlots), &
                  valueModel_tmp(GrapesPostvarIO%kdn, HorNumGridModel, GrapesPostvarIO%tdn))
        cellCntrModel(1, :) = RESHAPE(GrapesPostvarIO%lat2D, (/HorNumGridModel/))
        cellCntrModel(2, :) = RESHAPE(GrapesPostvarIO%lon2D, (/HorNumGridModel/))

      END IF
      CALL sg%mpddInfo_sg%bCast(HorNumGridModel)
      CALL sg%mpddInfo_sg%bCast(GrapesPostvarIO%tdn)
      CALL sg%mpddInfo_sg%bCast(GrapesPostvarIO%kdn)
      CALL sg%mpddInfo_sg%bCast(GrapesPostvarIO%idn)
      CALL sg%mpddInfo_sg%bCast(GrapesPostvarIO%jdn)

      ! Ingest the topo and update the height in MOTOR-DA
      ! CALL ModelCoupler%ingest_to_topo_and_update_zHght_d(GrapesPostvarIO%zs * 1.0D0, HorNumGridModel, sg_mp)
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

            ! if(latTopoFromNC(i) > 25.4 .and. latTopoFromNC(i) < 25.6 .and. lonTopoFromNC(j) > 100.5 .and. lonTopoFromNC(j) < 100.6) then
            !   valueTopo((j - 1)*lat_len + i) = 800.0
            ! ELSE
            !   valueTopo((j - 1)*lat_len + i) = 0.0
            ! end if

          END DO
        END DO

        PRINT *, 'MAXVAL(llTopo)', MAXVAL(llTopo(1, :)) * radian2degree, MAXVAL(llTopo(2, :)) * radian2degree
        PRINT *, 'MINVAL(llTopo)', MINVAL(llTopo(1, :)) * radian2degree, MINVAL(llTopo(2, :)) * radian2degree

        CALL ModelCouplerForTopo%Initialize(this%m_configFile, this%geometry)
        CALL ModelCouplerForTopo%gen_interp_coeffs(llTopo, lat_len * lon_len, sg)
        ! -----------------------------------------------------------
        CALL ModelCouplerForTopo%ingest_to_topo_and_update_zHght_d(valueTopo*1.0D0, lat_len*lon_len, sg_mp)
        ! -----------------------------------------------------------

        ! CALL ModelCouplerForTopo%interhpN%interp_singleLevel_d(valueTopo * 1.0D0, sg_mp%topo)
        ! CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
        ! CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)
        ! -----------------------------------------------------------

        PRINT *, 'here'
        ! STOP
      END BLOCK

      ! Initialize the model coupler.
      CALL ModelCoupler%Initialize(this%m_configFile, this%geometry)
      CALL ModelCoupler%gen_interp_coeffs(cellCntrModel, HorNumGridModel, sg_mp)

      IF (sg_mp%isBaseProc()) THEN
        timeModel = GrapesPostvarIO%time_unix
        timeDA = TimeSlots_target
      ELSE
        ALLOCATE (timeModel(GrapesPostvarIO%tdn), timeDA(sg%tSlots))
      END IF
      PRINT *, "ccc: ", GrapesPostvarIO%tdn, sg%tSlots, SHAPE(timeModel), SHAPE(timeDA)

      CALL sg%mpddInfo_sg%bCastReal1D(timeModel)
      CALL sg%mpddInfo_sg%bCastReal1D(timeDA)
      ! PRINT *, 'aaa: ', timeModel
      ! PRINT*, 'bbb: ', timeDA

      ! STOP

      ! postvar only work for surface analysis
      PRINT *, "sg_mp%vLevel", sg_mp%vLevel
      IF (sg_mp%vLevel == 1) THEN
        !! Ingest the 2D background field

        ! Ingest t2m-temperature at 2m
        IF (Xm%getVarIdx('temp') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesPostvarIO%t2m)) THEN
              PRINT *, "ERROR!!! There is no 't2m' in GrapesPostvar! "
              STOP
            ELSE
              valueModel_tmp = RESHAPE(GrapesPostvarIO%t2m, (/1, HorNumGridModel, GrapesPostvarIO%tdn/))
            END IF

            ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('temp'))%DATA, sg_mp)
          PRINT *, 'IOGrapesPostvar - temp has been loaded.', Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, 1, 1), &
            sg_mp%cell_cntr(:, 1), sg_mp%mpddInfo_sg%myrank
        END IF

        ! Ingest uwnd at 10m
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesPostvarIO%u10m)) THEN
              PRINT *, "ERROR!!! There is no 'u10m' in GrapesPostvar! "
              STOP
            ELSE
              valueModel_tmp = RESHAPE(GrapesPostvarIO%u10m, (/1, HorNumGridModel, GrapesPostvarIO%tdn/))
            END IF

            ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesPostvar - u10m has been loaded.'
        END IF

        ! Ingest vwnd at 10m
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesPostvarIO%v10m)) THEN
              PRINT *, "ERROR!!! There is no 'v10m' in GrapesPostvar! "
              STOP
            ELSE
              valueModel_tmp = RESHAPE(GrapesPostvarIO%v10m, (/1, HorNumGridModel, GrapesPostvarIO%tdn/))
            END IF

            ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapesPostvar - v10m has been loaded.'
        END IF

        ! Ingest q2m-qvapor at 2m
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            IF (.NOT. ALLOCATED(GrapesPostvarIO%q2m)) THEN
              PRINT *, "ERROR!!! There is no 'q2m' in GrapesPostvar! "
              STOP
            ELSE
              valueModel_tmp = RESHAPE(GrapesPostvarIO%q2m, (/1, HorNumGridModel, GrapesPostvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)
          PRINT *, 'IOGrapesPostvar - q2m has been loaded.'
        END IF

        ! Ingest pressure at ground
        IF (Xm%getVarIdx('pres') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            ! IF (getStrIndex(this%varNames, 'psfc') < 0) THEN
            IF (.NOT. ALLOCATED(GrapesPostvarIO%psfc)) THEN
              PRINT *, "ERROR!!! There is no 'psfc' in GrapesPostvar! "
              STOP
            ELSE
              valueModel_tmp = RESHAPE(GrapesPostvarIO%psfc, (/1, HorNumGridModel, GrapesPostvarIO%tdn/))
            END IF

            ! ! for time-scale interpolation
            CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                  timeDA, valueModel)
          END IF
          CALL ModelCoupler%ingest_to_field_data_MT(sg_mp%vLevel, valueModel, Xm%Fields(Xm%getVarIdx('pres'))%DATA, sg_mp)
          PRINT *, 'IOGrapesPostvar - psfc has been loaded.'
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
        PRINT *, "GrapesPostvar_t run in vertical level other than 1  vLevel: vLevel", sg_mp%vLevel

        ALLOCATE (zRHghtModelAtDA(GrapesPostvarIO%kdn, sg_mp%num_cell, sg%tSlots))
        BLOCK
          REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA_temp(:, :, :)
          REAL(r_kind), ALLOCATABLE :: valueTemp(:, :, :)

          ALLOCATE (zRHghtModelAtDA_temp(GrapesPostvarIO%kdn, sg_mp%num_cell, GrapesPostvarIO%tdn))
          ALLOCATE (valueTemp(GrapesPostvarIO%kdn, HorNumGridModel, GrapesPostvarIO%tdn))
          ! PRINT *, SHAPE(GrapesPostvarIO%H)

          PRINT *, 'GrapesPostvarIO%tdn: ', GrapesPostvarIO%tdn

          IF (sg%isBaseProc()) THEN
          DO k = 1, GrapesPostvarIO%tdn
            FORALL (i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
              valueTemp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%H(i, j, :, k)
            END FORALL
          END DO
          END IF
          CALL sg%mpddInfo_sg%barrier
          CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, &
                                                    valueTemp, zRHghtModelAtDA_temp, sg_mp, GrapesPostvarIO%tdn)          ! CALL sg%mpddInfo_sg%barrier
          CALL interp1d_3D_idx3(timeModel, zRHghtModelAtDA_temp, timeDA, zRHghtModelAtDA)

          PRINT *, "zRHghtModelAtDA", MAXVAL(zRHghtModelAtDA_temp), timeModel, timeDA
          ! STOP
          DEALLOCATE (zRHghtModelAtDA_temp, valueTemp)
        END BLOCK
        !

        BLOCK
          REAL(r_kind), ALLOCATABLE :: valueDA_mp(:, :, :)
          ALLOCATE (valueDA_mp(GrapesPostvarIO%kdn, sg_mp%num_cell, sg_mp%tSlots))

          ! Ingest uwnd
          IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesPostvarIO%uwnd)) THEN
                PRINT *, "ERROR!!! There is no 'uwnd' in GrapesPostvar! "
                STOP
              ELSE
                PRINT *, SHAPE(GrapesPostvarIO%uwnd)
                PRINT *, SHAPE(valueModel_tmp)
                FORALL (k=1:GrapesPostvarIO%tdn, i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
                  valueModel_tmp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%uwnd(i, j, :, k)
                END FORALL
              END IF
              PRINT *, SHAPE(valueModel_tmp)
              PRINT *, SHAPE(valueModel), timeModel, timeDA
              ! ! for time-scale interpolation

              CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                    timeDA, valueModel)
            END IF

            CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, valueModel, valueDA_mp, sg_mp)
            DO j = 1, sg_mp%tSlots
              DO i = 1, sg_mp%num_cell
                CALL interp1d(zRHghtModelAtDA(:, i, j), valueDA_mp(:, i, j), &
                              sg_mp%zHght(:, i), Xm%Fields(Xm%getVarIdx('uwnd'))%DATA(:, i, j))
              END DO
            END DO

            PRINT *, 'IOGrapesPostvar - uwnd has been loaded.'
          END IF

          ! Ingest uwnd
          IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesPostvarIO%vwnd)) THEN
                PRINT *, "ERROR!!! There is no 'vwnd' in GrapesPostvar! "
                STOP
              ELSE
                FORALL (k=1:GrapesPostvarIO%tdn, i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
                  valueModel_tmp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%vwnd(i, j, :, k)
                END FORALL
              END IF

              ! ! for time-scale interpolation
              CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                    timeDA, valueModel)
            END IF
            CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, valueModel, valueDA_mp, sg_mp)
            DO j = 1, sg_mp%tSlots
              DO i = 1, sg_mp%num_cell
                CALL interp1d(zRHghtModelAtDA(:, i, j), valueDA_mp(:, i, j), &
                              sg_mp%zHght(:, i), Xm%Fields(Xm%getVarIdx('vwnd'))%DATA(:, i, j))
              END DO
            END DO
            PRINT *, 'IOGrapesPostvar - vwnd has been loaded.'
          END IF

          ! Ingest temp
          IF (Xm%getVarIdx('temp') .NE. 0 .OR. Xm%getVarIdx('theta') .NE. 0) THEN

            CALL Xm%addVar('temp')

            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesPostvarIO%temp)) THEN
                PRINT *, "ERROR!!! There is no 'temp' in GrapesPostvar! "
                STOP
              ELSE
                FORALL (k=1:GrapesPostvarIO%tdn, i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
                  valueModel_tmp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%temp(i, j, :, k)
                END FORALL
              END IF

              ! ! for time-scale interpolation
              CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                    timeDA, valueModel)
            END IF
            CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, valueModel, valueDA_mp, sg_mp)
            ASSOCIATE (temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA)
            DO j = 1, sg_mp%tSlots
              DO i = 1, sg_mp%num_cell
                CALL interp1d(zRHghtModelAtDA(:, i, j), valueDA_mp(:, i, j), &
                              sg_mp%zHght(:, i), temp(:, i, j))

                ! Temp below the ground is set to 0.6 degree increase per 100m
                DO k = sg_mp%vLevel, 1, -1
                  IF (sg_mp%zHght(k, i) < zRHghtModelAtDA(1, i, j)) THEN
                    temp(k, i, j) = temp(k + 1, i, j) + 0.60D0 * (sg_mp%zHght(k + 1, i) - sg_mp%zHght(k, i)) / 100.0D0
                    ! PRINT*, 'temp(k, i, j) ', temp(k, i, j) , temp(k+1, i, j)
                  END IF
                END DO
                ! -----------------------------------------------------------

              END DO
            END DO
            END ASSOCIATE
            PRINT *, 'IOGrapesPostvar - temp has been loaded.'
          END IF

          ! Ingest pres
          IF (Xm%getVarIdx('pres') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesPostvarIO%pres)) THEN
                PRINT *, "ERROR!!! There is no 'pres' in GrapesPostvar! "
                STOP
              ELSE
                FORALL (k=1:GrapesPostvarIO%tdn, i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
                  valueModel_tmp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%pres(i, j, :, k)
                END FORALL
              END IF

              ! ! for time-scale interpolation
              CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                    timeDA, valueModel)
            END IF
            CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, valueModel, valueDA_mp, sg_mp)
            ASSOCIATE (pres => Xm%Fields(Xm%getVarIdx('pres'))%DATA, temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA)

              DO j = 1, sg_mp%tSlots
                DO i = 1, sg_mp%num_cell
                  CALL interp1d(zRHghtModelAtDA(:, i, j), valueDA_mp(:, i, j), &
                                sg_mp%zHght(:, i), Xm%Fields(Xm%getVarIdx('pres'))%DATA(:, i, j), 'LOG')
                  ! Temp below the ground is set to 0.6 degree increase per 100m
                  DO k = sg_mp%vLevel, 1, -1
                    IF (sg_mp%zHght(k, i) < zRHghtModelAtDA(1, i, j)) THEN
                      pres(k, i, j) = pres(k + 1, i, j) * EXP(-(sg_mp%zHght(k, i) - sg_mp%zHght(k + 1, i)) * g &
                                                              / (dry_air_gas_const * (temp(k, i, j) + temp(k + 1, i, j)) / 2.0D0))
                      ! PRINT*, 'pres(k, i, j) ', pres(k, i, j) , pres(k+1, i, j)
                    END IF
                  END DO
                  ! -----------------------------------------------------------
                END DO
              END DO
            END ASSOCIATE
            PRINT *, 'IOGrapesPostvar - pres has been loaded.'
          END IF

          ! Ingest qvapor
          IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
            IF (sg_mp%isBaseProc()) THEN
              IF (.NOT. ALLOCATED(GrapesPostvarIO%qvapor)) THEN
                PRINT *, "ERROR!!! There is no 'qvapor' in GrapesPostvar! "
                STOP
              ELSE
                FORALL (k=1:GrapesPostvarIO%tdn, i=1:GrapesPostvarIO%idn, j=1:GrapesPostvarIO%jdn)
                  valueModel_tmp(:, (j - 1) * GrapesPostvarIO%idn + i, k) = GrapesPostvarIO%qvapor(i, j, :, k)
                END FORALL
              END IF

              ! ! for time-scale interpolation
              CALL interp1d_3D_idx3(timeModel, valueModel_tmp, &
                                    timeDA, valueModel)
            END IF
            CALL ModelCoupler%ingest_to_field_data_MT(GrapesPostvarIO%kdn, valueModel, valueDA_mp, sg_mp)
            DO j = 1, sg_mp%tSlots
              DO i = 1, sg_mp%num_cell
                CALL interp1d(zRHghtModelAtDA(:, i, j), valueDA_mp(:, i, j), &
                              sg_mp%zHght(:, i), Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(:, i, j))
              END DO
            END DO
            PRINT *, 'IOGrapesPostvar - qvapor has been loaded.'
          END IF

        END BLOCK

        ! Xm%Fields(Xm%getVarIdx('uwnd'))%DATA = 10.0;
        ! Xm%Fields(Xm%getVarIdx('vwnd'))%DATA = 0.0;

        ! IF (Xm%getVarIdx('theta') .NE. 0 .AND. Xm%getVarIdx('theta_0') .NE. 0) THEN

        !   CALL Xm%addVar('theta')
        !   CALL Xm%addVar('theta_0')

        !   ASSOCIATE (theta => Xm%Fields(Xm%getVarIdx('theta'))%DATA, &
        !              theta_0 => Xm%Fields(Xm%getVarIdx('theta_0'))%DATA, &
        !              temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA, &
        !              pres => Xm%Fields(Xm%getVarIdx('pres'))%DATA)

        !     DO j = 1, sg_mp%tSlots
        !       DO i = 1, sg_mp%num_cell
        !         DO k = 1, sg_mp%vLevel
        !           theta(k, i, j) = temp(k, i, j) * (1000.0D0 / pres(k, i, j))**0.28557214D0
        !           theta_0(k, i, j) = theta(k, i, j) * (1000.0D0 / 1000.0D0)**0.28557214D0
        !         END DO
        !       END DO
        !     END DO

        !   END ASSOCIATE
        !   CALL Xm%rmVar('temp')
        ! END IF

        ! DO i = 1, SIZE(Xm%Fields)
        !   DO j = 1, sg%num_cell
        !     DO k = 1, sg%vLevel
        !       IF (sg%zHght(k, j) < sg%topo(j)) THEN
        !         Xm%Fields(i)%DATA(k, j, :) = 0.0D0
        !       END IF
        !     END DO
        !   END DO
        ! END DO

      END IF

      CALL Calc_sg_vars(this%m_configFile, Xm, sg_mp)
      CALL update_restrictionOfStatics(this%m_configFile, Xm, sg_mp)

      DEALLOCATE (zRHghtModelAtDA)
      CALL sg_mp%mpddInfo_sg%barrier

      

    END ASSOCIATE

    IF (ALLOCATED(cellCntrModel)) DEALLOCATE (cellCntrModel)
    IF (ALLOCATED(valueModel)) DEALLOCATE (valueModel)
    IF (ALLOCATED(valueDA)) DEALLOCATE (valueDA)
    IF (ALLOCATED(valueModel_tmp)) DEALLOCATE (valueModel_tmp)
    IF (ALLOCATED(TimeSlots_target)) DEALLOCATE (TimeSlots_target)
    IF (ALLOCATED(this%postvarFns)) DEALLOCATE (this%postvarFns)
    IF (ALLOCATED(this%varNames)) DEALLOCATE (this%varNames)
    IF (ALLOCATED(this%mdate)) DEALLOCATE (this%mdate)

    ! IF (ALLOCATED(this%inFileNames)) DEALLOCATE (this%inFileNames)

  END SUBROUTINE m_read_bcg_into_Xm

END MODULE IOGrapesPostvar_m
