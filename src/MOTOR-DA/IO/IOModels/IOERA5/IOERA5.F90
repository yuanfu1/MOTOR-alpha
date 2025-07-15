!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IOERA5
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/3/8, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/30, @GBA-MWF, Shenzhen
!     for
!         1. printing more info on saturated qvapor;
!         2. changing gravity acceleration to parameter;
!         3. modified the system call of cp inputfile to output, from -Lvr to -Lv as r requires other options
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE IOERA5_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE YAMLRead_m
  USE Interp1D_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE AdvanceTime_m
  USE MGOpts_m, ONLY: Calc_sg_vars, update_restrictionOfStatics
  USE conversions_m
  USE parameters_m, ONLY: g, degree2radian, rhov_ctl_limit, qv_limit, spec_heat_const_pres, dry_air_gas_const, surface_ref_pres
  USE Filter_m, ONLY: smoothField, guidedfilter

!#define HYDROSTATIC_CHECK

  IMPLICIT NONE

  TYPE, EXTENDS(IOModel_t) :: IOERA5_t
    CHARACTER(LEN=1024) :: m_modelFileName
    INTEGER(r_kind) :: nModelFiles
    TYPE(GrapesIO_t), ALLOCATABLE  :: GrapesIO(:)
    REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA(:, :) !< In dimension of (vLevel, numCell)
    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    REAL(r_kind), ALLOCATABLE :: piTotal(:, :, :)
    REAL(r_kind), ALLOCATABLE :: thTotal(:, :, :)
    REAL(r_kind), ALLOCATABLE :: pip(:, :, :)
    REAL(r_kind), ALLOCATABLE :: thp(:, :, :)
    REAL(r_kind), ALLOCATABLE :: unAmbigRange

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm
    PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
    FINAL :: destructor

  END TYPE IOERA5_t

CONTAINS
  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOERA5_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%initializeIOModel(configFile, geometry)
    CALL this%m_read_config_file
  END SUBROUTINE

  SUBROUTINE m_read_config_file(this)
    !USE NMLRead_m
    USE YAMLRead_m
    CLASS(IOERA5_t) :: this
    CHARACTER(LEN=1024) :: inputDir
    INTEGER :: i, ifile

    PRINT *, 'this%m_configFile ', TRIM(this%m_configFile)

    IF(yaml_get_var(this%m_configFile, 'Verify', 'ERA5_File', this%m_modelFileName) /= 0) THEN
      PRINT *, 'Error: Verify/ERA5_File not found in the configuration file.'
      STOP
    END IF
    

  END SUBROUTINE m_read_config_file

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    IMPLICIT NONE

    CLASS(IOERA5_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(ModelCoupler_t) :: ModelCoupler

    INTEGER(i_kind) :: vLevel_s, vLevel_u, istatus

    INTEGER(i_kind) :: i, j, k, nprocInterp, idn, jdn, kdn, varIdx, it
    REAL(r_kind), ALLOCATABLE :: cc(:, :), gph(:, :, :)

    CHARACTER(len=256) :: grapes_model_type
    INTEGER(i_kind) :: HorNumGridModel
    INTEGER(i_kind) :: nLat, nLon, nLev
    REAL(r_kind) :: scale_factor, add_offset
    REAL(r_kind), ALLOCATABLE :: presLevel(:)

    TYPE(NcDataset) :: nc
    TYPE(NcVariable)  :: var
    TYPE(NcDimension)  :: dim

    CALL this%m_read_config_file

    IF (sg%isBaseProc()) THEN
      BLOCK
        REAL(r_kind), ALLOCATABLE :: latitude(:), longitude(:), level(:)
        nc = NcDataset(this%m_modelFileName, "r")
        var = nc%getVariable("latitude"); CALL var%getData(latitude)
        var = nc%getVariable("longitude"); CALL var%getData(longitude)
        var = nc%getVariable("pressure_level"); CALL var%getData(level)
        var = nc%getVariable("z"); CALL var%getData(gph)
        gph = gph/9.806_R_KIND 
        PRINT *, "gph", MAXVAL(gph), MINVAL(gph)

        dim = nc%getDimension("longitude"); nLon = dim%getLength()
        dim = nc%getDimension("latitude"); nLat = dim%getLength()
        dim = nc%getDimension("pressure_level"); nLev = dim%getLength()

        latitude = latitude * degree2radian
        longitude = longitude * degree2radian

        HorNumGridModel = nLat * nLon
        ALLOCATE (this%cellCntrModel(2, HorNumGridModel))

        PRINT *, 'nLon, nLat, nLev: ', nLon, nLat, nLev
        PRINT *, 'HorNumGridModel: ', HorNumGridModel
        PRINT *, 'pressure_level: ', level, SHAPE(gph)

        FORALL (i=1:nLon, j=1:nLat)
          this%cellCntrModel(1, (j - 1) * nLon + i) = latitude(j)
          this%cellCntrModel(2, (j - 1) * nLon + i) = longitude(i)
        END FORALL

        presLevel = level
        DEALLOCATE (latitude, longitude, level)
      END BLOCK
    END IF

    CALL sg%mpddInfo_sg%bCast(nLat)
    CALL sg%mpddInfo_sg%bCast(nLon)
    CALL sg%mpddInfo_sg%bCast(nLev)
    IF (.NOT. ALLOCATED(presLevel)) ALLOCATE (presLevel(nLev))
    CALL sg%mpddInfo_sg%bCast(presLevel)
    PRINT *, 'presLevel: ', presLevel

    ! Initialize the model coupler.
    CALL ModelCoupler%initialize(this%m_configFile, this%geometry)
    PRINT*,'here'
    CALL ModelCoupler%gen_interp_coeffs(this%cellCntrModel, HorNumGridModel, sg)

    PRINT *, 'Finish initialize the ModelCoupler for ERA5.'

    ! ! Ingest surface parameter
    ! DO it = 1, sg%tSlots
    !   CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%snowc(:, :), sg%snowc(:, it), sg)
    !   CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%xice(:, :), sg%xice(:, it), sg)
    ! END DO
    ! CALL ModelCoupler%ingest_to_field_data_1D_nn(HorNumGridModel, GrapesIO(1)%hgrid%xland(:, :), sg%landmask, sg)
    ! CALL ModelCoupler%ingest_to_field_data_1D_nn(HorNumGridModel, GrapesIO(1)%hgrid%soil_type(:, :), &
    !                                              sg%soil_type, sg)
    ! CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%veg_fraction(:, :), sg%veg_frac, sg)

    ! FORALL (i=1:sg%num_cell, sg%landmask(i) == 2.0D0) sg%landmask(i) = 0.0D0

    BLOCK
      !   ! Generate the model height at location of MOTOR-DA cell
      REAL(r_single) :: zRHght_temp(nLev, HorNumGridModel)
      !   REAL(r_single) :: zRHght_u_temp(GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1, HorNumGridModel)

      !   ! PRINT *, 'Here generate the model height at location of MOTOR-DA cell'
      IF (sg%isBaseProc()) THEN
        !     vLevel_s = GrapesIO(1)%kde_s - GrapesIO(1)%kds_s + 1
        !     vLevel_u = GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1

        FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
          zRHght_temp(k, (j - 1) * nLon + i) = gph(i, j, k)
        END FORALL
      END IF

      ALLOCATE (this%zRHghtModelAtDA(nLev, sg%num_cell)); 
      CALL ModelCoupler%ingest_to_field_data_2D(nLev, zRHght_temp, this%zRHghtModelAtDA, sg)
    END BLOCK


    ! PRINT *, 'Start ingesting the model variables.'
    BLOCK
      REAL(r_single), ALLOCATABLE :: tmpModel(:, :)

      REAL(r_single), ALLOCATABLE :: tmpData_FromFile(:, :, :)
      REAL(r_kind), ALLOCATABLE :: tmpData_ModelAtDAvLevel(:, :)

      IF (sg%isBaseProc()) THEN
        ALLOCATE (tmpModel(nLev, HorNumGridModel))
      END IF

      IF (sg%isActiveProc()) THEN
        ALLOCATE (tmpData_ModelAtDAvLevel(nLev, sg%num_cell))
      END IF

      !   ! CALL sg%mpddInfo_sg%barrier
      !   ! Ingest uwnd
      IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("u"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)

          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('uwnd'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - uwnd has been loaded.'
      END IF

      IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("v"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)

          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('vwnd'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - vwnd has been loaded.'
      END IF

      IF (Xm%getVarIdx('temp') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("t"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)

          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('temp'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - temp has been loaded.'
      END IF

      IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("q"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)
          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - qvapor has been loaded.'
      END IF

      IF (Xm%getVarIdx('qcloud') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("clwc"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)
          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('qcloud'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - qcloud has been loaded.'
      END IF

      IF (Xm%getVarIdx('qice') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          var = nc%getVariable("ciwc"); CALL var%getData(tmpData_FromFile)
          tmpData_FromFile = tmpData_FromFile! * scale_factor + add_offset

          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = tmpData_FromFile(i, j, k)
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)
          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('qice'))%DATA(:, i, sg%tSlots))
          END DO
        END IF

        PRINT *, 'IOERA5 - qice has been loaded.'
      END IF

      IF (Xm%getVarIdx('pres') .NE. 0) THEN
        IF (sg%isBaseProc()) THEN
          FORALL (i=1:nLon, j=1:nLat, k=1:nLev)
            tmpModel(k, (j - 1) * nLon + i) = presLevel(k) * 100.0_R_KIND
          END FORALL
        END IF

        IF (sg%isActiveProc()) THEN
          CALL ModelCoupler%ingest_to_field_data_2D(nLev, tmpModel, tmpData_ModelAtDAvLevel, sg)

          DO i = 1, sg%num_cell
            CALL interp1d(this%zRHghtModelAtDA(:, i), tmpData_ModelAtDAvLevel(:, i), &
                          sg%zHght(:, i), Xm%Fields(Xm%getVarIdx('pres'))%DATA(:, i, sg%tSlots), 'LOG')
          END DO
        END IF

        PRINT *, 'IOERA5 - pres has been loaded.'
      END IF

    END BLOCK

    CALL sg%mpddInfo_sg%barrier
    ! Deallocate memories
  END SUBROUTINE m_read_bcg_into_Xm

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(IOERA5_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%GrapesIO)) DEALLOCATE (this%GrapesIO)
    IF (ALLOCATED(this%zRHghtModelAtDA)) DEALLOCATE (this%zRHghtModelAtDA)
    IF (ALLOCATED(this%cellCntrModel)) DEALLOCATE (this%cellCntrModel)
    PRINT *, 'Finish the destructor of IOERA5_t.'

  END SUBROUTINE destructor

END MODULE IOERA5_m
