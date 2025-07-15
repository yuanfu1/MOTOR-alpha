!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IOGrapes
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
MODULE IOGrapes_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE IOModel_m, ONLY: IOModel_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres, degree2radian, rhov_ctl_limit, qv_limit, g
  USE YAMLRead_m
  USE Interp1D_m
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE AdvanceTime_m
  USE MGOpts_m, ONLY: Calc_sg_vars, update_restrictionOfStatics
  USE conversions_m
  USE Filter_m, ONLY: smoothField, guidedfilter

!#define HYDROSTATIC_CHECK

  IMPLICIT NONE

  TYPE, EXTENDS(IOModel_t) :: IOGrapes_t
    CHARACTER(LEN=1024), ALLOCATABLE :: m_nmlFilename(:), m_modelFileName(:), m_qcqrFileName(:)
    INTEGER(r_kind) :: nModelFiles
    TYPE(GrapesIO_t), ALLOCATABLE  :: GrapesIO(:)
    REAL(r_kind), ALLOCATABLE :: zRHghtModelAtDA_s(:, :), & !< In dimension of (vLevel, numCell)
                                 zRHghtModelAtDA_u(:, :)!
    REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
    INTEGER(i_kind) :: HorNumGridModel
    REAL(r_kind), ALLOCATABLE :: piTotal(:, :, :)
    REAL(r_kind), ALLOCATABLE :: thTotal(:, :, :)
    REAL(r_kind), ALLOCATABLE :: pip(:, :, :)
    REAL(r_kind), ALLOCATABLE :: thp(:, :, :)
    REAL(r_kind), ALLOCATABLE :: unAmbigRange
    LOGICAL :: CloudyMode = .FALSE. ! If set to .FALSE., qcqr files will not be required. Default is .FALSE.

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PUBLIC, PASS(this) :: m_read_bcg_into_Xm
    PROCEDURE, PRIVATE, PASS(this) :: m_read_config_file
    PROCEDURE, PUBLIC, PASS(this) :: m_write_Xm_into_bcg
    FINAL :: destructor

  END TYPE IOGrapes_t

CONTAINS
  SUBROUTINE initialize(this, configFile, geometry)
    IMPLICIT NONE
    CLASS(IOGrapes_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%initializeIOModel(configFile, geometry)
  END SUBROUTINE

  SUBROUTINE m_read_config_file(this)
    !USE NMLRead_m
    USE YAMLRead_m
    CLASS(IOGrapes_t) :: this
    CHARACTER(LEN=1024) :: inputDir
    INTEGER :: i, ifile

    PRINT *, 'this%m_configFile ', TRIM(this%m_configFile)

    !CALL nml_read_string_array(this%m_configFile, 'NMLFileName', this%m_nmlFilename)
    !CALL nml_read_string_array(this%m_configFile, 'ModelFileName', this%m_modelFileName)
    ifile = yaml_get_var(this%m_configFile, 'IO', 'NMLFileName', this%m_nmlFilename)
    ifile = yaml_get_var(this%m_configFile, 'IO', 'ModelFileName', this%m_modelFileName)
    ifile = yaml_get_var(this%m_configFile, 'RTTOV', 'rttov_clouds', this%CloudyMode)
    IF (this%CloudyMode) ifile = yaml_get_var(this%m_configFile, 'IO', 'HydroFileName', this%m_qcqrFileName)
    ifile = yaml_get_var(this%m_configFile, 'IO', 'input_dir_model', inputDir)

    this%nModelFiles = size(this%m_nmlFilename)
    DO i = 1, this%nModelFiles
      this%m_nmlFilename(i) = TRIM(inputDir)//'/'//TRIM(this%m_nmlFilename(i))
      this%m_modelFileName(i) = TRIM(inputDir)//'/'//TRIM(this%m_modelFileName(i))
      IF (this%CloudyMode) this%m_qcqrFileName(i) = TRIM(inputDir)//'/'//TRIM(this%m_qcqrFileName(i))
      PRINT *, 'nmlFilename for grapes: ', TRIM(this%m_nmlFilename(i))
      PRINT *, 'modelFileName for grapes: ', TRIM(this%m_modelFileName(i))
      IF (this%CloudyMode) PRINT *, 'qcqrFileName for grapes: ', TRIM(this%m_qcqrFileName(i))
    END DO

  END SUBROUTINE m_read_config_file

!> @brief
!! Build the restriction matrix for the aggregation 1
! @see
! @note
! @warning
! @attention
  SUBROUTINE m_write_Xm_into_bcg(this, Xm, sg, Xb)
    USE NMLRead_m
    USE Interp1D_m, ONLY: interp1d
    USE RectangleGridUtility_m, ONLY: RectangleGridUtility_t

    CLASS(IOGrapes_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    CHARACTER(LEN=1024) :: outputFile, outputFileName, inputFile
    CHARACTER(LEN=1024) :: OutputQcqrName, inputQcqr, outputQcqr
    INTEGER(i_kind) :: i, j, k, err, fileStat
    TYPE(State_t), INTENT(IN) :: Xb

    IF (SIZE(this%m_modelFileName) == 0) CALL this%m_read_config_file
    PRINT *, 'Input file is: ', TRIM(this%m_modelFileName(SIZE(this%m_modelFileName)))
    inputFile = TRIM(this%m_modelFileName(SIZE(this%m_modelFileName)))
    IF (this%CloudyMode) inputQcqr = TRIM(this%m_qcqrFileName(SIZE(this%m_qcqrFileName)))

    fileStat = yaml_get_var(this%m_configFile, 'IO', 'output_dir', outputFile)
    fileStat = yaml_get_var(this%m_configFile, 'IO', 'OutputFileName', outputFileName)
    IF (this%CloudyMode) fileStat = yaml_get_var(this%m_configFile, 'IO', 'OutputQcqrName', OutputQcqrName)
    IF (this%CloudyMode) outputQcqr = TRIM(outputFile)//'/'//TRIM(OutputQcqrName)
    outputFile = TRIM(outputFile)//'/'//TRIM(outputFileName)
    PRINT *, 'Output file is: ', TRIM(outputFile)
    IF (this%CloudyMode) PRINT *, 'Output file is: ', TRIM(OutputQcqrName)

    ASSOCIATE (sg_mp => sg, HorNumGridModel => this%HorNumGridModel, GrapesIO => this%GrapesIO(this%nModelFiles))

      ! Copy the old file.
      IF (sg%isBaseProc()) THEN
        CALL SYSTEM('cp -Lv '//TRIM(inputFile)//' '//TRIM(outputFile))
        IF (this%CloudyMode) CALL SYSTEM('cp -Lv '//TRIM(inputQcqr)//' '//TRIM(outputQcqr))

        ! Check if the copy of the input file successful: Yuanfu Xie 2022/10/30
        IF (err .EQ. 0) THEN
          PRINT *, 'Done copy the output file: ', err
        ELSE
          WRITE (*, 3) TRIM(outputFileName)
3         FORMAT('Copying file: ', A, ' failed, please chmod your output directory to all can write and rerun!')
        END IF
        GrapesIO%giFileName = TRIM(outputFile)
        CALL GrapesIO%releaseHGrid
      END IF

      BLOCK
        INTEGER(i_kind) :: vLevel_s, vLevel_u
        TYPE(RectangleGridUtility_t) :: regular
        REAL(r_kind), ALLOCATABLE :: valueDA_mp_s(:, :)
        REAL(r_kind), ALLOCATABLE :: valueDA_sp_s(:, :)

        REAL(r_kind), ALLOCATABLE :: valueDA_mp_u(:, :)
        REAL(r_kind), ALLOCATABLE :: valueDA_sp_u(:, :)
        REAL(r_kind), ALLOCATABLE :: valueDA_Ca1(:, :)
        REAL(r_kind), ALLOCATABLE :: valueDA_Ca2(:, :)

        INTEGER(i_kind), ALLOCATABLE :: idxtgt(:, :)
        REAL(r_kind), ALLOCATABLE :: coetgt(:, :)
        ! REAL(r_kind) :: swap3d(GrapesIO%ids:GrapesIO%ide, &
        !                        GrapesIO%kds:GrapesIO%kde, &
        !                        GrapesIO%jds:GrapesIO%jde)

        vLevel_s = GrapesIO%kde_s - GrapesIO%kds_s + 1
        vLevel_u = GrapesIO%kde_u - GrapesIO%kds_u + 1

        CALL sg_mp%mpddInfo_sg%bCast(vLevel_s)
        CALL sg_mp%mpddInfo_sg%bCast(vLevel_u)

        ! Allocate memory and initialize the interpolation stencils and coefficients
        ALLOCATE (valueDA_mp_s(vLevel_s, sg%num_cell))
        ALLOCATE (valueDA_mp_u(vLevel_u, sg%num_cell))

        BLOCK
          REAL(r_kind), ALLOCATABLE :: cell_cntr(:, :)

          IF (sg_mp%isBaseProc()) ALLOCATE (cell_cntr(2, sg%num_icell_global))
          CALL sg_mp%aggrGridReal2D(sg_mp%cell_cntr, cell_cntr, [2, sg%num_icell_global])

          IF (sg%isBaseProc()) THEN
            PRINT *, 'vLevel_u', vLevel_u, 'sg%num_icell_global', sg%num_icell_global, HorNumGridModel, vLevel_s
            ALLOCATE (valueDA_sp_u(vLevel_u, sg%num_icell_global))
            ALLOCATE (valueDA_sp_s(vLevel_s, sg%num_icell_global))
            ALLOCATE (idxtgt(4, HorNumGridModel), coetgt(4, HorNumGridModel))

            CALL regular%RectangleHorizontalIntp(TRANSPOSE(cell_cntr), sg%num_icell_global, &
                                                 TRANSPOSE(this%cellCntrModel), this%HorNumGridModel, &
                                                 idxtgt, coetgt, sg%mpddInfo_sg%myrank)

            PRINT *, 'Finish the initialization', SHAPE(cell_cntr), SHAPE(this%cellCntrModel)
          END IF

          IF (sg_mp%isBaseProc()) DEALLOCATE (cell_cntr)
        END BLOCK

        ! If uwnd is available
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('uwnd'))%DATA(:, i, sg%tSlots) - Xb%fields(Xb%getVarIdx('uwnd'))%DATA(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_u(:, i), valueDA_mp_u(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_u, valueDA_sp_u, [vLevel_u, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%u(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))

            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%u(i, GrapesIO%kds_u:GrapesIO%kde_u, j) = &
                valueDA_sp_u(:, idxtgt(1, (j - 1) * GrapesIO%idn + i)) * coetgt(1, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(2, (j - 1) * GrapesIO%idn + i)) * coetgt(2, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(3, (j - 1) * GrapesIO%idn + i)) * coetgt(3, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(4, (j - 1) * GrapesIO%idn + i)) * coetgt(4, (j - 1) * GrapesIO%idn + i)
            END FORALL
          END IF

          PRINT *, 'uwnd has been interpolated to the grapes grid.'
        END IF

        ! If vwnd is available
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('vwnd'))%DATA(:, i, sg%tSlots) - Xb%fields(Xb%getVarIdx('vwnd'))%DATA(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_u(:, i), valueDA_mp_u(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_u, valueDA_sp_u, [vLevel_u, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%v(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%v(i, GrapesIO%kds_u:GrapesIO%kde_u, j) = &
                valueDA_sp_u(:, idxtgt(1, (j - 1) * GrapesIO%idn + i)) * coetgt(1, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(2, (j - 1) * GrapesIO%idn + i)) * coetgt(2, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(3, (j - 1) * GrapesIO%idn + i)) * coetgt(3, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_u(:, idxtgt(4, (j - 1) * GrapesIO%idn + i)) * coetgt(4, (j - 1) * GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'vwnd has been interpolated to the grapes grid.'

        END IF

        ! ! If pi
        ! IF (Xm%getVarIdx('rho') .ne. 0 .AND. Xm%getVarIdx('rhov') .ne. 0) THEN

        !   Xm%Fields(Xm%getVarIdx('rhov'))%data = Xm%Fields(Xm%getVarIdx('rhov'))%data/Xm%Fields(Xm%getVarIdx('rho'))%data

        !   ! Aggregate the matrix to the base process
        !   DO i = 1, sg%num_cell
        !     CALL interp1d(sg%zHght(:, i), Xm%fields(Xm%getVarIdx('rhov'))%data(:, i, sg%tSlots), &
        !                   this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
        !   END DO
        !   CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

        !   ! Interpolate the data from the MOTOR-DA grid to the model grid
        !   IF (sg%isBaseProc()) THEN
        !     ALLOCATE (GrapesIO%hgrid%q(GrapesIO%ids:GrapesIO%ide, &
        !                                GrapesIO%kds:GrapesIO%kde, &
        !                                GrapesIO%jds:GrapesIO%jde))
        !     FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
        !       GrapesIO%hgrid%q(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
        !         valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
        !         valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
        !         valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
        !         valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
        !     END FORALL

        !   END IF
        !   PRINT *, 'rhov has been interpolated to the grapes grid.'

        ! END IF

        ! If vwnd is available
        IF ((Xm%getVarIdx('rho') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) .OR. &
            (Xm%getVarIdx('pres') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) .OR. &
            (Xm%getVarIdx('lnp') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0)) THEN
          BLOCK
            REAL(r_kind), ALLOCATABLE :: swapDA(:, :, :)
            ALLOCATE (swapDA(sg_mp%vLevel, sg_mp%num_cell, sg_mp%tSlots))

            ! Convert the lnp to pres for grapesinput
            IF (Xm%getVarIdx('lnp') .NE. 0) THEN
              CALL Xm%Fields(Xm%getVarIdx('lnp'))%Set_Name('pres')
              Xm%Fields(Xm%getVarIdx('pres'))%DATA = EXP(Xm%Fields(Xm%getVarIdx('pres'))%DATA)
            ELSE IF (Xm%getVarIdx('psl') .NE. 0) THEN
              CALL Xm%Fields(Xm%getVarIdx('psl'))%Set_Name('pres')

              Xm%Fields(Xm%getVarIdx('pres'))%DATA(1, :, sg%tSlots) = psl_to_psfc( &
                                                                      Xm%Fields(Xm%getVarIdx('pres'))%DATA(1, :, sg%tSlots), &
                                                                      sg%zHght(1, :), Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, :, sg%tSlots), sg%tmsl)
            END IF

            IF (Xm%getVarIdx('rho') .NE. 0) THEN
              swapDA(1, :, :) = Xm%Fields(Xm%getVarIdx('rho'))%DATA(1, :, :) * &
                                Xm%Fields(Xm%getVarIdx('temp'))%DATA(1, :, :) * dry_air_gas_const
            ELSE IF (Xm%getVarIdx('pres') .NE. 0) THEN
              swapDA(1, :, :) = Xm%Fields(Xm%getVarIdx('pres'))%DATA(1, :, :)
            END IF

            ! Calculate the pressure at the model level with surface pressure and hydrosatic balance
            DO j = 2, sg%vLevel
              swapDA(j, :, sg_mp%tSlots) = swapDA(j - 1, :, sg_mp%tSlots) * EXP(-(sg%zHght(j, :) - sg%zHght(j - 1, :)) * g &
                                                                                / (dry_air_gas_const * &
                                                                                   (Xm%Fields(Xm%getVarIdx('temp'))%DATA(j, :, sg_mp%tSlots) * &
                                                                                    (Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(j, :, sg_mp%tSlots) * 0.608 + 1) + &
                                                                                    Xm%Fields(Xm%getVarIdx('temp'))%DATA(j - 1, :, sg_mp%tSlots) * &
                                                                                    (Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(j - 1, :, sg_mp%tSlots) * 0.608 + 1)) / 2.0D0))
            END DO
            Xm%Fields(Xm%getVarIdx('pres'))%DATA = swapDA

            ! Calculate the Saturated vapor at each level
            ASSOCIATE (temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA, pres => Xm%Fields(Xm%getVarIdx('pres'))%DATA, &
                       qvapor => Xm%fields(Xm%getVarIdx('qvapor'))%DATA)
              DO j = 1, sg%vLevel
                DO i = 1, sg%num_cell
                  IF (qvapor(j, i, sg%tSlots) > Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots))) THEN
! #ifdef DEBUG
!                     ! Yuanfu Xie modified this for more information:
!                     WRITE (*, 1) j, i, qvapor(j, i, sg%tSlots), Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots)), &
!                       sg%mpddInfo_sg%myrank
! 1                   FORMAT('Write_X - Saturated Vapor at vlvl/cell: ', I2, I8, ' q/svapor: ', 2D12.4, ' pc ', I4)
! #endif
                    qvapor(j, i, sg%tSlots) = Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots))
                  END IF
                END DO
              END DO
            END ASSOCIATE

            ! Recalculate the pressure at the model level with surface pressure and hydrosatic balance
            DO j = 2, sg%vLevel
              swapDA(j, :, sg_mp%tSlots) = swapDA(j - 1, :, sg_mp%tSlots) * EXP(-(sg%zHght(j, :) - sg%zHght(j - 1, :)) * g &
                                                                                / (dry_air_gas_const * &
                                                                                   (Xm%Fields(Xm%getVarIdx('temp'))%DATA(j, :, sg_mp%tSlots) * &
                                                                                    (Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(j, :, sg_mp%tSlots) * 0.608 + 1) + &
                                                                                    Xm%Fields(Xm%getVarIdx('temp'))%DATA(j - 1, :, sg_mp%tSlots) * &
                                                                                    (Xm%Fields(Xm%getVarIdx('qvapor'))%DATA(j - 1, :, sg_mp%tSlots) * 0.608 + 1)) / 2.0D0))
            END DO
            Xm%Fields(Xm%getVarIdx('pres'))%DATA = swapDA

            swapDA = (swapDA / surface_ref_pres / 100.0D0)**(dry_air_gas_const / spec_heat_const_pres)
            ! this%pip = swapDA - this%piTotal

            ! Aggregate the matrix to the base process
            ALLOCATE (valueDA_Ca1(vLevel_u, sg%num_cell), valueDA_Ca2(vLevel_u, sg%num_cell))
            DO i = 1, sg%num_cell
              ! CALL interp1d(sg%zHght(:, i), swapDA(:, i, sg%tSlots) - this%piTotal(:, i, sg%tSlots), &
              !               this%zRHghtModelAtDA_u(:, i), valueDA_mp_u(:, i))
              CALL interp1d(sg%zHght(:, i), swapDA(:, i, sg%tSlots), &
                            this%zRHghtModelAtDA_u(:, i), valueDA_Ca1(:, i), 'LOG')
              CALL interp1d(sg%zHght(:, i), this%piTotal(:, i, sg%tSlots), &
                            this%zRHghtModelAtDA_u(:, i), valueDA_Ca2(:, i), 'LOG')
              valueDA_mp_u(:, i) = valueDA_Ca1(:, i) - valueDA_Ca2(:, i)
            END DO
            DEALLOCATE (valueDA_Ca1, valueDA_Ca2)

            CALL sg%aggrGridReal2D(valueDA_mp_u, valueDA_sp_u, [vLevel_u, sg%num_icell_global])

            ! Interpolate the data from the MOTOR-DA grid to the model grid
            IF (sg%isBaseProc()) THEN

              ALLOCATE (GrapesIO%hgrid%pi(GrapesIO%ids:GrapesIO%ide, &
                                          GrapesIO%kds:GrapesIO%kde, &
                                          GrapesIO%jds:GrapesIO%jde))
              FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
                GrapesIO%hgrid%pi(i, GrapesIO%kds_u:GrapesIO%kde_u, j) = &
                  valueDA_sp_u(:, idxtgt(1, (j - 1) * GrapesIO%idn + i)) * coetgt(1, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_u(:, idxtgt(2, (j - 1) * GrapesIO%idn + i)) * coetgt(2, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_u(:, idxtgt(3, (j - 1) * GrapesIO%idn + i)) * coetgt(3, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_u(:, idxtgt(4, (j - 1) * GrapesIO%idn + i)) * coetgt(4, (j - 1) * GrapesIO%idn + i)

              END FORALL
            END IF
            PRINT *, 'pi has been interpolated to the grapes grid to pi.'
            !----------------------------------------------------------------

            !------------------------------------------------------------------------------------ theta
            swapDA = Xm%Fields(Xm%getVarIdx('temp'))%DATA / swapDA

            ! this%thp = swapDA - this%thTotal

            ! Aggregate the matrix to the base process
            DO i = 1, sg%num_cell
              CALL interp1d(sg%zHght(:, i), swapDA(:, i, sg%tSlots) - this%thTotal(:, i, sg%tSlots), &
                            this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
            END DO
            CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

            ! Interpolate the data from the MOTOR-DA grid to the model grid
            IF (sg%isBaseProc()) THEN

              ALLOCATE (GrapesIO%hgrid%th(GrapesIO%ids:GrapesIO%ide, &
                                          GrapesIO%kds:GrapesIO%kde, &
                                          GrapesIO%jds:GrapesIO%jde))
              FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
                GrapesIO%hgrid%th(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                  valueDA_sp_s(:, idxtgt(1, (j - 1) * GrapesIO%idn + i)) * coetgt(1, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_s(:, idxtgt(2, (j - 1) * GrapesIO%idn + i)) * coetgt(2, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_s(:, idxtgt(3, (j - 1) * GrapesIO%idn + i)) * coetgt(3, (j - 1) * GrapesIO%idn + i) + &
                  valueDA_sp_s(:, idxtgt(4, (j - 1) * GrapesIO%idn + i)) * coetgt(4, (j - 1) * GrapesIO%idn + i)
              END FORALL
            END IF
            PRINT *, 'th has been interpolated to the grapes grid to th.'

            ! ! ------------------------------------------------ pip
            ! ! Aggregate the matrix to the base process
            ! DO i = 1, sg%num_cell
            !   CALL interp1d(sg%zHght(:, i), this%pip(:, i, sg%tSlots), this%zRHghtModelAtDA_u(:, i), valueDA_mp_u(:, i))
            ! END DO
            ! CALL sg%aggrGridReal2D(valueDA_mp_u, valueDA_sp_u, [vLevel_u, sg%num_icell_global])

            ! ! Interpolate the data from the MOTOR-DA grid to the model grid
            ! IF (sg%isBaseProc()) THEN

            !   ALLOCATE (GrapesIO%hgrid%pip(GrapesIO%ids:GrapesIO%ide, &
            !                                GrapesIO%kds:GrapesIO%kde, &
            !                                GrapesIO%jds:GrapesIO%jde))

            !   FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
            !     GrapesIO%hgrid%pip(i, GrapesIO%kds_u:GrapesIO%kde_u, j) = &
            !       valueDA_sp_u(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_u(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_u(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_u(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)

            !   END FORALL
            ! END IF
            ! PRINT *, 'pip has been interpolated to the grapes grid to pi.'
            ! !----------------------------------------------------------------

            ! !------------------------------------------------------------------------------------ thetap

            ! ! Aggregate the matrix to the base process
            ! DO i = 1, sg%num_cell
            !   CALL interp1d(sg%zHght(:, i), this%thp(:, i, sg%tSlots), &
            !                 this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
            ! END DO
            ! CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

            ! ! Interpolate the data from the MOTOR-DA grid to the model grid
            ! IF (sg%isBaseProc()) THEN
            !   ALLOCATE (GrapesIO%hgrid%thp(GrapesIO%ids:GrapesIO%ide, &
            !                                GrapesIO%kds:GrapesIO%kde, &
            !                                GrapesIO%jds:GrapesIO%jde))
            !   FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
            !     GrapesIO%hgrid%thp(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
            !       valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
            !       valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            !   END FORALL
            ! END IF
            ! PRINT *, 'thp has been interpolated to the grapes grid to th.'

            DEALLOCATE (swapDA)
          END BLOCK

          ASSOCIATE (temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA, pres => Xm%Fields(Xm%getVarIdx('pres'))%DATA, &
                     qvapor => Xm%fields(Xm%getVarIdx('qvapor'))%DATA)
          DO j = 1, sg%vLevel
            DO i = 1, sg%num_cell
              IF (qvapor(j, i, sg%tSlots) > Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots))) THEN
#ifdef DEBUG
                ! Yuanfu Xie modified this for more information:
                WRITE (*, 2) j, i, qvapor(j, i, sg%tSlots), Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots)), &
                  sg%mpddInfo_sg%myrank
2               FORMAT('Write_X - 2nd Saturated Vapor at vlvl/cell: ', I2, I8, ' q/svapor: ', 2D12.4, ' pc ', I4)
#endif
                qvapor(j, i, sg%tSlots) = Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots))
              END IF
            END DO
          END DO

          END ASSOCIATE

          IF (Xm%getVarIdx('ref') .NE. 0) THEN
            ASSOCIATE (temp => Xm%Fields(Xm%getVarIdx('temp'))%DATA, pres => Xm%Fields(Xm%getVarIdx('pres'))%DATA, &
                       qvapor => Xm%fields(Xm%getVarIdx('qvapor'))%DATA, ref => Xm%fields(Xm%getVarIdx('ref'))%DATA)
              DO j = 1, sg%vLevel
                DO i = 1, sg%num_cell
                  IF (ref(j, i, sg%tSlots) > 5) qvapor(j, i, sg%tSlots) = Saturated_Vapor(temp(j, i, sg%tSlots), pres(j, i, sg%tSlots))
                END DO
              END DO
            END ASSOCIATE
          END IF

        END IF

        ! If qvapor is available
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qvapor'))%DATA(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qvapor'))%DATA(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%q(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%q(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1) * GrapesIO%idn + i)) * coetgt(1, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1) * GrapesIO%idn + i)) * coetgt(2, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1) * GrapesIO%idn + i)) * coetgt(3, (j - 1) * GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1) * GrapesIO%idn + i)) * coetgt(4, (j - 1) * GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qvapor has been interpolated to the grapes grid.'

        END IF

        ! If qcloud is available
        IF (Xm%getVarIdx('qcloud') .ne. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qcloud'))%data(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qcloud'))%data(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%qc(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%qc(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qcloud has been interpolated to the grapes grid.'

        END IF

        ! If qrain is available
        IF (Xm%getVarIdx('qrain') .ne. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qrain'))%data(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qrain'))%data(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%qr(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%qr(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qrain has been interpolated to the grapes grid.'

        END IF

        ! If qice is available
        IF (Xm%getVarIdx('qice') .ne. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qice'))%data(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qice'))%data(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%qi(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%qi(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qice has been interpolated to the grapes grid.'

        END IF

        ! If qsnow is available
        IF (Xm%getVarIdx('qsnow') .ne. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qsnow'))%data(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qsnow'))%data(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%qs(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%qs(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qsnow has been interpolated to the grapes grid.'

        END IF

        ! If qgraupel is available
        IF (Xm%getVarIdx('qgraupel') .ne. 0) THEN
          ! Aggregate the matrix to the base process
          DO i = 1, sg%num_cell
            CALL interp1d(sg%zHght(:, i), &
                          Xm%fields(Xm%getVarIdx('qgraupel'))%data(:, i, sg%tSlots) - &
                          Xb%fields(Xb%getVarIdx('qgraupel'))%data(:, i, sg%tSlots), &
                          this%zRHghtModelAtDA_s(:, i), valueDA_mp_s(:, i))
          END DO
          CALL sg%aggrGridReal2D(valueDA_mp_s, valueDA_sp_s, [vLevel_s, sg%num_icell_global])

          ! Interpolate the data from the MOTOR-DA grid to the model grid
          IF (sg%isBaseProc()) THEN
            ALLOCATE (GrapesIO%hgrid%qg(GrapesIO%ids:GrapesIO%ide, &
                                       GrapesIO%kds:GrapesIO%kde, &
                                       GrapesIO%jds:GrapesIO%jde))
            FORALL (i=1:GrapesIO%idn, j=1:GrapesIO%jdn)
              GrapesIO%hgrid%qg(i, GrapesIO%kds_s:GrapesIO%kde_s, j) = &
                valueDA_sp_s(:, idxtgt(1, (j - 1)*GrapesIO%idn + i))*coetgt(1, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(2, (j - 1)*GrapesIO%idn + i))*coetgt(2, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(3, (j - 1)*GrapesIO%idn + i))*coetgt(3, (j - 1)*GrapesIO%idn + i) + &
                valueDA_sp_s(:, idxtgt(4, (j - 1)*GrapesIO%idn + i))*coetgt(4, (j - 1)*GrapesIO%idn + i)
            END FORALL
          END IF
          PRINT *, 'qgraupel has been interpolated to the grapes grid.'

        END IF

        ! Flush the data to the file
        IF (sg%isBaseProc() .AND. this%CloudyMode) CALL GrapesIO%writeBackToQcqrInput(TRIM(inputQcqr), TRIM(outputQcqr))
        IF (sg%isBaseProc()) CALL GrapesIO%writeBackToGrapesInput
        IF (sg%isBaseProc()) CALL GrapesIO%releaseHGrid

        ! Release the memory
        IF (ALLOCATED(valueDA_mp_u)) DEALLOCATE (valueDA_mp_u)
        IF (ALLOCATED(valueDA_sp_u)) DEALLOCATE (valueDA_sp_u)

        IF (ALLOCATED(valueDA_mp_s)) DEALLOCATE (valueDA_mp_s)
        IF (ALLOCATED(valueDA_sp_s)) DEALLOCATE (valueDA_sp_s)

        IF (ALLOCATED(idxtgt)) DEALLOCATE (idxtgt)
        IF (ALLOCATED(coetgt)) DEALLOCATE (coetgt)

      END BLOCK
    END ASSOCIATE
  END SUBROUTINE m_write_Xm_into_bcg

  SUBROUTINE m_read_bcg_into_Xm(this, Xm, sg)
    IMPLICIT NONE

    CLASS(IOGrapes_t) :: this
    TYPE(State_t), INTENT(INOUT) :: Xm
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(ModelCoupler_t) :: ModelCoupler

    REAL(r_kind), ALLOCATABLE :: timeSlotsModel(:)
    INTEGER(i_kind) :: vLevel_s, vLevel_u, istatus

    INTEGER(i_kind) :: i, j, k, l, nprocInterp, idn, jdn, kdn, varIdx, it ! Yuanfu Xie added l for horizonSimilarity initialization
    REAL(r_kind), ALLOCATABLE :: cc(:, :)
    CHARACTER(len=256) :: grapes_model_type

    CALL this%m_read_config_file
    ALLOCATE (this%GrapesIO(this%nModelFiles))

    ASSOCIATE (sg_mp => sg, GrapesIO => this%GrapesIO, &
               HorNumGridModel => this%HorNumGridModel)

      ! Check the horizonSimilarity parameters:
      WRITE (*, 1) sg%gLevel, sg%mpddInfo_sg%myrank, sg%scales4Similarity
1     FORMAT('Similarity for landmask and topo at G', I3, ' pc: ', I3, ' scales land/topo: ', 2E14.6)

      IF (sg_mp%isBaseProc()) THEN
        ALLOCATE (timeSlotsModel(this%nModelFiles))

        istatus = yaml_get_var(this%m_configfile, 'IO', 'grapes_model_type', grapes_model_type)
        IF (istatus .NE. 0) grapes_model_type = 'CMA-GD'  ! If users do not specify fill option, it does not fill

        PRINT *, 'grapes_model_type is: ', TRIM(grapes_model_type), istatus

        ! Read the grapes fields.
        DO i = 1, size(this%m_nmlFilename)
          IF (this%CloudyMode) THEN
            GrapesIO(i) = GrapesIO_t(this%m_nmlFilename(i), this%m_modelFileName(i), grapes_model_type=grapes_model_type, &
            qcqrFileName=this%m_qcqrFileName(i))
          ELSE
            GrapesIO(i) = GrapesIO_t(this%m_nmlFilename(i), this%m_modelFileName(i), grapes_model_type=grapes_model_type)
          END IF
          timeSlotsModel(i) = Get_Time_GMT_to_Unix((/GrapesIO(i)%hgrid%config%start_year, &
                                                     GrapesIO(i)%hgrid%config%start_month, &
                                                     GrapesIO(i)%hgrid%config%start_day, &
                                                     GrapesIO(i)%hgrid%config%start_hour, &
                                                     GrapesIO(i)%hgrid%config%start_minute, &
                                                     GrapesIO(i)%hgrid%config%start_second/)); 
          PRINT *, 'Posix time is: ', timeSlotsModel(i)
        END DO

        HorNumGridModel = GrapesIO(1)%idn * GrapesIO(1)%jdn
        ALLOCATE (this%cellCntrModel(2, HorNumGridModel))

        this%cellCntrModel(1, :) = GrapesIO(1)%lat * degree2radian
        this%cellCntrModel(2, :) = GrapesIO(1)%lon * degree2radian
        PRINT *, 'Finish loading of Grapes input files.'

        idn = GrapesIO(1)%idn
        jdn = GrapesIO(1)%jdn
        kdn = GrapesIO(1)%kdn
      END IF

      CALL sg_mp%mpddInfo_sg%bCast(GrapesIO(1)%kdn)

      ! Initialize the model coupler.
      CALL ModelCoupler%initialize(this%m_configFile, this%geometry)
      print*,'Jilong test1'
      CALL ModelCoupler%gen_interp_coeffs(this%cellCntrModel, HorNumGridModel, sg_mp)
      PRINT *, 'Finish initialize the ModelCoupler.'
      PRINT *, 'Start loading the topo and grid height in MOTOR.'

      CALL sg_mp%mpddInfo_sg%barrier
      ! Ingest topo and update the zHght
      IF (sg_mp%isBaseProc()) THEN

        BLOCK
          REAL(r_kind), ALLOCATABLE :: sigIn(:)
          REAL(r_kind), ALLOCATABLE :: sigOut(:)
          REAL(r_kind), ALLOCATABLE :: hIn(:)

          ALLOCATE (sigIn(GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1), &
                    sigOut(sg%vLevel), &
                    hIn(GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1))

          sigIn = (/(i, i=0, GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1 - 1)/) * 1.0D0 / &
                  (GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1 - 1)
          sigOut = (/(i, i=0, sg%vLevel - 1)/) * 1.0D0 / (sg%vLevel - 1)
          PRINT *, 'sigIn: ', sigIn
          PRINT *, 'sigOut: ', sigOut
          PRINT *, 'GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1: ', GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1
          hIn = GrapesIO(1)%zSigma_u(GrapesIO(1)%kds_u:GrapesIO(1)%kde_u)
          hIn(1) = 0
          hIn(GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1) = GrapesIO(1)%zSigma_s(GrapesIO(1)%kde_s)
          CALL interp1d_1D(sigIn, hIn, sigOut, sg_mp%sigma)
          PRINT *, 'zh: ', GrapesIO(1)%zSigma_u(GrapesIO(1)%kds_u:GrapesIO(1)%kde_u)
          !PRINT *, 'sigma: ', sg_mp%sigma

          DEALLOCATE (sigIn, sigOut, hIn)
        END BLOCK

      END IF

      CALL sg_mp%mpddInfo_sg%bcast(sg_mp%sigma)
      PRINT *, 'sg_mp%sigma: ', sg_mp%sigma

      ! Guided filter
      ! IF (sg_mp%isBaseProc()) then
      !   ALLOCATE (cc(size(GrapesIO(1)%hgrid%ht, 1), size(GrapesIO(1)%hgrid%ht, 2)))

      !   BLOCK
      !     cc = GrapesIO(1)%hgrid%ht*1.0D0
      !     CALL guidedfilter(cc, cc, 3, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%ht = cc
      !   END BLOCK

      !   DO i = GrapesIO(1)%kds_u, GrapesIO(1)%kde_u
      !     cc = GrapesIO(1)%hgrid%u(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%u(:, i, :) = cc

      !     cc = GrapesIO(1)%hgrid%v(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%v(:, i, :) = cc

      !     cc = DLOG(GrapesIO(1)%hgrid%pi(:, i, :)*1.0D0)
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%pi(:, i, :) = EXP(cc)
      !   END DO

      !   DO i = GrapesIO(1)%kds_s, GrapesIO(1)%kde_s
      !     cc = GrapesIO(1)%hgrid%th(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%th(:, i, :) = cc

      !     cc = GrapesIO(1)%hgrid%q(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     GrapesIO(1)%hgrid%q(:, i, :) = cc
      !   END DO

      !   DEALLOCATE (cc)

      ! END IF

      ! Ingest topo
      CALL ModelCoupler%ingest_to_topo(GrapesIO(1)%hgrid%ht(:, :), HorNumGridModel, sg_mp)

      ! Smooth the topography:
      CALL smoothField(sg_mp, 1, 1, ModelCoupler%numSmoothTopo(sg_mp%glevel), sg_mp%topo, sg_mp%topo)

      ! Update Hght
      CALL sg_mp%update_zHght_from_topo_and_sigma

      ! Ingest surface parameter
      DO it = 1, sg_mp%tSlots
        ! CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%tsk(:, :), sg_mp%tskin(:,it), sg_mp)
        ! CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%psfc(:, :), sg_mp%psfc(:,it), sg_mp)
        ! CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%u(:, GrapesIO%kds_u, :), sg_mp%u10m(:,it), sg_mp)
        ! CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%v(:, GrapesIO%kds_u, :), sg_mp%v10m(:,it), sg_mp)
        CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%snowc(:, :), sg_mp%snowc(:, it), sg_mp)
        CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%xice(:, :), sg_mp%xice(:, it), sg_mp)
      END DO
      CALL ModelCoupler%ingest_to_field_data_1D_nn(HorNumGridModel, GrapesIO(1)%hgrid%xland(:, :), sg_mp%landmask, sg_mp)
      CALL ModelCoupler%ingest_to_field_data_1D_nn(HorNumGridModel, GrapesIO(1)%hgrid%soil_type(:, :), &
                                                   sg_mp%soil_type, sg_mp)
      CALL ModelCoupler%ingest_to_field_data_1D(HorNumGridModel, GrapesIO(1)%hgrid%veg_fraction(:, :), sg_mp%veg_frac, sg_mp)

      FORALL (i=1:sg_mp%num_cell, sg_mp%landmask(i) == 2.0D0) sg_mp%landmask(i) = 0.0D0

      ! Yuanfu Xie 2023/10/18 added assignment to initialize horizonSimilarity:
      DO i = 1, sg_mp%num_icell
        IF (sg_mp%cell_type(i) .NE. 0) CYCLE ! Apply to interior points only
        DO j = 1, sg_mp%num_edge
          DO k = 1, UBOUND(sg_mp%edge_stcl, 1)
            sg_mp%horizonSimilarity(i) = sg_mp%horizonSimilarity(i) * &
                                         EXP(-sg_mp%scales4Similarity(1) * (sg_mp%landmask(i) - sg_mp%landmask(sg_mp%edge_stcl(k, j, i)))**2) * &
                                         EXP(-sg_mp%scales4Similarity(2) * (sg_mp%topo(i) - sg_mp%topo(sg_mp%edge_stcl(k, j, i)))**2)
          END DO
        END DO
      END DO

      BLOCK
        ! Generate the model height at location of MOTOR-DA cell
        REAL(r_single) :: zRHght_s_temp(GrapesIO(1)%kde_s - GrapesIO(1)%kds_s + 1, HorNumGridModel)
        ! REAL(r_kind) :: sigma_u_bottom, sigma_u_top
        REAL(r_single) :: zRHght_u_temp(GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1, HorNumGridModel)

        ! PRINT *, 'Here generate the model height at location of MOTOR-DA cell'
        IF (sg_mp%isBaseProc()) THEN
          vLevel_s = GrapesIO(1)%kde_s - GrapesIO(1)%kds_s + 1
          vLevel_u = GrapesIO(1)%kde_u - GrapesIO(1)%kds_u + 1

          FORALL (i=1:idn, j=1:jdn)
            zRHght_s_temp(:, (j - 1) * GrapesIO(1)%idn + i) = GrapesIO(1)%zRHght_s(i, GrapesIO(1)%kds_s:GrapesIO(1)%kde_s, j)
            zRHght_u_temp(:, (j - 1) * GrapesIO(1)%idn + i) = GrapesIO(1)%zRHght_u(i, GrapesIO(1)%kds_u:GrapesIO(1)%kde_u, j)
          END FORALL
        END IF

        CALL sg_mp%mpddInfo_sg%bCast(vLevel_s)
        CALL sg_mp%mpddInfo_sg%bCast(vLevel_u)

        ALLOCATE (this%zRHghtModelAtDA_s(vLevel_s, sg_mp%num_cell)); 
        ALLOCATE (this%zRHghtModelAtDA_u(vLevel_u, sg_mp%num_cell)); 
        ! CALL sg_mp%mpddInfo_sg%barrier
        PRINT *, 'Start ingest_to_field_data_2D ', sg_mp%mpddInfo_sg%myrank
        ! CALL sg_mp%mpddInfo_sg%barrier

        CALL ModelCoupler%ingest_to_field_data_2D(vLevel_s, zRHght_s_temp, this%zRHghtModelAtDA_s, sg_mp)
        CALL ModelCoupler%ingest_to_field_data_2D(vLevel_u, zRHght_u_temp, this%zRHghtModelAtDA_u, sg_mp)
      END BLOCK

      PRINT *, 'Start ingesting the model variables with landmask: ', MAXVAL(sg_mp%landmask), MINVAL(sg_mp%landmask), sg_mp%mpddInfo_sg%myrank
      BLOCK
        REAL(r_single), ALLOCATABLE :: tmpModel_u(:, :, :)
        REAL(r_single), ALLOCATABLE :: tmp_u(:, :, :)

        REAL(r_single), ALLOCATABLE :: tmpModel_s(:, :, :)
        REAL(r_single), ALLOCATABLE:: tmp_s(:, :, :)

        REAL(r_kind) :: TMSL
        INTEGER(i_kind) :: vIdx_psl, vIdx_temp
        INTEGER :: it

        IF (sg_mp%isBaseProc()) THEN
          ALLOCATE (tmpModel_u(vLevel_u, HorNumGridModel, this%nModelFiles), &
                    tmp_u(vLevel_u, HorNumGridModel, sg_mp%tSlots), &
                    tmpModel_s(vLevel_s, HorNumGridModel, this%nModelFiles), &
                    tmp_s(vLevel_s, HorNumGridModel, sg_mp%tSlots))
        ELSE
          ALLOCATE (tmp_u(1, 1, sg_mp%tSlots), &
                    tmp_s(1, 1, sg_mp%tSlots))
        END IF

        ! CALL sg_mp%mpddInfo_sg%barrier
        ! Ingest uwnd
        IF (Xm%getVarIdx('uwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_u(:, (j - 1) * GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%u(i, GrapesIO(k)%kds_u:GrapesIO(k)%kde_u, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_u, sg_mp%tt, tmp_u)
            PRINT *, "End of interp1d_3D_idx3.",sg_mp%tSlots,ubound(tmp_u,3),sg_mp%mpddInfo_sg%myrank
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_u, this%zRHghtModelAtDA_u, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('uwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - uwnd has been loaded.'
        END IF

        ! CALL sg_mp%mpddInfo_sg%barrier
        ! Ingest vwnd
        IF (Xm%getVarIdx('vwnd') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_u(:, (j - 1) * GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%v(i, GrapesIO(k)%kds_u:GrapesIO(k)%kde_u, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_u, sg_mp%tt, tmp_u)
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_u, this%zRHghtModelAtDA_u, tmp_u, &
                                                             Xm%Fields(Xm%getVarIdx('vwnd'))%DATA, sg_mp)
          PRINT *, 'IOGrapes - vwnd has been loaded.'
        END IF

        ! CALL sg_mp%mpddInfo_sg%barrier
        ! Ingest qvapor
        IF (Xm%getVarIdx('qvapor') .NE. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1) * GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%q(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qvapor'))%DATA, sg_mp)

          varIdx = Xm%getVarIdx('qvapor')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%DATA(i, j, k) <= 1.0D-7) &
            Xm%Fields(varIdx)%DATA(i, j, k) = 1.0D-7
          PRINT *, 'IOGrapes - qvapor has been loaded.'
        END IF

        ! Ingest qcloud
        IF (Xm%getVarIdx('qcloud') .ne. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1)*GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%qc(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qcloud'))%data, sg_mp)

          varIdx = Xm%getVarIdx('qcloud')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%data(i, j, k) <= 1.0D-12) &
            Xm%Fields(varIdx)%data(i, j, k) = 0.0D0
          PRINT *, 'IOGrapes - qcloud has been loaded.'
        END IF

        ! Ingest qrain
        IF (Xm%getVarIdx('qrain') .ne. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1)*GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%qr(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qrain'))%data, sg_mp)

          varIdx = Xm%getVarIdx('qrain')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%data(i, j, k) <= 1.0D-12) &
            Xm%Fields(varIdx)%data(i, j, k) = 0.0D0
          PRINT *, 'IOGrapes - qrain has been loaded.'
        END IF

        ! Ingest qice
        IF (Xm%getVarIdx('qice') .ne. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1)*GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%qi(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qice'))%data, sg_mp)

          varIdx = Xm%getVarIdx('qice')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%data(i, j, k) <= 1.0D-12) &
            Xm%Fields(varIdx)%data(i, j, k) = 0.0D0
          PRINT *, 'IOGrapes - qice has been loaded.'
        END IF

        ! Ingest qsnow
        IF (Xm%getVarIdx('qsnow') .ne. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1)*GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%qs(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qsnow'))%data, sg_mp)

          varIdx = Xm%getVarIdx('qsnow')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%data(i, j, k) <= 1.0D-12) &
            Xm%Fields(varIdx)%data(i, j, k) = 0.0D0
          PRINT *, 'IOGrapes - qsnow has been loaded.'
        END IF

        ! Ingest qgraupel
        IF (Xm%getVarIdx('qgraupel') .ne. 0) THEN
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn)
              tmpModel_s(:, (j - 1)*GrapesIO(k)%idn + i, k) = GrapesIO(k)%hgrid%qg(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            END FORALL

            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s) ! , 'LOG')
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             Xm%Fields(Xm%getVarIdx('qgraupel'))%data, sg_mp)

          varIdx = Xm%getVarIdx('qgraupel')
          FORALL (i=1:sg%vLevel, j=1:sg%num_cell, k=1:sg%tSlots, Xm%Fields(varIdx)%data(i, j, k) <= 1.0D-12) &
            Xm%Fields(varIdx)%data(i, j, k) = 0.0D0
          PRINT *, 'IOGrapes - qgraupel has been loaded.'
        END IF

        ! CALL sg_mp%mpddInfo_sg%barrier
        ! Ingest temp
        IF ((Xm%getVarIdx('rho') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) .OR. &
            (Xm%getVarIdx('pres') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) .OR. &
            (Xm%getVarIdx('psl') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0) .OR. &
            (Xm%getVarIdx('lnp') .NE. 0 .AND. Xm%getVarIdx('temp') .NE. 0)) THEN

          ALLOCATE (this%piTotal(sg_mp%vLevel, sg_mp%num_cell, sg_mp%tSlots))
          ALLOCATE (this%thTotal(sg_mp%vLevel, sg_mp%num_cell, sg_mp%tSlots))
          ALLOCATE (this%pip(sg_mp%vLevel, sg_mp%num_cell, sg_mp%tSlots))
          ALLOCATE (this%thp(sg_mp%vLevel, sg_mp%num_cell, sg_mp%tSlots))

          ! Input pi
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn) tmpModel_u(:, (j - 1) * GrapesIO(k)%idn + i, k) = &
              GrapesIO(k)%hgrid%pi(i, GrapesIO(k)%kds_u:GrapesIO(k)%kde_u, j)
            CALL interp1d(timeSlotsModel, tmpModel_u, sg_mp%tt, tmp_u, 'LOG')
          END IF

          PRINT *, 'Start loading piTotal....'
          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_u, this%zRHghtModelAtDA_u, tmp_u, &
                                                             this%piTotal, sg_mp)
          PRINT *, 'IOGrapes - pi has been loaded.', MAXVAL(this%piTotal)

          ! Input theta
          IF (sg_mp%isBaseProc()) THEN
            FORALL (k=1:this%nModelFiles, i=1:idn, j=1:jdn) tmpModel_s(:, (j - 1) * GrapesIO(k)%idn + i, k) = &
              GrapesIO(k)%hgrid%th(i, GrapesIO(k)%kds_s:GrapesIO(k)%kde_s, j)
            CALL interp1d(timeSlotsModel, tmpModel_s, sg_mp%tt, tmp_s)
          END IF

          CALL ModelCoupler%ingest_to_field_data_3D_ZRHght_s(vLevel_s, this%zRHghtModelAtDA_s, tmp_s, &
                                                             this%thTotal, sg_mp)
          PRINT *, 'IOGrapes - theta has been loaded.'

          Xm%Fields(Xm%getVarIdx('temp'))%DATA = this%piTotal * this%thTotal

          ! Convert pres from theta and pi
          IF (Xm%getVarIdx('pres') .NE. 0) THEN
            Xm%Fields(Xm%getVarIdx('pres'))%DATA = surface_ref_pres * 100.0D0 * (this%piTotal) &
                                                   **(spec_heat_const_pres / dry_air_gas_const)
            PRINT *, 'IOGrapes - pres has been loaded.', MAXVAL(Xm%Fields(Xm%getVarIdx('pres'))%DATA)
            sg%pres = Xm%Fields(Xm%getVarIdx('pres'))%DATA(:, :, sg%tSlots)        ! Set a reference pressure

            ! ! Convert rho from theta and pi
            ! ELSE IF (Xm%getVarIdx('rho') .ne. 0) THEN
            !   Xm%Fields(Xm%getVarIdx('rho'))%data = surface_ref_pres*100.0D0/dry_air_gas_const/ &
            !                                         this%thTotal*this%piTotal**(spec_heat_const_pres/dry_air_gas_const - 1)
            !   PRINT *, 'IOGrapes - rho has been loaded.'
          ELSE IF (Xm%getVarIdx('lnp') .NE. 0) THEN
            Xm%Fields(Xm%getVarIdx('lnp'))%DATA = LOG(surface_ref_pres * 100.0D0 * (this%piTotal) &
                                                      **(spec_heat_const_pres / dry_air_gas_const))
            PRINT *, 'IOGrapes - lnp has been loaded.'
            sg%pres = EXP(Xm%Fields(Xm%getVarIdx('lnp'))%DATA(:, :, sg%tSlots))        ! Set a reference pressure

          ELSE IF (Xm%getVarIdx('psl') .NE. 0) THEN
            Xm%Fields(Xm%getVarIdx('psl'))%DATA = surface_ref_pres * 100.0D0 * (this%piTotal) &
                                                  **(spec_heat_const_pres / dry_air_gas_const)
            sg%pres = Xm%Fields(Xm%getVarIdx('psl'))%DATA(:, :, sg%tSlots)

            Xm%Fields(Xm%getVarIdx('psl'))%DATA(2:sg%vLevel, :, :) = 0.0D0

            TMSL = Xm%getMeanSeaLevelTempaure()
            sg%tmsl = TMSL

            vIdx_psl = Xm%getVarIdx('psl')
            vIdx_temp = Xm%getVarIdx('temp')
            ASSOCIATE (psl => Xm%Fields(vIdx_psl)%DATA, temp => Xm%Fields(vIdx_temp)%DATA)
              FORALL (i=1:sg%num_cell, j=1:sg%tSlots)
                psl(1, i, j) = psfc_to_psl(psl(1, i, j), sg%zHght(1, i), temp(1, i, j), TMSL)
              END FORALL

            END ASSOCIATE

            PRINT *, 'IOGrapes - psl has been loaded.', MAXVAL(Xm%Fields(Xm%getVarIdx('psl'))%DATA)
            PRINT *, 'Mean sea level temperature is: ', TMSL

          END IF
        END IF

        ! Set a reference pressure
        IF (Xm%getVarIdx('pres') /= 0) THEN
        ELSE IF (Xm%getVarIdx('lnp') /= 0) THEN
        END IF

        IF (sg_mp%isBaseProc()) THEN
          DEALLOCATE (tmpModel_u, &
                      tmp_u, &
                      tmpModel_s, &
                      tmp_s)
          ! IF (ALLOCATED(tmpModel_u)) DEALLOCATE (tmpModel_u)
          ! IF (ALLOCATED(tmp_u)) DEALLOCATE (tmp_u)
          ! IF (ALLOCATED(tmpModel_s)) DEALLOCATE (tmpModel_s)
          ! IF (ALLOCATED(tmp_s)) DEALLOCATE (tmp_s)
        END IF

        IF (ALLOCATED(timeSlotsModel)) DEALLOCATE (timeSlotsModel)

        DO i = 1, SIZE(this%m_nmlFilename)
          CALL this%GrapesIO(i)%releaseHGrid
        END DO

      END BLOCK
      CALL Calc_sg_vars(this%m_configFile, Xm, sg_mp)
      CALL update_restrictionOfStatics(this%m_configFile, Xm, sg_mp)

      CALL sg_mp%mpddInfo_sg%barrier

    END ASSOCIATE

    ! Deallocate memories
  END SUBROUTINE m_read_bcg_into_Xm

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(IOGrapes_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%GrapesIO)) DEALLOCATE (this%GrapesIO)
    IF (ALLOCATED(this%zRHghtModelAtDA_s)) DEALLOCATE (this%zRHghtModelAtDA_s)
    IF (ALLOCATED(this%zRHghtModelAtDA_u)) DEALLOCATE (this%zRHghtModelAtDA_u)
    IF (ALLOCATED(this%cellCntrModel)) DEALLOCATE (this%cellCntrModel)
    IF (ALLOCATED(this%piTotal)) DEALLOCATE (this%piTotal)
    IF (ALLOCATED(this%thTotal)) DEALLOCATE (this%thTotal)
    IF (ALLOCATED(this%thp)) DEALLOCATE (this%thp)
    IF (ALLOCATED(this%pip)) DEALLOCATE (this%pip)
    IF (ALLOCATED(this%m_nmlFilename)) DEALLOCATE (this%m_nmlFilename)
    IF (ALLOCATED(this%m_modelFileName)) DEALLOCATE (this%m_modelFileName)
    IF (ALLOCATED(this%m_qcqrFileName)) DEALLOCATE (this%m_qcqrFileName)

    PRINT *, 'Finish the destructor of IOGrapes_t.'

  END SUBROUTINE destructor

END MODULE IOGrapes_m
