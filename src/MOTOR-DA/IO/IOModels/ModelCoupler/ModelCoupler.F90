!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/20, @GBA-MWF, Shenzhen
!     For adding yaml read of multigrid topography smoothing info and smoothing capability
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE ModelCoupler_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Geometry_m, ONLY: Geometry_t
!   USE interhp_m, only: interhp_t
  USE InterpHPNew_m, ONLY: InterpHPNew_t !新的插值并行方法测试
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE YAMLRead_m

  TYPE ModelCoupler_t
    CHARACTER(len=1024) :: configFile
    INTEGER(i_kind), ALLOCATABLE :: numSmoothTopo(:)
    TYPE(Geometry_t), POINTER :: geometry
    !  TYPE(interhp_t) :: interhp
    TYPE(InterpHPNew_t) :: interhpN
    ! TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE, PASS(this) :: gen_interp_coeffs
    !PROCEDURE, PASS(this) :: ingest_to_field_data
    PROCEDURE, PASS(this) :: ingest_to_field_data_MT
    PROCEDURE, PASS(this) :: ingest_to_topo_and_update_zHght
    PROCEDURE, PASS(this) :: ingest_to_topo
    PROCEDURE, PASS(this) :: ingest_to_topo_and_update_zHght_d
    PROCEDURE, PASS(this) :: ingest_to_field_data_2D
    PROCEDURE, PASS(this) :: ingest_to_field_data_1D
    PROCEDURE, PASS(this) :: ingest_to_field_data_1D_nn

    ! PROCEDURE, PASS(this) :: ingest_to_field_data_3D_sigma
    PROCEDURE, PASS(this) :: ingest_to_field_data_3D_ZRHght_s
    PROCEDURE, PASS(this) :: ingest_to_field_data_3D_ZRHght_d
    PROCEDURE, PUBLIC :: initialize

    FINAL :: destructor
  END TYPE ModelCoupler_t

CONTAINS

  SUBROUTINE initialize(this, configFile, geometry)
    CLASS(ModelCoupler_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(Geometry_t), TARGET, INTENT(IN) :: geometry

    ! Local variables:
    INTEGER(i_kind) :: i, istatus
    INTEGER(i_kind), ALLOCATABLE :: numSmooth(:)

    this%configFile = configFile
    this%geometry => geometry

    ALLOCATE (this%numSmoothTopo(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
    this%numSmoothTopo = 0
    IF (yaml_get_var(this%configFile, 'geometry', 'topoSmooth', numSmooth) == 0) THEN
      DO i = this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
        this%numSmoothTopo(i) = numSmooth(i - this%geometry%mg%mg_coarsest + 1)
      END DO
    ELSE
      PRINT *, 'ModelCoupler: no topography smoothing yaml found. this%numSmoothTopo is all zeros.'
    END IF

    WRITE (*, 3) (this%numSmoothTopo(i), i=this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest)
3   FORMAT('ModelCoupler topography smooth: ', 20(1X, I2))
  END SUBROUTINE

  SUBROUTINE gen_interp_coeffs(this, cellCntrModel, HorNumGridModel, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_kind), INTENT(IN) :: cellCntrModel(2, HorNumGridModel)

    PRINT *, 'Before the constructor of InterpHPNew_t.'
    CALL this%interhpN%initializeInterpHPNew(this%geometry%mpdd, sg_mp, &
                                             cellCntrModel, HorNumGridModel)
  END SUBROUTINE gen_interp_coeffs

  ! SUBROUTINE ingest_to_field_data_3D_sigma(this, vLevel, zSigma, valueModel, fieldData, sg_sp, sg_mp)
  !   USE Interp1D_m, ONLY: interp1d
  !   CLASS(ModelCoupler_t) :: this
  !   TYPE(SingleGrid_t), INTENT(IN) :: sg_sp, sg_mp
  !   INTEGER(i_kind), INTENT(IN) :: vLevel
  !   REAL(r_kind), INTENT(IN) :: valueModel(:, :, :)
  !   REAL(r_kind), INTENT(IN) :: zSigma(:)
  !   REAL(r_kind), INTENT(INOUT) :: fieldData(:, :, :)
  !   REAL(r_kind), ALLOCATABLE :: valueDA_sp(:, :, :)
  !   REAL(r_kind), ALLOCATABLE :: valueDA_mp(:, :, :)
  !   INTEGER(i_kind) :: i

  !   ! Interpolat horizontally on single proc
  !   IF (sg_sp%isBaseProc()) ALLOCATE (valueDA_sp(vLevel, sg_sp%num_icell_global, sg_sp%tSlots))

  !   DO i = 1, sg_mp%tSlots
  !     CALL this%interhp%interp(vLevel, valueModel(:, :, i), valueDA_sp(:, :, i))
  !   END DO

  !   IF (sg_sp%isBaseProc()) FORALL (i=1:sg_sp%num_icell_global, sg_sp%cell_type(i) .gt. 0) valueDA_sp(:, i, :) = 0.0D0 ! Set the points on boundary to zero.

  !   ! Distribute to active procs
  !   ALLOCATE (valueDA_mp(vLevel, sg_mp%num_cell, sg_mp%tSlots))

  !   DO i = 1, sg_mp%tSlots
  !     CALL sg_sp%m_interp_points_on_bndy_linear(vLevel, valueDA_sp(:, :, i))
  !   END DO

  !   CALL sg_mp%DistGridRealWithHaloExForFieldGrid(valueDA_sp, valueDA_mp, &
  !                                                 [vLevel, sg_mp%num_icell_global, sg_mp%tSlots])

  !   ! Implement the vertical interpolation
  !   DO i = 1, sg_mp%tSlots
  !     CALL interp1d(zSigma, valueDA_mp(:, :, i), sg_mp%sigma, fieldData(:, :, i))
  !   END DO

  !   IF (sg_sp%isBaseProc()) DEALLOCATE (valueDA_sp)
  !   IF (ALLOCATED(valueDA_mp)) DEALLOCATE (valueDA_mp)
  ! END SUBROUTINE ingest_to_field_data_3D_sigma

  SUBROUTINE ingest_to_field_data_3D_ZRHght_d(this, vLevel, zRHghtModelAtDA, valueModel, fieldData, sg_mp)
    USE Interp1D_m, ONLY: interp1d
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), INTENT(IN) :: valueModel(:, :, :)
    REAL(r_kind), INTENT(IN) :: zRHghtModelAtDA(:, :)
    REAL(r_kind), INTENT(INOUT) :: fieldData(:, :, :)
    REAL(r_kind), ALLOCATABLE :: valueDA_mp(:, :, :)
    INTEGER(i_kind) :: i

    IF (sg_mp%isActiveProc()) THEN
      ! Interpolat horizontally on single proc
      ALLOCATE (valueDA_mp(vLevel, sg_mp%num_cell, sg_mp%tSlots))

      DO i = 1, sg_mp%tSlots
        CALL this%interhpN%interpN(vLevel, valueModel(:, :, i), valueDA_mp(:, :, i))
        CALL sg_mp%m_interp_points_on_bndy_linear(vLevel, valueDA_mp(:, :, i))
        CALL sg_mp%ExchangeMatOnHalo2D(vLevel, valueDA_mp(:, :, i))
      END DO

      ! Implement the vertical interpolation
      DO j = 1, sg_mp%tSlots
        ! Implement the vertical interpolation
        DO i = 1, sg_mp%num_cell
          CALL interp1d(zRHghtModelAtDA(:, i), valueDA_mp(:, i, j), &
                        sg_mp%zHght(:, i), fieldData(:, i, j))
        END DO
      END DO

      IF (ALLOCATED(valueDA_mp)) DEALLOCATE (valueDA_mp)
    END IF
  END SUBROUTINE

  SUBROUTINE ingest_to_field_data_3D_ZRHght_s(this, vLevel, zRHghtModelAtDA, valueModel, fieldData, sg_mp)
    USE Interp1D_m, ONLY: interp1d
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_single), INTENT(IN) :: valueModel(:, :, :)
    REAL(r_kind), INTENT(IN) :: zRHghtModelAtDA(:, :)
    REAL(r_kind), INTENT(INOUT) :: fieldData(:, :, :)
    REAL(r_kind), ALLOCATABLE :: valueDA_mp(:, :, :)
    INTEGER(i_kind) :: i

    IF (sg_mp%isActiveProc()) THEN
      ! Interpolat horizontally on single proc
      ALLOCATE (valueDA_mp(vLevel, sg_mp%num_cell, sg_mp%tSlots))

      DO i = 1, sg_mp%tSlots
        CALL this%interhpN%interpN(vLevel, valueModel(:, :, i), valueDA_mp(:, :, i))
        CALL sg_mp%m_interp_points_on_bndy_linear(vLevel, valueDA_mp(:, :, i))
        CALL sg_mp%ExchangeMatOnHalo2D(vLevel, valueDA_mp(:, :, i))
      END DO

      ! Implement the vertical interpolation
      DO j = 1, sg_mp%tSlots
        ! Implement the vertical interpolation
        DO i = 1, sg_mp%num_cell
          CALL interp1d(zRHghtModelAtDA(:, i), valueDA_mp(:, i, j), &
                        sg_mp%zHght(:, i), fieldData(:, i, j))
        END DO
      END DO

      IF (ALLOCATED(valueDA_mp)) DEALLOCATE (valueDA_mp)
    END IF
  END SUBROUTINE

  SUBROUTINE ingest_to_topo_and_update_zHght(this, topoModel, HorNumGridModel, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_single), INTENT(IN) :: topoModel(HorNumGridModel)

    CALL this%interhpN%interp_singleLevel_s(topoModel, sg_mp%topo)
    CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
    CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)

    CALL sg_mp%update_zHght_from_topo_and_sigma
  END SUBROUTINE ingest_to_topo_and_update_zHght

  SUBROUTINE ingest_to_topo(this, topoModel, HorNumGridModel, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_single), INTENT(IN) :: topoModel(HorNumGridModel)

    CALL this%interhpN%interp_singleLevel_s(topoModel, sg_mp%topo)
    CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
    CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)
  END SUBROUTINE

  SUBROUTINE ingest_to_topo_and_update_zHght_d(this, topoModel, HorNumGridModel, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_kind), INTENT(IN) :: topoModel(HorNumGridModel)

    CALL this%interhpN%interp_singleLevel_d(topoModel, sg_mp%topo)
    CALL sg_mp%m_interp_points_on_bndy_linear(1, sg_mp%topo)
    CALL sg_mp%ExchangeMatOnHalo2D(1, sg_mp%topo)

    CALL sg_mp%update_zHght_from_topo_and_sigma
  END SUBROUTINE ingest_to_topo_and_update_zHght_d

  SUBROUTINE ingest_to_field_data_2D(this, vLevel, input_sp, output_mp, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    REAL(r_single), INTENT(IN) :: input_sp(:, :)
    REAL(r_kind), INTENT(INOUT) :: output_mp(:, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel

    CALL this%interhpN%interp_multiLevel_s(vLevel, input_sp, output_mp)
    CALL sg_mp%m_interp_points_on_bndy_linear(vLevel, output_mp)
    CALL sg_mp%ExchangeMatOnHalo2D(vLevel, output_mp)

  END SUBROUTINE ingest_to_field_data_2D

  SUBROUTINE ingest_to_field_data_1D(this, HorNumGridModel, input_sp, output_mp, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_single), INTENT(IN) :: input_sp(HorNumGridModel)
    REAL(r_kind), INTENT(INOUT) :: output_mp(:)

    CALL this%interhpN%interp_singleLevel_s(input_sp, output_mp)
    CALL sg_mp%m_interp_points_on_bndy_linear(1, output_mp)
    CALL sg_mp%ExchangeMatOnHalo2D(1, output_mp)

  END SUBROUTINE ingest_to_field_data_1D

  SUBROUTINE ingest_to_field_data_1D_nn(this, HorNumGridModel, input_sp, output_mp, sg_mp)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_mp
    INTEGER(i_kind), INTENT(IN) :: HorNumGridModel
    REAL(r_single), INTENT(IN) :: input_sp(HorNumGridModel)
    REAL(r_kind), INTENT(INOUT) :: output_mp(:)

    CALL this%interhpN%interp_singleLevel_nearest_s(input_sp, output_mp)
    CALL sg_mp%m_interp_points_on_bndy_linear(1, output_mp)
    CALL sg_mp%ExchangeMatOnHalo2D(1, output_mp)

  END SUBROUTINE ingest_to_field_data_1D_nn

  SUBROUTINE ingest_to_field_data_MT(this, vLevel, valueModel, fieldData, sg_mp, nt)
    CLASS(ModelCoupler_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_mp
    REAL(r_kind), ALLOCATABLE, INTENT(IN) :: valueModel(:, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: fieldData(:, :, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel
    INTEGER(i_kind) :: t
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: nt

    IF (PRESENT(nt)) THEN
      DO t = 1, nt
        CALL this%interhpN%interpN(vLevel, valueModel(:, :, t), fieldData(:, :, t))
        CALL sg_mp%m_interp_points_on_bndy_linear(vLevel, fieldData(:, :, t))
        CALL sg_mp%ExchangeMatOnHalo2D(vLevel, fieldData(:, :, t))
      END DO
    ELSE
      DO t = 1, sg_mp%tSlots
        CALL this%interhpN%interpN(vLevel, valueModel(:, :, t), fieldData(:, :, t))
        CALL sg_mp%m_interp_points_on_bndy_linear(vLevel, fieldData(:, :, t))
        CALL sg_mp%ExchangeMatOnHalo2D(vLevel, fieldData(:, :, t))
      END DO
    END IF

  END SUBROUTINE ingest_to_field_data_MT

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(ModelCoupler_t), INTENT(INOUT) :: this

    PRINT *, 'In destructor of ModelCoupler_t'

  END SUBROUTINE destructor

END MODULE ModelCoupler_m
