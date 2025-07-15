PROGRAM unit_test_gzm_tlm_ad_grid
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  USE Verification_m, ONLY: Verify_Adjoint_AutomaticDifferentiation, Verify_TLM_AD
  IMPLICIT NONE

  ! 声明变量
  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  INTEGER(i_kind) :: i, j, k, ic
  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), d_oprand(:, :), adj_oprand(:, :), pert_oprand(:, :)
  REAL(r_kind) :: epsilon
  REAL(r_kind), ALLOCATABLE :: random_perturbation(:, :)

  CHARACTER(LEN=1024) :: configFile, ncOutputFile
  CHARACTER(LEN=1024) :: op_name

  !REAL(r_kind), ALLOCATABLE :: x(:), y_obs(:), H(:,:), R_inv(:,:), adj_output(:)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ! 初始化
  !epsilon = 0.005

  ! Associate SingleGrid with gzm and gzm_tlm
  sg => geometry%mg%sg(6)

  sg%bdy_type = sg%cell_type

  ! Initialize gzm, gzm_tlm and gzm_adj
  gzm%sg => sg
  gzm_tlm%sg => sg
  gzm_adj%sg => sg

  ! Ensure arrays are allocated and initialized properly
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))

  ALLOCATE (random_perturbation(SIZE(opr_left, 1), SIZE(opr_left, 2)))
!  CALL random_number(random_perturbation)
!   d_value = 0.005D0
!   value = 0.0005D0
!   adj_value = 0.0D0
!   DO i = 1, sg%vLevel
!     d_opr_left(i,:) = (DSIN(sg%cell_cntr(1,:))+DCOS(sg%cell_cntr(2,:)))/1.0D5/2.0D0
!     d_opr_right(i,:) = (DCOS(sg%cell_cntr(1,:))+DSIN(sg%cell_cntr(2,:)))/1.0D5/2.0D0
!     d_oprand(i,:) = DSIN(sg%cell_cntr(1,:)) / 1000D0
!     adj_opr_left(i,:) = 0.0D0
!     adj_opr_right(i,:) = 0.0D0
!     adj_oprand(i,:) = 0.0D0
!   END DO

!   opr_left = d_opr_left * 1.0D5
!   opr_right = d_opr_right * 1.0D5
!   oprand = d_oprand

  ! CALL random_number(opr_left)
  ! CALL random_number(opr_right)
  ! CALL random_number(d_opr_left)
  ! CALL random_number(d_opr_right)
  ! CALL random_number(adj_value)

  ! Adjust opr_left and opr_right values

!   ! 1. Multiply by a constant factor
!   opr_left = opr_left * 10.0_r_kind
!   opr_right = opr_right * 10.0_r_kind

!   ! 2. Add a random perturbation
!   opr_left = opr_left + random_perturbation * 0.1_r_kind
!   opr_right = opr_right + random_perturbation * 0.1_r_kind

!   ! 3. Apply exponential transformation
!   opr_left = EXP(opr_left)
!   opr_right = EXP(opr_right)

!   ! 4. Add a constant value
!   opr_left = opr_left + 5.0_r_kind
!   opr_right = opr_right + 5.0_r_kind

! opr_left = opr_left + random_perturbation * 0.5_r_kind
! opr_right = opr_right + random_perturbation * 0.5_r_kind
!   epsilon = 10.0E0_r_kind  ! Small perturbation

  ! Print sizes of arrays to ensure correct allocation
  PRINT *, "opr_left size:", SIZE(opr_left, 1), SIZE(opr_left, 2)
  PRINT *, "opr_right size:", SIZE(opr_right, 1), SIZE(opr_right, 2)
  PRINT *, "value size:", SIZE(VALUE, 1), SIZE(VALUE, 2)
  PRINT *, "sg%vLevel:", sg%vLevel
  PRINT *, "sg%num_icell:", sg%num_icell

  ! 调整为一维数组并分配内存
!   ALLOCATE(x(sg%vLevel * sg%num_icell))
!   ALLOCATE(y_obs(sg%vLevel * sg%num_icell))
!   ALLOCATE(H(SIZE(x), SIZE(x)))
!   ALLOCATE(R_inv(SIZE(x), SIZE(x)))
!   ALLOCATE(adj_output(SIZE(x)))

!   ! Print sizes of new arrays
!   PRINT *, "x size:", SIZE(x)
!   PRINT *, "y_obs size:", SIZE(y_obs)
!   PRINT *, "H size:", SIZE(H, 1), SIZE(H, 2)
!   PRINT *, "R_inv size:", SIZE(R_inv, 1), SIZE(R_inv, 2)
!   PRINT *, "adj_output size:", SIZE(adj_output)

  ! 随机初始化数据
  CALL RANDOM_NUMBER(opr_left)
  CALL RANDOM_NUMBER(opr_right)
  CALL RANDOM_NUMBER(d_opr_left)
  CALL RANDOM_NUMBER(d_opr_right)
  CALL RANDOM_NUMBER(adj_value)
  CALL RANDOM_NUMBER(random_perturbation)

  opr_left = opr_left * 10.0_R_KIND
  opr_right = opr_right * 10.0_R_KIND
  opr_left = opr_left + random_perturbation * 0.1_R_KIND
  opr_right = opr_right + random_perturbation * 0.1_R_KIND
  opr_left = EXP(opr_left)
  opr_right = EXP(opr_right)
  opr_left = opr_left + 5.0_R_KIND
  opr_right = opr_right + 5.0_R_KIND

  VALUE = 1.0_R_KIND
  pert_value = 1.1_R_KIND
  d_value = 0.1_R_KIND
  adj_value = 0.05_R_KIND
!   x = 1.0_r_kind
!   y_obs = 1.2_r_kind
!   H = 1.0_r_kind
!   R_inv = 1.0_r_kind
!   adj_output = 0.1_r_kind
  epsilon = 1.0E-5_R_KIND

  !PRINT *, "Before perturbation:"
  !PRINT *, "opr_left: ", opr_left
  !PRINT *, "d_opr_left: ", d_opr_left
  !PRINT *, "opr_right: ", opr_right
  !PRINT *, "d_opr_right: ", d_opr_right

  ! Test Divergen
  ! 计算原始散度
  !CALL gzm%Divergen(opr_left, opr_right, value)

  ! 计算切线模式
  CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)

!   ! Calculate perturbed opr_left and opr_right
! !   DO k = 1, sg%vLevel
! !     DO ic = 1, sg%num_icell
! !       IF (sg%bdy_type(ic) < 1) THEN
! !          pert_opr_left(k, ic) = opr_left(k, ic) + epsilon * d_opr_left(k, ic)
! !          pert_opr_right(k, ic) = opr_right(k, ic) + epsilon * d_opr_right(k, ic)
! !       ELSE
! !          pert_opr_left(k, ic) = opr_left(k, ic)
! !          pert_opr_right(k, ic) = opr_right(k, ic)
! !       END IF
! !     END DO
! !   END DO

!     ! 计算扰动后的值
!     CALL gzm%Divergen(opr_left + epsilon * d_opr_left, opr_right + epsilon * d_opr_right, pert_value)  ! 注意使用 pert_value
! print  *,  "Pert_value  ", pert_value(1,:)
!     ! 计算伴随模式
  CALL gzm_adj%Divergen_AD(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)
! PRINT *, "Opr_left:", opr_left(1,:)
! PRINT *, "Opr_right:", opr_right(1,:)
! PRINT *, "D_opr_left:", d_opr_left(1,:)
! PRINT *, "D_opr_right:", d_opr_right(1,:)
! PRINT *, "epsilon * d_opr_left:", epsilon * d_opr_left(1,:)
! PRINT *, "epsilon * d_opr_right:", epsilon * d_opr_right(1,:)

  ! 计算有限差分近似
  ! finite_diff = SUM((pert_value - value) * adj_value) / epsilon  ! 使用 pert_value 和 value
  ! dot_product_adj = SUM(d_opr_left * adj_opr_left) + SUM(d_opr_right * adj_opr_right)

  ! PRINT *, "Finite difference (FD) approximation: ", finite_diff
  ! PRINT *, "Dot product from Adjoint: ", dot_product_adj

  ! ! 检查有限差分与伴随模式计算结果是否接近
  ! IF (ABS(finite_diff - dot_product_adj) < 1.0e-6) THEN
  !     PRINT *, "The TLM and Adjoint models are consistent with FD!"
  ! ELSE
  !     PRINT *, "There is a discrepancy between the TLM/Adjoint models and FD."
  ! END IF

!    CALL Verify_TL_FiniteDifference(value, d_value, pert_value, epsilon, sg, "Divergen")

!   ! Initialize adjoint variables for Divergen
!   adj_value = d_value
!   adj_opr_left = 0.0_r_kind
!   adj_opr_right = 0.0_r_kind

!   ! Verify Divergen (TL-AD consistency)
  op_name = "Divergen"
  !CALL Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, "Divergen")
  CALL Verify_TLM_AD(d_value, adj_value, d_opr_right, adj_opr_right, sg, "Divergen")

  ! Call the Adjoint model verification using Automatic Differentiation
  !CALL Verify_Adjoint_AutomaticDifferentiation(x, y_obs, H, R_inv, adj_output, sg, op_name)

  PRINT *, 'Verification completed. Check verification_results.txt for details.'
!   ! Test Jacobian
!   CALL gzm%Jacobian(opr_left, opr_right, value)
!   CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, value, d_opr_left, d_opr_right, d_value)
!   pert_opr_left = opr_left + epsilon * d_opr_left
!   pert_opr_right = opr_right + epsilon * d_opr_right
!   CALL gzm%Jacobian(pert_opr_left, pert_opr_right, pert_value)
!   CALL Verify_TL_FiniteDifference(value, d_value, pert_value, epsilon, sg, "Jacobian")

!   ! Initialize adjoint variables for Jacobian
!   adj_value = d_value
!   adj_opr_left = 0.0_r_kind
!   adj_opr_right = 0.0_r_kind
!   CALL gzm_adj%Jacobian_AD(opr_left, opr_right, value, adj_opr_left, adj_opr_right, adj_value)

!   ! Verify Jacobian (TL-AD consistency)
  !CALL Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, "Jacobian")

!   ! Test Laplacia
!   CALL gzm%Laplacia(oprand, value)
!   CALL gzm_tlm%Laplacia_TLM(oprand, value, d_oprand, d_value)
!   pert_oprand = oprand + epsilon * d_oprand
!   CALL gzm%Laplacia(pert_oprand, pert_value)
!   CALL Verify_TL_FiniteDifference(value, d_value, pert_value, epsilon, sg, "Laplacia")

!   ! Initialize adjoint variables for Laplacia
!   adj_value = d_value
!   adj_oprand = 0.0_r_kind
!   CALL gzm_adj%Laplacia_AD(oprand, value, adj_oprand, adj_value)

!   ! Verify Laplacia (TL-AD consistency)
! ! Verify Laplacia (TL-AD consistency)
!   CALL Verify_TLM_AD(d_oprand, adj_oprand, d_oprand, adj_oprand, sg, "Laplacia")

!   ! 释放内存
!   DEALLOCATE(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, value, d_value, pert_value, adj_opr_left, adj_opr_right, adj_oprand, adj_value)

END PROGRAM unit_test_gzm_tlm_ad_grid
