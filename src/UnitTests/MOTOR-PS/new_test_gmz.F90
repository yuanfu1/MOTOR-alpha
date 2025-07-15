PROGRAM new_test_gmz
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  USE Verification_m
  USE gzm_tlm_Wrapper_m, ONLY: Divergen_TLM_Wrapper, Jacobian_TLM_Wrapper, Laplacia_TLM_Wrapper
  USE gzm_adj_Wrapper_m, ONLY: Divergen_AD_Wrapper, Jacobian_AD_Wrapper, Laplacia_AD_Wrapper
  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  INTEGER(i_kind) :: i
  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), d_oprand(:, :), adj_oprand(:, :), pert_oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: H(:, :), R_inv(:, :), grad_f(:), y_obs(:)
  REAL(r_kind) :: epsilon
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ! Associate SingleGrid with gzm and gzm_tlm
  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  ! Initialize gzm, gzm_tlm and gzm_adj
  gzm%sg => sg
  CALL gzm_tlm%initialize(sg)
  gzm_adj%sg => sg

  ! Ensure arrays are allocated and initialized properly
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))
  ALLOCATE (H(sg%vLevel, sg%num_icell), R_inv(sg%num_icell, sg%num_icell), y_obs(sg%vLevel), grad_f(sg%vLevel))

  ! Initialize input data
  opr_left = 0.05D0
  opr_right = 0.05D0
  oprand = 0.05D0
  d_value = 0.0D0
  VALUE = 0.0D0
  adj_value = 0.0D0
  H = 1.0_R_KIND
  R_inv = 1.0_R_KIND
  y_obs = 2.0_R_KIND
  DO i = 1, sg%vLevel
    d_opr_left(i, :) = (DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))) / 1.0D5 / 2.0D0
    d_opr_right(i, :) = (DCOS(sg%cell_cntr(1, :)) + DSIN(sg%cell_cntr(2, :))) / 1.0D5 / 2.0D0
    d_oprand(i, :) = DSIN(sg%cell_cntr(1, :)) / 1000D0
    adj_opr_left(i, :) = 0.0D0
    adj_opr_right(i, :) = 0.0D0
    adj_oprand(i, :) = 0.0D0
  END DO

  opr_left = d_opr_left * 1.0D5
  opr_right = d_opr_right * 1.0D5
  oprand = d_oprand

  epsilon = 1E-6_R_KIND  ! Small perturbation

  ! Test Divergen
  CALL gzm%Divergen(opr_left, opr_right, VALUE)
  CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
  pert_opr_left = opr_left + epsilon * d_opr_left
  pert_opr_right = opr_right + epsilon * d_opr_right
  !CALL gzm%Divergen(pert_opr_left, pert_opr_right, pert_value)
  !CALL Verify_TLM_FiniteDifference(gzm_tlm, opr_left, opr_right, d_value, Divergen_TLM_Wrapper, epsilon, sg, "Divergen")

  ! Initialize adjoint variables for Divergen
  adj_value = d_value
  adj_opr_left = 0.0_R_KIND
  adj_opr_right = 0.0_R_KIND
  CALL gzm_adj%Divergen_AD(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)

  ! Verify Divergen (TL-AD consistency)
  CALL Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, "Divergen")
  CALL Verify_AD_FiniteDifference(gzm_adj, opr_left, opr_right, adj_value, Divergen_AD_Wrapper, epsilon, sg, "Divergen")

  ! Test Jacobian
  CALL gzm%Jacobian(opr_left, opr_right, VALUE)
  CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
  pert_opr_left = opr_left + epsilon * d_opr_left
  pert_opr_right = opr_right + epsilon * d_opr_right
  CALL gzm%Jacobian(pert_opr_left, pert_opr_right, pert_value)
  CALL Verify_TLM_FiniteDifference(gzm_tlm, opr_left, opr_right, d_value, Jacobian_TLM_Wrapper, epsilon, sg, "Jacobian")

  ! Initialize adjoint variables for Jacobian
  adj_value = d_value
  adj_opr_left = 0.0_R_KIND
  adj_opr_right = 0.0_R_KIND
  CALL gzm_adj%Jacobian_AD(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)

  ! Verify Jacobian (TL-AD consistency)
  CALL Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, "Jacobian")
  CALL Verify_AD_FiniteDifference(gzm_adj, opr_left, opr_right, adj_value, Jacobian_AD_Wrapper, epsilon, sg, "Jacobian")

  ! Test Laplacia
  CALL gzm%Laplacia(oprand, VALUE)
  CALL gzm_tlm%Laplacia_TLM(oprand, VALUE, d_oprand, d_value)
  pert_oprand = oprand + epsilon * d_oprand
  CALL gzm%Laplacia(pert_oprand, pert_value)
  CALL Verify_TLM_FiniteDifference(gzm_tlm, oprand, oprand, d_value, Laplacia_TLM_Wrapper, epsilon, sg, "Laplacia")

  ! Initialize adjoint variables for Laplacia
  adj_value = d_value
  adj_oprand = 0.0_R_KIND
  CALL gzm_adj%Laplacia_AD(oprand, VALUE, adj_oprand, adj_value)

  ! Verify Laplacia (TL-AD consistency)
  CALL Verify_TLM_AD(d_oprand, adj_oprand, d_oprand, adj_oprand, sg, "Laplacia")
  CALL Verify_AD_FiniteDifference(gzm_adj, oprand, VALUE, adj_value, Laplacia_AD_Wrapper, epsilon, sg, "Laplacia")

  ! Clean up
  DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, pert_opr_left, pert_opr_right, pert_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_oprand, H, R_inv, y_obs, grad_f)

  CALL mpddGlob%finalize
END PROGRAM new_test_gmz
