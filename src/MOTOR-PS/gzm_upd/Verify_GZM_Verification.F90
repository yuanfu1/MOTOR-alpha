MODULE Verify_GZM_Verification_m
  USE kinds_m, ONLY: r_kind
  USE Verification_m
  USE gzm_m, ONLY: gzm_t
  USE gzm_tlm_m, ONLY: gzm_tlm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Test_And_Verify(gzm, gzm_tlm, gzm_adj, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value, adj_opr_left, adj_opr_right, adj_oprand, adj_value, pert_opr_left, pert_opr_right, pert_oprand, pert_value, epsilon, sg, op_name)
    CLASS(gzm_t), INTENT(INOUT) :: gzm
    CLASS(gzm_tlm_t), INTENT(INOUT) :: gzm_tlm
    CLASS(gzm_adj_t), INTENT(INOUT) :: gzm_adj
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :), d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), epsilon
    REAL(r_kind), INTENT(INOUT) :: VALUE(:, :), d_value(:, :), adj_value(:, :), pert_opr_left(:, :), pert_opr_right(:, :), pert_oprand(:, :), pert_value(:, :), adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
    REAL(r_kind), ALLOCATABLE :: gradient(:, :)
    TYPE(SingleGrid_t), POINTER :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name

    ! 调用前向模式的操作符
    SELECT CASE (op_name)
    CASE ("Divergen")
      CALL gzm%Divergen(opr_left, opr_right, VALUE)
      CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
      CALL gzm_adj%Divergen_AD(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)
    CASE ("Jacobian")
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
      CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
      CALL gzm_adj%Jacobian_AD(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)
    CASE ("Laplacia")
      CALL gzm%Laplacia(oprand, VALUE)
      CALL gzm_tlm%Laplacia_TLM(oprand, VALUE, d_oprand, d_value)
      CALL gzm_adj%Laplacia_AD(oprand, VALUE, adj_oprand, adj_value)
    END SELECT

    ! 使用有限差分法验证切线模式 (TLM)
    pert_opr_left = opr_left + epsilon * d_opr_left
    pert_opr_right = opr_right + epsilon * d_opr_right
    pert_oprand = oprand + epsilon * d_oprand

    SELECT CASE (op_name)
    CASE ("Divergen")
      CALL gzm%Divergen(pert_opr_left, pert_opr_right, pert_value)
    CASE ("Jacobian")
      CALL gzm%Jacobian(pert_opr_left, pert_opr_right, pert_value)
    CASE ("Laplacia")
      CALL gzm%Laplacia(pert_oprand, pert_value)
    END SELECT

    ! 调用所有验证方法
    CALL Verify_TL_FiniteDifference(VALUE, d_value, pert_value, epsilon, sg, op_name)
    CALL Verify_AD_FiniteDifference(VALUE, adj_value, pert_value, epsilon, sg, op_name)
    CALL Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, op_name)
    !CALL Verify_Kinetic_Energy(value, d_value, op_name)
    !CALL Verify_Total_Energy(value, adj_value, op_name)
    CALL Verify_Model_Equations(d_value, adj_value, op_name)
    !CALL Verify_Spectral_Radius(d_value, op_name)
    !CALL Perturbation_Growth_Analysis(value, d_value, epsilon, 100, op_name)
    !CALL Verify_Adjoint_Objective(opr_left(1,:), opr_right(1,:), oprand(1,:), d_oprand(1,:), adj_value(1,:), op_name)
    ! Allocate gradient array
    ALLOCATE (gradient(SIZE(d_value, 1), SIZE(d_value, 2)))
    ! Compute gradient before verifying the inner product
    CALL Compute_Gradient(adj_value, gradient, sg, op_name)
    CALL Verify_Inner_Product(adj_value, gradient, sg, op_name)

    !CALL Verify_Inner_Product(adj_value(1,:), d_value(1,:), op_name)
    CALL Verify_Adjoint(VALUE, adj_value, pert_value, epsilon, sg, op_name)
  END SUBROUTINE Test_And_Verify

END MODULE Verify_GZM_Verification_m
