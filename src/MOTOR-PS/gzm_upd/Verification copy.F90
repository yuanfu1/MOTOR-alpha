MODULE Verification_m
  USE kinds_m, ONLY: r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Objective_Function_m, ONLY: Compute_Objective_Function, Compute_Gradient_Objective_Function
  USE gzm_m, ONLY: gzm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  USE gzm_tlm_m, ONLY: gzm_tlm_t

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Verify_TLM_FiniteDifference(gzm_tlm, x, opr_right, tlm_output, func, epsilon, sg, op_name)
    USE gzm_tlm_m, ONLY: gzm_tlm_t
    IMPLICIT NONE
    TYPE(gzm_tlm_t) :: gzm_tlm
    REAL(r_kind), INTENT(IN) :: x(:, :), opr_right(:, :), tlm_output(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    INTERFACE
      SUBROUTINE func(this, x, y, z, tlm_input1, tlm_input2, tlm_output)
        USE kinds_m, ONLY: r_kind
        USE gzm_tlm_m, ONLY: gzm_tlm_t
        CLASS(gzm_tlm_t) :: this
        REAL(r_kind), INTENT(IN) :: x(:, :), y(:, :), z(:, :), tlm_input1(:, :), tlm_input2(:, :)
        REAL(r_kind), INTENT(OUT) :: tlm_output(:, :)
      END SUBROUTINE func
    END INTERFACE
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: obj_func, fd_approx
    INTEGER :: k, j, n
    REAL(r_kind), ALLOCATABLE :: x_pert(:, :), tlm_output_pert(:, :)
    REAL(r_kind) :: max_error, avg_error

    obj_func = 0.0_R_KIND
    ALLOCATE (x_pert(SIZE(x, 1), SIZE(x, 2)))
    ALLOCATE (tlm_output_pert(SIZE(tlm_output, 1), SIZE(tlm_output, 2)))
    max_error = 0.0_R_KIND
    avg_error = 0.0_R_KIND
    n = 0

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        x_pert = x
        x_pert(k, j) = x(k, j) + epsilon
        CALL func(gzm_tlm, x_pert, opr_right, opr_right, opr_right, opr_right, tlm_output_pert)
        fd_approx = (tlm_output_pert(k, j) - tlm_output(k, j)) / epsilon
        max_error = MAX(max_error, ABS(fd_approx - tlm_output(k, j)))
        avg_error = avg_error + ABS(fd_approx - tlm_output(k, j))
        n = n + 1
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (Finite Difference for TLM', op_name, '):', max_error
    PRINT *, 'Avg Error (Finite Difference for TLM', op_name, '):', avg_error

    IF (max_error < 1E-5) THEN
      PRINT *, 'Finite Difference for TLM', op_name, 'verification passed.'
    ELSE
      PRINT *, 'Finite Difference for TLM', op_name, 'verification failed.'
    END IF

    DEALLOCATE (x_pert, tlm_output_pert)
  END SUBROUTINE Verify_TLM_FiniteDifference

  SUBROUTINE Verify_AD_FiniteDifference(gzm_adj, x, opr_right, ad_output, grad_func, epsilon, sg, op_name)
    USE gzm_adj_m, ONLY: gzm_adj_t
    IMPLICIT NONE
    TYPE(gzm_adj_t) :: gzm_adj
    REAL(r_kind), INTENT(IN) :: x(:, :), opr_right(:, :), ad_output(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    INTERFACE
      SUBROUTINE grad_func(this, x, y, z, w, adj_output)
        USE kinds_m, ONLY: r_kind
        USE gzm_adj_m, ONLY: gzm_adj_t
        CLASS(gzm_adj_t) :: this
        REAL(r_kind), INTENT(IN) :: x(:, :), y(:, :), z(:, :), w(:, :)
        REAL(r_kind), INTENT(OUT) :: adj_output(:, :)
      END SUBROUTINE grad_func
    END INTERFACE
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: grad(:, :), grad_pert(:, :)
    REAL(r_kind), ALLOCATABLE :: x_temp(:, :)
    REAL(r_kind) :: fd_approx, max_error, avg_error
    INTEGER :: k, j, n

    ALLOCATE (grad(SIZE(ad_output, 1), SIZE(ad_output, 2)))
    ALLOCATE (grad_pert(SIZE(ad_output, 1), SIZE(ad_output, 2)))
    ALLOCATE (x_temp(SIZE(x, 1), SIZE(x, 2)))

    max_error = 0.0_R_KIND
    avg_error = 0.0_R_KIND
    n = 0

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        CALL grad_func(gzm_adj, x, opr_right, opr_right, opr_right, grad)
        x_temp = x
        x_temp(k, j) = x(k, j) + epsilon
        CALL grad_func(gzm_adj, x_temp, opr_right, opr_right, opr_right, grad_pert)
        fd_approx = (grad_pert(k, j) - grad(k, j)) / epsilon
        max_error = MAX(max_error, ABS(fd_approx - ad_output(k, j)))
        avg_error = avg_error + ABS(fd_approx - ad_output(k, j))
        n = n + 1
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (Finite Difference for AD', op_name, '):', max_error
    PRINT *, 'Avg Error (Finite Difference for AD', op_name, '):', avg_error

    IF (max_error < 1E-5) THEN
      PRINT *, 'Finite Difference for AD', op_name, 'verification passed.'
    ELSE
      PRINT *, 'Finite Difference for AD', op_name, 'verification failed.'
    END IF

    DEALLOCATE (grad, grad_pert, x_temp)
  END SUBROUTINE Verify_AD_FiniteDifference

  SUBROUTINE Verify_TLM(VALUE, d_value, pert_value, epsilon, sg, op_name)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), d_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
          avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (TLM', op_name, '):', max_error
    PRINT *, 'Avg Error (TLM', op_name, '):', avg_error

    IF (max_error < 1E-14) THEN
      PRINT *, 'TLM', op_name, 'verification passed.'
    ELSE
      PRINT *, 'TLM', op_name, 'verification failed.'
    END IF
  END SUBROUTINE Verify_TLM

  SUBROUTINE Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, op_name)
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), adj_opr_left(:, :)
    REAL(r_kind), INTENT(IN) :: d_opr_right(:, :), adj_opr_right(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: inner_product_tl, inner_product_ad, tolerance
    INTEGER :: k, j

    inner_product_tl = 0.0
    inner_product_ad = 0.0
    tolerance = 1E-14

    ! Compute the inner product <M delta x, delta y>
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_opr_left(k, j) * adj_opr_right(k, j)
      END DO
    END DO

    ! Compute the inner product <delta x, M* delta y>
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_ad = inner_product_ad + adj_opr_left(k, j) * d_opr_right(k, j)
      END DO
    END DO

    PRINT *, 'Inner Product (TL', op_name, '):', inner_product_tl
    PRINT *, 'Inner Product (AD', op_name, '):', inner_product_ad

    IF (ABS(inner_product_tl - inner_product_ad) < tolerance) THEN
      PRINT *, 'TL-AD', op_name, 'verification passed.'
    ELSE
      PRINT *, 'TL-AD', op_name, 'verification failed.'
    END IF
  END SUBROUTINE Verify_TLM_AD

  ! Verify the adjoint model (AD) with a unit perturbation method
  SUBROUTINE Verify_Adjoint(VALUE, adj_value, pert_value, epsilon, op_name)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), adj_value(:, :), pert_value(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, SIZE(VALUE, 1)
      DO j = 1, SIZE(VALUE, 2)
        max_error = MAX(max_error, ABS(adj_value(k, j) - pert_value(k, j) / epsilon))
        avg_error = avg_error + ABS(adj_value(k, j) - pert_value(k, j) / epsilon)
        n = n + 1
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (AD', op_name, '):', max_error
    PRINT *, 'Avg Error (AD', op_name, '):', avg_error

    IF (max_error < 1E-14) THEN
      PRINT *, 'AD', op_name, 'verification passed.'
    ELSE
      PRINT *, 'AD', op_name, 'verification failed.'
    END IF
  END SUBROUTINE Verify_Adjoint

  ! Verify kinetic energy balance
  SUBROUTINE Verify_Kinetic_Energy(initial_state, final_state, op_name)
    REAL(r_kind), INTENT(IN) :: initial_state(:, :), final_state(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: initial_kinetic_energy, final_kinetic_energy, diff

    initial_kinetic_energy = 0.5 * SUM(initial_state**2)
    final_kinetic_energy = 0.5 * SUM(final_state**2)
    diff = ABS(final_kinetic_energy - initial_kinetic_energy)

    PRINT *, "Verifying Kinetic Energy Balance for ", TRIM(op_name)
    PRINT *, "Initial Kinetic Energy: ", initial_kinetic_energy
    PRINT *, "Final Kinetic Energy: ", final_kinetic_energy
    PRINT *, "Difference: ", diff
  END SUBROUTINE Verify_Kinetic_Energy

  ! Verify total energy balance
  SUBROUTINE Verify_Total_Energy(initial_state, final_state, op_name)
    REAL(r_kind), INTENT(IN) :: initial_state(:, :), final_state(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: initial_total_energy, final_total_energy, diff

    initial_total_energy = 0.5 * SUM(initial_state**2)
    final_total_energy = 0.5 * SUM(final_state**2)
    diff = ABS(final_total_energy - initial_total_energy)

    PRINT *, "Verifying Total Energy Balance for ", TRIM(op_name)
    PRINT *, "Initial Total Energy: ", initial_total_energy
    PRINT *, "Final Total Energy: ", final_total_energy
    PRINT *, "Difference: ", diff
  END SUBROUTINE Verify_Total_Energy

  ! Verify the model equations using variational method
  SUBROUTINE Verify_Model_Equations(tlm_output, ad_output, op_name)
    REAL(r_kind), INTENT(IN) :: tlm_output(:), ad_output(:)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_diff

    max_diff = MAXVAL(ABS(tlm_output - ad_output))

    PRINT *, "Verifying Model Equations for ", TRIM(op_name)
    PRINT *, "Max difference: ", max_diff
  END SUBROUTINE Verify_Model_Equations

  ! Spectral radius analysis for the TLM
  SUBROUTINE Verify_Spectral_Radius(tlm_matrix, op_name)
    REAL(r_kind), INTENT(IN) :: tlm_matrix(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: eigenvalues(:)
    REAL(r_kind) :: spectral_radius

    ALLOCATE (eigenvalues(SIZE(tlm_matrix, 1)))
    CALL Eigenvalue_Analysis(tlm_matrix, eigenvalues)
    spectral_radius = MAXVAL(ABS(eigenvalues))

    PRINT *, "Spectral Radius Analysis for ", TRIM(op_name)
    PRINT *, "Spectral Radius: ", spectral_radius

    DEALLOCATE (eigenvalues)
  END SUBROUTINE Verify_Spectral_Radius

  ! Calculate eigenvalues of a matrix (simple implementation)
  SUBROUTINE Eigenvalue_Analysis(matrix, eigenvalues)
    REAL(r_kind), INTENT(IN) :: matrix(:, :)
    REAL(r_kind), INTENT(OUT) :: eigenvalues(:)
    INTEGER :: n, info
    REAL(r_kind), ALLOCATABLE :: work(:), wr(:), wi(:), vr(:, :), vl(:, :)

    n = SIZE(matrix, 1)
    ALLOCATE (work(4 * n), wr(n), wi(n), vr(n, n), vl(n, n))
    eigenvalues = 0.0_R_KIND

    CALL DGEEV('N', 'N', n, matrix, n, wr, wi, vl, n, vr, n, work, 4 * n, info)

    eigenvalues = CMPLX(wr, wi, kind=r_kind)

    DEALLOCATE (work, wr, wi, vr, vl)
  END SUBROUTINE Eigenvalue_Analysis

  ! Perturbation growth analysis
  SUBROUTINE Perturbation_Growth_Analysis(x, dx, dt, num_steps, op_name)
    REAL(r_kind), INTENT(IN) :: x(:, :), dx(:, :), dt
    INTEGER, INTENT(IN) :: num_steps
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: perturbed_x(:, :)
    REAL(r_kind) :: growth_rate, initial_norm, final_norm
    INTEGER :: i

    ALLOCATE (perturbed_x(SIZE(x, 1), SIZE(x, 2)))
    perturbed_x = x + dx
    initial_norm = SQRT(SUM(dx**2))

    ! Time integration loop
    DO i = 1, num_steps
      CALL TLM_Step(perturbed_x, dt)
    END DO

    final_norm = SQRT(SUM((perturbed_x - x)**2))
    growth_rate = LOG(final_norm / initial_norm) / (num_steps * dt)

    PRINT *, "Perturbation Growth Analysis for ", TRIM(op_name)
    PRINT *, "Initial Norm: ", initial_norm
    PRINT *, "Final Norm: ", final_norm
    PRINT *, "Growth Rate: ", growth_rate

    DEALLOCATE (perturbed_x)
  END SUBROUTINE Perturbation_Growth_Analysis

  ! Dummy subroutine for TLM step (to be implemented as per model specifics)
  SUBROUTINE TLM_Step(state, dt)
    REAL(r_kind), INTENT(INOUT) :: state(:, :)
    REAL(r_kind), INTENT(IN) :: dt
  END SUBROUTINE TLM_Step

  ! Verify the adjoint model using objective function and its gradient
  SUBROUTINE Verify_Adjoint_Objective(x, y_obs, H, R_inv, adj_output, op_name)
    REAL(r_kind), INTENT(IN) :: x(:), y_obs(:), H(:, :), R_inv(:, :), adj_output(:)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: grad_f(:)
    REAL(r_kind) :: inner_product_ad, inner_product_grad

    ALLOCATE (grad_f(SIZE(x)))

    ! Compute the gradient of the objective function
    CALL Compute_Gradient_Objective_Function(x, y_obs, H, R_inv, grad_f)
    ! Compute the inner products
    inner_product_ad = SUM(adj_output * grad_f)
    inner_product_grad = SUM(grad_f * grad_f)

    PRINT *, "Verifying Adjoint using Objective Function Gradient for ", TRIM(op_name)
    PRINT *, "Inner product (AD * Gradient): ", inner_product_ad
    PRINT *, "Inner product (Gradient * Gradient): ", inner_product_grad
    PRINT *, "Difference: ", ABS(inner_product_ad - inner_product_grad)

    DEALLOCATE (grad_f)
  END SUBROUTINE Verify_Adjoint_Objective

END MODULE Verification_m
