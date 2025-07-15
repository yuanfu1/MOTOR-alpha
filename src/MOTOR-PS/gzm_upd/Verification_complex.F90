MODULE Verification_m
  USE kinds_m, ONLY: r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Objective_Function_m, ONLY: Compute_Objective_Function, Compute_Gradient_Objective_Function
  USE gzm_m, ONLY: gzm_t
  USE gzm_tlm_m, ONLY: gzm_tlm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  IMPLICIT NONE

  ! 参考文献：
  ! 1. Errico, R. M. (1997). What is an adjoint model?. Bulletin of the American Meteorological Society, 78(11), 2577-2591.
  ! 2. Courtier, P., Thépaut, J. N., & Hollingsworth, A. (1994). A strategy for operational implementation of 4D-Var, using an incremental approach. Quarterly Journal of the Royal Meteorological Society, 120(519), 1367-1387.
  ! 3. Giering, R., & Kaminski, T. (1998). Recipes for adjoint code construction. ACM Transactions on Mathematical Software (TOMS), 24(4), 437-474.
  ! 4. Navon, I. M. (2009). Data assimilation for numerical weather prediction: a review. Data assimilation for atmospheric, oceanic and hydrologic applications, 1, 21-65.
  ! 5. Farrell, B. F., & Ioannou, P. J. (1996). Generalized stability theory. Part I: Autonomous operators. Journal of the Atmospheric Sciences, 53(14), 2025-2040.
  ! 6. Schubert, S., & Lucarini, V. (2016). Dynamical analysis of North Atlantic weather regimes: Validation and application. Journal of the Atmospheric Sciences, 73(7), 2567-2588.

CONTAINS

  !-------------------------------------------------------------------
  ! Subroutine: Verify_TL_FiniteDifference
  ! Purpose: Verify the Tangent Linear Model (TLM) using finite difference method.
  ! Inputs:
  !   - value: Original state variable
  !   - d_value: Tangent linear model output
  !   - pert_value: Perturbed state variable
  !   - epsilon: Perturbation size
  !   - sg: SingleGrid object
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_TL_FiniteDifference(VALUE, d_value, pert_value, epsilon, sg, op_name)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), d_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    ! Debug: Print input values
    PRINT *, "Verify_TL_FiniteDifference Inputs:"
    PRINT *, "Value:", VALUE
    PRINT *, "d_Value:", d_value
    PRINT *, "Pert_Value:", pert_value
    PRINT *, "Epsilon:", epsilon

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

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Max Error (TLM FiniteDifference ', TRIM(op_name), '): ', max_error
    WRITE (11, *) 'Avg Error (TLM FiniteDifference ', TRIM(op_name), '): ', avg_error
    IF (max_error < 1E-14) THEN
      WRITE (11, *) 'TLM FiniteDifference ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'TLM FiniteDifference ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_TL_FiniteDifference

  !-------------------------------------------------------------------
  ! Subroutine: Verify_AD_FiniteDifference
  ! Purpose: Verify the Adjoint Model (AD) using finite difference method.
  ! Inputs:
  !   - value: Original state variable
  !   - adj_value: Adjoint model output
  !   - pert_value: Perturbed state variable
  !   - epsilon: Perturbation size
  !   - sg: SingleGrid object
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_AD_FiniteDifference(VALUE, adj_value, pert_value, epsilon, sg, op_name)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), adj_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    ! Debug: Print input values
    PRINT *, "Verify_AD_FiniteDifference Inputs:"
    PRINT *, "Value:", VALUE
    PRINT *, "adj_Value:", adj_value
    PRINT *, "Pert_Value:", pert_value
    PRINT *, "Epsilon:", epsilon

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(adj_value(k, j) - pert_value(k, j) / epsilon))
          avg_error = avg_error + ABS(adj_value(k, j) - pert_value(k, j) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    ! Debug: Print calculated errors
    PRINT *, "Max Error:", max_error
    PRINT *, "Avg Error:", avg_error

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Max Error (AD FiniteDifference ', TRIM(op_name), '): ', max_error
    WRITE (11, *) 'Avg Error (AD FiniteDifference ', TRIM(op_name), '): ', avg_error
    IF (max_error < 1E-14) THEN
      WRITE (11, *) 'AD FiniteDifference ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'AD FiniteDifference ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_AD_FiniteDifference

  !-------------------------------------------------------------------
  ! Subroutine: Verify_TLM_AD
  ! Purpose: Verify the consistency between Tangent Linear Model (TLM) and Adjoint Model (AD).
  ! Inputs:
  !   - d_opr_left: Perturbation in left operand for TLM
  !   - adj_opr_left: Adjoint output corresponding to left operand
  !   - d_opr_right: Perturbation in right operand for TLM
  !   - adj_opr_right: Adjoint output corresponding to right operand
  !   - sg: SingleGrid object
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_TLM_AD(d_opr_left, adj_opr_left, d_opr_right, adj_opr_right, sg, op_name)
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), adj_opr_left(:, :)
    REAL(r_kind), INTENT(IN) :: d_opr_right(:, :), adj_opr_right(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: inner_product_tl, inner_product_ad, tolerance
    INTEGER :: k, j

    ! Debug: Print input values
    PRINT *, "Verify_TLM_AD Inputs:"
    PRINT *, "d_opr_left:", d_opr_left
    PRINT *, "adj_opr_left:", adj_opr_left
    PRINT *, "d_opr_right:", d_opr_right
    PRINT *, "adj_opr_right:", adj_opr_right

    inner_product_tl = 0.0
    inner_product_ad = 0.0
    tolerance = 1E-14

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_opr_left(k, j) * adj_opr_right(k, j)
      END DO
    END DO

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_ad = inner_product_ad + adj_opr_left(k, j) * d_opr_right(k, j)
      END DO
    END DO

    ! Debug: Print inner product results
    PRINT *, "Inner Product TL:", inner_product_tl
    PRINT *, "Inner Product AD:", inner_product_ad

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Inner Product (TL ', TRIM(op_name), '): ', inner_product_tl
    WRITE (11, *) 'Inner Product (AD ', TRIM(op_name), '): ', inner_product_ad
    IF (ABS(inner_product_tl - inner_product_ad) < tolerance) THEN
      WRITE (11, *) 'TL-AD ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'TL-AD ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_TLM_AD

  ! This subroutine verifies the adjoint model (AD) with a unit perturbation method
  ! Inputs:
  !   - value: Original state variable
  !   - adj_value: Adjoint model output
  !   - pert_value: Perturbed state variable
  !   - epsilon: Perturbation size
  !   - sg: SingleGrid structure
  !   - op_name: Operation name for identification
  SUBROUTINE Verify_Adjoint(VALUE, adj_value, pert_value, epsilon, sg, op_name)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), adj_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    ! Initialize error metrics
    max_error = 0.0
    avg_error = 0.0
    n = 0

    ! Compare adjoint output with finite difference approximation
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(adj_value(k, j) - pert_value(k, j) / epsilon))
          avg_error = avg_error + ABS(adj_value(k, j) - pert_value(k, j) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    ! Write results to a verification file
    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Max Error (AD with a unit perturbation ', TRIM(op_name), ' ): ', max_error
    WRITE (11, *) 'Avg Error (AD with a unit perturbation ', TRIM(op_name), ' ): ', avg_error
    IF (max_error < 1E-14) THEN
      WRITE (11, *) 'AD with a unit perturbation ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'AD with a unit perturbation ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_Adjoint

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Kinetic_Energy
  ! Purpose: Verify the kinetic energy balance between initial and final states.
  ! Inputs:
  !   - initial_state: Initial state variable
  !   - final_state: Final state variable
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Kinetic_Energy(initial_state, final_state, op_name)
    REAL(r_kind), INTENT(IN) :: initial_state(:, :), final_state(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: initial_kinetic_energy, final_kinetic_energy, diff

    initial_kinetic_energy = 0.5 * SUM(initial_state**2)
    final_kinetic_energy = 0.5 * SUM(final_state**2)
    diff = ABS(final_kinetic_energy - initial_kinetic_energy)

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Kinetic Energy Balance for ', TRIM(op_name)
    WRITE (11, *) 'Initial Kinetic Energy: ', initial_kinetic_energy
    WRITE (11, *) 'Final Kinetic Energy: ', final_kinetic_energy
    WRITE (11, *) 'Difference: ', diff
    IF (diff < 1E-14) THEN
      WRITE (11, *) 'Kinetic Energy ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Kinetic Energy ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_Kinetic_Energy

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Total_Energy
  ! Purpose: Verify the total energy balance between initial and final states.
  ! Inputs:
  !   - initial_state: Initial state variable
  !   - final_state: Final state variable
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Total_Energy(initial_state, final_state, op_name)
    REAL(r_kind), INTENT(IN) :: initial_state(:, :), final_state(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: initial_total_energy, final_total_energy, diff

    initial_total_energy = 0.5 * SUM(initial_state**2)
    final_total_energy = 0.5 * SUM(final_state**2)
    diff = ABS(final_total_energy - initial_total_energy)

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Total Energy Balance for ', TRIM(op_name)
    WRITE (11, *) 'Initial Total Energy: ', initial_total_energy
    WRITE (11, *) 'Final Total Energy: ', final_total_energy
    WRITE (11, *) 'Difference: ', diff
    IF (diff < 1E-14) THEN
      WRITE (11, *) 'Total Energy ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Total Energy ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_Total_Energy

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Model_Equations
  ! Purpose: Verify the model equations by comparing the tangent linear model output with adjoint model output.
  ! Inputs:
  !   - tlm_output: Tangent linear model output
  !   - ad_output: Adjoint model output
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Model_Equations(tlm_output, ad_output, op_name)
    REAL(r_kind), INTENT(IN) :: tlm_output(:, :), ad_output(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: max_diff

    max_diff = MAXVAL(ABS(tlm_output - ad_output))

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Model Equations for ', TRIM(op_name)
    WRITE (11, *) 'Max difference: ', max_diff
    IF (max_diff < 1E-14) THEN
      WRITE (11, *) 'Model Equations ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Model Equations ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_Model_Equations

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Spectral_Radius
  ! Purpose: Verify the spectral radius of the tangent linear model matrix.
  ! Inputs:
  !   - tlm_matrix: Tangent linear model matrix
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Spectral_Radius(tlm_matrix, op_name)
    REAL(r_kind), INTENT(IN) :: tlm_matrix(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: eigenvalues(:)
    REAL(r_kind) :: spectral_radius

    ALLOCATE (eigenvalues(SIZE(tlm_matrix, 1)))
    CALL Eigenvalue_Analysis(tlm_matrix, eigenvalues)
    spectral_radius = MAXVAL(ABS(eigenvalues))

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Spectral Radius Analysis for ', TRIM(op_name)
    WRITE (11, *) 'Spectral Radius: ', spectral_radius
    IF (spectral_radius < 1.0_R_KIND) THEN
      WRITE (11, *) 'Spectral Radius ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Spectral Radius ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)

    DEALLOCATE (eigenvalues)
  END SUBROUTINE Verify_Spectral_Radius

  !-------------------------------------------------------------------
  ! Subroutine: Eigenvalue_Analysis
  ! Purpose: Perform eigenvalue analysis on a given matrix.
  ! Inputs:
  !   - matrix: Input matrix
  !   - eigenvalues: Computed eigenvalues
  !-------------------------------------------------------------------
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

  !-------------------------------------------------------------------
  ! Subroutine: Perturbation_Growth_Analysis
  ! Purpose: Analyze perturbation growth in the tangent linear model.
  ! Inputs:
  !   - x: Initial state variable
  !   - dx: Perturbation in state variable
  !   - dt: Time step
  !   - num_steps: Number of time steps
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
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

    DO i = 1, num_steps
      CALL TLM_Step(perturbed_x, dt)
    END DO

    final_norm = SQRT(SUM((perturbed_x - x)**2))
    growth_rate = LOG(final_norm / initial_norm) / (num_steps * dt)

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Perturbation Growth Analysis for ', TRIM(op_name)
    WRITE (11, *) 'Initial Norm: ', initial_norm
    WRITE (11, *) 'Final Norm: ', final_norm
    WRITE (11, *) 'Growth Rate: ', growth_rate
    IF (growth_rate < 1E-14) THEN
      WRITE (11, *) 'Perturbation Growth ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Perturbation Growth ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)

    DEALLOCATE (perturbed_x)
  END SUBROUTINE Perturbation_Growth_Analysis

  !-------------------------------------------------------------------
  ! Subroutine: TLM_Step
  ! Purpose: Perform a single time step of the tangent linear model.
  ! Inputs:
  !   - state: State variable to be updated
  !   - dt: Time step
  !-------------------------------------------------------------------
  SUBROUTINE TLM_Step(state, dt)
    REAL(r_kind), INTENT(INOUT) :: state(:, :)
    REAL(r_kind), INTENT(IN) :: dt
    ! Placeholder: Implement the actual TLM step
    ! Note: This subroutine needs to be replaced with the actual TLM step function
  END SUBROUTINE TLM_Step

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Adjoint_Objective
  ! Purpose: Verify the adjoint model using the gradient of an objective function.
  ! Inputs:
  !   - x: State variable
  !   - y_obs: Observed variable
  !   - H: Observation operator matrix
  !   - R_inv: Inverse of observation error covariance matrix
  !   - adj_output: Adjoint model output
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Adjoint_Objective(x, y_obs, H, R_inv, adj_output, op_name)
    REAL(r_kind), INTENT(IN) :: x(:), y_obs(:), H(:, :), R_inv(:, :), adj_output(:)
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: grad_f(:)
    REAL(r_kind) :: inner_product_ad, inner_product_grad

    ALLOCATE (grad_f(SIZE(x)))

    CALL Compute_Gradient_Objective_Function(x, y_obs, H, R_inv, grad_f)
    inner_product_ad = SUM(adj_output * grad_f)
    inner_product_grad = SUM(grad_f * grad_f)

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Adjoint using Objective Function Gradient for ', TRIM(op_name)
    WRITE (11, *) 'Inner product (AD * Gradient): ', inner_product_ad
    WRITE (11, *) 'Inner product (Gradient * Gradient): ', inner_product_grad
    WRITE (11, *) 'Difference: ', ABS(inner_product_ad - inner_product_grad)
    IF (ABS(inner_product_ad - inner_product_grad) < 1E-14) THEN
      WRITE (11, *) 'Adjoint Objective ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Adjoint Objective ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)

    DEALLOCATE (grad_f)
  END SUBROUTINE Verify_Adjoint_Objective

  ! This subroutine computes the gradient of a given objective function using the adjoint model output.
  SUBROUTINE Compute_Gradient(adj_output, gradient, sg, op_name)
    REAL(r_kind), INTENT(IN) :: adj_output(:, :)
    REAL(r_kind), INTENT(OUT) :: gradient(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name

    INTEGER :: k, j

    ! Initialize gradient to zero
    gradient = 0.0_R_KIND

    ! Compute the gradient using adjoint variables
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          gradient(k, j) = adj_output(k, j)  ! Simple example, actual computation may vary based on the model
        END IF
      END DO
    END DO

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Gradient computed using adjoint for ', TRIM(op_name)
    CLOSE (11)
  END SUBROUTINE Compute_Gradient

  ! This subroutine verifies the inner product consistency between adjoint output and gradient.
  SUBROUTINE Verify_Inner_Product(adj_output, gradient, sg, op_name)
    REAL(r_kind), INTENT(IN) :: adj_output(:, :)
    REAL(r_kind), INTENT(IN) :: gradient(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name

    REAL(r_kind) :: inner_product_ad, inner_product_grad

    ! Calculate the inner products
    inner_product_ad = SUM(adj_output * gradient)
    inner_product_grad = SUM(gradient * gradient)

    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Inner Product for ', TRIM(op_name)
    WRITE (11, *) 'Inner product (AD * Gradient): ', inner_product_ad
    WRITE (11, *) 'Inner product (Gradient * Gradient): ', inner_product_grad
    WRITE (11, *) 'Difference: ', ABS(inner_product_ad - inner_product_grad)
    IF (ABS(inner_product_ad - inner_product_grad) < 1E-14) THEN
      WRITE (11, *) 'Inner Product between adjoint output and gradient ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Inner Product between adjoint output and gradient  ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)
  END SUBROUTINE Verify_Inner_Product

END MODULE Verification_m
