MODULE Verification_m
  USE kinds_m, ONLY: r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Objective_Function_m, ONLY: Compute_Objective_Function, Compute_Gradient_Objective_Function
  USE gzm_m, ONLY: gzm_t
  USE gzm_tlm_m, ONLY: gzm_tlm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  IMPLICIT NONE

CONTAINS

  !-------------------------------------------------------------------
  ! Subroutine: Verify_TLM_AD
  ! Purpose: Verify the consistency between Tangent Linear Model (TLM) and Adjoint Model (AD).
  ! Inputs:
  !   - tlm_output: Tangent linear model output
  !   - adj_output: Adjoint model output
  !   - d_opr_right: Perturbation in right operand for TLM
  !   - sg: SingleGrid object
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_TLM_AD(tlm_output, adj_output, tlm_perturbation, adj_perturbation, sg, op_name)
    REAL(r_kind), INTENT(IN) :: tlm_output(:, :), adj_output(:, :)
    REAL(r_kind), INTENT(IN) :: tlm_perturbation(:, :), adj_perturbation(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind) :: inner_product_tl, inner_product_ad, tolerance
    INTEGER :: k, j

    ! Initialize the inner products
    inner_product_tl = 0.0_R_KIND
    inner_product_ad = 0.0_R_KIND
    tolerance = 1E-14_R_KIND

    ! Calculate the inner product for the TLM output with the perturbation
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + tlm_output(k, j) * adj_perturbation(k, j)
      END DO
    END DO

    ! Calculate the inner product for the adjoint output with the perturbation
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_ad = inner_product_ad + adj_output(k, j) * tlm_perturbation(k, j)
      END DO
    END DO

    ! Print the inner product results for debugging
    PRINT *, "Inner Product TL:", inner_product_tl
    PRINT *, "Inner Product AD:", inner_product_ad

    ! Write the results to a verification file
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

  !-------------------------------------------------------------------
  ! Subroutine: Verify_Adjoint_AutomaticDifferentiation
  ! Purpose: Verify the adjoint model using Automatic Differentiation (AD).
  ! Inputs:
  !   - x: State variable
  !   - y_obs: Observed variable
  !   - H: Observation operator matrix
  !   - R_inv: Inverse of observation error covariance matrix
  !   - adj_output: Adjoint model output
  !   - sg: SingleGrid object
  !   - op_name: Operation name (for output labeling)
  !-------------------------------------------------------------------
  SUBROUTINE Verify_Adjoint_AutomaticDifferentiation(x, y_obs, H, R_inv, adj_output, sg, op_name)
    REAL(r_kind), INTENT(IN) :: x(:), y_obs(:), H(:, :), R_inv(:, :), adj_output(:)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: op_name
    REAL(r_kind), ALLOCATABLE :: grad_f(:)
    REAL(r_kind) :: inner_product_ad, inner_product_grad

    ALLOCATE (grad_f(SIZE(x)))

    ! Compute gradient using Automatic Differentiation
    CALL Compute_Gradient_Objective_Function(x, y_obs, H, R_inv, grad_f)

    inner_product_ad = SUM(adj_output * grad_f)
    inner_product_grad = SUM(grad_f * grad_f)

    ! Debug: Print results
    PRINT *, "Inner Product (AD * Gradient): ", inner_product_ad
    PRINT *, "Inner Product (Gradient * Gradient): ", inner_product_grad

    ! Write results to verification file
    OPEN (unit=11, file='verification_results.txt', status='unknown', position='append')
    WRITE (11, *) 'Verifying Adjoint using Automatic Differentiation for ', TRIM(op_name)
    WRITE (11, *) 'Inner product (AD * Gradient): ', inner_product_ad
    WRITE (11, *) 'Inner product (Gradient * Gradient): ', inner_product_grad
    WRITE (11, *) 'Difference: ', ABS(inner_product_ad - inner_product_grad)
    IF (ABS(inner_product_ad - inner_product_grad) < 1E-14_R_KIND) THEN
      WRITE (11, *) 'Adjoint Automatic Differentiation ', TRIM(op_name), ' verification passed.'
    ELSE
      WRITE (11, *) 'Adjoint Automatic Differentiation ', TRIM(op_name), ' verification failed.'
    END IF
    CLOSE (11)

    DEALLOCATE (grad_f)
  END SUBROUTINE Verify_Adjoint_AutomaticDifferentiation

END MODULE Verification_m
