SUBROUTINE Compute_Gradient_Objective_Function(x, y_obs, H, R_inv, grad_f)
  USE kinds_m, ONLY: r_kind
  REAL(r_kind), INTENT(IN) :: x(:), y_obs(:), H(:, :), R_inv(:, :)
  REAL(r_kind), INTENT(OUT) :: grad_f(:)

  ! Implementation for Compute_Gradient_Objective_Function
  grad_f = MATMUL(H, x - y_obs)  ! Placeholder for actual implementation
END SUBROUTINE Compute_Gradient_Objective_Function
