MODULE Objective_Function_m
  USE kinds_m, ONLY: r_kind
  IMPLICIT NONE

  ! This module provides subroutines for computing the objective function and its gradient,
  ! which are essential components in data assimilation systems. The objective function measures
  ! the difference between model predictions and observations, while the gradient provides the
  ! direction to adjust the model state to minimize this difference.

  ! References:
  ! 1. Errico, R. M. (1997). What is an adjoint model?. Bulletin of the American Meteorological Society, 78(11), 2577-2591.
  ! 2. Courtier, P., Thépaut, J. N., & Hollingsworth, A. (1994). A strategy for operational implementation of 4D-Var, using an incremental approach. Quarterly Journal of the Royal Meteorological Society, 120(519), 1367-1387.
  ! 3. Giering, R., & Kaminski, T. (1998). Recipes for adjoint code construction. ACM Transactions on Mathematical Software (TOMS), 24(4), 437-474.
  ! 4. Navon, I. M. (2009). Data assimilation for numerical weather prediction: a review. Data assimilation for atmospheric, oceanic and hydrologic applications, 1, 21-65.

CONTAINS

  ! Compute the objective function value
  ! x: state variable
  ! y_obs: observed variable
  ! H: observation operator matrix
  ! R_inv: inverse of observation error covariance matrix
  ! J: objective function value (output)
  SUBROUTINE Compute_Objective_Function(x, y_obs, H, R_inv, J)
    REAL(r_kind), INTENT(IN) :: x(:)
    REAL(r_kind), INTENT(IN) :: y_obs(:)
    REAL(r_kind), INTENT(IN) :: H(:, :)
    REAL(r_kind), INTENT(IN) :: R_inv(:, :)
    REAL(r_kind), INTENT(OUT) :: J
    REAL(r_kind), ALLOCATABLE :: y_model(:)
    REAL(r_kind), ALLOCATABLE :: innovation(:)
    REAL(r_kind), ALLOCATABLE :: temp(:)

    ! Allocate temporary arrays
    ALLOCATE (y_model(SIZE(y_obs)))
    ALLOCATE (innovation(SIZE(y_obs)))
    ALLOCATE (temp(SIZE(y_obs)))

    ! Compute the model equivalent of observations
    y_model = MATMUL(H, x)

    ! Compute the innovation (observation minus model equivalent)
    innovation = y_obs - y_model

    ! Compute the objective function value
    temp = MATMUL(R_inv, innovation)
    J = 0.5 * SUM(innovation * temp)

    ! Deallocate temporary arrays
    DEALLOCATE (y_model)
    DEALLOCATE (innovation)
    DEALLOCATE (temp)
  END SUBROUTINE Compute_Objective_Function

  ! Compute the gradient of the objective function
  ! x: state variable
  ! y_obs: observed variable
  ! H: observation operator matrix
  ! R_inv: inverse of observation error covariance matrix
  ! grad_J: gradient of the objective function (output)
  SUBROUTINE Compute_Gradient_Objective_Function(x, y_obs, H, R_inv, grad_J)
    REAL(r_kind), INTENT(IN) :: x(:)
    REAL(r_kind), INTENT(IN) :: y_obs(:)
    REAL(r_kind), INTENT(IN) :: H(:, :)
    REAL(r_kind), INTENT(IN) :: R_inv(:, :)
    REAL(r_kind), INTENT(OUT) :: grad_J(:)
    REAL(r_kind), ALLOCATABLE :: y_model(:)
    REAL(r_kind), ALLOCATABLE :: innovation(:)
    REAL(r_kind), ALLOCATABLE :: adjoint_innovation(:)

    ! Allocate temporary arrays
    ALLOCATE (y_model(SIZE(y_obs)))
    ALLOCATE (innovation(SIZE(y_obs)))
    ALLOCATE (adjoint_innovation(SIZE(y_obs)))

    ! Compute the model equivalent of observations
    y_model = MATMUL(H, x)

    ! Compute the innovation (observation minus model equivalent)
    innovation = y_obs - y_model

    ! Compute the adjoint of the innovation
    adjoint_innovation = MATMUL(TRANSPOSE(R_inv), innovation)

    ! Compute the gradient of the objective function
    grad_J = -MATMUL(TRANSPOSE(H), adjoint_innovation)

    ! Deallocate temporary arrays
    DEALLOCATE (y_model)
    DEALLOCATE (innovation)
    DEALLOCATE (adjoint_innovation)
  END SUBROUTINE Compute_Gradient_Objective_Function

END MODULE Objective_Function_m
