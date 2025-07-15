MODULE TimeIntegrationAB3Adjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE rhsAdjointBase_m, ONLY: rhsAdjointBase_t
  USE TimeIntegrationAdjointBase_m, ONLY: TimeIntegrationAdjointBase_t
  IMPLICIT NONE

  TYPE, EXTENDS(TimeIntegrationAdjointBase_t) :: TimeIntegrationAB3Adjoint_t
    REAL(r_kind) :: coeff(3)
  CONTAINS
    PROCEDURE, PUBLIC :: initialization_s => initializationAB3Adj_s
    PROCEDURE, PUBLIC :: TimeIntegrationAdjointBase_s => TimeIntegrationAB3Adj_s
  END TYPE TimeIntegrationAB3Adjoint_t

CONTAINS

  !> Initialize the adjoint AB3 time integration scheme
  SUBROUTINE initializationAB3Adj_s(this, rhsAdj, adj_X)
    CLASS(TimeIntegrationAB3Adjoint_t) :: this
    CLASS(rhsAdjointBase_t), TARGET :: rhsAdj
    TYPE(State_t), INTENT(IN) :: adj_X

    ! Initialize base class components
    this%rhsAdj => rhsAdj

    ! Initialize AB3 coefficients (same as forward)
    this%coeff = (/23.0_R_KIND / 12.0_R_KIND, 16.0_R_KIND / 12.0_R_KIND, 5.0_R_KIND / 12.0_R_KIND/)
  END SUBROUTINE initializationAB3Adj_s

  !> Perform adjoint time integration using AB3 scheme
  SUBROUTINE TimeIntegrationAB3Adj_s(this, it, dt, adj_X)
    CLASS(TimeIntegrationAB3Adjoint_t) :: this
    INTEGER(i_kind), INTENT(IN) :: it
    REAL(r_kind), INTENT(IN) :: dt
    TYPE(State_t), INTENT(INOUT) :: adj_X

    ! Local variables
    INTEGER(i_kind) :: i, j, num_vars
    REAL(r_kind), ALLOCATABLE :: adj_rhsSaved(:, :, :, :), rhsSaved(:, :, :, :)
    REAL(r_kind), ALLOCATABLE :: adj_yk(:, :, :)

    ! Number of variables in the state
    num_vars = UBOUND(adj_X%fields, 1)

    ! Allocate arrays for adjoint computations and forward-mode saved variables
    ALLOCATE (adj_rhsSaved(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 3))
    ALLOCATE (rhsSaved(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 3))
    ALLOCATE (adj_yk(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars))

    ! Retrieve the forward-saved variables (rhsSaved) for adjoint computation
    CALL retrieve_forward_variables(it, rhsSaved)

    ! Initialize adjoint variables (only accumulate for the last stage)
    adj_yk = 0.0_R_KIND

    ! Loop over the three stages of AB3 in reverse order
    DO i = 3, 1, -1
      ! Update adjoint variables for each stage
      IF (i == 3) THEN
        adj_yk(:, :, :) = adj_yk(:, :, :) + dt * this%coeff(1) * adj_rhsSaved(:, :, :, i)
      ELSE IF (i == 2) THEN
        adj_yk(:, :, :) = adj_yk(:, :, :) + dt * this%coeff(2) * adj_rhsSaved(:, :, :, i)
      ELSE IF (i == 1) THEN
        adj_yk(:, :, :) = adj_yk(:, :, :) + dt * this%coeff(3) * adj_rhsSaved(:, :, :, i)
      END IF

      ! Compute the adjoint RHS for this stage
      CALL this%rhsAdj%rightHandSideAdjoint(adj_X, it, adj_yk, adj_rhsSaved(:, :, :, i), rhsSaved(:, :, :, i))
    END DO

    ! Accumulate into adj_X
    DO j = 1, num_vars
      adj_X%fields(j)%DATA(:, :, it) = adj_X%fields(j)%DATA(:, :, it) + adj_yk(:, :, j)
    END DO

    ! Deallocate temporary arrays
    DEALLOCATE (adj_rhsSaved, rhsSaved, adj_yk)
  END SUBROUTINE TimeIntegrationAB3Adj_s

  !> Retrieve the forward-saved variables for the adjoint AB3 method
  SUBROUTINE retrieve_forward_variables(it, rhsSaved)
    INTEGER(i_kind), INTENT(IN) :: it
    REAL(r_kind), INTENT(OUT) :: rhsSaved(:, :, :, :)

    ! Retrieve the forward variables from storage (implementation depends on the storage mechanism)
    ! This subroutine assumes that forward rhsSaved values for each stage are saved during the forward pass
    ! Implement the logic to fetch those values based on the forward integration structure
    ! ...
  END SUBROUTINE retrieve_forward_variables

END MODULE TimeIntegrationAB3Adjoint_m
