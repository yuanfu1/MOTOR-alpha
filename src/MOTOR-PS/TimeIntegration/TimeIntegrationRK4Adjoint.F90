MODULE TimeIntegrationRK4Adjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE rhsAdjointBase_m, ONLY: rhsAdjointBase_t
  USE TimeIntegrationAdjointBase_m, ONLY: TimeIntegrationAdjointBase_t
  IMPLICIT NONE

  TYPE, EXTENDS(TimeIntegrationAdjointBase_t) :: TimeIntegrationRK4Adjoint_t
    REAL(r_kind) :: steps(4), coeff(4)
  CONTAINS
    PROCEDURE, PUBLIC :: initialization_s => initializationRK4Adj_s
    PROCEDURE, PUBLIC :: TimeIntegrationAdjointBase_s => TimeIntegrationRK4Adj_s
  END TYPE TimeIntegrationRK4Adjoint_t

CONTAINS

  !> Initialize the adjoint RK4 time integration scheme
  SUBROUTINE initializationRK4Adj_s(this, rhsAdj, adj_X)
    CLASS(TimeIntegrationRK4Adjoint_t) :: this
    CLASS(rhsAdjointBase_t), TARGET :: rhsAdj
    TYPE(State_t), INTENT(IN) :: adj_X

    ! Initialize base class components
    this%rhsAdj => rhsAdj

    ! Initialize RK4 coefficients (same as forward)
    this%steps = (/0.0_R_KIND, 0.5_R_KIND, 0.5_R_KIND, 1.0_R_KIND/)
    this%coeff = (/1.0_R_KIND / 6.0_R_KIND, 1.0_R_KIND / 3.0_R_KIND, 1.0_R_KIND / 3.0_R_KIND, 1.0_R_KIND / 6.0_R_KIND/)
  END SUBROUTINE initializationRK4Adj_s

  !> Perform adjoint time integration using RK4 scheme
  SUBROUTINE TimeIntegrationRK4Adj_s(this, it, dt, adj_X)
    CLASS(TimeIntegrationRK4Adjoint_t) :: this
    INTEGER(i_kind), INTENT(IN) :: it
    REAL(r_kind), INTENT(IN) :: dt
    TYPE(State_t), INTENT(INOUT) :: adj_X

    ! Local variables
    INTEGER(i_kind) :: i, j, num_vars
    INTEGER :: stage
    REAL(r_kind), ALLOCATABLE :: adj_rk4(:, :, :, :), adj_yk4(:, :, :, :)
    REAL(r_kind), ALLOCATABLE :: yk4(:, :, :, :), rk4(:, :, :, :)

    num_vars = UBOUND(adj_X%fields, 1)
    ALLOCATE (adj_rk4(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 4))
    ALLOCATE (adj_yk4(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 4))
    ALLOCATE (yk4(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 4))
    ALLOCATE (rk4(adj_X%sg%vLevel, adj_X%sg%num_cell, num_vars, 4))

    ! The adjoint integration proceeds backward in time
    ! Loop over RK4 stages in reverse order
    DO stage = 4, 1, -1
      ! Retrieve stored forward variables for this stage
      CALL retrieve_forward_variables(it, stage, yk4(:, :, :, stage), rk4(:, :, :, stage))

      ! Initialize adjoint variables for this stage
      IF (stage == 4) THEN
        adj_yk4(:, :, :, stage) = 0.0_R_KIND
      ELSE
        adj_yk4(:, :, :, stage) = adj_yk4(:, :, :, stage) + dt * this%coeff(stage + 1) * adj_rk4(:, :, :, stage + 1)
      END IF

      ! Compute the adjoint RHS for this stage
      CALL this%rhsAdj%rightHandSideAdjoint(adj_X, it, adj_yk4(:, :, :, stage), adj_rk4(:, :, :, stage), &
                                            yk4(:, :, :, stage), rk4(:, :, :, stage))

      ! Update adjoint state variables
      IF (stage > 1) THEN
        adj_yk4(:, :, :, stage - 1) = adj_yk4(:, :, :, stage - 1) + this%steps(stage) * adj_rk4(:, :, :, stage)
      ELSE
        ! Accumulate into adj_X
        DO j = 1, num_vars
          adj_X%fields(j)%DATA(:, :, it) = adj_X%fields(j)%DATA(:, :, it) + adj_yk4(:, :, j, stage)
        END DO
      END IF
    END DO

    ! Deallocate temporary arrays
    DEALLOCATE (adj_rk4, adj_yk4, yk4, rk4)
  END SUBROUTINE TimeIntegrationRK4Adj_s

END MODULE TimeIntegrationRK4Adjoint_m
