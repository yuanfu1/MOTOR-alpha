MODULE rhsShallowWaterAdjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE rhsBaseAdjoint_m, ONLY: rhsBaseAdjoint_t
  USE State_m, ONLY: State_t
  USE ShallowWaterCommonAdjoint_m, ONLY: adjoint_rhs_vorticity, adjoint_rhs_divergence, adjoint_rhs_height, adjoint_kinetic_energy
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  IMPLICIT NONE

  TYPE, EXTENDS(rhsBaseAdjoint_t) :: rhsShallowWaterAdjoint_t
    TYPE(poissonSolver_adj_t) :: poissonSolverAdj
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s => initialize_shallowWaterAdjoint
    PROCEDURE, PUBLIC :: rightHandSideAdjoint => rhsShallowWaterAdjoint_s
  END TYPE rhsShallowWaterAdjoint_t

CONTAINS

  ! Initialize adjoint variable indices
  SUBROUTINE initialize_shallowWaterAdjoint(this)
    CLASS(rhsShallowWaterAdjoint_t), INTENT(INOUT) :: this

    ! Allocate space for adjoint variable indices (vorticity, divergence, height, velocity_u, velocity_v)
    ALLOCATE (this%variableIdx(5))
    this%variableIdx(1) = 1  ! Vorticity
    this%variableIdx(2) = 2  ! Divergence
    this%variableIdx(3) = 3  ! Height
    this%variableIdx(4) = 4  ! Velocity_u
    this%variableIdx(5) = 5  ! Velocity_v
  END SUBROUTINE initialize_shallowWaterAdjoint

  ! Adjoint RHS calculation for the shallow water model
  SUBROUTINE rhsShallowWaterAdjoint_s(this, X, it, adj_current, adj_forward)
    CLASS(rhsShallowWaterAdjoint_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: it
    REAL(r_kind), TARGET, INTENT(IN) :: adj_current(X%sg%vLevel, X%sg%num_cell, this%num_eqns)
    REAL(r_kind), INTENT(OUT) :: adj_forward(X%sg%vLevel, X%sg%num_cell, this%num_eqns)

    ! Local variables for adjoint vorticity, divergence, height, and velocities
    REAL(r_kind), POINTER :: adj_vorticity(:, :), adj_divergence(:, :), adj_height_star(:, :)
    REAL(r_kind), POINTER :: adj_velocity_u(:, :), adj_velocity_v(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_rhs_vorticity(:, :), adj_rhs_divergence(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_rhs_height(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_streamfunc(:, :), adj_velocity_pot(:, :), adj_K(:, :)
    REAL(r_kind), POINTER :: adj_rhs_velocity_u(:, :), adj_rhs_velocity_v(:, :)
    REAL(r_kind), POINTER :: adj_value(:, :)

    ! Step 1: Allocate memory for adjoint RHS calculations
    ALLOCATE (adj_rhs_vorticity(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_rhs_divergence(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_rhs_height(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_rhs_velocity_u(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_rhs_velocity_v(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_streamfunc(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_velocity_pot(SIZE(adj_current, 1), SIZE(adj_current, 2)))
    ALLOCATE (adj_K(SIZE(adj_current, 1), SIZE(adj_current, 2)))

    ! Extract adjoint variables from the adj_current state using pointer assignment
    adj_vorticity => adj_current(:, :, this%variableIdx(1))
    adj_divergence => adj_current(:, :, this%variableIdx(2))
    adj_height_star => adj_current(:, :, this%variableIdx(3))
    adj_velocity_u => adj_current(:, :, this%variableIdx(4))
    adj_velocity_v => adj_current(:, :, this%variableIdx(5))

    ! Step 2: Adjoint RHS calculations
    CALL adjoint_rhs_height(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_height_star, adj_rhs_height, adj_value)

    CALL adjoint_rhs_divergence(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_K, adj_height_star, adj_rhs_divergence, adj_value)

    CALL adjoint_rhs_vorticity(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_rhs_vorticity, adj_value)

    ! Step 3: Solve adjoint Poisson equations for streamfunction and velocity potential
    CALL this%poissonSolverAdj%PoissonSol_sphere_adjoint(X%sg%gLevel, X%sg%vLevel, adj_vorticity, adj_divergence, adj_streamfunc, adj_velocity_pot)

    ! Step 4: Compute adjoint kinetic energy RHS
    CALL adjoint_kinetic_energy(adj_streamfunc, adj_velocity_pot, adj_K, adj_rhs_velocity_u, adj_rhs_velocity_v)

    ! Step 5: Assign computed adjoint RHS values to adj_forward array
    adj_forward(:, :, this%variableIdx(1)) = adj_rhs_vorticity
    adj_forward(:, :, this%variableIdx(2)) = adj_rhs_divergence
    adj_forward(:, :, this%variableIdx(3)) = adj_rhs_height
    adj_forward(:, :, this%variableIdx(4)) = adj_rhs_velocity_u
    adj_forward(:, :, this%variableIdx(5)) = adj_rhs_velocity_v

    ! Deallocate temporary arrays
    DEALLOCATE (adj_rhs_vorticity, adj_rhs_divergence, adj_rhs_height, adj_streamfunc, adj_velocity_pot, adj_K, adj_rhs_velocity_u, adj_rhs_velocity_v)
  END SUBROUTINE rhsShallowWaterAdjoint_s

END MODULE rhsShallowWaterAdjoint_m
