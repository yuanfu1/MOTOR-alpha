MODULE ShallowWaterAdjointModel_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE gzm_m, ONLY: gzm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE CalVerDer_adj_m, ONLY: CalVerDer_adj_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  USE CoriolisForce_m, ONLY: compute_coriolis_force_adjoint
  USE TimeIntegrationRK4Adjoint_m, ONLY: TimeIntegrationRK4Adjoint_t
  USE rhsAdjointBase_m, ONLY: rhsAdjointBase_t
  IMPLICIT NONE

  ! Global variables or objects
  TYPE(gzm_t) :: gzmOps
  TYPE(gzm_adj_t) :: gzmAdjOps
  TYPE(CalVerDer_t) :: calVerOps
  TYPE(CalVerDer_adj_t) :: calVerAdjOps
  TYPE(poissonSolver_t) :: psOps
  TYPE(poissonSolver_adj_t) :: psAdjOps
  TYPE(SingleGrid_t) :: sg

CONTAINS

  ! Define the adjoint RHS of the shallow water equations by extending rhsAdjointBase_t
  TYPE, EXTENDS(rhsAdjointBase_t) :: ShallowWaterAdjointRHS_t
    TYPE(poissonSolver_adj_t) :: psAdjOps
    TYPE(gzm_adj_t) :: gzmAdjOps
    TYPE(CalVerDer_adj_t) :: calVerAdjOps
  CONTAINS
    PROCEDURE :: rightHandSideAdjoint
    PROCEDURE :: initialize_s => initialize_sw_adj_rhs
  END TYPE ShallowWaterAdjointRHS_t

  ! Initialize the ShallowWaterAdjointRHS_t object
  SUBROUTINE initialize_sw_adj_rhs(this)
    CLASS(ShallowWaterAdjointRHS_t), INTENT(INOUT) :: this
    ! Initialize base class components
    CALL this%rhsAdjointBase_t%initialize_s()
    ! Set the number of equations and variable indices
    this%num_eqns = 3  ! Vorticity, Divergence, Height
    ALLOCATE (this%variableIdx(3))
    this%variableIdx = (/1, 2, 3/)
  END SUBROUTINE initialize_sw_adj_rhs

  ! Implement the rightHandSideAdjoint procedure
!  SUBROUTINE rightHandSideAdjoint(this, adj_X, it, adj_current, adj_forward, yk4_stage, rk4_stage)
!   CLASS(ShallowWaterAdjointRHS_t), INTENT(IN) :: this
!   TYPE(State_t), INTENT(IN) :: adj_X
!   INTEGER(i_kind), INTENT(IN) :: it
!   REAL(r_kind), INTENT(IN) :: adj_current(adj_X%sg%vLevel, adj_X%sg%num_cell, UBOUND(adj_X%fields,1))
!   REAL(r_kind), INTENT(OUT) :: adj_forward(adj_X%sg%vLevel, adj_X%sg%num_cell, UBOUND(adj_X%fields,1))
!   REAL(r_kind), INTENT(IN) :: yk4_stage(adj_X%sg%vLevel, adj_X%sg%num_cell, UBOUND(adj_X%fields,1))
!   INTEGER(i_kind), INTENT(IN) :: rk4_stage

!     ! Extract adjoint variables from 'adj_current' array
!   REAL(r_kind), POINTER :: adj_vorticity(:,:), adj_divergence(:,:), adj_height_star(:,:)

!     ! Declare necessary variables for adjoint computations
!   REAL(r_kind), ALLOCATABLE :: vorticity(:,:), divergence(:,:), height_star(:,:)
!   REAL(r_kind), ALLOCATABLE :: adj_streamfunc(:,:), adj_velocity_pot(:,:), adj_K(:,:)
!   REAL(r_kind), ALLOCATABLE :: streamfunc(:,:), velocity_pot(:,:), K(:,:), u(:,:), v(:,:)
!   REAL(r_kind), ALLOCATABLE :: adj_u(:,:), adj_v(:,:)

!   adj_vorticity => adj_current(:,:,1)
!   adj_divergence => adj_current(:,:,2)
!   adj_height_star => adj_current(:,:,3)

!   ! Initialize adjoint RHS arrays
!   adj_forward(:,:,1) = 0.0_r_kind  ! Adjoint RHS for vorticity
!   adj_forward(:,:,2) = 0.0_r_kind  ! Adjoint RHS for divergence
!   adj_forward(:,:,3) = 0.0_r_kind  ! Adjoint RHS for height_star

!   ! Allocate arrays
!   ALLOCATE(vorticity(adj_X%sg%vLevel, adj_X%sg%num_cell))
!   ALLOCATE(divergence(adj_X%sg%vLevel, adj_X%sg%num_cell))
!   ALLOCATE(height_star(adj_X%sg%vLevel, adj_X%sg%num_cell))

!   ALLOCATE(adj_streamfunc(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(adj_velocity_pot(SIZE(divergence,1), SIZE(divergence,2)))
!   ALLOCATE(adj_K(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(adj_u(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(adj_v(SIZE(vorticity,1), SIZE(vorticity,2)))

!   ALLOCATE(streamfunc(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(velocity_pot(SIZE(divergence,1), SIZE(divergence,2)))
!   ALLOCATE(K(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(u(SIZE(vorticity,1), SIZE(vorticity,2)))
!   ALLOCATE(v(SIZE(vorticity,1), SIZE(vorticity,2)))

!   ! Initialize adjoint variables
!   adj_streamfunc = 0.0_r_kind
!   adj_velocity_pot = 0.0_r_kind
!   adj_K = 0.0_r_kind
!   adj_u = 0.0_r_kind
!   adj_v = 0.0_r_kind

!   ! Retrieve forward variables for this time step and stage
!   CALL retrieve_forward_variables(it, rk4_stage, vorticity, divergence, height_star, K, streamfunc, velocity_pot)

!   ! Compute kinetic energy and velocities (if not stored, recompute)
!   ! If K, streamfunc, velocity_pot are retrieved, you can skip recomputing them
!   ! If necessary, recompute K, u, v
!   CALL kinetic_energy(streamfunc, velocity_pot, K, u, v)

!   ! Adjoint of divergence equation RHS
!   CALL adjoint_rhs_divergence(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_K, adj_height_star, &
!                               adj_forward(:,:,2), vorticity, velocity_pot, streamfunc, K, height_star)

!   ! Adjoint of vorticity equation RHS
!   CALL adjoint_rhs_vorticity(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_forward(:,:,1), &
!                              vorticity, velocity_pot, streamfunc)

!   ! Adjoint of height equation RHS
!   CALL adjoint_rhs_height(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_height_star, adj_forward(:,:,3), &
!                           vorticity, velocity_pot, streamfunc, height_star)

!   ! Adjoint of kinetic energy computation
!   CALL adjoint_kinetic_energy(adj_streamfunc, adj_velocity_pot, adj_K, adj_u, adj_v, streamfunc, velocity_pot, u, v)

!   ! Solve adjoint Poisson equations
!   CALL psAdjOps%PoissonSol_sphere_adjoint(1, 1, adj_streamfunc, adj_velocity_pot, adj_vorticity, adj_divergence, sg)

!   ! Deallocate arrays
!   DEALLOCATE(vorticity, divergence, height_star)
!   DEALLOCATE(adj_streamfunc, adj_velocity_pot, adj_K, adj_u, adj_v)
!   DEALLOCATE(streamfunc, velocity_pot, K, u, v)
! END SUBROUTINE rightHandSideAdjoint

  ! Run the adjoint model using the adjoint RK4 time integration scheme
  ! SUBROUTINE run_adjoint_model(adj_X, adj_Y, dt, num_steps)
  !   TYPE(State_t), INTENT(INOUT) :: adj_X   ! Input: gradient of the cost function w.r.t. state variables
  !   TYPE(State_t), INTENT(OUT) :: adj_Y     ! Output: gradient of the cost function w.r.t. control variables
  !   REAL(r_kind), INTENT(IN) :: dt
  !   INTEGER, INTENT(IN) :: num_steps

  !   TYPE(ShallowWaterAdjointRHS_t) :: adj_rhs
  !   TYPE(TimeIntegrationRK4Adjoint_t) :: rk4AdjIntegrator

  !       ! Time-stepping loop (backward in time)
  !   INTEGER :: it

  !   ! Initialize adjoint RHS object
  !   CALL adj_rhs%initialize_s()
  !   adj_rhs%psAdjOps = psAdjOps
  !   adj_rhs%gzmAdjOps = gzmAdjOps
  !   adj_rhs%calVerAdjOps = calVerAdjOps

  !   ! Initialize the adjoint time integrator
  !   CALL rk4AdjIntegrator%initialization_s(adj_rhs, adj_X)

  !   ! Initialize the adjoint state variable adj_Y
  !   adj_Y = adj_X

  !   DO it = num_steps, 1, -1
  !     CALL rk4AdjIntegrator%TimeIntegrationAdjointBase_s(it, dt, adj_Y)
  !   END DO
  ! END SUBROUTINE run_adjoint_model

  ! Implement the adjoint RHS subroutines

  ! Adjoint subroutine for kinetic energy computation
  SUBROUTINE adjoint_kinetic_energy(adj_streamfunc, adj_velocity_pot, adj_K, adj_u, adj_v, streamfunc, velocity_pot, u, v)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: adj_streamfunc(:, :), adj_velocity_pot(:, :), adj_K(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_u(:, :), adj_v(:, :)
    REAL(r_kind), INTENT(IN) :: streamfunc(:, :), velocity_pot(:, :), u(:, :), v(:, :)

    REAL(r_kind), ALLOCATABLE :: grad_psi_lat(:, :), grad_psi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: grad_chi_lat(:, :), grad_chi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_grad_psi_lat(:, :), adj_grad_psi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_grad_chi_lat(:, :), adj_grad_chi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_kinetic_energy_field(:, :)

    ! Allocate necessary arrays
    ALLOCATE (grad_psi_lat(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_psi_lon(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_chi_lat(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (grad_chi_lon(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (adj_grad_psi_lat(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (adj_grad_psi_lon(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (adj_grad_chi_lat(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (adj_grad_chi_lon(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (adj_kinetic_energy_field(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))

    ! Initialize adjoint variables
    adj_grad_psi_lat = 0.0_R_KIND
    adj_grad_psi_lon = 0.0_R_KIND
    adj_grad_chi_lat = 0.0_R_KIND
    adj_grad_chi_lon = 0.0_R_KIND
    adj_kinetic_energy_field = adj_K

    ! Forward calculations (to get necessary intermediate variables)
    CALL gzmOps%Gradient(streamfunc, grad_psi_lat, grad_psi_lon)
    CALL gzmOps%Gradient(velocity_pot, grad_chi_lat, grad_chi_lon)

    ! Compute kinetic energy field
    ! kinetic_energy_field = 0.5 * ((grad_psi_lat)^2 + (grad_psi_lon)^2 + (grad_chi_lat)^2 + (grad_chi_lon)^2)
    ! u = grad_psi_lon + grad_chi_lat
    ! v = -grad_psi_lat + grad_chi_lon

    ! Adjoint calculations

    ! Adjoint of kinetic energy field
    adj_grad_psi_lat = adj_grad_psi_lat + adj_kinetic_energy_field * grad_psi_lat
    adj_grad_psi_lon = adj_grad_psi_lon + adj_kinetic_energy_field * grad_psi_lon
    adj_grad_chi_lat = adj_grad_chi_lat + adj_kinetic_energy_field * grad_chi_lat
    adj_grad_chi_lon = adj_grad_chi_lon + adj_kinetic_energy_field * grad_chi_lon

    ! Adjoint of velocities
    adj_grad_psi_lon = adj_grad_psi_lon + adj_u
    adj_grad_chi_lat = adj_grad_chi_lat + adj_u
    adj_grad_psi_lat = adj_grad_psi_lat - adj_v
    adj_grad_chi_lon = adj_grad_chi_lon + adj_v

    ! Adjoint of gradients
    CALL gzmAdjOps%Gradient_adjoint(streamfunc, adj_streamfunc, adj_grad_psi_lat, adj_grad_psi_lon)
    CALL gzmAdjOps%Gradient_adjoint(velocity_pot, adj_velocity_pot, adj_grad_chi_lat, adj_grad_chi_lon)

    ! Deallocate temporary arrays
    DEALLOCATE (grad_psi_lat, grad_psi_lon, grad_chi_lat, grad_chi_lon)
    DEALLOCATE (adj_grad_psi_lat, adj_grad_psi_lon, adj_grad_chi_lat, adj_grad_chi_lon)
    DEALLOCATE (adj_kinetic_energy_field)
  END SUBROUTINE adjoint_kinetic_energy

  ! Adjoint subroutine for vorticity equation RHS
  SUBROUTINE adjoint_rhs_vorticity(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_rhs_vorticity, vorticity, velocity_pot, streamfunc)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :)
    REAL(r_kind), INTENT(IN) :: adj_rhs_vorticity(:, :)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_div_eta_grad_chi(:, :), adj_jacobian_eta_psi(:, :), adj_vorticity_2nd(:, :)
    REAL(r_kind), ALLOCATABLE :: div_eta_grad_chi(:, :), jacobian_eta_psi(:, :), vorticity_2nd(:, :)

    ! Allocate necessary arrays
    ALLOCATE (adj_div_eta_grad_chi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (adj_jacobian_eta_psi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (adj_vorticity_2nd(SIZE(vorticity, 1), SIZE(vorticity, 2)))

    ! Initialize adjoint variables
    adj_div_eta_grad_chi = -adj_rhs_vorticity
    adj_jacobian_eta_psi = adj_rhs_vorticity
    adj_vorticity_2nd = adj_rhs_vorticity

    ! Forward calculations
    CALL gzmOps%Divergen(vorticity, velocity_pot, div_eta_grad_chi)
    CALL gzmOps%Jacobian(vorticity, streamfunc, jacobian_eta_psi)
    CALL calVerOps%SecondOrder(vorticity, vorticity_2nd)

    ! Adjoint of second-order derivative
    CALL calVerAdjOps%SecondOrder_adjoint(vorticity, adj_vorticity, adj_vorticity_2nd)

    ! Adjoint of Jacobian
    CALL gzmAdjOps%Jacobian_adjoint(vorticity, streamfunc, adj_vorticity, adj_streamfunc, adj_jacobian_eta_psi)

    ! Adjoint of divergence
    CALL gzmAdjOps%Divergen_adjoint(vorticity, velocity_pot, adj_vorticity, adj_velocity_pot, adj_div_eta_grad_chi)

    ! Deallocate temporary arrays
    DEALLOCATE (adj_div_eta_grad_chi, adj_jacobian_eta_psi, adj_vorticity_2nd)
  END SUBROUTINE adjoint_rhs_vorticity

  ! Adjoint subroutine for divergence equation RHS
  SUBROUTINE adjoint_rhs_divergence(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_K, adj_height_star, adj_rhs_divergence, &
                                    vorticity, velocity_pot, streamfunc, K, height_star)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :), adj_K(:, :), adj_height_star(:, :)
    REAL(r_kind), INTENT(IN) :: adj_rhs_divergence(:, :)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :), K(:, :), height_star(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_div_eta_grad_psi(:, :), adj_jacobian_eta_chi(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_laplacian_term(:, :), adj_divergence_2nd(:, :)
    REAL(r_kind), ALLOCATABLE :: temp_adj(:, :)

    ! Allocate necessary arrays
    ALLOCATE (adj_div_eta_grad_psi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (adj_jacobian_eta_chi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (adj_laplacian_term(SIZE(K, 1), SIZE(K, 2)))
    ALLOCATE (adj_divergence_2nd(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (temp_adj(SIZE(K, 1), SIZE(K, 2)))

    ! Initialize adjoint variables
    adj_div_eta_grad_psi = -adj_rhs_divergence
    adj_jacobian_eta_chi = -adj_rhs_divergence
    adj_laplacian_term = adj_rhs_divergence
    adj_divergence_2nd = adj_rhs_divergence

    ! Forward calculations
    CALL gzmOps%Divergen(vorticity, streamfunc, temp_adj)
    CALL gzmOps%Jacobian(vorticity, velocity_pot, temp_adj)
    CALL gzmOps%Laplacia(K + 9.81_R_KIND * height_star, temp_adj)
    CALL calVerOps%SecondOrder(divergence, temp_adj)

    ! Adjoint of second-order derivative
    CALL calVerAdjOps%SecondOrder_adjoint(divergence, adj_divergence, adj_divergence_2nd)

    ! Adjoint of Laplacian term
    CALL gzmAdjOps%Laplacia_adjoint(K + 9.81_R_KIND * height_star, temp_adj, adj_laplacian_term)
    adj_K = adj_K + temp_adj
    adj_height_star = adj_height_star + temp_adj * 9.81_R_KIND

    ! Adjoint of Jacobian
    CALL gzmAdjOps%Jacobian_adjoint(vorticity, velocity_pot, adj_vorticity, adj_velocity_pot, adj_jacobian_eta_chi)

    ! Adjoint of divergence
    CALL gzmAdjOps%Divergen_adjoint(vorticity, streamfunc, adj_vorticity, adj_streamfunc, adj_div_eta_grad_psi)

    ! Deallocate temporary arrays
    DEALLOCATE (adj_div_eta_grad_psi, adj_jacobian_eta_chi, adj_laplacian_term, adj_divergence_2nd, temp_adj)
  END SUBROUTINE adjoint_rhs_divergence

  ! Adjoint subroutine for height equation RHS
  SUBROUTINE adjoint_rhs_height(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_height_star, adj_rhs_height, &
                                vorticity, velocity_pot, streamfunc, height_star)
    IMPLICIT NONE
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :), adj_height_star(:, :)
    REAL(r_kind), INTENT(IN) :: adj_rhs_height(:, :)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :), height_star(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_div_h_star_grad_chi(:, :), adj_jacobian_h_star_psi(:, :), adj_height_star_2nd(:, :)

    ! Allocate necessary arrays
    ALLOCATE (adj_div_h_star_grad_chi(SIZE(height_star, 1), SIZE(height_star, 2)))
    ALLOCATE (adj_jacobian_h_star_psi(SIZE(height_star, 1), SIZE(height_star, 2)))
    ALLOCATE (adj_height_star_2nd(SIZE(height_star, 1), SIZE(height_star, 2)))

    ! Initialize adjoint variables
    adj_div_h_star_grad_chi = -adj_rhs_height
    adj_jacobian_h_star_psi = adj_rhs_height
    adj_height_star_2nd = adj_rhs_height

    ! Forward calculations
    CALL gzmOps%Divergen(height_star, velocity_pot, adj_div_h_star_grad_chi)
    CALL gzmOps%Jacobian(height_star, streamfunc, adj_jacobian_h_star_psi)
    CALL calVerOps%SecondOrder(height_star, adj_height_star_2nd)

    ! Adjoint of second-order derivative
    CALL calVerAdjOps%SecondOrder_adjoint(height_star, adj_height_star, adj_height_star_2nd)

    ! Adjoint of Jacobian
    CALL gzmAdjOps%Jacobian_adjoint(height_star, streamfunc, adj_height_star, adj_streamfunc, adj_jacobian_h_star_psi)

    ! Adjoint of divergence
    CALL gzmAdjOps%Divergen_adjoint(height_star, velocity_pot, adj_height_star, adj_velocity_pot, adj_div_h_star_grad_chi)

    ! Deallocate temporary arrays
    DEALLOCATE (adj_div_h_star_grad_chi, adj_jacobian_h_star_psi, adj_height_star_2nd)
  END SUBROUTINE adjoint_rhs_height

  SUBROUTINE retrieve_forward_variables(it, stage, vorticity, divergence, height_star, K, streamfunc, velocity_pot)
    INTEGER, INTENT(IN) :: it, stage
    REAL(r_kind), INTENT(OUT) :: vorticity(:, :), divergence(:, :), height_star(:, :)
    REAL(r_kind), INTENT(OUT) :: K(:, :), streamfunc(:, :), velocity_pot(:, :)

    CHARACTER(len=1024) :: filename
    TYPE(NCOutput_t) :: NCOutput
    INTEGER :: ncid, varid

    ! Construct filename based on time step and stage
    WRITE (filename, '(A,I0,A,I0,A)') 'forward_vars_it', it, '_stage', stage, '.nc'

    ! Initialize NCOutput with the filename
    CALL NCOutput%initialize(filename)

    ! Retrieve each variable from the NetCDF file
    CALL NCOutput%getVar('vorticity', vorticity)
    CALL NCOutput%getVar('divergence', divergence)
    CALL NCOutput%getVar('height_star', height_star)
    CALL NCOutput%getVar('K', K)
    CALL NCOutput%getVar('streamfunc', streamfunc)
    CALL NCOutput%getVar('velocity_pot', velocity_pot)

    ! Close the NetCDF file
    CALL NCOutput%CLOSE()

  END SUBROUTINE retrieve_forward_variables

END MODULE ShallowWaterAdjointModel_m
