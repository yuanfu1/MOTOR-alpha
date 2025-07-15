MODULE ShallowWaterCommonAdjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  USE CalVerDer_adj_m, ONLY: CalVerDer_adj_t
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  USE parameters_m, ONLY: g
  USE CoriolisForce_m, ONLY: compute_coriolis_force
  IMPLICIT NONE

  TYPE :: ForwardVarsAdjoint_t
    REAL(r_kind), ALLOCATABLE :: adj_vorticity(:, :, :)
    REAL(r_kind), ALLOCATABLE :: adj_divergence(:, :, :)
    REAL(r_kind), ALLOCATABLE :: adj_height_star(:, :, :)
  CONTAINS
    PROCEDURE :: initializeAdj
    PROCEDURE :: storeAdj
    PROCEDURE :: retrieveAdj
  END TYPE ForwardVarsAdjoint_t

  TYPE(gzm_adj_t) :: gzmAdjOps
  TYPE(CalVerDer_adj_t) :: calVerOpsAdj
  TYPE(poissonSolver_adj_t) :: psAdjOps
  TYPE(SingleGrid_t) :: sg

CONTAINS

  ! Initialization for adjoint variables
  SUBROUTINE initializeAdj(this, num_steps, num_stages, grid_size)
    CLASS(ForwardVarsAdjoint_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: num_steps, num_stages, grid_size

    ! Allocate adjoint variables for vorticity, divergence, and height_star
    ALLOCATE (this%adj_vorticity(grid_size, num_steps, num_stages))
    ALLOCATE (this%adj_divergence(grid_size, num_steps, num_stages))
    ALLOCATE (this%adj_height_star(grid_size, num_steps, num_stages))
  END SUBROUTINE initializeAdj

  ! Store adjoint variables
  SUBROUTINE storeAdj(this, it, stage, adj_Y)
    CLASS(ForwardVarsAdjoint_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(IN) :: adj_Y

    ! Store the adjoint variables into adjoint variables arrays
    this%adj_vorticity(:, :, it) = RESHAPE(adj_Y%fields(1)%DATA, SHAPE(this%adj_vorticity(:, :, it)))
    this%adj_divergence(:, :, it) = RESHAPE(adj_Y%fields(2)%DATA, SHAPE(this%adj_divergence(:, :, it)))
    this%adj_height_star(:, :, it) = RESHAPE(adj_Y%fields(3)%DATA, SHAPE(this%adj_height_star(:, :, it)))
  END SUBROUTINE storeAdj

  ! Retrieve adjoint variables
  SUBROUTINE retrieveAdj(this, it, stage, adj_Y)
    CLASS(ForwardVarsAdjoint_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(INOUT) :: adj_Y

    ! Retrieve the stored adjoint variables
    adj_Y%fields(1)%DATA(:, :, it) = RESHAPE(this%adj_vorticity(:, stage, it), SHAPE(adj_Y%fields(1)%DATA(:, :, it)))
    adj_Y%fields(2)%DATA(:, :, it) = RESHAPE(this%adj_divergence(:, stage, it), SHAPE(adj_Y%fields(2)%DATA(:, :, it)))
    adj_Y%fields(3)%DATA(:, :, it) = RESHAPE(this%adj_height_star(:, stage, it), SHAPE(adj_Y%fields(3)%DATA(:, :, it)))
  END SUBROUTINE retrieveAdj

  ! Adjoint of kinetic energy computation
  SUBROUTINE adjoint_kinetic_energy(adj_streamfunc, adj_velocity_pot, adj_K, adj_u, adj_v)
    REAL(r_kind), INTENT(INOUT) :: adj_streamfunc(:, :), adj_velocity_pot(:, :)
    REAL(r_kind), POINTER :: adj_K(:, :), adj_u(:, :), adj_v(:, :)

    REAL(r_kind), POINTER :: adj_grad_psi_lat(:, :), adj_grad_psi_lon(:, :)
    REAL(r_kind), POINTER :: adj_grad_chi_lat(:, :), adj_grad_chi_lon(:, :)

    ALLOCATE (adj_grad_psi_lat(SIZE(adj_streamfunc, 1), SIZE(adj_streamfunc, 2)))
    ALLOCATE (adj_grad_psi_lon(SIZE(adj_streamfunc, 1), SIZE(adj_streamfunc, 2)))
    ALLOCATE (adj_grad_chi_lat(SIZE(adj_velocity_pot, 1), SIZE(adj_velocity_pot, 2)))
    ALLOCATE (adj_grad_chi_lon(SIZE(adj_velocity_pot, 1), SIZE(adj_velocity_pot, 2)))

    ! Step 1: Compute adjoint gradients of ψ and χ
    CALL gzmAdjOps%Grad_Lat_AD(adj_streamfunc, adj_grad_psi_lat, adj_u)
    CALL gzmAdjOps%Grad_Lon_AD(adj_streamfunc, adj_grad_psi_lon, adj_u)
    CALL gzmAdjOps%Grad_Lat_AD(adj_velocity_pot, adj_grad_chi_lat, adj_v)
    CALL gzmAdjOps%Grad_Lon_AD(adj_velocity_pot, adj_grad_chi_lon, adj_v)

    ! Step 2: Compute adjoint of kinetic energy field and vertical second-order derivatives
    CALL calVerOpsAdj%SecondOrder_AD(adj_K, adj_grad_psi_lat, adj_grad_psi_lon)

    ! Deallocate temporary arrays
    DEALLOCATE (adj_grad_psi_lat, adj_grad_psi_lon, adj_grad_chi_lat, adj_grad_chi_lon)
  END SUBROUTINE adjoint_kinetic_energy

  ! Adjoint of RHS for vorticity equation
  SUBROUTINE adjoint_rhs_vorticity(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_rhs_vorticity, adj_value)
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :)
    REAL(r_kind), INTENT(OUT) :: adj_rhs_vorticity(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)  ! Add adjoint of the output

    REAL(r_kind), ALLOCATABLE :: adj_div_eta_grad_chi(:, :), adj_jacobian_eta_psi(:, :)

    ALLOCATE (adj_div_eta_grad_chi(SIZE(adj_vorticity, 1), SIZE(adj_vorticity, 2)))
    ALLOCATE (adj_jacobian_eta_psi(SIZE(adj_vorticity, 1), SIZE(adj_vorticity, 2)))

    ! Compute adjoint of the divergence and Jacobian terms
    CALL gzmAdjOps%Divergen_AD(adj_vorticity, adj_velocity_pot, adj_div_eta_grad_chi, adj_rhs_vorticity, adj_value)
    CALL gzmAdjOps%Jacobian_AD(adj_vorticity, adj_streamfunc, adj_jacobian_eta_psi, adj_rhs_vorticity, adj_value)

    ! Deallocate temporary arrays
    DEALLOCATE (adj_div_eta_grad_chi, adj_jacobian_eta_psi)
  END SUBROUTINE adjoint_rhs_vorticity

  ! Adjoint of RHS for divergence equation
  SUBROUTINE adjoint_rhs_divergence(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_K, adj_height_star, adj_rhs_divergence, adj_value)
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :), adj_K(:, :), adj_height_star(:, :)
    REAL(r_kind), INTENT(OUT) :: adj_rhs_divergence(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_div_eta_grad_psi(:, :), adj_jacobian_eta_chi(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_laplacian_term(:, :)

    ALLOCATE (adj_div_eta_grad_psi(SIZE(adj_vorticity, 1), SIZE(adj_vorticity, 2)))
    ALLOCATE (adj_jacobian_eta_chi(SIZE(adj_vorticity, 1), SIZE(adj_vorticity, 2)))
    ALLOCATE (adj_laplacian_term(SIZE(adj_K, 1), SIZE(adj_K, 2)))

    ! Compute adjoint of the divergence and Jacobian terms
    CALL gzmAdjOps%Divergen_AD(adj_vorticity, adj_streamfunc, adj_div_eta_grad_psi, adj_rhs_divergence, adj_value)
    CALL gzmAdjOps%Jacobian_AD(adj_vorticity, adj_velocity_pot, adj_jacobian_eta_chi, adj_rhs_divergence, adj_value)

    ! Compute adjoint of Laplacian
    CALL gzmAdjOps%Laplacia_AD(adj_K, adj_laplacian_term, adj_rhs_divergence)

    DEALLOCATE (adj_div_eta_grad_psi, adj_jacobian_eta_chi, adj_laplacian_term)
  END SUBROUTINE adjoint_rhs_divergence

  ! Adjoint of RHS for height equation
  SUBROUTINE adjoint_rhs_height(adj_vorticity, adj_velocity_pot, adj_streamfunc, adj_height_star, adj_rhs_height, adj_value)
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_velocity_pot(:, :), adj_streamfunc(:, :), adj_height_star(:, :)
    REAL(r_kind), INTENT(OUT) :: adj_rhs_height(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_div_h_star_grad_chi(:, :), adj_jacobian_h_star_psi(:, :)

    ALLOCATE (adj_div_h_star_grad_chi(SIZE(adj_height_star, 1), SIZE(adj_height_star, 2)))
    ALLOCATE (adj_jacobian_h_star_psi(SIZE(adj_height_star, 1), SIZE(adj_height_star, 2)))

    ! Compute adjoint of the divergence and Jacobian terms
    CALL gzmAdjOps%Divergen_AD(adj_height_star, adj_velocity_pot, adj_div_h_star_grad_chi, adj_rhs_height, adj_value)
    CALL gzmAdjOps%Jacobian_AD(adj_height_star, adj_streamfunc, adj_jacobian_h_star_psi, adj_rhs_height, adj_value)

    DEALLOCATE (adj_div_h_star_grad_chi, adj_jacobian_h_star_psi)
  END SUBROUTINE adjoint_rhs_height

  ! Adjoint of Poisson solver
  SUBROUTINE solve_poisson_adjoint(poissonSolverAdj, vLevel, gLevel, adj_vorticity, adj_divergence, adj_streamfunc, adj_velocity_pot, sg)
    USE kinds_m, ONLY: i_kind, r_kind
    USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    IMPLICIT NONE

    CLASS(poissonSolver_adj_t), INTENT(INOUT) :: poissonSolverAdj
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(INOUT) :: adj_vorticity(:, :), adj_divergence(:, :)
    REAL(r_kind), INTENT(OUT) :: adj_streamfunc(:, :), adj_velocity_pot(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_rhs_psi(:, :), adj_rhs_chi(:, :)

    ALLOCATE (adj_rhs_psi(SIZE(adj_vorticity, 1), SIZE(adj_vorticity, 2)))
    ALLOCATE (adj_rhs_chi(SIZE(adj_divergence, 1), SIZE(adj_divergence, 2)))

    ! Solve adjoint Poisson equations for adjoint stream function and velocity potential
    CALL poissonSolverAdj%PoissonSol_sphere_adjoint(gLevel, vLevel, adj_vorticity, adj_streamfunc)
    CALL poissonSolverAdj%PoissonSol_sphere_adjoint(gLevel, vLevel, adj_divergence, adj_velocity_pot)

    DEALLOCATE (adj_rhs_psi, adj_rhs_chi)
  END SUBROUTINE solve_poisson_adjoint

  ! Adjoint of velocity computation
  SUBROUTINE compute_velocity_adjoint(adj_streamfunc, adj_velocity_pot, adj_u, adj_v)
    REAL(r_kind), INTENT(IN) :: adj_streamfunc(:, :), adj_velocity_pot(:, :)
    REAL(r_kind), INTENT(OUT) :: adj_u(:, :), adj_v(:, :)

    REAL(r_kind), ALLOCATABLE :: adj_grad_psi_lat(:, :), adj_grad_psi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_grad_chi_lat(:, :), adj_grad_chi_lon(:, :)

    ALLOCATE (adj_grad_psi_lat(SIZE(adj_streamfunc, 1), SIZE(adj_streamfunc, 2)))
    ALLOCATE (adj_grad_psi_lon(SIZE(adj_streamfunc, 1), SIZE(adj_streamfunc, 2)))
    ALLOCATE (adj_grad_chi_lat(SIZE(adj_velocity_pot, 1), SIZE(adj_velocity_pot, 2)))
    ALLOCATE (adj_grad_chi_lon(SIZE(adj_velocity_pot, 1), SIZE(adj_velocity_pot, 2)))

    ! Compute adjoint gradients of ψ and χ
    CALL gzmAdjOps%Grad_Lat_AD(adj_streamfunc, adj_grad_psi_lat, adj_u)
    CALL gzmAdjOps%Grad_Lon_AD(adj_streamfunc, adj_grad_psi_lon, adj_u)
    CALL gzmAdjOps%Grad_Lat_AD(adj_velocity_pot, adj_grad_chi_lat, adj_v)
    CALL gzmAdjOps%Grad_Lon_AD(adj_velocity_pot, adj_grad_chi_lon, adj_v)

    DEALLOCATE (adj_grad_psi_lat, adj_grad_psi_lon, adj_grad_chi_lat, adj_grad_chi_lon)
  END SUBROUTINE compute_velocity_adjoint

END MODULE ShallowWaterCommonAdjoint_m
