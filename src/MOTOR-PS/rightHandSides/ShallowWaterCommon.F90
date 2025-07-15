MODULE ShallowWaterCommon_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE parameters_m, ONLY: g
  USE CoriolisForce_m, ONLY: compute_coriolis_force
  IMPLICIT NONE

  TYPE :: ForwardVars_t
    REAL(r_kind), ALLOCATABLE :: vorticity(:, :, :)
    REAL(r_kind), ALLOCATABLE :: divergence(:, :, :)
    REAL(r_kind), ALLOCATABLE :: height_star(:, :, :)
  CONTAINS
    PROCEDURE :: initialize
    PROCEDURE :: store
    PROCEDURE :: retrieve
  END TYPE ForwardVars_t

  TYPE(gzm_t) :: gzmOps
  TYPE(CalVerDer_t) :: calVerOps
  TYPE(poissonSolver_t) :: psOps    ! Poisson solver object
  TYPE(SingleGrid_t) :: sg          ! SingleGrid object
  !TYPE(CoriolisParam_t) :: coriolisOps   ! Coriolis parameter operations

CONTAINS

  SUBROUTINE initialize(this, num_steps, num_stages, grid_size)
    CLASS(ForwardVars_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: num_steps, num_stages, grid_size

    ALLOCATE (this%vorticity(grid_size, num_steps, num_stages))
    ALLOCATE (this%divergence(grid_size, num_steps, num_stages))
    ALLOCATE (this%height_star(grid_size, num_steps, num_stages))

    ! Yuanfu Xie initializes gzmOps:
    PRINT*,'ForwardVars_t - SG: ',sg%num_cell
    !CALL gzmOps%initialize(sg)
  END SUBROUTINE initialize

  SUBROUTINE store(this, it, stage, Y)
    CLASS(ForwardVars_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(IN) :: Y

    ! Store the state variables into forwardVars
    this%vorticity(:, :, it) = RESHAPE(Y%fields(1)%DATA, SHAPE(this%vorticity(:, :, it)))
    this%divergence(:, :, it) = RESHAPE(Y%fields(2)%DATA, SHAPE(this%divergence(:, :, it)))
    this%height_star(:, :, it) = RESHAPE(Y%fields(3)%DATA, SHAPE(this%height_star(:, :, it)))
  END SUBROUTINE store

  SUBROUTINE retrieve(this, it, stage, Y)
    CLASS(ForwardVars_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(INOUT) :: Y
    INTEGER :: vLevel, num_cell, num_steps

    vLevel = Y%sg%vLevel
    num_cell = Y%sg%num_cell
    num_steps = SIZE(this%vorticity, 3)
    ALLOCATE (Y%fields(1)%DATA(vLevel, num_cell, num_steps))
    ALLOCATE (Y%fields(2)%DATA(vLevel, num_cell, num_steps))
    ALLOCATE (Y%fields(3)%DATA(vLevel, num_cell, num_steps))

    ! Retrieve the stored forward variables
    Y%fields(1)%DATA(:, :, it) = RESHAPE(this%vorticity(:, stage, it), SHAPE(Y%fields(1)%DATA(:, :, it)))
    Y%fields(2)%DATA(:, :, it) = RESHAPE(this%divergence(:, stage, it), SHAPE(Y%fields(2)%DATA(:, :, it)))
    Y%fields(3)%DATA(:, :, it) = RESHAPE(this%height_star(:, stage, it), SHAPE(Y%fields(3)%DATA(:, :, it)))

  END SUBROUTINE retrieve

  SUBROUTINE kinetic_energy(streamfunc, velocity_pot, K, u, v)
    !------------------------------------------------------------------
    ! This subroutine computes the kinetic energy based on the stream
    ! function (ψ, psi) and velocity potential (χ, chi), following the
    ! method given in the referenced literature.
    !
    ! The kinetic energy is given by:
    !
    !     K = 0.5 * (|∇ψ|^2 + |∇χ|^2) + J(ψ, χ)
    !
    ! where:
    ! - ψ (psi) is the stream function.
    ! - χ (chi) is the velocity potential.
    ! - ∇ψ, ∇χ are the gradients of ψ and χ respectively.
    ! - J(ψ, χ) is the Jacobian determinant of ψ and χ, which represents
    !   the non-linear interaction between these fields.
    !
    ! Additionally, it computes the wind components u and v:
    !
    !     u = ∂ψ/∂λ + ∂χ/∂φ
    !     v = -∂ψ/∂φ + ∂χ/∂λ
    !
    !------------------------------------------------------------------
    REAL(r_kind), INTENT(IN) :: streamfunc(:, :), velocity_pot(:, :)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: K(:, :)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: u(:, :), v(:, :)

    REAL(r_kind), ALLOCATABLE :: grad_psi_lat(:, :), grad_psi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: grad_chi_lat(:, :), grad_chi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: grad_psi_sq(:, :), grad_chi_sq(:, :)
    REAL(r_kind), ALLOCATABLE :: jacobian_psi_chi(:, :)
    REAL(r_kind), ALLOCATABLE :: kinetic_energy_field(:, :), kinetic_energy_2nd(:, :)

PRINT*,'Entering kinetic energy... ',SIZE(streamfunc, 1), SIZE(streamfunc, 2)
    ! Allocate memory for intermediate calculations
    ALLOCATE (grad_psi_lat(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_psi_lon(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_chi_lat(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (grad_chi_lon(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (grad_psi_sq(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_chi_sq(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (jacobian_psi_chi(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (kinetic_energy_field(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (kinetic_energy_2nd(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (u(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (v(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))

    ! Step 1: Compute the gradients of ψ (stream function)
    PRINT*,'Calling gradient 1...'
    CALL gzmOps%Gradient(streamfunc, grad_psi_lat, grad_psi_lon)

    ! Step 2: Compute the gradients of χ (velocity potential)
    PRINT*,'Calling gradient 2...'
    CALL gzmOps%Gradient(velocity_pot, grad_chi_lat, grad_chi_lon)

    ! Step 3: Compute the wind components u and v
    ! u = ∂ψ/∂λ + ∂χ/∂φ
    ! v = -∂ψ/∂φ + ∂χ/∂λ
    u = grad_psi_lon + grad_chi_lat
    v = -grad_psi_lat + grad_chi_lon

    ! Step 4: Compute the magnitude squared of the gradients
    ! |∇ψ|^2 = (∂ψ/∂φ)^2 + (∂ψ/∂λ)^2
    ! |∇χ|^2 = (∂χ/∂φ)^2 + (∂χ/∂λ)^2
    grad_psi_sq = grad_psi_lat**2 + grad_psi_lon**2
    grad_chi_sq = grad_chi_lat**2 + grad_chi_lon**2

    ! Step 5: Compute the Jacobian J(ψ, χ)
    PRINT*,'Calculating the Jacobian...'
    CALL gzmOps%Jacobian(streamfunc, velocity_pot, jacobian_psi_chi)

    ! Step 6: Compute the kinetic energy field
    kinetic_energy_field = 0.5_R_KIND * (grad_psi_sq + grad_chi_sq) + jacobian_psi_chi

    ! Step 7: Compute the vertical second-order derivative of kinetic energy
    PRINT*,'DD K /DD...'
    CALL calVerOps%SecondOrder(kinetic_energy_field, kinetic_energy_2nd)

    ! Step 8: Combine kinetic energy and its vertical derivative
    K = kinetic_energy_field + kinetic_energy_2nd

    ! Deallocate temporary arrays to free memory
    DEALLOCATE (grad_psi_lat, grad_psi_lon, grad_chi_lat, grad_chi_lon)
    DEALLOCATE (grad_psi_sq, grad_chi_sq)
    DEALLOCATE (jacobian_psi_chi)
    DEALLOCATE (kinetic_energy_field, kinetic_energy_2nd)
    DEALLOCATE (u, v)

  END SUBROUTINE kinetic_energy

  ! Compute RHS for the vorticity equation based on the reference formula
  SUBROUTINE compute_rhs_vorticity(vorticity, velocity_pot, streamfunc, rhs_vorticity)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :)
    REAL(r_kind), INTENT(OUT) :: rhs_vorticity(:, :)   ! Output: Right-hand side for vorticity equation

    REAL(r_kind), ALLOCATABLE :: div_eta_grad_chi(:, :), jacobian_eta_psi(:, :)
    REAL(r_kind), ALLOCATABLE :: vorticity_2nd(:, :)

    ! Allocate temporary arrays for intermediate calculations
    ALLOCATE (div_eta_grad_chi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (jacobian_eta_psi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (vorticity_2nd(SIZE(vorticity, 1), SIZE(vorticity, 2)))  ! Allocate memory for second-order derivative

    ! Compute the divergence of (vorticity * gradient of velocity potential), i.e., ∇ · (η ∇χ)
    CALL gzmOps%Divergen(vorticity, velocity_pot, div_eta_grad_chi)

    ! Compute the Jacobian of (vorticity, stream function), i.e., J(η, ψ)
    CALL gzmOps%Jacobian(vorticity, streamfunc, jacobian_eta_psi)

    ! Compute the second-order vertical derivative of vorticity
    CALL calVerOps%SecondOrder(vorticity, vorticity_2nd)

    ! The RHS of the vorticity equation is given by the sum of these terms
    rhs_vorticity = -div_eta_grad_chi + jacobian_eta_psi + vorticity_2nd

    ! Deallocate temporary arrays
    DEALLOCATE (div_eta_grad_chi, jacobian_eta_psi, vorticity_2nd)
  END SUBROUTINE compute_rhs_vorticity

  ! Compute RHS for the divergence equation based on the reference formula
  SUBROUTINE compute_rhs_divergence(vorticity, velocity_pot, streamfunc, K, height_star, rhs_divergence)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :), K(:, :), height_star(:, :)
    REAL(r_kind), INTENT(OUT) :: rhs_divergence(:, :)   ! Output: Right-hand side for divergence equation

    REAL(r_kind), ALLOCATABLE :: div_eta_grad_psi(:, :), jacobian_eta_chi(:, :)
    REAL(r_kind), ALLOCATABLE :: laplacian_term(:, :)
    REAL(r_kind), ALLOCATABLE :: divergence_2nd(:, :)    ! For second-order vertical derivative of divergence

    ! Allocate temporary arrays for intermediate calculations
    ALLOCATE (div_eta_grad_psi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (jacobian_eta_chi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (laplacian_term(SIZE(K, 1), SIZE(K, 2)))
    ALLOCATE (divergence_2nd(SIZE(vorticity, 1), SIZE(vorticity, 2)))  ! Allocate for second-order derivative

    ! Compute the divergence of (vorticity * gradient of stream function), i.e., ∇ · (η ∇ψ)
    CALL gzmOps%Divergen(vorticity, streamfunc, div_eta_grad_psi)

    ! Compute the Jacobian of (vorticity, velocity potential), i.e., J(η, χ)
    CALL gzmOps%Jacobian(vorticity, velocity_pot, jacobian_eta_chi)

    ! Compute Laplacian of (K + gh), i.e., ∇²(K + gh)
    CALL gzmOps%Laplacia(K + g * height_star, laplacian_term)

    ! Compute the second-order vertical derivative of divergence
    CALL calVerOps%SecondOrder(vorticity, divergence_2nd)

    ! The RHS of the divergence equation is given by the sum of these terms
    rhs_divergence = -div_eta_grad_psi - jacobian_eta_chi + laplacian_term + divergence_2nd

    ! Deallocate temporary arrays
    DEALLOCATE (div_eta_grad_psi, jacobian_eta_chi, laplacian_term, divergence_2nd)
  END SUBROUTINE compute_rhs_divergence

  ! Compute RHS for the height equation based on the reference formula
  SUBROUTINE compute_rhs_height(vorticity, velocity_pot, streamfunc, height_star, rhs_height)
    REAL(r_kind), INTENT(IN) :: vorticity(:, :), velocity_pot(:, :), streamfunc(:, :), height_star(:, :)
    REAL(r_kind), INTENT(OUT) :: rhs_height(:, :)   ! Output: Right-hand side for height equation
    REAL(r_kind), ALLOCATABLE :: div_h_star_grad_chi(:, :), jacobian_h_star_psi(:, :)
    REAL(r_kind), ALLOCATABLE :: height_star_2nd(:, :)   ! For vertical second-order derivative

    ! Allocate temporary arrays for intermediate calculations
    ALLOCATE (div_h_star_grad_chi(SIZE(height_star, 1), SIZE(height_star, 2)))
    ALLOCATE (jacobian_h_star_psi(SIZE(height_star, 1), SIZE(height_star, 2)))
    ALLOCATE (height_star_2nd(SIZE(height_star, 1), SIZE(height_star, 2)))   ! Allocate for second-order derivative

    ! Compute the divergence of (height_star * gradient of velocity potential), i.e., ∇ · (h* ∇χ)
    CALL gzmOps%Divergen(height_star, velocity_pot, div_h_star_grad_chi)

    ! Compute the Jacobian of (height_star, stream function), i.e., J(h*, ψ)
    CALL gzmOps%Jacobian(height_star, streamfunc, jacobian_h_star_psi)

    ! Compute the second-order vertical derivative of height_star
    CALL calVerOps%SecondOrder(height_star, height_star_2nd)

    ! The RHS of the height equation is given by the combination of these terms
    rhs_height = -div_h_star_grad_chi + jacobian_h_star_psi + height_star_2nd   ! Adding the vertical second-order derivative

    ! Deallocate temporary arrays
    DEALLOCATE (div_h_star_grad_chi, jacobian_h_star_psi, height_star_2nd)
  END SUBROUTINE compute_rhs_height

  SUBROUTINE solve_poisson(poissonSolver, vLevel, gLevel, vorticity, divergence, streamfunc, velocity_pot, sg)
    USE kinds_m, ONLY: i_kind, r_kind
    USE poissonSolver_m, ONLY: poissonSolver_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    IMPLICIT NONE

    ! Input/output parameters
    CLASS(poissonSolver_t), INTENT(INOUT) :: poissonSolver   ! Poisson solver object
    INTEGER(i_kind), INTENT(IN) :: vLevel                    ! Vertical level
    INTEGER(i_kind), INTENT(IN) :: gLevel                    ! Horizontal grid level
    TYPE(SingleGrid_t), INTENT(IN) :: sg                     ! SingleGrid object containing latitudes

    ! 2D Arrays for vorticity, divergence, stream function, and velocity potential
    REAL(r_kind), INTENT(IN) :: vorticity(:, :)            ! Vorticity input
    REAL(r_kind), INTENT(IN) :: divergence(:, :)           ! Divergence input
    REAL(r_kind), INTENT(OUT) :: streamfunc(:, :)             ! Output stream function (ψ)
    REAL(r_kind), INTENT(OUT) :: velocity_pot(:, :)           ! Output velocity potential (χ)

    REAL(r_kind), ALLOCATABLE :: rhs_psi(:, :), rhs_chi(:, :)
    REAL(r_kind), ALLOCATABLE :: coriolis_force(:, :)         ! Coriolis force array

    ! Step 0: Compute the Coriolis force based on grid latitude
    CALL compute_coriolis_force(coriolis_force, sg)  ! Retrieve the computed Coriolis parameter array

    ! Allocate RHS arrays for the Poisson equations
    ALLOCATE (rhs_psi(SIZE(vorticity, 1), SIZE(vorticity, 2)))
    ALLOCATE (rhs_chi(SIZE(divergence, 1), SIZE(divergence, 2)))

    ! Step 1: Solve for the stream function (ψ) using the vorticity field (η - coriolis_force)
    rhs_psi = vorticity - coriolis_force  ! RHS for stream function equation
    CALL poissonSolver%PoissonSol_sphere(gLevel, vLevel, rhs_psi, streamfunc)

    ! Step 2: Solve for the velocity potential (χ) using the divergence field (δ)
    rhs_chi = divergence  ! RHS for velocity potential equation
    CALL poissonSolver%PoissonSol_sphere(gLevel, vLevel, rhs_chi, velocity_pot)

    ! Deallocate RHS arrays
    DEALLOCATE (rhs_psi, rhs_chi)
    DEALLOCATE (coriolis_force)

  END SUBROUTINE solve_poisson

  SUBROUTINE compute_velocity(streamfunc, velocity_pot, u, v)
    REAL(r_kind), INTENT(IN) :: streamfunc(:, :), velocity_pot(:, :)
    REAL(r_kind), INTENT(OUT) :: u(:, :), v(:, :)

    REAL(r_kind), ALLOCATABLE :: grad_psi_lat(:, :), grad_psi_lon(:, :)
    REAL(r_kind), ALLOCATABLE :: grad_chi_lat(:, :), grad_chi_lon(:, :)

    ! 分配数组
    ALLOCATE (grad_psi_lat(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_psi_lon(SIZE(streamfunc, 1), SIZE(streamfunc, 2)))
    ALLOCATE (grad_chi_lat(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))
    ALLOCATE (grad_chi_lon(SIZE(velocity_pot, 1), SIZE(velocity_pot, 2)))

    CALL gzmOps%Gradient(streamfunc, grad_psi_lat, grad_psi_lon)
    CALL gzmOps%Gradient(velocity_pot, grad_chi_lat, grad_chi_lon)

    u = grad_psi_lon + grad_chi_lat
    v = -grad_psi_lat + grad_chi_lon

    DEALLOCATE (grad_psi_lat, grad_psi_lon, grad_chi_lat, grad_chi_lon)
  END SUBROUTINE compute_velocity

END MODULE ShallowWaterCommon_m
