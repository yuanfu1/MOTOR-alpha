MODULE psMatrix_adjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE singleGrid_m, ONLY: singleGrid_t
  USE psMatrix_m, ONLY: psMatrix_t
  USE gzm_adj_m, ONLY: gzm_adj_t  ! Include the gzm_adj_m module

  IMPLICIT NONE

  TYPE, EXTENDS(psMatrix_t) :: psMatrix_adj_t
    ! Adjoint variables
    REAL(r_kind), ALLOCATABLE :: adj_row_elem(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_jacobian(:, :)  ! Adjoint Jacobian for nonlinear problems
    LOGICAL :: is_nonlinear = .FALSE.  ! Flag to indicate if the problem is nonlinear
    TYPE(gzm_adj_t) :: gzmOps  ! Define the gzm_adj_t instance for adjoint operations
  CONTAINS
    FINAL :: destructor_adj

    PROCEDURE, PUBLIC :: loadMatrix_adj
    PROCEDURE, PUBLIC :: CalResidus_adj
    PROCEDURE, PUBLIC :: CalJacobian_adj   ! New procedure for nonlinear problems
  END TYPE psMatrix_adj_t

  INTERFACE psMatrix_adj_t
    MODULE PROCEDURE :: constructor_adj
  END INTERFACE psMatrix_adj_t

CONTAINS

  !------------------------------------------------------------------------------
  ! Constructor for psMatrix_adj_t
  !------------------------------------------------------------------------------
  FUNCTION constructor_adj(max_band, sg, nonlinear_flag) RESULT(this)
    IMPLICIT NONE
    TYPE(psMatrix_adj_t) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(IN) :: max_band
    LOGICAL, INTENT(IN), OPTIONAL :: nonlinear_flag  ! Optional flag for nonlinearity

    ! Set the nonlinearity flag if provided
    IF (PRESENT(nonlinear_flag)) THEN
      this%is_nonlinear = nonlinear_flag
    END IF

    ! Call the base constructor
    this%num_eqns = sg%num_icell
    this%max_band = max_band

    ALLOCATE (this%num_elem(this%num_eqns))
    ALLOCATE (this%idx_cols(max_band, this%num_eqns))
    ALLOCATE (this%row_elem(max_band, this%num_eqns))

    ! Initialize adjoint variables
    ALLOCATE (this%adj_row_elem(max_band, this%num_eqns))
    ALLOCATE (this%adj_jacobian(max_band, this%num_eqns))  ! Jacobian for nonlinear adjoint
    this%adj_row_elem = 0.0_R_KIND
    this%adj_jacobian = 0.0_R_KIND

    ! Initialize the gzmOps object for adjoint calculations
    CALL this%gzmOps%constructor(sg)

    ! Call the forward loadMatrix
    CALL this%loadMatrix_adj(sg)

  END FUNCTION constructor_adj

  !------------------------------------------------------------------------------
  ! Adjoint destructor for psMatrix_adj_t
  !------------------------------------------------------------------------------
  IMPURE ELEMENTAL SUBROUTINE destructor_adj(this)
    TYPE(psMatrix_adj_t), INTENT(INOUT) :: this

    ! Deallocate adjoint variables
    IF (ALLOCATED(this%adj_row_elem)) DEALLOCATE (this%adj_row_elem)
    IF (ALLOCATED(this%adj_jacobian)) DEALLOCATE (this%adj_jacobian)

  END SUBROUTINE destructor_adj

  !------------------------------------------------------------------------------
  ! Adjoint subroutine for loadMatrix
  !------------------------------------------------------------------------------
  SUBROUTINE loadMatrix_adj(this, sg)
    CLASS(psMatrix_adj_t), INTENT(INOUT) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg

    ! Reset adjoint variables associated with row_elem and Jacobian
    this%adj_row_elem = 0.0_R_KIND
    this%adj_jacobian = 0.0_R_KIND

    ! Logic from the forward mode can be reused here if needed
    CALL this%loadMatrix(sg)

  END SUBROUTINE loadMatrix_adj

  !------------------------------------------------------------------------------
  ! Adjoint subroutine for CalResidus (linear and nonlinear)
  !------------------------------------------------------------------------------
  SUBROUTINE CalResidus_adj(this, vlevel, rhs, sol, rsd, sg, adj_rhs, adj_sol, adj_rsd)
    CLASS(psMatrix_adj_t), INTENT(INOUT) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(IN) :: vlevel ! Number of vertical levels
    REAL(r_kind), INTENT(IN) :: rhs(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(IN) :: sol(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(IN) :: rsd(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(INOUT) :: adj_rhs(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(INOUT) :: adj_sol(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(INOUT) :: adj_rsd(vlevel, sg%num_cell)

    INTEGER(i_kind) :: i, j, col_idx
    REAL(r_kind) :: mat_elem

    ! Ensure adj_row_elem and adj_jacobian are allocated
    IF (.NOT. ALLOCATED(this%adj_row_elem)) THEN
      ALLOCATE (this%adj_row_elem(SIZE(this%row_elem, 1), SIZE(this%row_elem, 2)))
      this%adj_row_elem = 0.0_R_KIND
    END IF

    IF (.NOT. ALLOCATED(this%adj_jacobian)) THEN
      ALLOCATE (this%adj_jacobian(SIZE(this%row_elem, 1), SIZE(this%row_elem, 2)))
      this%adj_jacobian = 0.0_R_KIND
    END IF

    ! Exchange adj_rsd across halo regions to ensure consistency
    CALL sg%ExchangeMatOnHalo2D(vlevel, adj_rsd)

    ! Loop over each equation in reverse order (for reverse mode adjoint)
    DO i = sg%num_icell, 1, -1
      ! Accumulate adjoint of rhs
      adj_rhs(:, i) = adj_rhs(:, i) + adj_rsd(:, i)

      ! Loop over non-zero elements in the row
      DO j = this%num_elem(i), 1, -1
        col_idx = this%idx_cols(j, i)
        mat_elem = this%row_elem(j, i)

        ! Accumulate adjoint of sol (linear term)
        adj_sol(:, col_idx) = adj_sol(:, col_idx) - mat_elem * adj_rsd(:, i)

        ! Accumulate adjoint of row_elem (linear term)
        this%adj_row_elem(j, i) = this%adj_row_elem(j, i) - SUM(adj_rsd(:, i) * sol(:, col_idx))

        ! Handle nonlinear terms: if nonlinear, call Jacobian adjoint accumulation
        IF (this%is_nonlinear) THEN
          CALL this%CalJacobian_adj(vlevel, adj_sol, adj_rsd, this%adj_jacobian, sol)
        END IF
      END DO

      ! Reset adj_rsd after use
      adj_rsd(:, i) = 0.0_R_KIND
    END DO
  END SUBROUTINE CalResidus_adj

  !------------------------------------------------------------------------------
  ! Adjoint subroutine for Jacobian in nonlinear problems
  !------------------------------------------------------------------------------
  SUBROUTINE CalJacobian_adj(this, vlevel, adj_sol, adj_rsd, adj_jacobian, sol)
    !-------------------------------------------------------------------
    ! This subroutine computes the adjoint Jacobian using reverse-mode differentiation
    ! It makes use of the adjoint procedures in the gzm_adj_m module.
    !-------------------------------------------------------------------
    IMPLICIT NONE
    CLASS(psMatrix_adj_t), INTENT(INOUT) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel  ! Number of vertical levels
    REAL(r_kind), INTENT(INOUT) :: adj_sol(:, :)  ! Adjoint solution (e.g., adjoint stream function, velocity potential)
    REAL(r_kind), INTENT(INOUT) :: adj_rsd(:, :)  ! Adjoint residuals
    REAL(r_kind), INTENT(INOUT) :: adj_jacobian(:, :)  ! Adjoint Jacobian matrix
    REAL(r_kind), INTENT(IN) :: sol(:, :)  ! Forward-mode solution (stream function, velocity potential, etc.)

    INTEGER(i_kind) :: i, j
    REAL(r_kind), ALLOCATABLE :: jacobian_term(:, :)
    REAL(r_kind), ALLOCATABLE :: psi(:, :), chi(:, :), eta(:, :), h_star(:, :)

    ! Allocate memory for intermediate terms used in Jacobian computation
    ALLOCATE (jacobian_term(SIZE(sol, 1), SIZE(sol, 2)))
    ALLOCATE (psi(SIZE(sol, 1), SIZE(sol, 2)))
    ALLOCATE (chi(SIZE(sol, 1), SIZE(sol, 2)))
    ALLOCATE (eta(SIZE(sol, 1), SIZE(sol, 2)))
    ALLOCATE (h_star(SIZE(sol, 1), SIZE(sol, 2)))

    ! Example: Assuming adj_sol contains stream function (ψ) and velocity potential (χ)
    psi = sol(:, :)  ! Stream function (ψ)
    chi = sol(:, :)  ! Velocity potential (χ)
    eta = sol(:, :)  ! Vorticity (η)
    h_star = sol(:, :)  ! Height (h*)

    ! Use adjoint Jacobian procedure from gzm_adj_m to compute adjoints
    ! First adjoint Jacobian of (ψ, χ)
    CALL this%gzmOps%Jacobian_AD(psi, chi, adj_sol, adj_rsd, jacobian_term)
    adj_jacobian = adj_jacobian + jacobian_term  ! Accumulate into adjoint Jacobian

    ! Next, adjoint Jacobian of (η, ψ)
    CALL this%gzmOps%Jacobian_AD(eta, psi, adj_sol, adj_rsd, jacobian_term)
    adj_jacobian = adj_jacobian + jacobian_term  ! Accumulate

    ! Adjoint Jacobian of (η, χ)
    CALL this%gzmOps%Jacobian_AD(eta, chi, adj_sol, adj_rsd, jacobian_term)
    adj_jacobian = adj_jacobian + jacobian_term  ! Accumulate

    ! Adjoint Jacobian of (h*, ψ)
    CALL this%gzmOps%Jacobian_AD(h_star, psi, adj_sol, adj_rsd, jacobian_term)
    adj_jacobian = adj_jacobian + jacobian_term  ! Accumulate

    ! Deallocate memory used for temporary arrays
    DEALLOCATE (jacobian_term, psi, chi, eta, h_star)

  END SUBROUTINE CalJacobian_adj

END MODULE psMatrix_adjoint_m
