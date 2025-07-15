MODULE CoriolisForce_m
  USE kinds_m, ONLY: r_kind, i_kind
  USE SingleGrid_m, ONLY: SingleGrid_t   ! Use SingleGrid type for accessing grid structure
  USE parameters_m, ONLY: Omega          ! Use Omega defined in parameters_m
  IMPLICIT NONE

CONTAINS

  !------------------------------------------------------------------
  !> Initialization subroutine for Coriolis force. Allocates and zeroes the array.
  !! @param coriolis_force  Output array for the computed Coriolis force (2D array)
  !! @param grid            Input grid structure containing latitude information
  SUBROUTINE initialize_coriolis_force(coriolis_force, grid)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: coriolis_force(:, :)   ! Coriolis force 2D array
    TYPE(SingleGrid_t), INTENT(IN) :: grid                          ! Input grid structure containing latitude
    INTEGER(i_kind) :: i, j                                         ! Loop indices

    ! Allocate the Coriolis force array based on the grid's structure.
    IF (.NOT. ALLOCATED(coriolis_force)) THEN
      ! Coriolis is an array of vertical levels and horizontal grid:
      ! ALLOCATE(coriolis_force(grid%num_icell, UBOUND(grid%cell_stcl, 1))) ! This is incorrect
      ALLOCATE (coriolis_force(grid%vLevel, grid%num_cell))
    END IF

    ! Initialize the Coriolis force array to zero.
    coriolis_force = 0.0_R_KIND

  END SUBROUTINE initialize_coriolis_force

  !------------------------------------------------------------------
  !> Forward mode: Compute the Coriolis force (f = 2 * Omega * DSIN(latitude)).
  !! @param coriolis_force  Output array for the computed Coriolis force (2D array)
  !! @param grid            Input grid structure containing latitude information
  SUBROUTINE compute_coriolis_force(coriolis_force, grid)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: coriolis_force(:, :)   ! Coriolis force 2D array
    TYPE(SingleGrid_t), INTENT(IN) :: grid                          ! Input grid structure containing latitude
    INTEGER(i_kind) :: i, j                                         ! Loop indices
    REAL(r_kind) :: latitude                                        ! Latitude for each cell (in radians)

    ! Check if the Coriolis force array is allocated; if not, allocate it.
    IF (.NOT. ALLOCATED(coriolis_force)) THEN
      CALL initialize_coriolis_force(coriolis_force, grid)
    END IF

    ! Loop over each inner cell and compute Coriolis force based on latitude
    ! Yuanfu Xie has changed the do loop from num_icell to num_cell 2024-10-08:
    DO i = 1, grid%num_cell    ! Loop over all inner cells
      !DO j = 1, UBOUND(grid%cell_stcl, 1)   ! Loop over neighbors (stencils) for each cell
      latitude = grid%cell_cntr(1, i)  ! Retrieve the latitude (in radians) for the cell
      coriolis_force(:, i) = 2.0_R_KIND * Omega * DSIN(latitude)  ! Calculate Coriolis force using DSIN for precision
      !END DO
    END DO

  END SUBROUTINE compute_coriolis_force

  !------------------------------------------------------------------
  !> Adjoint mode: Compute the adjoint (derivatives) of the Coriolis force
  !!  with respect to the input variables. This is used for sensitivity analysis.
  !! @param adj_coriolis_force  Input adjoint of Coriolis force (2D array)
  !! @param adj_latitude        Output adjoint of the latitude (derivatives with respect to latitude)
  !! @param grid                Input grid structure containing latitude information
  SUBROUTINE compute_coriolis_force_adjoint(adj_coriolis_force, adj_latitude, grid)
    REAL(r_kind), INTENT(IN) :: adj_coriolis_force(:, :)   ! Adjoint input from Coriolis force
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: adj_latitude(:)  ! Adjoint output (derivatives with respect to latitude)
    TYPE(SingleGrid_t), INTENT(IN) :: grid                ! Grid structure containing latitude information
    INTEGER(i_kind) :: i, j                               ! Loop indices
    REAL(r_kind) :: latitude                              ! Latitude of the cell

    ! Allocate adjoint latitude array
    IF (.NOT. ALLOCATED(adj_latitude)) THEN
      ALLOCATE (adj_latitude(grid%num_icell))  ! Adjoint latitude is a 1D array
    END IF

    ! Initialize adjoint latitude to zero
    adj_latitude = 0.0_R_KIND

    ! Loop over all inner cells and neighbors
    DO i = 1, grid%num_icell
      DO j = 1, UBOUND(grid%cell_stcl, 1)
        latitude = grid%cell_cntr(1, i)  ! Retrieve the latitude for this cell
        ! Derivative of the Coriolis force with respect to latitude is (2 * Omega * DCOS(latitude))
        adj_latitude(i) = adj_latitude(i) + adj_coriolis_force(i, j) * 2.0_R_KIND * Omega * DCOS(latitude)
      END DO
    END DO

  END SUBROUTINE compute_coriolis_force_adjoint

END MODULE CoriolisForce_m
