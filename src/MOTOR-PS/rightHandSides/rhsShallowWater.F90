!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.rightHandSide.rhsShallowWater
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-11   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a shallow water equation system using a Z-grid model.
MODULE rhsShallowWater_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE rhsBase_m, ONLY: rhsBase_t
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE geometry_m, ONLY: geometry_t

  USE gzm_m, ONLY: gzm_t

  IMPLICIT NONE

  TYPE, EXTENDS(rhsBase_t) :: rhsShallowWater_t
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s => initialize_shallowWater
    PROCEDURE, PUBLIC :: rightHandSide => rhsShallowWater_s
  END TYPE rhsShallowWater_t

CONTAINS

  ! Initialization routine for rhsShallowWater_t
  SUBROUTINE initialize_shallowWater(this, sg, psolver)
    CLASS(rhsShallowWater_t) :: this
    TYPE(poissonSolver_t), TARGET :: psolver
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    ! Allocate space for variable indices (for vorticity, divergence, height, velocity_u, velocity_v)
    ALLOCATE (this%variableIdx(5))
    this%variableIdx(1) = 1  ! Vorticity
    this%variableIdx(2) = 2  ! Divergence
    this%variableIdx(3) = 3  ! Height
    this%variableIdx(4) = 4  ! Velocity_u
    this%variableIdx(5) = 5  ! Velocity_v

    this%poissonSolver => psolver  ! Assign the Poisson solver

  END SUBROUTINE initialize_shallowWater

  ! Right-hand side calculation for shallow water model
  SUBROUTINE rhsShallowWater_s(this, it, current, forward)
    CLASS(rhsShallowWater_t) :: this
    INTEGER(i_kind), INTENT(IN) :: it
    TYPE(State_t), INTENT(IN) :: current
    TYPE(State_t), INTENT(INOUT) :: forward

    ! Local variables:
    ! REAL(r_kind), ALLOCATABLE :: rhs_vorticity(:, :), rhs_divergence(:, :)
    ! REAL(r_kind), ALLOCATABLE :: rhs_height(:, :)
    ! REAL(r_kind), ALLOCATABLE :: streamfunc(:, :), velocity_pot(:, :), K(:, :)
    ! REAL(r_kind), ALLOCATABLE :: rhs_velocity_u(:, :), rhs_velocity_v(:, :)

    ! Extract variables from the current state using pointer assignment
    PRINT*,'Entering rhsShallowWater_s...'
    forward = forward%zeroCopy()  ! initialize forward state

    ! Vorticity equation:
    ! CALL this%forward%fields(this%variableIdx(1))%DATA(:, :, it) = &
    !   current%fields(this%variableIdx(1))%DATA(:, :, it) +
    ! ASSOCIATE ( &
    !   vlevel => current%sg%vLevel, num_cell => current%sg%num_cell, &
    !   vorticity => current%fields(this%variableIdx(1))%DATA, &
    !   divergence => current%fields(this%variableIdx(2))%DATA, &
    !   heightStar => current%fields(this%variableIdx(3))%DATA, &
    !   velocity_u => current%fields(this%variableIdx(4))%DATA, &
    !   velocity_v => current%fields(this%variableIdx(5))%DATA)

    !   ! Allocate memory for right-hand side (RHS) calculations and Poisson solver outputs
    !   ALLOCATE (rhs_vorticity(vlevel, num_cell))
    !   ALLOCATE (rhs_divergence(vlevel, num_cell))
    !   ALLOCATE (rhs_height(vlevel, num_cell))
    !   ALLOCATE (rhs_velocity_u(vlevel, num_cell))
    !   ALLOCATE (rhs_velocity_v(vlevel, num_cell))
    !   ALLOCATE (streamfunc(vlevel, num_cell))
    !   ALLOCATE (velocity_pot(vlevel, num_cell))
    !   ALLOCATE (K(vlevel, num_cell))

    !   ! Step 1: Solve Poisson equation for streamfunction (ψ) and velocity potential (χ)
    !   PRINT *, 'rhsShallowWater - Solving a Poisson...'
    !   CALL solve_poisson(this%poissonSolver, current%sg%vLevel, current%sg%gLevel, &
    !                      vorticity(:, :, it), divergence(:, :, it), streamfunc, velocity_pot, current%sg)
    !   PRINT *, 'rhsShallowWater - The Poisson solver is done'

    !   ! Step 2: Compute RHS for kinetic energy
    !   CALL kinetic_energy(streamfunc, velocity_pot, K, rhs_velocity_u, rhs_velocity_v)

    !   ! Step 3: Compute RHS for vorticity, divergence, and height using the updated velocity fields
    !   PRINT *, 'rhsShallowWater - Calculating RHS for vorticity...'
    !   CALL compute_rhs_vorticity(vorticity(:, :, it), velocity_pot, streamfunc, rhs_vorticity)
    !   CALL compute_rhs_divergence(vorticity(:, :, it), velocity_pot, streamfunc, K, heightStar(:, :, it), rhs_divergence)
    !   CALL compute_rhs_height(vorticity(:, :, it), velocity_pot, streamfunc, heightStar(:, :, it), rhs_height)

    !   ! Assign computed RHS values to the forward array
    !   forward%fields(this%variableIdx(1))%DATA(:, :, it) = rhs_vorticity
    !   forward%fields(this%variableIdx(2))%DATA(:, :, it) = rhs_divergence
    !   forward%fields(this%variableIdx(3))%DATA(:, :, it) = rhs_height

    ! END ASSOCIATE

    ! Deallocate temporary arrays
    ! DEALLOCATE (rhs_vorticity, rhs_divergence, rhs_height, streamfunc, velocity_pot, rhs_velocity_u, rhs_velocity_v)
  END SUBROUTINE rhsShallowWater_s

END MODULE rhsShallowWater_m
