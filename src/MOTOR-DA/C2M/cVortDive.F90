!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/4/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines a control variable case using vorticity and divergence as control variables.
!! Control variables: $\zeta$, $\delta$, $\ln(P)$, $temp$ and $rh$, where rh is scaled mixing ratio
!! by an exponential function profile.
!!
MODULE cVortDive_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE C2MBase_m, ONLY: C2MBase_t
  USE State_m, ONLY: State_t
  USE ShallowWaterCommon_m, ONLY: solve_poisson, kinetic_energy  ! Poisson solver and energy computation
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t  ! RHS computation for vorticity and divergence
  USE SingleGrid_m, ONLY: SingleGrid_t  ! Multigrid model
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE SingleGrid_m, ONLY: SingleGrid_t 

  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t):: cVortDive_t
    TYPE(rhsShallowWater_t) :: rhs  ! RHS solver for vorticity/divergence
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => cVortDiveForward
    PROCEDURE, PUBLIC :: bckwardOpr_s => cVortDiveBckward
    PROCEDURE, PUBLIC :: tangentOpr_s => cVortDiveTangent
    PROCEDURE, PUBLIC :: adjointOpr_s => cVortDiveAdjoint
  END TYPE cVortDive_t

CONTAINS

  ! Initialization routine
  SUBROUTINE initialize_s(this, sg, X, dyc)

    IMPLICIT NONE

    CLASS(cVortDive_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), INTENT(IN) :: X
    CLASS(dyCoresBase_t), TARGET, OPTIONAL :: dyc

    ! Local variables:
    INTEGER(i_kind) :: i

    ! Set analysis to the multigrid model states as default:
    this%analysis = X
    this%sg => sg

    ! Loop over the fields and check their names
    DO i = 1, UBOUND(this%analysis%fields, 1)
      PRINT *, 'Module state: ', i, this%analysis%fields(i)%Get_Name()

      SELECT CASE (TRIM(this%analysis%fields(i)%Get_Name()))
      CASE ('vorticity', 'divergence', 'ln(P)', 'temp', 'rh')
        ! Set as control variables
        this%analysis%fields(i)%numTotalMask = sg%num_icell_global
        this%analysis%fields(i)%maskHorizontal = 1
        this%analysis%fields(i)%maskTemporal = 0
        this%analysis%fields(i)%maskTemporal(1) = 1
        this%analysis%fields(i)%maskVertical = 0
      CASE DEFAULT
        ! Non-control variables
        this%analysis%fields(i)%numTotalMask = 0
        this%analysis%fields(i)%maskHorizontal = 0
        this%analysis%fields(i)%maskTemporal = 0
        this%analysis%fields(i)%maskVertical = 0
      END SELECT

    END DO

  END SUBROUTINE initialize_s

  ! Forward operator (computes the state evolution based on vorticity/divergence control variables)
  FUNCTION cVortDiveForward(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables
    INTEGER(i_kind) :: i, num_steps
    REAL(r_kind), ALLOCATABLE :: current(:, :, :), forward(:, :, :)
    REAL(r_kind), ALLOCATABLE :: streamfunc(:, :), velocity_pot(:, :), height_star(:, :)
    REAL(r_kind) :: dt, T

    !T = 288.0_r_kind  ! Constant temperature (can be replaced by another model)

    ! Set the time step and number of steps
    dt = 600.0_R_KIND  ! 600 seconds time step
    num_steps = 24  ! 24 time steps

    ! Initialize Y with the input state X
    Y = X

    ! Allocate space for current and forward fields
    ALLOCATE (current(X%sg%vLevel, X%sg%num_cell, 3))  ! vorticity, divergence, lnP
    ALLOCATE (forward(X%sg%vLevel, X%sg%num_cell, 3))

    ! Allocate space for streamfunction, velocity potential, and height_star
    ALLOCATE (streamfunc(X%sg%vLevel, X%sg%num_cell))
    ALLOCATE (velocity_pot(X%sg%vLevel, X%sg%num_cell))
    ALLOCATE (height_star(X%sg%vLevel, X%sg%num_cell))

    ! Time-stepping loop using the Runge-Kutta method
    DO i = 1, num_steps
      PRINT *, 'Time step: ', i

      ! Step 1: Extract current state fields
      current(:, :, 1) = X%fields(1)%DATA(:, :, i)  ! Vorticity
      current(:, :, 2) = X%fields(2)%DATA(:, :, i)  ! Divergence
      current(:, :, 3) = X%fields(3)%DATA(:, :, i)  ! lnP (log of pressure)

      ! Step 2: Solve Poisson equation for streamfunction and velocity potential
      CALL solve_poisson(this%rhs%poissonSolver, X%sg%vLevel, X%sg%gLevel, current(:, :, 1), current(:, :, 2), streamfunc, velocity_pot, X%sg)

      ! Step 3: Compute RHS using the shallow water equations
      CALL this%rhs%rightHandSide(i, X, Y)

      ! ! Step 4: Calculate height_star using ln(P)
      ! height_star = (dry_air_gas_const * T / g) * (LOG(surface_ref_pres) - current(:, :, 3))

      ! Step 5: Update the model state
      Y%fields(1)%DATA(:, :, i) = forward(:, :, 1)  ! Vorticity
      Y%fields(2)%DATA(:, :, i) = forward(:, :, 2)  ! Divergence
      Y%fields(3)%DATA(:, :, i) = forward(:, :, 3)  ! lnP (log of pressure)

    END DO

    ! Deallocate temporary arrays
    DEALLOCATE (current, forward, streamfunc, velocity_pot, height_star)

  END FUNCTION cVortDiveForward

  ! Placeholder backward operator (extend as needed)
  FUNCTION cVortDiveBckward(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
    Y = X
  END FUNCTION cVortDiveBckward

  ! Placeholder tangent operator (extend as needed)
  FUNCTION cVortDiveTangent(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
    Y = X
  END FUNCTION cVortDiveTangent

  ! Placeholder adjoint operator (extend as needed)
  FUNCTION cVortDiveAdjoint(this, X, dt) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(IN), OPTIONAL :: dt
    TYPE(State_t) :: Y
    Y = X
  END FUNCTION cVortDiveAdjoint

END MODULE cVortDive_m
