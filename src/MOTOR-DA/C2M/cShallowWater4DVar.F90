!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/09/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines a control variable case using a Z-grid shallow water equation model to
!! determine model states at other time frames.
!! Control variables: $\zeta$, $\delta$, $thick$ where thick is thickness perturbed from bottom
!! topography.
!! Reference:
!! "Numerical Integration of the Shallow-Water Equations on a Twisted Icosahedral Grid. Part I:
!!  Basic Design and Results of Tests" by Ross Heikes and David A. Randall, (1995) Mon. Wea. Rev.
!!  pp 1862-1880.
!!
MODULE cShallowWater4DVar_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE C2MBase_m, ONLY: C2MBase_t
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ShallowWaterCommon_m, ONLY: solve_poisson, kinetic_energy
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t) :: cShallowWater4DVar_t
    REAL(r_kind), ALLOCATABLE :: verticalProfiles(:, :, :)  ! For scaling pressure and water vapor
    TYPE(rhsShallowWater_t) :: rhs   ! RHS solver for shallow water equations
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => cShallowWaterForward
    PROCEDURE, PUBLIC :: adjointOpr_s => cShallowWaterAdjoint  ! Placeholder adjoint operator
    PROCEDURE, PUBLIC :: tangentOpr_s => cShallowWaterTangent  ! Placeholder tangent operator
    PROCEDURE, PUBLIC :: bckwardOpr_s => cShallowWaterBackward ! Placeholder backward operator
  END TYPE cShallowWater4DVar_t

CONTAINS

  !---------------------------------------------------------------------------
  !> Initialization routine for the control variables
  !! This routine sets up the control variables and initializes the RHS solver.
  SUBROUTINE initialize_s(this, mgStart, mgEnd, X, dyc)

    IMPLICIT NONE

    CLASS(cShallowWater4DVar_t) :: this
    INTEGER(i_kind) :: mgStart, mgEnd
    TYPE(State_t), INTENT(IN) :: X(mgStart:mgEnd)
    CLASS(dyCoresBase_t), TARGET, OPTIONAL :: dyc

    ! Local variables
    INTEGER(i_kind) :: i, img
    TYPE(SingleGrid_t) :: sg

    ! Initialize the rhs of the shallow water equations:
    CALL this%rhs%initialize_s()

    ! Allocate and assign model states to analysis
    ALLOCATE (this%analysis(mgStart:mgEnd))
    this%analysis = X

    ! Loop over the grid points and set control variables
    DO img = mgStart, mgEnd
      sg = this%analysis(img)%sg

      DO i = LBOUND(this%analysis(img)%fields, 1), UBOUND(this%analysis(img)%fields, 1)
        PRINT *, 'Module state: ', i, this%analysis(img)%fields(i)%Get_Name()

        SELECT CASE (TRIM(this%analysis(img)%fields(i)%Get_Name()))
        CASE ('vorticity', 'divergence', 'thick')
          ! Set control variables related to initial conditions
          this%analysis(img)%fields(i)%numTotalMask = sg%num_icell_global
          this%analysis(img)%fields(i)%maskHorizontal = 1
          this%analysis(img)%fields(i)%maskTemporal = 0
          this%analysis(img)%fields(i)%maskTemporal(1) = 1
          this%analysis(img)%fields(i)%maskVertical = 0
        CASE DEFAULT
          ! Non-control variables
          this%analysis(img)%fields(i)%numTotalMask = 0
          this%analysis(img)%fields(i)%maskHorizontal = 0
          this%analysis(img)%fields(i)%maskTemporal = 0
          this%analysis(img)%fields(i)%maskVertical = 0
        END SELECT
      END DO
    END DO

  END SUBROUTINE initialize_s

  !---------------------------------------------------------------------------
  !> Forward operator for the Z-grid shallow water model
  !! This function computes the forward evolution of the state over time, using
  !! the shallow water equations. It applies the Runge-Kutta time integration
  !! method and stores the results at each time step.
  FUNCTION cShallowWaterForward(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables
    INTEGER(i_kind) :: i, num_steps
    REAL(r_kind), ALLOCATABLE :: current(:, :, :), forward(:, :, :)
    REAL(r_kind), ALLOCATABLE :: streamfunc(:, :), velocity_pot(:, :)
    REAL(r_kind) :: dt

    ! Set the time step and number of steps
    dt = 600.0_R_KIND
    num_steps = 24

    ! Initialize Y with the input state X
    Y = X

    ! Allocate space for current and forward fields
    ALLOCATE (current(X%sg%vLevel, X%sg%num_cell, 3))  ! vorticity, divergence, height
    ALLOCATE (forward(X%sg%vLevel, X%sg%num_cell, 3))

    ! Allocate space for streamfunction and velocity potential
    ALLOCATE (streamfunc(X%sg%vLevel, X%sg%num_cell))
    ALLOCATE (velocity_pot(X%sg%vLevel, X%sg%num_cell))

    ! Time-stepping loop using the Runge-Kutta method
    DO i = 1, num_steps
      PRINT *, 'Time step: ', i

      ! Extract current state fields: vorticity, divergence, height
      current(:, :, 1) = X%fields(1)%DATA(:, :, i)  ! Vorticity for time step i
      current(:, :, 2) = X%fields(2)%DATA(:, :, i)  ! Divergence for time step i
      current(:, :, 3) = X%fields(3)%DATA(:, :, i)  ! Height for time step i

      ! Step 1: Solve Poisson equation to get streamfunction and velocity potential
      CALL solve_poisson(this%rhs%poissonSolver, X%sg%vLevel, X%sg%gLevel, current(:, :, 1), current(:, :, 2), streamfunc, velocity_pot, X%sg)

      ! Step 2: Compute the right-hand side (RHS) of the shallow water equations
      CALL this%rhs%rightHandSide(i, X, Y)

      ! Step 3: Update the state fields for the next time step
      Y%fields(1)%DATA(:, :, i) = forward(:, :, 1)  ! Vorticity
      Y%fields(2)%DATA(:, :, i) = forward(:, :, 2)  ! Divergence
      Y%fields(3)%DATA(:, :, i) = forward(:, :, 3)  ! Height

    END DO

    ! Deallocate temporary arrays
    DEALLOCATE (current, forward, streamfunc, velocity_pot)

  END FUNCTION cShallowWaterForward

  !---------------------------------------------------------------------------
  !> Placeholder for adjoint operator
  !! This is a placeholder function for the adjoint operator. It simply returns
  !! the input state. A full implementation will be needed for adjoint mode.
  FUNCTION cShallowWaterAdjoint(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Placeholder: simply return the input state for now
    Y = X
  END FUNCTION cShallowWaterAdjoint

  !---------------------------------------------------------------------------
  !> Placeholder for tangent operator
  !! This is a placeholder function for the tangent operator. It simply returns
  !! the input state. A full implementation will be needed for tangent mode.
  FUNCTION cShallowWaterTangent(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Placeholder: simply return the input state for now
    Y = X
  END FUNCTION cShallowWaterTangent

  !---------------------------------------------------------------------------
  !> Placeholder for backward operator
  !! This is a placeholder function for the backward operator. It simply returns
  !! the input state. A full implementation will be needed for backward mode.
  FUNCTION cShallowWaterBackward(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Placeholder: simply return the input state for now
    Y = X
  END FUNCTION cShallowWaterBackward

END MODULE cShallowWater4DVar_m
