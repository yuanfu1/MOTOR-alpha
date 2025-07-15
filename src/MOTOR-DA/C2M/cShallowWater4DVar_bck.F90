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
  !USE rhsBase_m, ONLY: rhsBase_t
  ! USE rhsShallowWater_m,  ONLY: rhsShallowWater_t
  USE ShallowWaterForwardModel_m, ONLY: ShallowWaterRHS_t
  USE ShallowWaterAdjointModel_m, ONLY: ShallowWaterAdjointRHS_t
  USE TimeIntegrationRK4_m, ONLY: TimeIntegrationRK4_t
  USE TimeIntegrationAdjointRK4_m, ONLY: TimeIntegrationAdjointRK4_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE gzm_m, ONLY: gzm_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  USE CalVerDer_adj_m, ONLY: CalVerDer_adj_t
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  USE CoriolisForce_m, ONLY: compute_coriolis_force, compute_coriolis_force_adjoint
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE rhsAdjointBase_m, ONLY: rhsAdjointBase_t

  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t) :: cShallowWater4DVar_t
    ! Forward model components
    TYPE(TimeIntegrationRK4_t) :: rk4Integrator            ! Forward time integrator
    TYPE(ShallowWaterRHS_t) :: rhs                         ! Forward RHS calculator

    ! Adjoint model components
    TYPE(TimeIntegrationAdjointRK4_t) :: adjRk4Integrator  ! Adjoint time integrator
    TYPE(ShallowWaterAdjointRHS_t) :: adjRhs               ! Adjoint RHS calculator

    ! Operators (shared or separate)
    TYPE(gzm_t) :: gzmOps                                  ! Gradient and divergence operator
    TYPE(CalVerDer_t) :: calVerOps                         ! Vertical derivative operator
    TYPE(poissonSolver_t) :: psOps                         ! Poisson solver

    ! Adjoint operators (if separate operators are needed)
    TYPE(gzmAdjoint_t) :: gzmAdjOps                        ! Adjoint gradient/divergence operator
    TYPE(CalVerDerAdjoint_t) :: calVerAdjOps               ! Adjoint vertical derivative operator
    TYPE(poissonSolverAdjoint_t) :: psAdjOps               ! Adjoint Poisson solver

    ! Grid information
    TYPE(SingleGrid_t) :: sg                               ! Grid information

    ! Storage for forward variables needed in the adjoint model
    TYPE(ForwardVariableStorage_t) :: forwardVars          ! Stores forward model variables

    TYPE(State_t), ALLOCATABLE :: analysis(:)              ! Analysis fields

  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => cShallowWaterForward
    PROCEDURE, PUBLIC :: tangentOpr_s => cShallowWaterTangent
    PROCEDURE, PUBLIC :: adjointOpr_s => cShallowWaterAdjoint
    PROCEDURE, PUBLIC :: bckwardOpr_s => cShallowWaterBackward
  END TYPE cShallowWater4DVar_t

  ! Define ForwardVariableStorage_t within the module
  TYPE :: ForwardVars_t
    ! Define storage variables here
    REAL(r_kind), ALLOCATABLE :: vorticity(:, :, :)
    REAL(r_kind), ALLOCATABLE :: divergence(:, :, :)
    REAL(r_kind), ALLOCATABLE :: height_star(:, :, :)
    ! Additional variables as needed
  CONTAINS
    PROCEDURE :: initialize
    PROCEDURE :: store
    PROCEDURE :: retrieve
  END TYPE ForwardVars_t

CONTAINS

! Implement the methods for ForwardVariableStorage_t
  SUBROUTINE initialize(this, num_steps, num_stages, grid_size)
    CLASS(ForwardVars_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: num_steps, num_stages, grid_size

    ALLOCATE (this%vorticity(grid_size, num_steps, num_stages))
    ALLOCATE (this%divergence(grid_size, num_steps, num_stages))
    ALLOCATE (this%height_star(grid_size, num_steps, num_stages))
    ! Initialize additional variables as needed
  END SUBROUTINE initialize

  SUBROUTINE store(this, it, stage, Y)
    CLASS(ForwardVars_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(IN) :: Y

    ! Extract data from Y and store it
    this%vorticity(:, it, stage) = Y%fields(1)%DATA
    this%divergence(:, it, stage) = Y%fields(2)%DATA
    this%height_star(:, it, stage) = Y%fields(3)%DATA
    ! Store additional variables as needed
  END SUBROUTINE store

  SUBROUTINE retrieve(this, it, stage, Y)
    CLASS(ForwardVars_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: it, stage
    TYPE(State_t), INTENT(INOUT) :: Y

    ! Retrieve data and assign to Y
    Y%fields(1)%DATA = this%vorticity(:, it, stage)
    Y%fields(2)%DATA = this%divergence(:, it, stage)
    Y%fields(3)%DATA = this%height_star(:, it, stage)
    ! Retrieve additional variables as needed
  END SUBROUTINE retrieve

  SUBROUTINE initialize_s(this, mgStart, mgEnd, X)
    CLASS(cShallowWater4DVar_t) :: this
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd
    TYPE(State_t), INTENT(IN) :: X(mgStart:mgEnd)
    INTEGER(i_kind) :: img, i
    INTEGER :: num_steps
    INTEGER :: num_stages, grid_size

    ! Initialize the analysis fields
    ALLOCATE (this%analysis(mgStart:mgEnd))
    this%analysis = X

    ! Initialize grid and operators
    this%sg = X(mgStart)%sg
    this%gzmOps = gzm_t(this%sg)
    this%calVerOps = CalVerDer_t(this%sg)
    this%psOps = poissonSolver_t(this%sg)

    ! Initialize the forward RHS calculator
    CALL this%rhs%initialize_s()
    this%rhs%psOps => this%psOps
    this%rhs%gzmOps => this%gzmOps
    this%rhs%calVerOps => this%c alVerOps

    ! Initialize the adjoint RHS calculator
    CALL this%adjRhs%initialize_s()
    this%adjRhs%psOps => this%psOps
    this%adjRhs%gzmOps => this%gzmOps
    this%adjRhs%calVerOps => this%calVerOps

    ! Initialize forward and adjoint time integrators
    CALL this%rk4Integrator%initialize_s(this%rhs, X)
    CALL this%adjRk4Integrator%initialize_s(this%adjRhs, X)

    ! Set the number of time steps (should match the forward model)
    num_steps = 24  ! For example
    num_stages = 4  ! For RK4
    grid_size = this%sg%num_icell_global  ! Assuming this is the grid size

    ! Initialize forward variable storage
    CALL this%forwardVars%initialize(num_steps, num_stages, grid_size)

    ! Setting the control variables (vorticity, divergence, thickness)
    DO img = mgStart, mgEnd
      ASSOCIATE (sg => this%analysis(img)%sg)
        DO i = LBOUND(this%analysis(img)%fields, 1), UBOUND(this%analysis(img)%fields, 1)
          SELECT CASE (TRIM(this%analysis(img)%fields(i)%Get_Name()))
          CASE ('vorticity', 'divergence', 'thick')
            ! Only initial condition as controls:
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
      END ASSOCIATE
    END DO

  END SUBROUTINE initialize_s

  FUNCTION cShallowWaterForward(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
    INTEGER :: it, stage
    REAL(r_kind) :: dt
    INTEGER :: num_steps

    ! Set time step and number of steps
    dt = 600.0_R_KIND
    num_steps = 24

    ! Initialize state variable Y
    Y = X

    ! Time-stepping loop
    DO it = 1, num_steps
      DO stage = 1, 4  ! Assuming 4 stages for RK4
        ! Forward integration step
        CALL this%rk4Integrator%TimeIntegrationStage_s(it, stage, dt, Y)

        ! Store forward variables for adjoint model
        CALL this%forwardVars%store(it, stage, Y)
      END DO
    END DO
  END FUNCTION cShallowWaterForward

  FUNCTION cShallowWaterBckward(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = this%analysis(X%sg%gLevel)
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) / X%sg%SVapor
        END DO
      CASE DEFAULT
        Y%fields(j)%DATA = X%fields(j)%DATA
      END SELECT
    END DO
  END FUNCTION cShallowWaterBckward

  FUNCTION cShallowWaterTangent(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    Y = X ! Identical tangent
  END FUNCTION cShallowWaterTangent

  FUNCTION cShallowWaterAdjoint(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X   ! Adjoint variables
    TYPE(State_t) :: Y
    INTEGER :: it, stage
    REAL(r_kind) :: dt
    INTEGER :: num_steps
    TYPE(State_t) :: forward_state  ! Declaration of forward_state

    ! Set time step and number of steps
    dt = 600.0_R_KIND
    num_steps = 24

    ! Initialize adjoint state variable Y
    Y = X

    ! Backward time-stepping loop
    DO it = num_steps, 1, -1
      DO stage = 4, 1, -1  ! Reverse stages
        ! Retrieve forward variables
        CALL this%forwardVars%retrieve(it, stage, forward_state)

        ! Adjoint integration step
        CALL this%adjRk4Integrator%TimeIntegrationStage_s(it, stage, dt, Y, forward_state)
      END DO
    END DO
  END FUNCTION cShallowWaterAdjoint

  FUNCTION cShallowWaterBckAdjoint(this, X) RESULT(Y)
    CLASS(cShallowWater4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = X ! Identity matrix transpose
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) / X%sg%SVapor
        END DO
      END SELECT
    END DO
  END FUNCTION cShallowWaterBckAdjoint
END MODULE cShallowWater4DVar_m
