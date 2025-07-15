MODULE ShallowWaterForwardModel_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE ShallowWaterCommon_m, ONLY: ForwardVars_t, solve_poisson, kinetic_energy
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t
  USE TimeIntegrationRK4_m, ONLY: TimeIntegrationRK4_t
  USE gzm_m, ONLY: gzm_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State2NC_m, ONLY: Output_NC_State_SV  ! To store the data in NetCDF format
  IMPLICIT NONE

  PUBLIC :: ShallowWaterForwardModel_t

  TYPE :: ShallowWaterForwardModel_t
    TYPE(TimeIntegrationRK4_t) :: rk4Integrator  ! Use the RK4 integrator type
    TYPE(ForwardVars_t) :: forwardVars
    TYPE(SingleGrid_t) :: sg
    TYPE(rhsShallowWater_t) :: rhs   ! Use the rhsShallowWater type for RHS computation
  CONTAINS
    PROCEDURE :: initialize_s => initializeShallowWaterForwardModel
    PROCEDURE :: runModel_s => runShallowWaterForwardModel
  END TYPE ShallowWaterForwardModel_t

CONTAINS

  ! Initialize the shallow water forward model
  SUBROUTINE initializeShallowWaterForwardModel(this, X)
    CLASS(ShallowWaterForwardModel_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X

    ! Initialize the forward variable storage
    CALL this%forwardVars%initialize(num_steps=24, num_stages=4, grid_size=X%sg%num_icell_global)

    ! Initialize the time integrator (RK4)
    CALL this%rk4Integrator%initialization_s(this%rhs, X, this%rk4Integrator%yamlFile)

    ! Initialize the grid and the RHS calculator
    this%sg = X%sg
    CALL this%rhs%initialize_s()   ! Initialize RHS model
  END SUBROUTINE initializeShallowWaterForwardModel

  ! Run the shallow water forward model
  SUBROUTINE runShallowWaterForwardModel(this, X, Y)
    CLASS(ShallowWaterForwardModel_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(State_t), INTENT(OUT) :: Y
    INTEGER :: it, num_steps, stage
    REAL(r_kind) :: dt
    CHARACTER(LEN=1024) :: output_path

    ! Set time step and number of steps
    dt = 600.0_R_KIND
    num_steps = 24
    output_path = "output_folder"  ! Specify the output folder

    ! Initialize state Y from input state X
    Y = X

    ! Time-stepping loop
    DO it = 1, num_steps
      ! Perform RK4 time integration using the `TimeIntegrationBase_s` method from RK4
      CALL this%rk4Integrator%TimeIntegrationBase_fwd(dt, it, X, Y)

      ! Store the forward variables (vorticity, divergence, height, etc.) into NetCDF
      CALL store_forward_variables(it, Y, output_path, "ShallowWaterOutput")
    END DO
  END SUBROUTINE runShallowWaterForwardModel

  SUBROUTINE store_forward_variables(it, Y, output_path, filename_prefix)
    INTEGER, INTENT(IN) :: it
    TYPE(State_t), INTENT(IN) :: Y
    CHARACTER(LEN=*), INTENT(IN) :: output_path, filename_prefix
    CHARACTER(LEN=1024) :: filename

    ! Construct the filename with the current time step
    WRITE (filename, '(A,A,I0,A)') TRIM(output_path), "/", TRIM(filename_prefix), "_it", it, ".nc"

    ! Call the NetCDF output subroutine to write variables to file
    CALL Output_NC_State_SV(Y, output_path, filename_prefix, "vorticity")
    CALL Output_NC_State_SV(Y, output_path, filename_prefix, "divergence")
    CALL Output_NC_State_SV(Y, output_path, filename_prefix, "height_star")
  END SUBROUTINE store_forward_variables

END MODULE ShallowWaterForwardModel_m
