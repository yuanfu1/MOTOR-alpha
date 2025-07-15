MODULE ShallowWaterAdjointModel_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE ShallowWaterCommonAdjoint_m, ONLY: solve_poisson_adjoint, adjoint_kinetic_energy
  USE rhsShallowWaterAdjoint_m, ONLY: rhsShallowWaterAdjoint_t
  USE TimeIntegrationRK4Adjoint_m, ONLY: TimeIntegrationRK4Adjoint_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  USE CalVerDer_adj_m, ONLY: CalVerDer_adj_t
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State2NC_m, ONLY: Output_NC_State_SV_Adj
  IMPLICIT NONE

  PUBLIC :: ShallowWaterAdjointModel_t

  TYPE :: ShallowWaterAdjointModel_t
    TYPE(TimeIntegrationRK4Adjoint_t) :: rk4AdjointIntegrator
    TYPE(gzm_adj_t) :: gzmAdjOps
    TYPE(CalVerDer_adj_t) :: calVerOpsAdj
    TYPE(ForwardVars_t) :: forwardVars
    TYPE(SingleGrid_t) :: sg
    TYPE(rhsShallowWaterAdjoint_t) :: rhsAdjoint
  CONTAINS
    PROCEDURE :: initialize_s => initializeShallowWaterAdjointModel
    PROCEDURE :: runAdjointModel_s => runShallowWaterAdjointModel
  END TYPE ShallowWaterAdjointModel_t

CONTAINS

  ! Initialize the shallow water adjoint model
  SUBROUTINE initializeShallowWaterAdjointModel(this, X)
    CLASS(ShallowWaterAdjointModel_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X

    ! Initialize the forward variable storage
    CALL this%forwardVars%initialize(num_steps=24, num_stages=4, grid_size=X%sg%num_icell_global)

    ! Initialize the time integrator (RK4 adjoint)
    CALL this%rk4AdjointIntegrator%initialization_s(this%rhsAdjoint, X)

    ! Initialize the grid and the adjoint RHS calculator
    this%sg = X%sg
    CALL this%rhsAdjoint%initialize_s()  ! Initialize adjoint RHS model
  END SUBROUTINE initializeShallowWaterAdjointModel

  ! Run the adjoint model for the shallow water system
  SUBROUTINE runShallowWaterAdjointModel(this, adj_X, adj_Y)
    CLASS(ShallowWaterAdjointModel_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: adj_X
    TYPE(State_t), INTENT(OUT) :: adj_Y
    INTEGER :: it, num_steps, stage
    REAL(r_kind) :: dt
    CHARACTER(LEN=1024) :: output_path

    ! Set time step and number of steps (same as forward model)
    dt = 600.0_R_KIND
    num_steps = 24
    output_path = "output_folder_adj"  ! Specify the output folder for adjoint outputs

    ! Initialize adjoint state Y from input state X
    adj_Y = adj_X

    ! Adjoint time-stepping loop (reverse order from forward mode)
    DO it = num_steps, 1, -1
      ! Perform adjoint RK4 time integration using the `TimeIntegrationAdjointBase_s` method from RK4 adjoint
      CALL this%rk4AdjointIntegrator%TimeIntegrationAdjointBase_s(it, dt, adj_Y)

      ! Store the adjoint variables (vorticity, divergence, height, etc.) into NetCDF
      CALL store_adjoint_variables(it, adj_Y, output_path, "ShallowWaterAdjointOutput")
    END DO
  END SUBROUTINE runShallowWaterAdjointModel

  SUBROUTINE store_adjoint_variables(it, adj_Y, output_path, filename_prefix)
    INTEGER, INTENT(IN) :: it
    TYPE(State_t), INTENT(IN) :: adj_Y
    CHARACTER(LEN=*), INTENT(IN) :: output_path, filename_prefix
    CHARACTER(LEN=1024) :: filename

    ! Construct the filename with the current time step for adjoint output
    WRITE (filename, '(A,A,I0,A)') TRIM(output_path), "/", TRIM(filename_prefix), "_it", it, ".nc"

    ! Call the NetCDF output subroutine to write adjoint variables to file
    CALL Output_NC_State_SV_Adj(adj_Y, output_path, filename_prefix, "adjoint_vorticity")
    CALL Output_NC_State_SV_Adj(adj_Y, output_path, filename_prefix, "adjoint_divergence")
    CALL Output_NC_State_SV_Adj(adj_Y, output_path, filename_prefix, "adjoint_height_star")
  END SUBROUTINE store_adjoint_variables

END MODULE ShallowWaterAdjointModel_m
