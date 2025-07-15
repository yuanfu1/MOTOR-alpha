PROGRAM TestForwardModel
  USE kinds_m, ONLY: r_kind, i_kind
  USE ShallowWaterForwardModel_m, ONLY: initialize_objects, time_step_rk4
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE PoissonSolver_m, ONLY: poissonSolver_t
  USE cDefault_m, ONLY: cDefault_t
  IMPLICIT NONE

  ! Declare variables for the shallow water model
  TYPE(State_t), ALLOCATABLE :: state(:)      ! Holds the state variables
  TYPE(SingleGrid_t) :: sg                    ! Grid structure object
  TYPE(poissonSolver_t) :: poissonSolver      ! Poisson solver for the model
  TYPE(mpddGlob_t) :: mpddGlob                ! Global model parameters
  TYPE(cDefault_t) :: control                 ! Control object for initialization

  INTEGER(i_kind) :: mgStart, mgEnd           ! Multigrid start and end levels
  REAL(r_kind) :: dt                          ! Time step size
  REAL(r_kind) :: total_time, current_time    ! Total and current simulation time
  INTEGER(i_kind) :: nSteps, step             ! Number of steps and current step
  INTEGER :: output_interval                  ! How often to output results
  CHARACTER(LEN=256) :: yaml_file             ! YAML configuration file path
  CHARACTER(LEN=256) :: line
  INTEGER :: iostat                          ! Status for reading the YAML file

  ! Initialize default values in case they are missing in YAML file
  mgStart = 1
  mgEnd = 10
  dt = 300.0_R_KIND
  total_time = 86400.0_R_KIND  ! Default to 1 day
  output_interval = 50

  ! Set YAML configuration file path
  yaml_file = 'path/to/config.yaml'  ! Modify to your file path

  ! Open the YAML configuration file and parse it
  OPEN (UNIT=10, FILE=yaml_file, STATUS='OLD', IOSTAT=iostat)
  IF (iostat /= 0) THEN
    PRINT *, "Error opening YAML file!"
    STOP
  END IF

  DO
    READ (10, '(A)', IOSTAT=iostat) line
    IF (iostat /= 0) EXIT  ! Exit the loop when the file ends

    ! Parse mgStart and mgEnd
    IF (INDEX(line, 'mgStart') > 0) READ (line, '(A, I5)') mgStart
    IF (INDEX(line, 'mgEnd') > 0) READ (line, '(A, I5)') mgEnd

    ! Parse time step (dt)
    IF (INDEX(line, 'time_steps_g14') > 0) READ (line, '(A, F10.5)') dt

    ! Parse total simulation time (in seconds)
    IF (INDEX(line, 'end_time') > 0) THEN
      total_time = 86400.0_R_KIND  ! Example: setting total time as 1 day
    END IF

    ! Parse output interval
    IF (INDEX(line, 'nCycle') > 0) READ (line, '(A, I5)') output_interval

  END DO

  ! Close the YAML configuration file
  CLOSE (UNIT=10)

  ! Initialize objects (grid, state variables, etc.)
  ALLOCATE (state(mgStart:mgEnd))
  CALL control%initialize_s(mgStart, mgEnd, state)
  CALL initialize_objects(mpddGlob, mgStart, mgEnd, state, 1)

  ! Main time-stepping simulation loop
  current_time = 0.0_R_KIND
  nSteps = INT(total_time / dt)

  DO step = 1, nSteps
    ! Perform one time step using RK4 method
    ! Access 2D slices from the 3D arrays
    CALL time_step_rk4(dt, state(1)%fields(1)%DATA(:, :, 1), state(1)%fields(2)%DATA(:, :, 1), state(1)%fields(3)%DATA(:, :, 1), poissonSolver, sg)

    ! Update current simulation time
    current_time = current_time + dt

    ! Output results at specified intervals
    IF (MOD(step, output_interval) == 0) THEN
      PRINT *, 'Step: ', step, ' Time: ', current_time
      CALL output_results(state, sg)
    END IF
  END DO

  PRINT *, 'Simulation complete. Final time: ', current_time

CONTAINS

  ! Subroutine to output state variables at each output step
  SUBROUTINE output_results(state, sg)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    TYPE(State_t), INTENT(IN) :: state(:)
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    ! Access 3D arrays but print a specific 2D slice
    PRINT *, 'Vorticity: ', state(1)%fields(1)%DATA(1, 1, 1)
    PRINT *, 'Divergence: ', state(1)%fields(2)%DATA(1, 1, 1)
    PRINT *, 'Height: ', state(1)%fields(3)%DATA(1, 1, 1)
  END SUBROUTINE output_results

END PROGRAM TestForwardModel
