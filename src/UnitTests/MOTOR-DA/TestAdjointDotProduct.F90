PROGRAM TestAdjointDotProduct
  USE kinds_m, ONLY: r_kind, i_kind
  USE ShallowWaterForwardModel_m, ONLY: initialize_objects, time_step_rk4
  USE ShallowWaterAdjointModel_m, ONLY: time_step_rk4_AD
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE cDefault_m, ONLY: cDefault_t
  IMPLICIT NONE

  ! Declare variables
  TYPE(State_t), ALLOCATABLE :: state_forward(:), state_adjoint(:)
  TYPE(SingleGrid_t) :: sg                    ! Grid structure object
  TYPE(mpddGlob_t) :: mpddGlob                ! Global model parameters
  TYPE(cDefault_t) :: control                 ! Control object for initialization

  REAL(r_kind), ALLOCATABLE :: perturbation(:, :)  ! Perturbation array
  REAL(r_kind) :: dot_forward, dot_adjoint
  INTEGER(i_kind) :: mgStart, mgEnd           ! Multigrid start and end levels
  REAL(r_kind) :: dt                          ! Time step size
  REAL(r_kind) :: total_time                  ! Total simulation time
  INTEGER(i_kind) :: nSteps                   ! Number of steps
  INTEGER :: nx, ny                           ! Grid dimensions
  REAL(r_kind) :: epsilon                     ! Perturbation size
  CHARACTER(LEN=256) :: yaml_file             ! YAML configuration file path
  INTEGER :: i, j, step, iostat

  ! Set YAML configuration file path (modify this path to your YAML file location)
  yaml_file = '~/project/motor_case_ggf3km_2024040500/App_3DVarVerification.yaml'

  ! Parse YAML file to set parameters
  CALL ParseYAML(yaml_file, dt, total_time, nx, ny, epsilon)

  ! Initialize control objects based on parsed parameters
  mgStart = 1
  mgEnd = 1  ! Assuming single-level grid for simplicity
  ALLOCATE (state_forward(mgStart:mgEnd), state_adjoint(mgStart:mgEnd))

  CALL control%initialize_s(mgStart, mgEnd, state_forward)
  CALL initialize_objects(mpddGlob, mgStart, mgEnd, state_forward, 1)

  ! Allocate arrays for perturbation and initialize it
  ALLOCATE (perturbation(nx, ny))
  perturbation = epsilon  ! Initialize with the perturbation value from YAML

  ! Number of steps for the forward and adjoint time integration
  nSteps = INT(total_time / dt)

  ! ---- Forward Model ----
  ! Run the forward model and store the state
  DO step = 1, nSteps
    CALL time_step_rk4(dt, state_forward(1)%fields(1)%DATA(:, :, 1), &
                       state_forward(1)%fields(2)%DATA(:, :, 1), &
                       state_forward(1)%fields(3)%DATA(:, :, 1), poissonSolver, sg)
  END DO

  ! ---- Perturbation ----
  ! Apply perturbation to the forward state
  state_adjoint(1)%fields(1)%DATA(:, :, 1) = state_forward(1)%fields(1)%DATA(:, :, 1) + perturbation
  state_adjoint(1)%fields(2)%DATA(:, :, 1) = state_forward(1)%fields(2)%DATA(:, :, 1) + perturbation
  state_adjoint(1)%fields(3)%DATA(:, :, 1) = state_forward(1)%fields(3)%DATA(:, :, 1) + perturbation

  ! ---- Adjoint Model ----
  ! Run the adjoint model with the perturbed state
  DO step = nSteps, 1, -1
    CALL time_step_rk4_AD(dt, state_adjoint(1)%fields(1)%DATA(:, :, 1), &
                          state_adjoint(1)%fields(2)%DATA(:, :, 1), &
                          state_adjoint(1)%fields(3)%DATA(:, :, 1), poissonSolver, sg)
  END DO

  ! ---- Dot Product Calculation ----
  ! Calculate the dot product between forward state and perturbation
  dot_forward = 0.0_R_KIND
  dot_adjoint = 0.0_R_KIND
  DO i = 1, nx
    DO j = 1, ny
      dot_forward = dot_forward + state_forward(1)%fields(1)%DATA(i, j, 1) * perturbation(i, j) &
                    + state_forward(1)%fields(2)%DATA(i, j, 1) * perturbation(i, j) &
                    + state_forward(1)%fields(3)%DATA(i, j, 1) * perturbation(i, j)

      dot_adjoint = dot_adjoint + state_adjoint(1)%fields(1)%DATA(i, j, 1) * perturbation(i, j) &
                    + state_adjoint(1)%fields(2)%DATA(i, j, 1) * perturbation(i, j) &
                    + state_adjoint(1)%fields(3)%DATA(i, j, 1) * perturbation(i, j)
    END DO
  END DO

  ! Output the results
  PRINT *, "Dot product (Forward): ", dot_forward
  PRINT *, "Dot product (Adjoint): ", dot_adjoint

  ! Compare the results
  IF (ABS(dot_forward - dot_adjoint) < 1.0E-6_R_KIND) THEN
    PRINT *, "Test Passed: Forward and Adjoint dot products are consistent."
  ELSE
    PRINT *, "Test Failed: Forward and Adjoint dot products are not consistent."
  END IF

CONTAINS

  ! Subroutine to parse the YAML configuration file
  SUBROUTINE ParseYAML(file, dt, total_time, nx, ny, epsilon)
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(r_kind), INTENT(OUT) :: dt, total_time, epsilon
    INTEGER, INTENT(OUT) :: nx, ny
    CHARACTER(LEN=256) :: line
    INTEGER :: iostat

    ! Open the YAML file
    OPEN (UNIT=10, FILE=file, STATUS='OLD', IOSTAT=iostat)
    IF (iostat /= 0) THEN
      PRINT *, "Error opening YAML file!"
      STOP
    END IF

    ! Read the YAML file line by line and extract parameters
    DO
      READ (10, '(A)', IOSTAT=iostat) line
      IF (iostat /= 0) EXIT

      ! Parse time step (dt)
      IF (INDEX(line, 'dt') > 0) READ (line, '(A, F10.5)') dt

      ! Parse total simulation time (total_time)
      IF (INDEX(line, 'total_time') > 0) READ (line, '(A, F10.5)') total_time

      ! Parse grid dimensions (nx, ny)
      IF (INDEX(line, 'nx') > 0) READ (line, '(A, I5)') nx
      IF (INDEX(line, 'ny') > 0) READ (line, '(A, I5)') ny

      ! Parse perturbation size (epsilon)
      IF (INDEX(line, 'epsilon') > 0) READ (line, '(A, F10.5)') epsilon
    END DO

    ! Close the YAML file
    CLOSE (UNIT=10)
  END SUBROUTINE ParseYAML

END PROGRAM TestAdjointDotProduct
