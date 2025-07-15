MODULE namelist_ps_adj_m

  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

CONTAINS

  SUBROUTINE namelist_ps_adj(configFile, solver, nCycle, nIterPre, nIterPost, nRelax, omegas, max_band)
    ! Model namelist variables:
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(LEN=3), INTENT(OUT) :: solver
    INTEGER(i_kind), INTENT(OUT) :: &
      nCycle, &
      nIterPre, &
      nIterPost, &
      nRelax, &
      max_band
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: omegas(:)
    INTEGER(i_kind) :: istatus

    ! Initialize output values in case of failure
    solver = '   '
    nCycle = -1
    nIterPre = -1
    nIterPost = -1
    nRelax = -1
    max_band = -1

    ! Read the solver configurations from the YAML file
    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'solver', solver)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading solver from config file.'
      ERROR STOP 'Failed to read solver.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'nCycle', nCycle)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading nCycle from config file.'
      ERROR STOP 'Failed to read nCycle.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'nRelax', nRelax)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading nRelax from config file.'
      ERROR STOP 'Failed to read nRelax.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'omegas', omegas)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading omegas from config file.'
      ERROR STOP 'Failed to read omegas.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'max_band', max_band)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading max_band from config file.'
      ERROR STOP 'Failed to read max_band.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'nIterPre', nIterPre)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading nIterPre from config file.'
      ERROR STOP 'Failed to read nIterPre.'
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'poissonSolver', 'nIterPost', nIterPost)
    IF (istatus /= 0) THEN
      PRINT *, 'Error reading nIterPost from config file.'
      ERROR STOP 'Failed to read nIterPost.'
    END IF

  END SUBROUTINE namelist_ps_adj

END MODULE namelist_ps_adj_m
