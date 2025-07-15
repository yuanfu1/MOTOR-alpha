PROGRAM test_PS_Adjoint
  USE kinds_m, ONLY: i_kind, r_kind
  USE poissonSolver_adjoint_m, ONLY: poissonSolver_adj_t
  USE YAMLRead_m
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: geometry_t

  INCLUDE "mpif.h"

  TYPE(geometry_t), TARGET :: geometry
  TYPE(GenContainers_t) :: GenContainers
  TYPE(poissonSolver_adj_t) :: ps_adj

  REAL(r_kind), ALLOCATABLE :: adjoint_rhs(:, :)      ! Container for adjoint RHS (right-hand side)
  REAL(r_kind), ALLOCATABLE :: adjoint_solution(:, :) ! Container for adjoint solution
  REAL(r_kind), ALLOCATABLE :: sigma_3d(:, :)         ! Sigma values for spherical geometry

  INTEGER(i_kind) :: vLevel
  INTEGER(i_kind) :: i
  REAL(r_kind) :: tempMax, tempMin, mx, mn, t1, t2
  CHARACTER(LEN=1024) :: configFile

  ! Get the configuration file
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"

  ! Read the vLevel variable from the configuration file
  IF (yaml_get_var(configFile, 'modelState', 'vLevel', vLevel, i) /= 0) STOP

  ! Initialize the GenContainers and geometry objects
  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry)
  ps_adj = poissonSolver_adj_t(configFile, geometry)   ! Initialize the adjoint Poisson solver

  ! Associate with the finest grid level in the geometry
  ASSOCIATE (sgFinest => geometry%mg%sg(6))

    ! Allocate memory for the adjoint right-hand side (RHS) and solution
    CALL sgFinest%allocateMat(vLevel, adjoint_rhs)
    CALL sgFinest%allocateMat(vLevel, adjoint_solution)
    CALL sgFinest%allocateMat(vLevel, sigma_3d)

    ! Initialize the sigma_3d array based on the grid sigma values
    DO i = 1, sgFinest%num_cell
      sigma_3d(:, i) = sgFinest%sigma
    END DO

    ! Set up the adjoint right-hand side based on the initial conditions
    ! In practice, this would come from the model's adjoint source terms or observations.
    adjoint_rhs = 0.0_R_KIND
    DO i = 1, sgFinest%num_icell
      ! Example initialization: setting RHS to some function (adjust as needed)
      IF (sgFinest%cell_type(i) .EQ. 0) THEN
        adjoint_rhs(:, i) = SIN(sgFinest%cell_cntr(1, i)) * COS(sgFinest%cell_cntr(2, i))
      ELSE
        adjoint_rhs(:, i) = 0.0_R_KIND  ! Boundary cells
      END IF
    END DO

    ! Solve the adjoint Poisson equation on the sphere
    CALL GenContainers%mpddGlob%barrier
    CALL CPU_TIME(t1)

    adjoint_solution = 0.0_R_KIND  ! Initialize the adjoint solution
    CALL ps_adj%PoissonSol_sphere_adjoint(sgFinest%gLevel, vLevel, adjoint_rhs, adjoint_solution)

    CALL GenContainers%mpddGlob%barrier
    CALL CPU_TIME(t2)

    IF (GenContainers%mpddGlob%isBaseProc()) PRINT *, '+----------------------------------+'
    IF (GenContainers%mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'
    IF (GenContainers%mpddGlob%isBaseProc()) PRINT *, '+----------------------------------+'

    ! Calculate and print max and min values of the adjoint solution
    tempMax = MAXVAL(adjoint_solution(1, 1:sgFinest%num_icell))
    tempMin = MINVAL(adjoint_solution(1, 1:sgFinest%num_icell))

    CALL MPI_REDUCE(tempMax, mx, 1, mpi_double_precision, mpi_max, GenContainers%mpddGlob%rankBase, GenContainers%mpddGlob%comm, GenContainers%mpddGlob%ierr)
    CALL MPI_REDUCE(tempMin, mn, 1, mpi_double_precision, mpi_min, GenContainers%mpddGlob%rankBase, GenContainers%mpddGlob%comm, GenContainers%mpddGlob%ierr)

    IF (GenContainers%mpddGlob%isBaseProc()) WRITE (*, 10) mx, mn
10  FORMAT('Max/min adjoint solution values: ', 2E30.22)

  END ASSOCIATE

  ! Deallocate memory
  DEALLOCATE (adjoint_rhs, adjoint_solution, sigma_3d)

  ! Test results
  IF (GenContainers%mpddGlob%isBaseProc()) THEN
    IF ((mx < 0.73E1) .AND. (mn > -0.73E1)) THEN
      PRINT *, 'Adjoint test passed!'
    ELSE
      PRINT *, 'Adjoint test failed!', mx, mn
    END IF
  END IF

  CALL GenContainers%mpddGlob%finalize

END PROGRAM test_PS_Adjoint
