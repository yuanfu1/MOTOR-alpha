MODULE TimeIntegral_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE nwp_m, ONLY: nwp_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE YAMLRead_m
  USE CalPres_m, ONLY: CalPres_t

  TYPE, ABSTRACT:: rk4_t
    TYPE(State_t), ALLOCATABLE :: stepsHalf(:), stepsFull(:)
    TYPE(geometry_t), POINTER :: geometryHalf, geometryFull
    ! a Poisson solver:
    TYPE(poissonSolver_t) :: poisson
    TYPE(nwp_t) :: nwp
    INTEGER(i_kind) :: gLevel, numSteps, kt
    REAL(r_kind) :: dt
    CHARACTER(LEN=1024) :: configFile
  CONTAINS
    FINAL :: destructor
    PROCEDURE :: timeIntegrals => timeIntegrals_rk4
  END TYPE rk4_t

  TYPE() rkr3
