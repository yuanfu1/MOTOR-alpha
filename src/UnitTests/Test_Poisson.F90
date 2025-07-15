!!--------------------------------------------------------------------------------------------------
! PROJECT           : Unitest.testMixedSolver.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie on 2024-12-09 for showing a pointer initialization by a constructor failed
!!--------------------------------------------------------------------------------------------------
PROGRAM testMixedSolver
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE mixedPsiChiSolver_Dirichlet_m, ONLY: mixedPsiChiSolver_Dirichlet_t

  IMPLICIT NONE

  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  TYPE(poissonSolver_t) :: ps
  TYPE(mixedPsiChiSolver_Dirichlet_t) :: mix

  CHARACTER(LEN=1024) :: yamlFile
  INTEGER(c_int) :: mgStart, mgEnd   ! mgEnd: Glevel solution is solved

  ! Get the configuration file
  CALL getarg(1, yamlFile)

  ! Set up MPDD and Geometry:
  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile,mpddGlob)
  mgStart = 2; mgEnd = 2
  PRINT*,'initializing mix...'
  ! CALL mix%initialize_s(mgStart,mgEnd,geo)
  CALL initializePoissonSolver(mgStart,mgEnd,geo)
  ! ps = poissonSolver_t(yamlFile, geo) 

  CALL mpddGlob%finalize()

  PRINT*,'Normal exit'
END PROGRAM testMixedSolver

SUBROUTINE initializePoissonSolver(mgStart,mgEnd,geo)
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  INTEGER(c_int), INTENT(IN) :: mgStart, mgEnd
  TYPE(geometry_t), INTENT(IN) :: geo

  CHARACTER(LEN=1024) :: yamlFile
  TYPE(poissonSolver_t), POINTER :: ps  ! POINTER made this initialization failed but removing it works

  ! Get the configuration file
  CALL getarg(1, yamlFile)

  ! Initialize the poissonSolver
  ps = poissonSolver_t(yamlFile, geo) 
  PRINT*,'Poisson solver is initiated: ',TRIM(yamlFile)
END SUBROUTINE initializePoissonSolver