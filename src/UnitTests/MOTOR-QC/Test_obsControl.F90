! Test observation control subroutine
! Created by Ting Shu @ 20230208
PROGRAM test_obsControl
  USE obsTools_m
  USE ObsSet_m, ONLY: ObsSet_t
  USE State_m, ONLY: State_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE mpObs_m, ONLY: mpObs_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: staticDir, configFile
  TYPE(ObsSet_t) :: my_thinObs
  TYPE(State_t) :: state
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(mpObs_t), TARGET, ALLOCATABLE :: mpObs(:)
  TYPE(SingleGrid_t), POINTER :: sgFine

  INTEGER(i_kind) :: mgStart = 1, mgEnd = 5, istatus, i

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest/motor-qc_test.yaml"

  ! Initializer
  CALL mpddGlob%initialize()

  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)

  ! Allocate mpObs:
  ALLOCATE (mpObs(mgStart:mgEnd))

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)

  DO i = mgStart, mgEnd
    CALL mpObs(i)%initializeMPObs(geometry%mg%sg(i))
  END DO

  sgFine => geometry%mg%sg(mgEnd - 1)
  CALL state%initialize(configfile, sgFine)

  CALL ObsControl(configFile, my_thinObs, state)
END PROGRAM test_obsControl
