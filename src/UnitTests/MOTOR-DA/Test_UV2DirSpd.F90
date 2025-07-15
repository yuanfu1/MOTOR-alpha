PROGRAM Test_UV2DirSpd
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE State_m, ONLY: State_t
  USE UV2DirSpd_m, ONLY: UV2DirSpd_t

  IMPLICIT NONE

  TYPE(ObsSet_t) :: Y
  TYPE(State_t) :: X, dX
  TYPE(UV2DirSpd_t) :: UV2DirSpd
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geometry

  CHARACTER(LEN=1024) :: configFile

  ! Get the configuration file
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/App_CMA_GD.exe configFile', /, &
           ' Check your configure file and rerun!')
    STOP
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  ! Auxtypes
  CALL mpddGlob%initialize()

  CALL geometry%initialize(configFile, mpddGlob)

  CALL mpddGlob%finalize()

END PROGRAM Test_UV2DirSpd
