!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application Test for BEC Generated through Ensembles online.
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-04-04, created by Jiongming Pang.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-04-04, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of MOTOR-DA with surfaces obs analysis for EnLoc based on CMA-GD.
!
PROGRAM App_CMA_GD_EnLoc_BEC

  USE Applications_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile
  ! CHARACTER(LEN=200), PARAMETER :: namelist = 'Application/App_CMA_GD_EnLoc.yaml'
  ! INTEGER(i_kind),    PARAMETER :: mgStart = 3, mgEnd = 6

  CHARACTER(LEN=8), PARAMETER :: &
    ! obsList(3) = (/'surfaces', 'sounding', 'profiler'/) ! Individual test of cons
    ! obsList(4) = (/'surfaces', 'sounding', 'profiler','radarRef'/)
    ! obsList(1) = (/'radarRef'/)     ! Individual test of reflectivity
    obsList(1) = (/'surfaces'/) !surfaces, shipRpts, sounding
  INTEGER(i_kind) :: istatus
  LOGICAL :: unitTest
  TYPE(Applications_t) :: app
  TYPE(Applications_t), ALLOCATABLE :: appB

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

  unitTest = .FALSE. !.TRUE.

  ! Calculate BEC
  ALLOCATE (appB)
  CALL appB%initial(configFile, 2, istatus)
  IF (istatus .EQ. 0) THEN
    CALL appB%calEnsB()
    PRINT *, 'APP - Calculation of BEC has been processed successively.'

    ! CALL appB%destroy
  END IF
  DEALLOCATE (appB)

END PROGRAM App_CMA_GD_EnLoc_BEC
