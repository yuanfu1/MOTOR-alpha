!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application Test for CMA_GD
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-04-22, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/04/22, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of MOTOR-DA with full observations analysis for initializing CMA_GD
!
PROGRAM App_CMA_GD
  USE Applications_m
  USE YAMLRead_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile

  CHARACTER(LEN=8), PARAMETER :: &
    ! obsList(1) = (/'radarRef'/) ! Individual test of cons
    ! obsList(3) = (/'surfaces', 'sounding', 'profiler'/) ! Individual test of cons
    obsList(4) = (/'surfaces', 'sounding', 'profiler', 'cloudWnd'/)
  ! obsList(1) = (/'profiler'/)     ! Individual test of reflectivity
  INTEGER(i_kind) :: istatus
  LOGICAL :: unitTest
  TYPE(Applications_t) :: app
  REAL(r_kind) :: t1, t2

  ! Calculate CPU time:
  CALL CPU_TIME(t1)

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

  unitTest = .FALSE.

  CALL app%initial(configFile, 2, istatus)
  IF (istatus .EQ. 0) THEN

    ! Get background:
    CALL app%backgrd(app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 3)
3   FORMAT('Background fields have been processed successively')

    ! Get observations:
    CALL app%readObs(obsList, app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 4)
4   FORMAT('Observations have been processed successively')

    ! MOTOR-DA analysis:
    ! Note: (mgStart, mgEnd) can be different from backgrd
    !   but must be in between.
    CALL app%analyss(obsList, app%mgStart, app%mgEnd)
    ! CALL app%analyze(obsList, app%mgStart, app%mgEnd)
    WRITE (*, 5)
5   FORMAT('Analysis has been successive')
    WRITE (*, 16)
16  FORMAT('Test passed')

    ! Unit test:
    ! IF (unitTest) CALL app%verifyA(obsList, app%mgEnd, unitTest)

    ! Dump Background Fields
    CALL app%dumpAna(unitTest)
  END IF

  CALL CPU_TIME(t2)
  WRITE (*, 6) t2 - t1, app%mpddGlob%myrank
6 FORMAT('App_CMA_GD: time used in total: ', D12.4, ' at proc: ', I3)

  CALL app%destroy
END PROGRAM App_CMA_GD
