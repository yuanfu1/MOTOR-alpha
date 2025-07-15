!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application Test for BEC Generated through Ensembles online.
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie, Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2022-04-22, created by Yuanfu Xie;
!                     2022-05-16, modified by Jiongming Pang.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/04/22, @GBA-MWF, Shenzhen
!   Midified by Jiongming Pang (pang.j.m@hotmail.com), 2022/05/16, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of MOTOR-DA with surfaces obs analysis for EnLoc based on CMA-GD.
!
PROGRAM App_CMA_GD_HybInc

  USE Applications_HybInc_m

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
  TYPE(Applications_HybInc_t) :: app
  TYPE(Applications_HybInc_t), ALLOCATABLE :: appB

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

  CALL app%initial(configFile, 2, istatus)
  IF (istatus .EQ. 0) THEN

    ! Get background:
    CALL app%backgrd(app%mgStart, app%mgEnd, unitTest)
    PRINT *, 'APP - Background fields have been processed successively.'

    ! Get ensemble background:
    CALL app%backgrdEns()
    PRINT *, 'APP - Ensemble background fields have been processed successively.'

    ! Get observations:
    CALL app%readObs(obsList, app%mgStart, app%mgEnd, unitTest)
    PRINT *, 'APP - Observations have been processed successively.'

    ! MOTOR-DA analysis:
    ! Note: (mgStart, mgEnd) can be different from backgrd
    !       but must be in between.
    CALL app%analyss(obsList, app%mgStart, app%mgEnd)
    PRINT *, 'Analysis has been successive.'

    ! ! Unit test:
    ! CALL app%verifyA(mgEnd)
  END IF

  CALL app%destroy

END PROGRAM App_CMA_GD_HybInc
