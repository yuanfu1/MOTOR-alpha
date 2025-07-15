!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application Test for BEC Generated through Ensembles online.
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2022-04-22, created by Yuanfu Xie;
!                     2022-05-16, modified by Jiongming Pang.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/04/22, @GBA-MWF, Shenzhen
!   Midified by Jiongming Pang (pang.j.m@hotmail.com), 2022/05/16, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is an APP for the calculation of BEC through ensemble method.
!
PROGRAM App_CalculateBEC_EnLoc

  USE Applications_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile
  ! CHARACTER(LEN=200), PARAMETER :: namelist = 'Application/App_CMA_GD_EnLoc.nl'
  ! CHARACTER(LEN=200), PARAMETER :: namelist = 'Application/App_CMA_GD_EnLoc.yaml'
  ! INTEGER(i_kind),    PARAMETER :: mgStart = 3, mgEnd = 6

  CHARACTER(LEN=8), PARAMETER :: &
    ! obsList(3) = (/'surfaces', 'sounding', 'profiler'/) ! Individual test of cons
    ! obsList(4) = (/'surfaces', 'sounding', 'profiler','radarRef'/)
    ! obsList(1) = (/'radarRef'/)     ! Individual test of reflectivity
    obsList(1) = (/'surfaces'/)
  INTEGER(i_kind) :: istatus
  LOGICAL :: unitTest
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
  PRINT *, 'BECsolver: ', TRIM(appB%BECsolver)
  IF (istatus .EQ. 0 .AND. TRIM(appB%BECsolver) == 'Ensemble') THEN
    CALL appB%calEnsB()
    PRINT *, 'APP - Calculation of BEC has been processed successively.'
  ELSE
    PRINT *, ">>>>>>> ERROR! 'CalculateBEC' should use 'Ensemble' as the BECsolver."
    PRINT *, ">>>>>>> Please check the setting of BMat/BECsolver in your yaml file!"
    DEALLOCATE (appB)
    STOP
  END IF

  DEALLOCATE (appB)

END PROGRAM App_CalculateBEC_EnLoc
