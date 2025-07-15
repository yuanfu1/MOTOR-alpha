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
PROGRAM App_SurfaceAnalysis
  USE Applications_m
  USE YAMLRead_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile

  CHARACTER(LEN=8), PARAMETER :: &
    obsList(1) = (/'surfaces'/) ! Individual test of cons
  ! obsList(4) = (/'surfaces', 'sounding', 'profiler','radarRef'/)
  ! obsList(1) = (/'radarRef'/)     ! Individual test of reflectivity
  INTEGER(i_kind) :: istatus
  LOGICAL :: unitTest
  TYPE(Applications_t) :: app
  REAL(r_kind) :: t1, t2

  ! Get the configuration file
  CALL getarg(1, configFile)
  ! print *, 'arguement: ', configFile
  IF (TRIM(configFile) .EQ. '') THEN
    WRITE (*, 4)
4   FORMAT('Usage of this driver: mpirun -n <n> Debug/App_CMA_GD.exe configFile', /, &
           ' Check your configure file and rerun!')
    STOP
  ELSE
    WRITE (*, 5) TRIM(configFile)
5   FORMAT('ConfigFile is: ', A)
  END IF

  unitTest = .FALSE.

  CALL CPU_TIME(t1)

  CALL app%initial(configFile, 2, istatus)
  IF (istatus .EQ. 0) THEN

    ! Get background:
    CALL app%backgrd(app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 1)
1   FORMAT('Background fields have been processed successively')

    ! Get observations:
    CALL app%readObs(obsList, app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 2)
2   FORMAT('Observations have been processed successively')

    ! MOTOR-DA analysis:
    ! Note: (mgStart, mgEnd) can be different from backgrd
    !   but must be in between.
    CALL app%analyss(obsList, app%mgStart, app%mgEnd)
    WRITE (*, 3)
3   FORMAT('Analysis has been successive')

    ! Dump Background Fields
    CALL app%dumpAna(unitTest)

    ! Analysis verification:
    ! PRINT *, 'Verifying??? ', MAXVAL(ABS(app%verifyLatlon)), MAXVAL(ABS(app%verifyTime))
    ! IF (MAXVAL(ABS(app%verifyLatlon)) .NE. 0.0D0 .AND. MAXVAL(ABS(app%verifyTime)) .NE. 0.0D0) &
      ! CALL app%verifyA(obsList, app%mgEnd, unitTest)
  END IF

  CALL CPU_TIME(t2)
  IF (app%mpddGlob%isBaseProc()) PRINT *, 'Total time cose:', t2 - t1, 's'

  CALL app%destroy
END PROGRAM App_SurfaceAnalysis
