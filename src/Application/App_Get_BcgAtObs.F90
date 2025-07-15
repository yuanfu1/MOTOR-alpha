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
PROGRAM App_Get_BcgAtObs
  USE Applications_m
  USE YAMLRead_m
  USE ObsRaw_m, ONLY: ObsRaw_t

  IMPLICIT NONE

  ! BLOCK
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  CHARACTER(LEN=8), PARAMETER :: &
    ! obsList(3) = (/'surfaces','sounding','cloudWnd'/) ! Individual test of cons
    ! obsList(2) = (/'surfaces','cloudWnd'/) ! Individual test of cons
    ! obsList(3) = (/'radarVel', 'sounding', 'profiler'/) ! Individual test of cons
    !  obsList(5) = (/'surfaces', 'sounding', 'profiler', 'radarVel', 'radarRef'/)
    obsList(4) = (/'surfaces', 'sounding', 'profiler', 'radarVel'/)
  ! obsList(7) = (/'surfaces', 'sounding', 'profiler', 'radarVel', 'shipRpts','buoyance','pilotRpt'/)

  !obsList(1) = (/'sounding'/)     ! Individual test of reflectivity
  INTEGER(i_kind) :: istatus, i
  LOGICAL :: unitTest
  TYPE(Applications_t) :: app
  REAL(r_kind) :: t1, t2

  ! Calculate CPU time:
  CALL CPU_TIME(t1)

  ! Get the configuration file
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    PRINT *, TRIM(configFile)

    CALL getarg(0, configFile)
    PRINT *, TRIM(configFile)
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/App_CMA_GD.exe configFile', /, &
           ' Check your configure file and rerun!')
    configFile = './App_3DVarVerification.yaml'
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  unitTest = .FALSE.

  CALL app%initial(configFile, 2, istatus)
  ! Testing by Yuanfu Xie 2022-10-02: turning the analysis off
  ! IF (.FALSE. .AND. istatus .EQ. 0) THEN
  IF (istatus .EQ. 0) THEN

    ! Get background:
    CALL app%backgrd(app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 3)
3   FORMAT('Background fields have been processed successively')

    ! Get observations:
    CALL app%readObs(obsList, app%mgStart, app%mgEnd, unitTest)
    WRITE (*, 4)
4   FORMAT('Observations have been processed successively')

    DO i = 1, 11
      BLOCK
        TYPE(State_t) :: X, XbRef
        TYPE(ObsRaw_t) :: rawObs, bcgObs
        INTEGER :: sgNum

        sgNum = 7

        XbRef = app%Ctl2State%fwdNL_opr(app%XbMG(sgNum))

        PRINT *, 'Before GetBcgAtObs.'
        CALL app%vel%radars(i)%ObsPrepareForSg(XbRef)
        CALL app%vel%radars(i)%GetBcgAtObs(XbRef, rawObs, bcgObs)
        PRINT *, 'After GetBcgAtObs.'
        istatus = yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile)

        CALL rawObs%outputNCForTest(TRIM(ncOutputFile)//'/Radar_RAW')
        CALL bcgObs%outputNCForTest(TRIM(ncOutputFile)//'/Radar_BCG')
      END BLOCK
    END DO
  END IF

  CALL CPU_TIME(t2)
  WRITE (*, 6) t2 - t1, app%mpddGlob%myrank
6 FORMAT('Application: time used in total: ', D12.4, ' at proc: ', I3)
  CALL app%destroy
  ! END BLOCK
END PROGRAM
