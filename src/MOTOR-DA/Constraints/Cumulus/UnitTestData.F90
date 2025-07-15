!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

MODULE UnitTestData_m
  USE State_m, ONLY: State_t

CONTAINS

  SUBROUTINE UnitTestData(X)
    TYPE(State_t) :: X
    CHARACTER(LEN=1024) :: InDataFile, &
                           sigmainput, uwndinput, vwndinput, wwndinput, &
                           thetainput, qvaporinput, presinput, gphinput, &
                           topoinput !< Config and output file name.

    CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", InDataFile)
    sigmainput = TRIM(InDataFile)//"/CumUniTest/sigma.txt"
    uwndinput = TRIM(InDataFile)//"/CumUniTest/uwnd.txt"
    vwndinput = TRIM(InDataFile)//"/CumUniTest/vwnd.txt"
    wwndinput = TRIM(InDataFile)//"/CumUniTest/wwnd.txt"
    thetainput = TRIM(InDataFile)//"/CumUniTest/theta.txt"
    qvaporinput = TRIM(InDataFile)//"/CumUniTest/qvapor.txt"
    presinput = TRIM(InDataFile)//"/CumUniTest/pres.txt"
    gphinput = TRIM(InDataFile)//"/CumUniTest/gph.txt"
    topoinput = TRIM(InDataFile)//"/CumUniTest/topo.txt"

    OPEN (unit=10, file=sigmainput, status="old")
    DO i = 1, 57
      READ (10, *) X%sg%sigma(i)
    END DO
    CLOSE (10)

    OPEN (unit=11, file=uwndinput, status="old")
    OPEN (unit=12, file=vwndinput, status="old")
    OPEN (unit=13, file=wwndinput, status="old")
    OPEN (unit=14, file=thetainput, status="old")
    OPEN (unit=15, file=qvaporinput, status="old")
    OPEN (unit=16, file=presinput, status="old")
    OPEN (unit=17, file=gphinput, status="old")

    DO j = 1, 9
      ! print *, j
      DO i = 1, 57
!  print *, i
        READ (11, *) X%Fields(1)%DATA(i, j, 1)
        READ (12, *) X%Fields(2)%DATA(i, j, 1)
        READ (13, *) X%Fields(3)%DATA(i, j, 1)
        READ (14, *) X%Fields(4)%DATA(i, j, 1)
        READ (15, *) X%Fields(5)%DATA(i, j, 1)
        READ (16, *) X%Fields(6)%DATA(i, j, 1)
        READ (17, *) X%sg%zHght(i, j)
      END DO
    END DO
    CLOSE (11)
    CLOSE (12)
    CLOSE (13)
    CLOSE (14)
    CLOSE (15)
    CLOSE (16)
    CLOSE (17)

    OPEN (unit=18, file=topoinput, status="old")
    DO j = 1, 9
      READ (18, *) X%sg%topo(j)
    END DO
    CLOSE (18)

  END SUBROUTINE

END MODULE
