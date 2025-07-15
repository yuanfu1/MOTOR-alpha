!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application Unit Test
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-03-19, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/03/19, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a unit test of MOTOR-DA conventional observations analysis.
!!   Uses a module of a multiscale analytic function to test the multiscale capabity of MOTOR-DA
!! This is a unit test using real conventional observation stations replacing their observation values with
!!   the analytic test function values.
!
PROGRAM unitTest_Conventions

  USE Applications_m

  IMPLICIT NONE

  CHARACTER(LEN=200), PARAMETER :: NAMELIST = 'UnitTest/testConvention_MultiObs.yaml'
  INTEGER(i_kind), PARAMETER :: mgStart = 3, mgEnd = 6

  CHARACTER(LEN=8), PARAMETER :: obsList(3) = (/'surfaces', 'sounding', 'profiler'/)
  INTEGER(i_kind) :: istatus
  TYPE(Applications_t) :: app

  CALL app%initial(NAMELIST, 2, istatus)
  IF (istatus .EQ. 0) THEN

    ! Get background:
    CALL app%backgrd(mgStart, mgEnd, .TRUE.)

    ! Get observations:
    CALL app%readObs(obsList, mgStart, mgEnd, .TRUE.)

    ! MOTOR-DA analysis:
    ! Note: (mgStart, mgEnd) can be different from backgrd
    !   but must be in between.
    CALL app%analyss(obsList, mgStart, mgEnd)

    ! Unit test:
    CALL app%verifyA(obsList, mgEnd, .TRUE.)
  END IF

  CALL app%destroy
END PROGRAM unitTest_Conventions
