!!--------------------------------------------------------------------------------------------------
! PROJECT           : Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/05/29, @GBA-MWF, Shenzhen
!   Modified by ..... (???@gmail.com), YYYY/MM/DD, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! # Grid Generation Test Driver
!!
!!  *This main driver tests grid generation*
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!! @test This is a beta testing for generating lat-lon grid.
PROGRAM test_gridGen
  USE mgGen_m, ONLY: mgGen_t, mgGenLatlon_t, mgGenIcosTr_t

  IMPLICIT NONE

  TYPE(mgGenLatlon_t) :: grid

  ! Get the configFile
  LOGICAL :: PASS
  CHARACTER(LEN=1024) :: configFile

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/test.yaml"
  PASS = .TRUE.

  PRINT *, '+--------------------+'
  PRINT *, '|  Start testing...  |'
  PRINT *, '+--------------------+'
  WRITE (*, 10)
10 FORMAT('Calling geoGrid...')
  CALL grid%mgGen_geoGrid(configFile)

  WRITE (*, 11)
11 FORMAT('Calling precals...')
  CALL grid%mgGen_precals

  WRITE (*, 12)
12 FORMAT('Calling a unit test...')
  CALL unitTestLatLon(configFile, PASS)

  IF (PASS) THEN
    WRITE (*, 13)
13  FORMAT('Test passed')
  ELSE
    WRITE (*, 14)
14  FORMAT('Test failed')
  END IF

END PROGRAM test_gridGen
