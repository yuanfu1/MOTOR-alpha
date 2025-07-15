!!--------------------------------------------------------------------------------------------------
! PROJECT           : Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/05/29, @GBA-MWF, Shenzhen
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2020/12/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Grid Generation Test Driver
!!
!!  *This main driver tests icos triangle grid generation*
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!! @test This is a beta testing for converting the ml_itmesh file to new format.
!!
!!===================================================================
PROGRAM test_gridGen

  USE mgGen_m, ONLY: mgGen_t, mgGenLatlon_t, mgGenIcosTr_t

  IMPLICIT NONE

  TYPE(mgGenIcosTr_t) :: icos
  LOGICAL :: PASS
  CHARACTER(LEN=1024) :: configFile

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/template.yaml"
  PASS = .TRUE.

  PRINT *, ''
  PRINT *, 'Calling icos geoGrid...'
  CALL icos%mgGen_geoGrid(configFile)

  IF (PASS) THEN
    WRITE (*, *) 'Test passed'
  ELSE
    WRITE (*, *) 'Test failed'
  END IF

END PROGRAM test_gridGen
