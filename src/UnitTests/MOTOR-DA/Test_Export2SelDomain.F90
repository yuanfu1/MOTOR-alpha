!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
PROGRAM Test_Export2SelDomain
  USE Export2SelDomain_m, ONLY: Export2SelDomain_t
  IMPLICIT NONE

  TYPE(Export2SelDomain_t) :: Export2SelDomain
  CHARACTER(len=1024) :: configFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/Application/App_SurfaceAnalysis.yaml"

  ! Export2SelDomain = Export2SelDomain_t(configFile)
END PROGRAM
