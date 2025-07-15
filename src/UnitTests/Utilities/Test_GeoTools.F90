!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/12, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM Test_GeoTools
  USE geoTools_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m
  USE conversions_m

  REAL(r_kind) :: locSrc(2), &        !> in lat/lon
                  locDst(2)           !> in lat/lon

  locSrc(1) = 20.0D0 * degree2radian
  locSrc(2) = 100.0D0 * degree2radian

  locDst(1) = 20.0D0 * degree2radian
  locDst(2) = 101.0D0 * degree2radian

  PRINT *, 'Distance between the two is: ', distance(locSrc, locDst)

  BLOCK
    REAL(r_kind) :: u, v, windSp, windDir

    u = 10
    v = 0

    CALL uv_to_wind(u, v, windSp, windDir)
    PRINT *, 'windSp windDir', windSp, windDir
  END BLOCK

  ! 1.111987663689131e+05
END PROGRAM Test_GeoTools
