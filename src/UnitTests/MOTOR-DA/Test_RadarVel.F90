!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
! Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/14, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_RadarVel
  USE OprRadarVel_m, ONLY: OprRadarVel_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE parameters_m

  REAL(r_kind) :: uWnd, vWnd, wWnd, locRadar(3), locGrid(3), factors(3)
  TYPE(OprRadarVel_t) :: OprRadarVel

  uWnd = 1.0D0
  vWnd = 1.0D0
  wWnd = 1.0D0

  !
  locRadar = (/20.0D0 * degree2radian, 100.0D0 * degree2radian, 0.0D0/)
  locGrid = (/20.01D0 * degree2radian, 100.00D0 * degree2radian, 0.0D0/)

  !
  factors = OprRadarVel%calRadarVelFactorElem(locRadar, locGrid)
  PRINT *, 'factors: ', factors, SUM(factors * factors)
END PROGRAM Test_RadarVel
