!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/27, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention

MODULE configuration_m
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE configuration_t
    CHARACTER(LEN=256) :: path2static

  CONTAINS

  END TYPE configuration_t
CONTAINS

END MODULE configuration_m
