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
PROGRAM Test_AdvTime
  USE AdvanceTime_m
  USE kinds_m, ONLY: i_kind, r_kind

  INTEGER uTime_sec
  INTEGER gTime(6)

  uTime_sec = 100000

  CALL Time_Unix_to_GMT(uTime_sec, gTime)

  PRINT *, 'gTime', gTime

  PRINT *, 'uTime_sec input: ', uTime_sec

  CALL Time_GMT_to_Unix(gTime, uTime_sec)

  PRINT *, 'uTime_sec ouput: ', uTime_sec

END PROGRAM
