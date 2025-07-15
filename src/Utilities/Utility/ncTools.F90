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
MODULE ncTools_m
  USE NETCDF

CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check

END MODULE ncTools_m
