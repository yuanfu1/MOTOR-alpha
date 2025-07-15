!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   2020-4 created by Yuanfu Xie
!   This file is reforged from fsw_namelist_m by Yuanfu Xie
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a module reading needed configs in namelist file.
MODULE NML_BMatrix_m
  USE kinds_m, ONLY: i_kind, r_kind
  !USE NMLRead_m, ONLY: namelist_read
  USE YAMLRead_m

CONTAINS
  SUBROUTINE NML_BMatrix(configFile, varList)
    ! Model namelist variables:
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(*), ALLOCATABLE :: varList(:)
    INTEGER(i_kind) :: ifile

    !CALL namelist_read(configFile, 'varList', varList)
    ifile = yaml_get_var(configFile, 'modelState', 'varList', varList)

  END SUBROUTINE NML_BMatrix

END MODULE NML_BMatrix_m
