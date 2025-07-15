!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

PROGRAM abc
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE kinds_m, ONLY: i_kind, r_kind, r_double, i_short, i_long, i_llong, i_byte

  ! type(NcDataset) :: nc
  ! type(NcVariable)  :: var
  ! INTEGER(i_byte), ALLOCATABLE :: data1(:,:)

  ! nc = NcDataset('/public/home/simi/optest/3DVarVerification/srf/220527/220527_0000_3km/input/obs/synop/20220526_2320', "r")
  ! var = nc%getVariable("stationId"); CALL var%getData(data1)

  ! print*, data1

  ! CALL nc%close()
  CHARACTER(len=10) :: ccc = 'abcd'

  PRINT *, 'Is: ', ccc .EQ. 'abc'

END PROGRAM
