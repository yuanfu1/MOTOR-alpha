!!--------------------------------------------------------------------------------------------------
! PROJECT           : Utilities/Utility Test for netCDFReader.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2023-08-15, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/08/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of netCDFReader.F90 extracting data from a netCDF file using a ncTable file
!
PROGRAM Test_netCDFReader
  USE kinds_m, ONLY: i_kind, r_kind
  USE netCDFTable_m, ONLY: netCDFTable_t
  USE netCDFReader_m, ONLY: netCDFReader_t
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup

  CHARACTER(LEN=1024) :: filename
  TYPE(netCDFTable_t) :: table
  TYPE(netCDFReader_t) :: reader

  ! Data:
  REAL(r_kind) :: domain(2, 2), missingValue
  REAL(r_kind), ALLOCATABLE :: ncData(:, :), latlonhgttim(:, :)

  domain(1, 1) = 16.0D0; domain(2, 1) = 31.36
  domain(1, 2) = 96.0D0; domain(2, 2) = 123.36
  missingValue = 99999999.0D0

  filename = "../../cdw20220527/FY4A-_AGRI--_N_DISK_1047E_L2-_AMV-_C009_NUL_20220527000000_20220527001459_064KM_V0001.NC"
  CALL table%getTable("template_ncTable.yaml")
  CALL reader%ncReader(table, TRIM(filename), domain, missingValue, ncData, latlonhgttim)
END PROGRAM Test_netCDFReader
