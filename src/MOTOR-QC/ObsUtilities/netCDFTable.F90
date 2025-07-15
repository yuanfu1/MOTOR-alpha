!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-Utilities/Utility.netCDFTable
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
! This is a netCDF table module designed for various data types defined in a yaml file.
MODULE netCDFTable_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

  INTEGER(i_kind), PARAMETER :: VarNameLen = 20

  TYPE :: netCDFTable_t
    CHARACTER(LEN=VarNameLen) :: Instrument, DataSource
    CHARACTER(LEN=VarNameLen), ALLOCATABLE :: VarNames(:), &
                                              latlonNames(:), timeName(:), heightName(:)
    INTEGER(i_kind), ALLOCATABLE :: required(:)
    REAL(r_kind) :: missingValue
  CONTAINS
    PROCEDURE :: getTable
  END TYPE netCDFTable_t

CONTAINS

  ! To read the table content
  SUBROUTINE getTable(this, tableFile)

    IMPLICIT NONE

    CLASS(netCDFTable_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: tableFile

    ! Local variables:
    CHARACTER(LEN=1024) :: path
    INTEGER(i_kind) :: i, istatus

    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", path)
    path = TRIM(path)//'/ncTables/'//TRIM(tableFile)

    istatus = yaml_get_var(TRIM(path), 'data_source', 'instrument', this%Instrument)
    istatus = yaml_get_var(TRIM(path), 'data_source', 'dataSource', this%DataSource)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'latlon', this%latlonNames)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'variables', this%VarNames)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'required', this%required)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'missingValue', this%missingValue)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'time', this%timeName)
    istatus = yaml_get_var(TRIM(path), 'data_content', 'height', this%heightName)
    WRITE (*, 5)
5   FORMAT('+----------------------------- netCDF Table -----------------------------+')
    WRITE (*, 1) TRIM(this%Instrument), TRIM(this%DataSource)
1   FORMAT('| Table read - Instrument: ', A, ' | DataSource: ', A)
    WRITE (*, 2) (TRIM(this%VarNames(i)), this%required(i), i=1, UBOUND(this%VarNames, 1))
2   FORMAT('| Table read - VarNames:', 50(' | ', A, ',', I1))
    WRITE (*, 3) this%missingValue
3   FORMAT('| Table read - missingValue: ', E14.6)
    WRITE (*, 4) TRIM(this%timeName(1)), TRIM(this%heightName(1))
4   FORMAT('| Table read - time: ', A, ' | height: ', A)
    WRITE (*, 6)
6   FORMAT('+----------------------------- End of Table -----------------------------+')
  END SUBROUTINE getTable
END MODULE netCDFTable_m
