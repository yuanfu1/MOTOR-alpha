!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2023/6/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
MODULE module_io_grapes_m
  USE module_domain, ONLY: domain_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE NameValueMap_m
  IMPLICIT NONE

  TYPE :: module_io_grapes_t
    CHARACTER(LEN=1024)       :: giFileName
    TYPE(NameValueMap_t)      :: valName2FileDisp

  CONTAINS
    PROCEDURE, PUBLIC, NOPASS :: init_grapes_map_modes
    PROCEDURE, PUBLIC, NOPASS :: input_grapes_sequential
    PROCEDURE, PUBLIC, NOPASS :: input_grapes_hydrometeors
    PROCEDURE, PUBLIC :: read_value_1d
    PROCEDURE, PUBLIC :: read_value_2d
    PROCEDURE, PUBLIC :: read_value_3d
    PROCEDURE, PUBLIC :: write_value_1d
    PROCEDURE, PUBLIC :: write_value_2d
    PROCEDURE, PUBLIC :: write_value_3d
    PROCEDURE, PUBLIC :: initialize
  END TYPE

CONTAINS

  SUBROUTINE initialize(this, hgrid, nlFileName, giFileName)
    CLASS(module_io_grapes_t)   :: this
    TYPE(domain_t), TARGET      :: hgrid
    CHARACTER(len=*)            :: giFileName, nlFileName

    this%giFileName = giFileName
  END SUBROUTINE

  SUBROUTINE input_grapes_hydrometeors(HydroFileName, grid)
    USE module_configure
    IMPLICIT NONE
    CHARACTER(*)            :: HydroFileName
    TYPE(domain_t)          :: grid

  END SUBROUTINE input_grapes_hydrometeors
  
  SUBROUTINE input_grapes_sequential(nlFileName, inputfn, grid)
    CHARACTER*(*)           :: inputfn, nlFileName
    TYPE(domain_t)          :: grid
  END SUBROUTINE

  FUNCTION init_grapes_map_modes(nlFileName, inputfn, grid, mode)
    USE module_configure
    USE NameValueMap_m
    IMPLICIT NONE
    CHARACTER*(*)           :: inputfn, nlFileName, mode
    TYPE(domain_t)         :: grid
    TYPE(NameValueMap_t)  :: init_grapes_map_modes
  END FUNCTION

  SUBROUTINE read_value_3d(this, vname, res, hgrid)
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

  SUBROUTINE read_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

  SUBROUTINE read_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

  SUBROUTINE write_value_3d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :, :) :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

  SUBROUTINE write_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :) :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

  SUBROUTINE write_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_t) :: this
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:) :: res
    TYPE(domain_t)          :: hgrid
  END SUBROUTINE

END MODULE
