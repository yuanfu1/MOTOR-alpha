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
MODULE ObsRaw_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
  USE parameters_m, ONLY: degree2radian

  TYPE ObsRaw_t
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &
                                 obsErrs(:, :), &
                                 olatlon(:, :), &
                                 height(:), &
                                 tempor(:)
    CHARACTER(LEN=20), ALLOCATABLE :: varNames(:)      ! Obs variable names
    CHARACTER(LEN=20) :: obsType
  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: outputNCForTest
  END TYPE

CONTAINS

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(ObsRaw_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%obsData)) DEALLOCATE (this%obsData)
    IF (ALLOCATED(this%obsErrs)) DEALLOCATE (this%obsErrs)
    IF (ALLOCATED(this%olatlon)) DEALLOCATE (this%olatlon)
    IF (ALLOCATED(this%height)) DEALLOCATE (this%height)
    IF (ALLOCATED(this%tempor)) DEALLOCATE (this%tempor)
    IF (ALLOCATED(this%varNames)) DEALLOCATE (this%varNames)

  END SUBROUTINE destructor

  SUBROUTINE outputNCForTest(this, filename)
    IMPLICIT NONE
    CLASS(ObsRaw_t) :: this
    CHARACTER(*) :: filename

    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1
    TYPE(NcVariable)  :: var, scalar
    TYPE(NcGroup)     :: grp

    INTEGER(i_kind) :: recNum, i, varNum
! open a dataset in write mode
! args:
!     filename
!     mode ("w": write, "r": read-only, "a": read-write)
    nc = NcDataset(filename//TRIM(this%obsType)//'.nc', "w")

    recNum = SIZE(this%height)
    varNum = SIZE(this%varNames)
    dim1 = nc%setDimension("recNum", recNum)

    var = nc%setVariable("lat", "f32", (/dim1/)); CALL var%setData(this%olatlon(1, :) / degree2radian)
    var = nc%setVariable("lon", "f32", (/dim1/)); CALL var%setData(this%olatlon(2, :) / degree2radian)
    var = nc%setVariable("time", "f32", (/dim1/)); CALL var%setData(this%tempor)
    var = nc%setVariable("alt", "f32", (/dim1/)); CALL var%setData(this%height)

    DO i = 1, varNum
      var = nc%setVariable(TRIM(this%varNames(i)), "f32", (/dim1/)); CALL var%setData(this%obsData(:, i))
    END DO
    PRINT *, 'Output Obs files.', filename//TRIM(this%obsType)//'.nc'

    ! close the dataset
    CALL nc%CLOSE()
  END SUBROUTINE

END MODULE
