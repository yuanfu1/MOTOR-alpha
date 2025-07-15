!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/23, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/19, @GBA-MWF, Shenzhen
!     for adding an option currentField preparing to add a conversion of pressure obs to log(pressure).
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BFieldBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE ObsSet_m, ONLY: ObsSet_t

  TYPE, ABSTRACT, EXTENDS(AuxTypeSG_t):: BFieldBase_t
    REAL(r_kind), ALLOCATABLE :: sigma(:, :, :)
    CHARACTER(len=10) :: name
    REAL(r_kind) :: scaleParaX, scaleParaY, scaleParaZ, scaleParaT
    REAL(r_kind) :: ts_DataAggr
    REAL(r_kind) :: ts_DataDist
    REAL(r_kind) :: ts_DataCalc
    REAL(r_kind) :: ts_DataCalc_FWD
    REAL(r_kind) :: ts_DataCalc_ADJ
    REAL(r_kind) :: ts_DataLoad

  CONTAINS
    PROCEDURE, PRIVATE :: loadBMatFiles

    PROCEDURE, PUBLIC :: inverse_multiply
    PROCEDURE, PUBLIC :: initialize

    PROCEDURE(sqrt_inverse_multiply), DEFERRED :: sqrtInvMul
    PROCEDURE(sqrt_inverse_multiply_adjoint), DEFERRED :: sqrtInvMulAdj

    PROCEDURE, PUBLIC :: Get_Name
  END TYPE BFieldBase_t

  ABSTRACT INTERFACE
    ! Adding an option currentField preparing to add a conversion of pressure obs to log(pressure) obs
    SUBROUTINE sqrt_inverse_multiply_adjoint(this, field, currentField)
      IMPORT :: BFieldBase_t, Field_t
      CLASS(BFieldBase_t) :: this
      TYPE(Field_t), INTENT(INOUT) :: field
      TYPE(Field_t), OPTIONAL, INTENT(IN) :: currentField
    END SUBROUTINE

    SUBROUTINE sqrt_inverse_multiply(this, field)
      IMPORT :: BFieldBase_t, Field_t
      IMPLICIT NONE
      CLASS(BFieldBase_t) :: this
      TYPE(Field_t), INTENT(INOUT) :: field
    END SUBROUTINE sqrt_inverse_multiply
  END INTERFACE

CONTAINS

  SUBROUTINE initialize(this, configFile, sg, varName, Y)
    CLASS(BFieldBase_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN)  :: configFile
    CHARACTER(LEN=10), INTENT(IN)  :: varName
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    ! Yuanfu Xie added this comment for users to check the right routine on 2025-02-13:
    ! This is deferred routine, please to see its implementation !

  END SUBROUTINE initialize

  SUBROUTINE loadBMatFiles(this, varName)
    CLASS(BFieldBase_t) :: this
    CHARACTER(LEN=10), INTENT(IN) :: varName

    ! ! Here is the code of monk for simulating the fileinput of B mat files
    ! ALLOCATE (this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots));
    ! ! Set all the sigma to 1.0D0
    ! this%sigma = 1.0D0
  END SUBROUTINE loadBMatFiles

  SUBROUTINE inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldBase_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field

  END SUBROUTINE inverse_multiply

  !> @brief
!! Return the type of this field.
  FUNCTION Get_Name(this) RESULT(name)
    CLASS(BFieldBase_t) :: this
    CHARACTER(LEN=10) :: name

    name = this%name
  END FUNCTION
END MODULE BFieldBase_m
