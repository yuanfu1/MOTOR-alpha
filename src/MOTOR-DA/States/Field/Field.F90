!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.Field
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/3/29, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/09/03, @GBA-MWF, Shenzhen
!     Added a new function of sum up of grid function values of the neighbors of a given
!     grid cell in nr nested levels, e.g. for a grid cell, j, it sums up the grid values
!     of recursively nr depth of neighbors. For nr = 2, it sums up j's neighbors and the
!     neighbor's neighbors. Note this involves grid cell cross multiple processors.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/04/18, @GBA-MWF, Shenzhen
!     Added integer variables for determining the control variables in space, and time,
!     replacing the large array of mask.
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE Field_m
  USE kinds_m, ONLY: i_kind, r_kind, i_byte
  USE singleGrid_m, ONLY: singleGrid_t, GridIdx_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t

  TYPE, EXTENDS(AuxTypeSG_t) :: Field_t
    PRIVATE

    ! CHARACTER(len=10) :: type
    CHARACTER(len=20) :: name                     !< Var names
    CHARACTER(len=20) :: TYPE                     !< The type of this var field, 3D surface(3DS)/4D stagggered(4DS)/4D unstagggered(4DU)
    REAL(r_kind), ALLOCATABLE, PUBLIC :: DATA(:, :, :)    !< In dimension of (vLevel, num_cell, tSlots)
    ! INTEGER(i_byte), ALLOCATABLE, PUBLIC :: mask(:, :, :)
    ! Use small arrays to replace mask: Yuanfu Xie 2024-04-18
    INTEGER(i_byte), PUBLIC :: maskHorizontal
    INTEGER(i_byte), ALLOCATABLE, PUBLIC :: maskVertical(:), maskTemporal(:)
    INTEGER(i_kind), PUBLIC :: numTotalMask = 0
    INTEGER(i_kind), PUBLIC :: lowBound = 0       !< 0 default no lower bound constraint; 1 with a constraint

  CONTAINS
    PROCEDURE, PUBLIC :: initialize

    PROCEDURE, PUBLIC :: add
    PROCEDURE, PUBLIC :: multiply
    PROCEDURE, PUBLIC :: subtract

    PROCEDURE, PUBLIC :: dot_multiply, divideByValue, divideByVec
    PROCEDURE, PUBLIC :: Get_Name
    PROCEDURE, PUBLIC :: Set_Name
    PROCEDURE, PUBLIC :: Get_Type
    PROCEDURE, PUBLIC :: Set_Type

    PROCEDURE, PUBLIC :: Get_Value
    PROCEDURE, PUBLIC :: Set_Value
    PROCEDURE, PUBLIC :: Add_Value

    ! Yuanfu Xie Added a new function of sum up of grid function values of the neighbors of a given
    ! grid cell in nr nested levels, e.g. for a grid cell, j, it sums up the grid values
    ! of recursively nr depth of neighbors. For nr = 2, it sums up j's neighbors and the
    ! neighbor's neighbors. Note this involves grid cell cross multiple processors.
    PROCEDURE, PUBLIC :: Sum_Neighbors

    FINAL :: destructor
    GENERIC :: OPERATOR(.DOT.) => dot_multiply
    GENERIC :: OPERATOR(*) => multiply
    GENERIC :: OPERATOR(+) => add
    GENERIC :: OPERATOR(-) => subtract
    GENERIC :: OPERATOR(/) => divideByValue, divideByVec

  END TYPE Field_t

  ! ! Yuanfu Xie adds this 2-D field for holding topography, landFactor etc.
  ! TYPE, EXTENDS(AuxTypeSG_t) :: Field2D_t
  !   CHARACTER(len=10) :: type
  !   REAL(r_kind), ALLOCATABLE :: data(:, :)
  ! END TYPE Field2D_t

  INTERFACE Field_t
    PROCEDURE :: constructor
  END INTERFACE Field_t

CONTAINS

  FUNCTION constructor(sg, varName, TYPE) RESULT(this)
    TYPE(Field_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(*), INTENT(IN) :: varName
    CHARACTER(*), INTENT(IN) :: TYPE

    CALL this%initialize(sg, varName, TYPE)
  END FUNCTION

  SUBROUTINE initialize(this, sg, varName, TYPE, numTimes)
    CLASS(Field_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(*), INTENT(IN) :: varName
    CHARACTER(*), INTENT(IN) :: TYPE
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: numTimes !< Yuanfu Xie added this for model IC

    !this%AuxTypeSG_t = AuxTypeSG_t(sg)
    CALL this%AuxTypeSG_t%aux_initialize(sg)

    ! Yuanfu Xie added this option for saving space for model IC
    IF (PRESENT(numTimes)) THEN
      ALLOCATE (this%DATA(this%sg%vlevel, this%sg%num_cell, numTimes))
      ALLOCATE (this%maskVertical(this%sg%vlevel), this%maskTemporal(numTimes))
    ELSE
      IF (TRIM(TYPE) .EQ. '4DU') THEN
        ALLOCATE (this%DATA(this%sg%vlevel, this%sg%num_cell, this%sg%tSlots))
      ELSE IF (TRIM(TYPE) .EQ. '3DS') THEN
        ALLOCATE (this%DATA(1, this%sg%num_cell, this%sg%tSlots))
      END IF
      ALLOCATE (this%maskVertical(this%sg%vlevel), this%maskTemporal(this%sg%tSlots))
    END IF

    this%DATA = 0.0D0
    this%maskVertical = 1   ! Default value, all dimensions are controls
    this%maskHorizontal = 1
    this%maskTemporal = 1
    this%name = TRIM(varName)
    this%TYPE = TRIM(TYPE)
  END SUBROUTINE

  FUNCTION add(X1, X2) RESULT(X3)
    CLASS(Field_t), INTENT(IN) :: X1
    TYPE(Field_t), INTENT(IN) :: X2
    TYPE(Field_t) :: X3

    X3 = X1
    X3%DATA = X1%DATA + X2%DATA

  END FUNCTION add

  FUNCTION multiply(X1, mul) RESULT(X2)
    CLASS(Field_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: mul
    TYPE(Field_t) :: X2

    X2 = X1
    X2%DATA = X1%DATA * mul
  END FUNCTION

  FUNCTION divideByVec(X1, X0) RESULT(X2)
    CLASS(Field_t), INTENT(IN) :: X1
    TYPE(Field_t), INTENT(IN) :: X0
    TYPE(Field_t) :: X2
    X2 = X1
    X2%DATA = X1%DATA / X0%DATA
  END FUNCTION

  FUNCTION divideByValue(X1, divider) RESULT(X2)
    CLASS(Field_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: divider
    TYPE(Field_t) :: X2
    X2 = X1
    X2%DATA = X1%DATA / divider
  END FUNCTION

  FUNCTION subtract(X1, X2) RESULT(X3)
    CLASS(Field_t), INTENT(IN) :: X1
    TYPE(Field_t), INTENT(IN) :: X2
    TYPE(Field_t) ::X3

    X3 = X2
    X3%DATA = X2%DATA - X1%DATA
  END FUNCTION subtract

!> @brief
!! Get the specific value on idx of this field.
  FUNCTION Get_Value(this, idx) RESULT(rsl)
    CLASS(Field_t) :: this
    TYPE(GridIdx_t) :: idx
    REAL(r_kind) :: rsl

    rsl = this%DATA(idx%vIdx, idx%hIdx, idx%tIdx)
  END FUNCTION

!> @brief
!! Get the specific value on idx of this field.
  SUBROUTINE Set_Value(this, idx, rsl)
    CLASS(Field_t) :: this
    TYPE(GridIdx_t) :: idx
    REAL(r_kind) :: rsl

    this%DATA(idx%vIdx, idx%hIdx, idx%tIdx) = rsl
  END SUBROUTINE

!> @brief
!! Get the specific value on idx of this field.
  SUBROUTINE Add_Value(this, idx, rsl)
    CLASS(Field_t) :: this
    TYPE(GridIdx_t) :: idx
    REAL(r_kind) :: rsl

    this%DATA(idx%vIdx, idx%hIdx, idx%tIdx) = this%DATA(idx%vIdx, idx%hIdx, idx%tIdx) + rsl
  END SUBROUTINE

  FUNCTION dot_multiply(X1, X2)
    CLASS(Field_t), INTENT(IN) :: X1
    TYPE(Field_t), INTENT(IN) :: X2
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind)    :: rs
    REAL(r_kind) :: dot_multiply

    ! rs = 0.0D0
    dot_multiply = 0.0D0

    ! DO j = 1, X1%sg%tSlots
    !   DO i = 1, X1%sg%num_icell
    !     rs = rs + SUM(X1%data(:, i, j)*X2%data(:, i, j))
    !   END DO
    ! END DO

    rs = SUM(X1%DATA(:, 1:X1%sg%num_icell, :) * X2%DATA(:, 1:X1%sg%num_icell, :))
    !print*,'Field DOT: rs',rs
    CALL X1%mpddSub%AllReduceSumReal(rs, dot_multiply)
  END FUNCTION

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(Field_t), INTENT(inout) :: this

    IF (ALLOCATED(this%DATA)) DEALLOCATE (this%DATA)
    IF (ALLOCATED(this%maskVertical)) DEALLOCATE (this%maskVertical)
    IF (ALLOCATED(this%maskTemporal)) DEALLOCATE (this%maskTemporal)
  END SUBROUTINE destructor

  !> @brief
!! Return the name of this field.
  FUNCTION Get_Name(this) RESULT(name)
    IMPLICIT NONE
    CLASS(Field_t) :: this
    CHARACTER(LEN=20) :: name

    name = this%name
  END FUNCTION

  SUBROUTINE Set_Name(this, name)
    CLASS(Field_t) :: this
    CHARACTER(LEN=*) :: name

    this%name = name
  END SUBROUTINE

  FUNCTION Get_Type(this) RESULT(TYPE)
    IMPLICIT NONE
    CLASS(Field_t) :: this
    CHARACTER(LEN=20) :: TYPE

    TYPE = TRIM(this%TYPE)
  END FUNCTION

  SUBROUTINE Set_Type(this, TYPE)
    IMPLICIT NONE
    CLASS(Field_t) :: this
    CHARACTER(LEN=*) :: TYPE

    this%TYPE = TRIM(TYPE)
  END SUBROUTINE

  SUBROUTINE Sum_Neighbors(this)
    IMPLICIT NONE
    CLASS(Field_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: ic, in
    REAL(r_kind), ALLOCATABLE :: buf(:, :, :)

    ALLOCATE (buf(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    buf = 0.0D0
    buf(:, 1:this%sg%num_icell, :) = this%DATA(:, 1:this%sg%num_icell, :)

    ! Sum the neighbors over icells:
    DO ic = 1, this%sg%num_icell
      DO in = 1, this%sg%num_edge
        IF (this%sg%cell_adnb(in, ic) .GT. 0) &
          buf(:, ic, :) = buf(:, ic, :) + this%DATA(:, this%sg%cell_adnb(in, ic), :)
      END DO
    END DO

    ! Modify the data by the sum:
    this%DATA(:, :, :) = buf(:, :, :)

    DEALLOCATE (buf)
  END SUBROUTINE Sum_Neighbors

END MODULE Field_m
