!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (), 2022/8/14, @GBA-MWF, Shenzhen
!     enhanced obsField_t to have get_valType and set_valType routines for checking as valType is private
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE ObsField_m
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE AuxTypeObs_m, ONLY: AuxTypeObs_t
  USE MPObs_m, ONLY: MPObs_t
  USE ObsAttr_m, ONLY: ObsAttr_t, ObsAttrSAT_t

  TYPE, EXTENDS(AuxTypeObs_t) :: ObsField_t
    PRIVATE

    TYPE(GridIdx_t), ALLOCATABLE, PUBLIC :: idx(:)
    ! CHARACTER(len=10) :: type
    REAL(r_kind), ALLOCATABLE, PUBLIC :: values(:)
    REAL(r_kind), ALLOCATABLE, PUBLIC :: valueArray(:, :)
    REAL(r_kind), PUBLIC :: locObs(3)  ! lat, lon, alt

    CHARACTER(len=20), PUBLIC :: name
    CHARACTER(len=20), PUBLIC :: valType       ! Type can be SP(Single point)/MC(Multi Channel, satellite)
    CHARACTER(len=20), PUBLIC :: obsType

    CLASS(ObsAttr_t), ALLOCATABLE, PUBLIC :: attr
    TYPE(ObsAttrSAT_t), PUBLIC :: ObsAttrSat

    REAL(r_kind), ALLOCATABLE, PUBLIC :: errors(:)
  CONTAINS
    PROCEDURE, PUBLIC :: dot_multiply, divideByValue, divideByVec
    PROCEDURE, PUBLIC :: subtraction
    PROCEDURE, PUBLIC :: addition
    PROCEDURE, PUBLIC :: multiply

    FINAL :: destructor
    PROCEDURE, PUBLIC :: Get_Name
    PROCEDURE, PUBLIC :: Set_Name
    PROCEDURE, PUBLIC :: Get_ObsType
    PROCEDURE, PUBLIC :: Set_ObsType
    PROCEDURE, PUBLIC :: Get_ValType
    PROCEDURE, PUBLIC :: Set_ValType
    PROCEDURE, PUBLIC :: Get_Id_Name

    ! PROCEDURE, PUBLIC :: dot_multiply
    GENERIC :: OPERATOR(.DOT.) => dot_multiply
    GENERIC :: OPERATOR(-) => subtraction
    GENERIC :: OPERATOR(+) => addition
    GENERIC :: OPERATOR(/) => divideByValue, divideByVec
    GENERIC :: OPERATOR(*) => multiply

  END TYPE ObsField_t

  INTERFACE ObsField_t
    PROCEDURE :: constructor
  END INTERFACE ObsField_t

CONTAINS

  FUNCTION addition(Y1, Y2) RESULT(X3)
    CLASS(ObsField_t), INTENT(IN) :: Y1
    TYPE(ObsField_t), INTENT(IN) :: Y2
    TYPE(ObsField_t) ::X3

    X3 = Y2
    IF (TRIM(X3%valType) .EQ. 'SP') THEN
      X3%values = Y2%values + Y1%values
    ELSE
      X3%valueArray = Y2%valueArray + Y1%valueArray
    END IF
  END FUNCTION addition

  FUNCTION subtraction(Y1, Y2) RESULT(X3)
    CLASS(ObsField_t), INTENT(IN) :: Y1
    TYPE(ObsField_t), INTENT(IN) :: Y2
    TYPE(ObsField_t) ::X3

    X3 = Y2
    IF (TRIM(X3%valType) .EQ. 'SP') THEN
      X3%values = Y2%values - Y1%values
    ELSE
      IF (.NOT. ALLOCATED(X3%valueArray)) THEN
        PRINT *, '+========================================================================================='
        PRINT *, 'ObsField subtraction - This obsField valueArray is not allocated for this obsField type!'
        PRINT *, 'Check your obsField initialization and rerun!!!'
        PRINT *, '+========================================================================================='
        STOP
      END IF
      X3%valueArray = Y2%valueArray - Y1%valueArray
    END IF
  END FUNCTION subtraction

  FUNCTION constructor(configFile, mpObs, TYPE) RESULT(this)
    IMPLICIT NONE
    TYPE(ObsField_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: TYPE

    this%valType = 'SP'
    IF (PRESENT(TYPE)) this%valType = TRIM(TYPE)
    this%AuxTypeObs_t = AuxTypeObs_t(mpObs)
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(ObsField_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%idx)) DEALLOCATE (this%idx)
    IF (ALLOCATED(this%values)) DEALLOCATE (this%values)
    IF (ALLOCATED(this%errors)) DEALLOCATE (this%errors)
    IF (ALLOCATED(this%attr)) DEALLOCATE (this%attr)
    IF (ALLOCATED(this%valueArray)) DEALLOCATE (this%valueArray)
  END SUBROUTINE destructor

  !> @brief
!! Return the name of this field.
  FUNCTION Get_Name(this) RESULT(name)
    IMPLICIT NONE
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=20) :: name

    name = this%name
  END FUNCTION

  SUBROUTINE Set_Name(this, name)
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=*) :: name

    this%name = name
  END SUBROUTINE

  FUNCTION Get_ObsType(this) RESULT(obsType)
    IMPLICIT NONE
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=20) :: obsType

    obsType = TRIM(this%obsType)
  END FUNCTION

  SUBROUTINE Set_ObsType(this, obsType)
    IMPLICIT NONE
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=*) :: obsType

    this%obsType = TRIM(obsType)
  END SUBROUTINE

  FUNCTION Get_ValType(this) RESULT(ValType)
    IMPLICIT NONE
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=20) :: ValType

    ValType = TRIM(this%ValType)
  END FUNCTION Get_ValType

  SUBROUTINE Set_ValType(this, ValType)
    IMPLICIT NONE
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=*) :: ValType

    this%ValType = TRIM(ValType)
  END SUBROUTINE Set_ValType

  FUNCTION Get_Id_Name(this) RESULT(id_name)
    CLASS(ObsField_t) :: this
    CHARACTER(LEN=20) :: id_name

    id_name = TRIM(this%Get_ObsType())//"_"//TRIM(this%Get_Name())
  END FUNCTION

  FUNCTION dot_multiply(Y1, Y2) RESULT(d0)
    CLASS(ObsField_t), INTENT(IN) :: Y1
    TYPE(ObsField_t), INTENT(IN) :: Y2
    INTEGER(i_kind) :: i
    REAL(r_kind)    :: rs
    REAL(r_kind) :: d0

    rs = 0.0D0
    d0 = 0.0D0

    IF (TRIM(Y1%valType) .EQ. 'SP') THEN
      rs = rs + SUM(Y1%values * Y2%values)
    ELSE
      rs = rs + SUM(Y1%valueArray * Y2%valueArray)
    END IF

    CALL Y1%mpObs%AllReduceSumReal(rs, d0)

  END FUNCTION

  FUNCTION divideByVec(Y1, Y0) RESULT(Y2)
    CLASS(ObsField_t), INTENT(IN) :: Y1
    TYPE(ObsField_t), INTENT(IN) :: Y0
    TYPE(ObsField_t) :: Y2
    Y2 = Y1
    Y2%values = Y1%values / Y0%values
  END FUNCTION

  FUNCTION divideByValue(Y1, divider) RESULT(Y2)
    CLASS(ObsField_t), INTENT(IN) :: Y1
    REAL(r_kind), INTENT(IN) :: divider
    TYPE(ObsField_t) :: Y2
    Y2 = Y1
    Y2%values = Y1%values / divider
  END FUNCTION

  FUNCTION multiply(X1, mul) RESULT(X2)
    CLASS(ObsField_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: mul
    TYPE(ObsField_t) :: X2

    X2 = X1
    X2%values = X1%values * mul

  END FUNCTION
END MODULE ObsField_m
