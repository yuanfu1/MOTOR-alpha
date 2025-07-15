!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/08/04, @GBA-MWF, Shenzhen
!     Added addition operator to ObsSet_t
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/13, @GBA-MWF, Shenzhen
!     Added obsMap and obsInv changes obs values and inverse it back.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/12/03, @GBA-MWF, Shenzhen
!     Added a new function of empty ObsSet.
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE ObsSet_m
  USE ObsField_m, ONLY: ObsField_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE MPObs_m, ONLY: MPObs_t
  USE AuxTypeObs_m, ONLY: AuxTypeObs_t
! #define TRACK_DEBUG_INFO

  TYPE, EXTENDS(AuxTypeObs_t) :: ObsSet_t
    TYPE(ObsField_t), ALLOCATABLE :: ObsFields(:)
  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: dot_multiply, divideByValue, divideByVec
    PROCEDURE, PUBLIC :: subtraction
    PROCEDURE, PUBLIC :: addition
    PROCEDURE, PUBLIC :: getObsIdx
    PROCEDURE, PUBLIC :: getObsTypeIdx
    PROCEDURE, PUBLIC :: zeroCopy
    PROCEDURE, PUBLIC :: multiply
    PROCEDURE, PUBLIC :: getObs
    PROCEDURE, PUBLIC :: rmVar
    PROCEDURE, PUBLIC :: getObsIdxWithIDName

    ! Yuanfu Xie added an empty obsSet function:
    PROCEDURE, PUBLIC :: empty

    PROCEDURE, PUBLIC :: obsMap
    PROCEDURE, PUBLIC :: obsFun

    GENERIC :: OPERATOR(.DOT.) => dot_multiply
    GENERIC :: OPERATOR(+) => addition
    GENERIC :: OPERATOR(-) => subtraction
    GENERIC :: OPERATOR(/) => divideByValue, divideByVec
    GENERIC :: OPERATOR(*) => multiply

  END TYPE ObsSet_t

  INTERFACE ObsSet_t
    PROCEDURE :: constructor
  END INTERFACE ObsSet_t

CONTAINS

  SUBROUTINE rmVar(this, varName)
    CLASS(ObsSet_t) :: this
    CHARACTER(*) :: varName
    INTEGER(i_kind) :: i, j, varIdx
    TYPE(ObsField_t), ALLOCATABLE :: temp(:)

    varIdx = this%getObsIdxWithIDName(varName)
    DO WHILE (varIdx /= 0)
      temp = this%ObsFields
      DEALLOCATE (this%ObsFields)
      ALLOCATE (this%ObsFields(SIZE(temp) - 1))
      j = 1
      DO i = 1, SIZE(temp)
        IF (i .NE. varIdx) THEN
          this%ObsFields(j) = temp(i)
          j = j + 1
        END IF
      END DO
      DEALLOCATE (temp)
      varIdx = this%getObsIdxWithIDName(varName)
    END DO

    ! &
    ! this%fields = pack(this%fields, (/(i, i=1, SIZE(this%fields))/) /= this%getVarIdx(varName))
  END SUBROUTINE

  FUNCTION getObs(this, obsName) RESULT(obs)
    CLASS(ObsSet_t), TARGET :: this
    CHARACTER(*) :: obsName
    TYPE(ObsField_t), POINTER :: obs

    DO i = LBOUND(this%ObsFields, 1), UBOUND(this%ObsFields, 1)
      IF (TRIM(this%ObsFields(i)%Get_Id_Name()) .EQ. TRIM(obsName)) THEN
        obs => this%ObsFields(i)
        RETURN
      END IF
    END DO
  END FUNCTION getObs

!> @brief
! @see
! @note
! @warning There can be multiple obsFields has the some obsName, but this function
!! returns the idx of the first only.
  FUNCTION getObsIdx(this, obsName) RESULT(idx)
    CLASS(ObsSet_t) :: this
    CHARACTER(*) :: obsName
    INTEGER(i_kind) :: idx
    INTEGER(i_kind) :: i

    idx = 0
    IF (ALLOCATED(this%ObsFields)) THEN
      DO i = LBOUND(this%ObsFields, 1), UBOUND(this%ObsFields, 1)
        IF (TRIM(this%ObsFields(i)%Get_Name()) .EQ. TRIM(obsName)) THEN
          idx = i
          RETURN
        END IF
      END DO
    END IF
  END FUNCTION

  FUNCTION getObsIdxWithIDName(this, obsName) RESULT(idx)
    CLASS(ObsSet_t) :: this
    CHARACTER(*) :: obsName
    INTEGER(i_kind) :: idx
    INTEGER(i_kind) :: i

    idx = 0
    IF (ALLOCATED(this%ObsFields)) THEN
      DO i = LBOUND(this%ObsFields, 1), UBOUND(this%ObsFields, 1)
        IF (INDEX(TRIM(this%ObsFields(i)%Get_Id_Name()), TRIM(obsName)) .NE. 0) THEN
          idx = i
          RETURN
        END IF
      END DO
    END IF
  END FUNCTION

  !@brief
  ! Note that this type addition does not meet the symmetry requirement, i.e. X+Y /= Y+X
  FUNCTION addition(X1, X2) RESULT(X3)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    TYPE(ObsSet_t), INTENT(IN) :: X2
    TYPE(ObsSet_t) ::X3

    X3 = X1
    DO i = LBOUND(X3%ObsFields, 1), UBOUND(X3%ObsFields, 1)
      X3%ObsFields(i) = X2%ObsFields(i) + X1%ObsFields(i)
    END DO

  END FUNCTION addition

  FUNCTION getObsTypeIdx(this, TypeName) RESULT(idx)
    CLASS(ObsSet_t) :: this
    CHARACTER(*) :: TypeName
    INTEGER(i_kind) :: idx
    INTEGER(i_kind) :: i

    idx = 0
    IF (ALLOCATED(this%ObsFields)) THEN
    DO i = LBOUND(this%ObsFields, 1), UBOUND(this%ObsFields, 1)
      IF (INDEX(TRIM(this%ObsFields(i)%Get_ObsType()), TRIM(TypeName)) .NE. 0) THEN
        idx = i
        RETURN
      END IF
    END DO
    END IF
  END FUNCTION

  SUBROUTINE showInfo(this)
    CLASS(ObsSet_t) :: this
    INTEGER(i_kind) :: i

    PRINT *, 'ObsSet var List:'
    DO i = LBOUND(this%ObsFields, 1), UBOUND(this%ObsFields, 1)
      PRINT *, this%ObsFields(i)%Get_Id_Name(), ' in OBS: ', this%ObsFields(i)%Get_ObsType()
    END DO
  END SUBROUTINE

  FUNCTION subtraction(X1, X2) RESULT(X3)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    TYPE(ObsSet_t), INTENT(IN) :: X2
    TYPE(ObsSet_t) ::X3

    X3 = X1
    DO i = 1, UBOUND(X3%ObsFields, 1)
      X3%ObsFields(i) = X2%ObsFields(i) - X1%ObsFields(i)
    END DO

  END FUNCTION subtraction

  REAL(r_kind) FUNCTION dot_multiply(X1, X2)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    TYPE(ObsSet_t), INTENT(IN) :: X2
    INTEGER(i_kind) :: i, j

    ! Local variables:
    CHARACTER(LEN=6000) :: names
    REAL(r_kind) :: prod(300)

    dot_multiply = 0.0D0
    prod = 0.0D0
    ! j = 0
    DO i = LBOUND(X1%ObsFields, 1), UBOUND(X1%ObsFields, 1)
      prod(i) = (X1%ObsFields(i) .DOT.X2%ObsFields(i))
!       WRITE(names(j+1:j+16),2) TRIM(X1%ObsFields(i)%Get_Name()),prod(i)
! 2     FORMAT(A4,1X,D11.4,1X)
!       j=j+16
      dot_multiply = prod(i) + dot_multiply
    END DO

#ifdef TRACK_DEBUG_INFO
    WRITE (*, 1) names(1:16 * UBOUND(X1%ObsFields, 1)), X1%mpObs%sg%gLevel, X1%mpObs%myrank
1   FORMAT('ObsSet-norms: ', A, ' at G:', I2, ' proc: ', I2)
#endif

  END FUNCTION

  FUNCTION multiply(X1, mul) RESULT(X2)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: mul
    TYPE(ObsSet_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%ObsFields, 1), UBOUND(X1%ObsFields, 1)
      X2%ObsFields(i) = X1%ObsFields(i) * mul
    END DO
  END FUNCTION

  FUNCTION divideByVec(X1, X0) RESULT(X2)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    TYPE(ObsSet_t), INTENT(IN) :: X0
    TYPE(ObsSet_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%ObsFields, 1), UBOUND(X1%ObsFields, 1)
      X2%ObsFields(i) = X1%ObsFields(i) / X0%ObsFields(i)
    END DO
  END FUNCTION

  FUNCTION divideByValue(X1, divider) RESULT(X2)
    CLASS(ObsSet_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: divider
    TYPE(ObsSet_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%ObsFields, 1), UBOUND(X1%ObsFields, 1)
      X2%ObsFields(i) = X1%ObsFields(i) / divider
    END DO
  END FUNCTION

  SUBROUTINE obsMap(this, nameFrom, nameTo, fwdFunc)
    IMPLICIT NONE
    CLASS(ObsSet_t), INTENT(INOUT) :: this
    CHARACTER(*) nameFrom, nameTo
    REAL(r_kind) :: fwdFunc

    ! Local variables:
    INTEGER(i_kind) :: iv, i

    ! Check if the name exists:
    iv = this%getObsIdx(TRIM(nameFrom))
    IF (iv .GT. 0) THEN
      DO i = LBOUND(this%ObsFields(iv)%values, 1), UBOUND(this%ObsFields(iv)%values, 1)
        this%ObsFields(iv)%values(i) = &
          fwdFunc(this%ObsFields(iv)%values(i), TRIM(nameFrom))
      END DO
    END IF
    CALL this%ObsFields(iv)%set_name(TRIM(nameTo))
  END SUBROUTINE obsMap

  REAL FUNCTION obsFun(this, VALUE, name) RESULT(toVal)
    IMPLICIT NONE
    CLASS(ObsSet_t), INTENT(IN) :: this
    CHARACTER(*), INTENT(IN) :: name
    REAL(r_kind), INTENT(IN) :: VALUE

    toVal = 0.0D0
  END FUNCTION obsFun

  FUNCTION constructor(configFile, mpObs) RESULT(this)
    IMPLICIT NONE
    TYPE(ObsSet_t) :: this
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    this%AuxTypeObs_t = AuxTypeObs_t(mpObs)
  END FUNCTION constructor

  FUNCTION zeroCopy(this) RESULT(Yout)
    CLASS(ObsSet_t) :: this
    TYPE(ObsSet_t) :: Yout
    INTEGER :: i

    Yout = this
    Yout%AuxTypeObs_t = AuxTypeObs_t(this%mpObs)

    ! DO i = LBOUND(Yout%ObsFields, 1), UBOUND(Yout%ObsFields, 1)
    DO i = 1, UBOUND(Yout%ObsFields, 1) ! exclude zero
      IF (SIZE(Yout%ObsFields(i)%values, 1) .GT. 0) &
        Yout%ObsFields(i)%values = 0.0D0
    END DO
  END FUNCTION

  ! Yuanfu Xie added this on 2022-12-03: it is an empty obsSet
  FUNCTION empty(this) RESULT(Yout)
    CLASS(ObsSet_t) :: this
    TYPE(ObsSet_t) :: Yout
    INTEGER :: i

    Yout = this
    Yout%AuxTypeObs_t = AuxTypeObs_t(this%mpObs)

    ALLOCATE (Yout%ObsFields(0))
  END FUNCTION empty

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%ObsFields)) DEALLOCATE (this%ObsFields)
  END SUBROUTINE destructor

END MODULE ObsSet_m
