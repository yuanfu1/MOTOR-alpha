!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/20, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the HydroStatic Object which is tranform the variables from control variable space to
!! model space.
MODULE HydroStatic_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: HydroStatic_t
    TYPE(State_t), POINTER :: X

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    PROCEDURE, PUBLIC, NOPASS :: transElemForward
    PROCEDURE, PUBLIC, NOPASS :: transElemAdjointMultiply

    FINAL :: destructor
  END TYPE HydroStatic_t

  INTERFACE HydroStatic_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(HydroStatic_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(HydroStatic_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k, hStcl, vStcl

    IF ((X%getVarIdx(TRIM('temp')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('qvapor')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('pres')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in HydroStatic_t, STOP! '
      STOP
    END IF

    ! ASSOCIATE (temp => X%Fields(X%getVarIdx(TRIM('temp')))%data, &
    !            qvapor => X%Fields(X%getVarIdx(TRIM('qvapor')))%data, &
    !            pres => X%Fields(X%getVarIdx(TRIM('pres')))%data)
    !   DO j = 2, X%sg%vLevel
    !     swapDA(j, :, sg_mp%tSlots) = swapDA(j - 1, :, sg_mp%tSlots)*exp(-(sg%zHght(j, :) - sg%zHght(j - 1, :))*g &
    !                                                                     /(dry_air_gas_const* &
    !                                                                     (Xm%Fields(Xm%getVarIdx('temp'))%data(j, :, sg_mp%tSlots)* &
    !                                                       (Xm%Fields(Xm%getVarIdx('qvapor'))%data(j, :, sg_mp%tSlots)*0.608 + 1) + &
    !                                                                  Xm%Fields(Xm%getVarIdx('temp'))%data(j - 1, :, sg_mp%tSlots)* &
    !                                              (Xm%Fields(Xm%getVarIdx('qvapor'))%data(j - 1, :, sg_mp%tSlots)*0.608 + 1))/2.0D0))
    !   END DO

    ! END ASSOCIATE

    ! CALL X%exHalo()
    CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, X%Fields(X%getVarIdx(TRIM('pres')))%DATA)

  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind)  :: i, j, k, hStcl, vStcl
    TYPE(State_t) :: dXi

    IF ((dX%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('wwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('vwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transAdjMultiply in HydroStatic_t, STOP! '
      STOP
    END IF

    ! Construct an empty state to contain increment
    dXi = dX%zeroCopy()
    dXi%fields(dXi%getVarIdx('wwnd')) = dX%fields(dX%getVarIdx('wwnd'))

    ASSOCIATE (uwnd => dXi%Fields(dXi%getVarIdx(TRIM('uwnd')))%DATA, &
               vwnd => dXi%Fields(dXi%getVarIdx(TRIM('vwnd')))%DATA, &
               wwnd => dXi%Fields(dXi%getVarIdx(TRIM('wwnd')))%DATA)

    END ASSOCIATE

    CALL dXi%rmVar('wwnd')
    CALL dXi%exHaloRevSum()
    CALL dXi%exHalo()

    ! Add the increment with dX
    CALL dx%rmVar('wwnd')
    dX = dX + dXi
  END SUBROUTINE transAdjMultiply

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX
    CALL this%transAdjMultiply(dX1)
  END FUNCTION

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X
    CALL this%transFwdNonLinear(X1)
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(HydroStatic_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  FUNCTION transElemForward(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION

  FUNCTION transElemAdjointMultiply(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION
END MODULE HydroStatic_m
