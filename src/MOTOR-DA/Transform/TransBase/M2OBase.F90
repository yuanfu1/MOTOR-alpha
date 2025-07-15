!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE M2OBase_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t

  TYPE, ABSTRACT :: M2OBase_t

  CONTAINS
    PROCEDURE(transFwdNonLinear), DEFERRED :: fwdNL
    PROCEDURE(transFwdTanLinear), DEFERRED :: fwdTL
    PROCEDURE(transAdjMultiply), DEFERRED :: adjMul

    PROCEDURE(transFwdNonLinear_opr), DEFERRED :: fwdNL_opr
    PROCEDURE(transFwdTanLinear_opr), DEFERRED :: fwdTL_opr
    PROCEDURE(transAdjMultiply_opr), DEFERRED :: adjMul_opr
  END TYPE

  ABSTRACT INTERFACE
    SUBROUTINE transFwdNonLinear(this, X, Y)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(State_t), INTENT(in) :: X
      TYPE(ObsSet_t), INTENT(INOUT) :: Y
    END SUBROUTINE

    SUBROUTINE transFwdTanLinear(this, dX, dY, X)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(State_t), INTENT(in) :: dX
      TYPE(ObsSet_t), INTENT(INOUT) :: dY
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END SUBROUTINE

    SUBROUTINE transAdjMultiply(this, D, dX, X)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(ObsSet_t), INTENT(in) :: D
      TYPE(State_t), INTENT(INOUT) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END SUBROUTINE

    FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(State_t), INTENT(in) :: X
      TYPE(ObsSet_t) :: Y
    END FUNCTION

    FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(State_t), INTENT(in) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
      TYPE(ObsSet_t) :: dY
    END FUNCTION

    FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
      IMPORT :: State_t, ObsSet_t, M2OBase_t
      CLASS(M2OBase_t) :: this
      TYPE(ObsSet_t), INTENT(IN) :: D
      TYPE(State_t) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END FUNCTION

  END INTERFACE

CONTAINS

END MODULE M2OBase_m
