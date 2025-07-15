!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief

MODULE M2MBase_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t

  TYPE, ABSTRACT :: M2MBase_t

  CONTAINS
    PROCEDURE(transFwdNonLinear), DEFERRED :: fwdNL
    PROCEDURE(transFwdTanLinear), DEFERRED :: fwdTL
    PROCEDURE(transAdjMultiply), DEFERRED :: adjMul

    PROCEDURE(transFwdNonLinear_opr), DEFERRED :: fwdNL_opr
    PROCEDURE(transFwdTanLinear_opr), DEFERRED :: fwdTL_opr
    PROCEDURE(transAdjMultiply_opr), DEFERRED :: adjMul_opr
  END TYPE

  ABSTRACT INTERFACE
    SUBROUTINE transFwdNonLinear(this, X)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(INOUT) :: X
    END SUBROUTINE

    SUBROUTINE transFwdTanLinear(this, dX, X)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(INOUT) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END SUBROUTINE

    SUBROUTINE transAdjMultiply(this, dX, X)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(INOUT) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END SUBROUTINE

    FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(in) :: X
      TYPE(State_t):: X1
    END FUNCTION

    FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(in) :: dX
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
      TYPE(State_t):: dX1
    END FUNCTION

    FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
      IMPORT :: State_t, ObsSet_t, M2MBase_t
      CLASS(M2MBase_t) :: this
      TYPE(State_t), INTENT(IN) :: dX
      TYPE(State_t) :: dX1
      TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    END FUNCTION

  END INTERFACE

CONTAINS

END MODULE M2MBase_m
