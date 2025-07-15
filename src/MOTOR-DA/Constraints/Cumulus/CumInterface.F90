!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

MODULE CumInterface_m
  USE PhyCon_m, ONLY: PhyCon_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t

  TYPE, EXTENDS(PhyCon_t) :: CumInterface_t

    TYPE(State_t), POINTER :: X
    TYPE(State_t), POINTER :: dX

  CONTAINS
    ! PROCEDURE, PUBLIC, PASS(this) :: Initialize
    PROCEDURE, PUBLIC, PASS(this) :: CumForward
    PROCEDURE, PUBLIC, PASS(this) :: CumAdjoint

    FINAL :: destructor
  END TYPE CumInterface_t

  INTERFACE CumInterface_t
    PROCEDURE :: constructor
  END INTERFACE CumInterface_t

CONTAINS

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(CumInterface_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(CumInterface_t), INTENT(INOUT) :: this
  END SUBROUTINE destructor

  ! SUBROUTINE Initialize(this, X, dX)
  !   IMPLICIT NONE
  !   CLASS(CumInterface_t) :: this
  !   TYPE(State_t), INTENT(IN) :: X
  !   TYPE(State_t), INTENT(IN) :: dX

  !   CALL this%DataPre(X, dX)

  ! END SUBROUTINE

  SUBROUTINE CumForward(this, X)

    IMPLICIT NONE
    CLASS(CumInterface_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i_this, i_x, size_dim_time
    CALL this%FwdDataPre(X)
    CALL this%RK4(X)

    size_dim_time = SIZE(this%qvapor_out, dim=3)

    i_x = 1
    DO i_this = 1, size_dim_time
      IF (MOD(i_this, this%kt) .EQ. 1) THEN
        X%Fields(X%getVarIdx('uwnd'))%DATA(:, :, i_x) = this%uwnd_out(:, :, i_this)
        X%Fields(X%getVarIdx('vwnd'))%DATA(:, :, i_x) = this%vwnd_out(:, :, i_this)
        X%Fields(X%getVarIdx('wwnd'))%DATA(:, :, i_x) = this%wwnd_out(:, :, i_this)
        X%Fields(X%getVarIdx('theta'))%DATA(:, :, i_x) = this%theta_out(:, :, i_this)
        X%Fields(X%getVarIdx('qvapor'))%DATA(:, :, i_x) = this%qvapor_out(:, :, i_this)
        X%Fields(X%getVarIdx('pres'))%DATA(:, :, i_x) = this%pres_out(:, :, i_this)
        i_x = i_x + 1
      END IF
    END DO

  END SUBROUTINE

  SUBROUTINE CumAdjoint(this, X, dX)

    IMPLICIT NONE
    CLASS(CumInterface_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: dX
    INTEGER(i_kind) :: i_this, i_x, size_dim_time

    CALL this%AdjDataPre(X, dX)

    CALL this%RK4_AD(X)

    size_dim_time = SIZE(this%qvapor_out, dim=3)

    i_x = 1
    DO i_this = 1, size_dim_time
      IF (MOD(i_this, this%kt) .EQ. 1) THEN
        dX%Fields(dX%getVarIdx('uwnd'))%DATA(:, :, i_x) = this%uwnd_outb(:, :, i_this)
        dX%Fields(dX%getVarIdx('vwnd'))%DATA(:, :, i_x) = this%vwnd_outb(:, :, i_this)
        dX%Fields(dX%getVarIdx('wwnd'))%DATA(:, :, i_x) = this%wwnd_outb(:, :, i_this)
        dX%Fields(dX%getVarIdx('theta'))%DATA(:, :, i_x) = this%theta_outb(:, :, i_this)
        dX%Fields(dX%getVarIdx('qvapor'))%DATA(:, :, i_x) = this%qvapor_outb(:, :, i_this)
        dX%Fields(dX%getVarIdx('pres'))%DATA(:, :, i_x) = this%pres_outb(:, :, i_this)
        i_x = i_x + 1
      END IF
    END DO

  END SUBROUTINE

END MODULE CumInterface_m
