!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/7/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines a control variable case using vorticity and divergence at initial time as
!! control variables and Z-grid model to calculate model states at other time frames.
!! Control variables: $\zeta$, $\delta$, $\ln(P)$, $temp$ and $rh$, where rh is scaled mixing ratio
!! by an exponential function profile.
!!
MODULE c4DVar_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE C2MBase_m, ONLY: C2MBase_t
  USE State_m, ONLY: State_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t) :: c4DVar_t
    CLASS(dyCoresBase_t), POINTER :: dycore
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => c4DVarForward
    PROCEDURE, PUBLIC :: adjointOpr_s => c4DVarAdjoint
    PROCEDURE, PUBLIC :: tangentOpr_s => c4DVarTangent
    PROCEDURE, PUBLIC :: bckwardOpr_s => c4DVarBckward
  END TYPE c4DVar_t

CONTAINS
  SUBROUTINE initialize_s(this, sg, X, dyc)
    CLASS(c4DVar_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), INTENT(IN) :: X
    CLASS(dyCoresBase_t), TARGET, OPTIONAL :: dyc

    ! Local variables:
    INTEGER(i_kind) :: i

    ! Check the options:
    IF (.NOT. PRESENT(dyc)) THEN
      PRINT *, 'c4DVarAdjoint -- A 4DVar requires a dynamic core, please check ans rerun...'
      STOP
    END IF

    ! Pass the 4DVar elements:
    this%dycore => dyc
    this%analysis = X
    this%sg => sg

    ! Set control mask and lower bounds:
    DO i = LBOUND(this%analysis%fields, 1), UBOUND(this%analysis%fields, 1)

      ! The total controls:
      this%analysis%fields(i)%numTotalMask = &
        sg%num_icell_global * sg%vLevel   ! a 4DVar uses initial conditions as controls

      ! The mask of controls:
      this%analysis%fields(i)%maskHorizontal = 1
      this%analysis%fields(i)%maskVertical = 1
      this%analysis%fields(i)%maskTemporal = 0
      this%analysis%fields(i)%maskTemporal(1) = 1  ! Only the initial time as controls

      ! Lower bounds:
      IF (this%analysis%fields(i)%Get_Name() .EQ. 'temp') &
          this%analysis%fields(i)%lowBound = 1
      IF (this%analysis%fields(i)%Get_Name() .EQ. 'qvapor' .OR. &
          this%analysis%fields(i)%Get_Name() .EQ. 'qvapor_ctl') &
          this%analysis%fields(i)%lowBound = 1
    END DO
  END SUBROUTINE initialize_s

  ! Forward forecast:
  FUNCTION c4DVarForward(this, X) RESULT(Y)
    CLASS(c4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    REAL(r_kind) :: dt
    TYPE(State_t) :: iC

    dt = 900.0D0

    iC = X
    Y = X   ! Initialize Y by Yuanfu Xie 2025-01-11
    WRITE(*,1) UBOUND(iC%fields, 1), X%sg%gLevel, X%sg%tSlots
1   FORMAT('c4DVarForward - dycore fwd...numVar: ',I1, ' gLevel: ',I2,' tSlots: ',I2)
    CALL this%dycore%dycore_s(iC, Y)
    this%analysis = Y
    PRINT*,'c4DVarForward - done fwd...',(Y .DOT. Y)
  END FUNCTION c4DVarForward

  !< Adjoint:
  !< Note: Y is allocated according the time integration schemes, in the return
  !< it contains the gradient of initial conditions in terms of the cost func.
  !< For a 4DVar, X holds the H' R^(-1) (H(X)-O), thinned obs gradient. The current
  !< state is held in the this%analysis
  FUNCTION c4DVarAdjoint(this, X, dt) RESULT(Y)
    CLASS(c4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(IN), OPTIONAL :: dt
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: it, iv

    ! Initialize Y:
    Y = X

    IF (.NOT. PRESENT(dt)) THEN
      WRITE (*, 1)
1     FORMAT('c4DVarAdjoint -- A 4DVar requires a model integration time step length: dt. Please check and rerun')
      STOP
    END IF

    ! Adjoint model integration:
    PRINT*,'c4DVar adjoint integrating...'
    DO it = X%sg%tSlots, 1, -1
      PRINT*,'Adjoint a step backward...'
      ! CALL this%dycore%tim%TimeIntegrationBase_adj(dt, it, this%analysis(X%sg%gLevel), Y)
      ! PRINT*,'Adjoint a step backward completed'

      ! Add the 4DVar forcing: thinned obs saved in X time frames
      DO iv = 1, UBOUND(this%analysis%fields, 1)
        Y%fields(iv)%DATA(:, :, it) = Y%fields(iv)%DATA(:, :, it) + X%fields(iv)%DATA(:, :, it)
      END DO
    END DO
  END FUNCTION c4DVarAdjoint

  ! Tangent linear:
  FUNCTION c4DVarTangent(this, X) RESULT(Y)
    CLASS(c4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION c4DVarTangent

  ! Backward mapping: this is not a part of a 4DVar
  FUNCTION c4DVarBckward(this, X) RESULT(Y)
    CLASS(c4DVar_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    Y = X
  END FUNCTION c4DVarBckward
END MODULE c4DVar_m
