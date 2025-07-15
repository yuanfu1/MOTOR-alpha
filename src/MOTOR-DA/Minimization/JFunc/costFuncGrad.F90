!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/11, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a new version of MOTOR-DA variational analysis cost function and gradient. It improves
!! the original JFunc.F90 in several aspects for simplicity and unified routines.
MODULE costFuncGrad_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE C2MBase_m, ONLY: C2MBase_t
  USE c4DVar_m, ONLY: c4DVar_t
  USE C2O_m, ONLY: C2O_t
  USE ObsSetArray_m, ONLY: ObsSetArray_t
  USE obsMG_m, ONLY: obsMG_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsUtilities_m

  TYPE costFuncGrad_t
    TYPE(BMatrix_t), POINTER :: L
    TYPE(BMatrix_t), POINTER :: B
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(State_t), POINTER :: X
    CLASS(C2MBase_t), POINTER :: C
    TYPE(ObsMG_t), POINTER :: Y
    TYPE(C2O_t), POINTER :: H
    TYPE(State_t), POINTER :: Xb
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    FINAL :: destructor_s

    PROCEDURE, PUBLIC :: costFuncGrad_s
  END TYPE costFuncGrad_t

CONTAINS
  SUBROUTINE initialize_s(this, configFile, X, C, Y, H, L, B, sg, Xb)
    IMPLICIT NONE
    CLASS(costFuncGrad_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET :: X
    CLASS(C2MBase_t), TARGET :: C
    TYPE(obsMG_t), TARGET :: Y
    TYPE(C2O_t), TARGET :: H
    TYPE(BMatrix_t), TARGET :: L  ! Temporarily hold for Laplace covariance
    TYPE(BMatrix_t), TARGET :: B
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), TARGET :: Xb

    ! Local variables:
    INTEGER(i_kind) :: istatus

    this%X => X
    this%C => C
    this%Y => Y
    this%H => H
    this%L => L
    this%B => B
    this%sg => sg
    this%Xb => Xb

  END SUBROUTINE initialize_s

  ! Destructor:
  IMPURE ELEMENTAL SUBROUTINE destructor_s(this)
    IMPLICIT NONE
    TYPE(costFuncGrad_t), INTENT(INOUT) :: this

    IF (this%sg%mpddGlob%isActiveProc()) WRITE (*, 1)
1   FORMAT('End of destructor of costFuncGrad_t')
  END SUBROUTINE destructor_s

  ! A unified cost function and gradient calculation:
  ! flag is choice of computing function and/or gradient
  ! The cost function formulation can be seen in
  ! The "MOTOR-DA General Form of Analysis" doc under
  ! $(HOME)/etc/doc
  ! The cost function is in form of J(W)
  SUBROUTINE costFuncGrad_s(this, ctlBkgd, weights, W, cost, grad, iflag)
    IMPLICIT NONE
    CLASS(costFuncGrad_t) :: this
    TYPE(State_t), TARGET :: ctlBkgd
    REAL(r_kind), INTENT(IN) :: weights(3)  ! Laplacian, ensemble and static background covariances
    TYPE(State_t), INTENT(IN) :: W          ! Control increment W = C - C_b
    REAL(r_kind), INTENT(OUT) :: cost       ! Cost function value
    TYPE(State_t), INTENT(OUT) :: grad ! B^{-1/2} (C-C_b) and gradient
    INTEGER(i_kind), INTENT(IN) :: iflag     ! choice of calculation
    ! flag = 0 for computing function value only;
    ! flag = 1 for computing gradient value only;
    ! flag = 2 for both function and gradient;
    ! flag = 3 for computing conjugate parameter;

    ! Local variables:
    INTEGER(i_kind) :: iotype
    TYPE(State_t) :: BHalfW, XM, workState              ! Laplacian times control increment L(C-C_b)
    TYPE(ObsSetArray_t) :: workObs, obsTemps

    ! Initialize the cost output:
    cost = 0.0D0

    ! Note: Historical reason: all the operators SQRTINVMUL are essentially SQRTMUL, no inverse
    ! This may be modified in the future noted by Yuanfu Xie

    ! Compute B^{1/2} W:
    IF (weights(2) .GT. 0.0D0) THEN
      BHalfW = this%B.SQRTINVMUL.W
    ELSE
      BHalfW = W
    END IF

    !> There are total three terms in the general form of variational cost function,
    !> Jo, Jb and Jl, where Jl is the Laplace terms. There is possible weak constraint term
    !> that will be considered later. Yuanfu Xie, 2024-05-28

    ! 1. Jo term:
    IF (iflag .LE. 2) THEN
    print*,'Starting Jo0',cost, (BHalfW .DOT. BHalfW), (ctlBkgd .DOT. ctlBkgd), weights
      XM = this%C%forwardOpr_s(ctlBkgd + BHalfW)         ! X_m(C_b + B^{1/2} W) back to model space

      ! H(X) - O:
    print*,'Starting Jo1',cost, iflag, (this%Y%thinnedObs.DOT.this%Y%thinnedObs)
      obsTemps = this%Y%forwardOpr_s(XM) !- this%Y%thinnedObs
      PRINT*,'obsTemps', (obsTemps .DOT.obsTemps), (XM .DOT.XM), (this%Y%thinnedObs .DOT.this%Y%thinnedObs)
      ! R^(-1/2) (H(x) - O):
      CALL this%Y%rSqrtInverseMultiply(obsTemps, workObs)
    print*,'Starting Jo2',cost, (obsTemps .DOT. obsTemps), (workObs.DOT.workObs)
    ELSE IF (iflag .EQ. 3) THEN
      XM = this%C%forwardOpr_s(BHalfW)            ! Perturbation and Hessian product approximation
      ! Used by the FRCG for conjudate coefficient gAg
      CALL this%Y%rSqrtInverseMultiply(this%Y%forwardOpr_s(XM), workObs)
    ELSE
      WRITE (*, 2) iflag
2     FORMAT('costFuncGrad - iflag parameter is not supported, please check and rerun! ', I2)
      STOP
    END IF

    print*,'Jo 1',cost,iflag

    ! Here: workObs = R^(-1/2) (H(x) - O) holds at this step
    IF (iflag .NE. 1) THEN
      cost = (workObs.DOT.workObs)
    END IF

    IF (iflag .EQ. 1 .OR. iflag .EQ. 2) THEN
      ! H^T (H(X) - O):
      obsTemps = this%Y.SQRTINVMULADJ.workObs

    print*,'Jo 2',cost

      grad = this%Y%adjointOpr_s(obsTemps, XM)

    print*,'Jo 3',cost

      ! \nable Xm^T (grad):
      SELECT TYPE (p => this%C)  ! Check the type of C
      TYPE IS (c4DVar_t)
        ! For c4DVar, we need to multiply by dt
        PRINT*,'Calculating adjoint for c4DVar_t', p%dycore%dt, grad%sg%gLevel
        grad = this%C%adjointOpr_s(grad, p%dyCore%dt)
      CLASS DEFAULT
        ! For other types, no need to multiply by dt
        grad = this%C%adjointOpr_s(grad)
      END SELECT

    print*,'Jo 4',cost

    END IF

    ! 2. Adding W^T W term: Only when weights(2) > 0
    IF (weights(2) .GT. 0.0D0) THEN
      IF (iflag .EQ. 0 .OR. iflag .GE. 2) &
        cost = cost + (W.DOT.W) * weights(2)
      IF (iflag .EQ. 1 .OR. iflag .EQ. 2) &
        grad = grad + ((this%B.SQRTINVMUL.grad) + W) * weights(2)
      ! Note this needs to be replaced
      ! by a .SQRTMUL. operator instead of SQRTINVMUL. Zilong is adding them in and will test
    END IF

    print*,'Jo 5',cost

    ! 3. Laplace term: Only when weights(1) > 0
    IF (weights(1) .GT. 0.0D0) THEN
      ! Note the current implementation applies Laplacian operators in X space instead of control
      ! Therefore, it requires an additional call of the this%C%forwardOpr_s(BHalfW):
      BHalfW = this%L.SQRTINVMUL. (this%C%forwardOpr_s(BHalfW))  ! L B^{1/2} W

      IF (iflag .EQ. 1 .OR. iflag .EQ. 2) THEN
        workState = this%L%sqrt_inverse_multiply_adjoint(BHalfW, XM)
        IF (weights(1) .GT. 0.0D0) &
          grad = grad + workState * weights(1)
        IF (weights(2) .GT. 0.0D0) THEN
          grad = grad + this%B%sqrt_inverse_multiply_adjoint(workState, XM) * weights(2)
          ! Note this needs to be replaced
          ! by a B%sqrt_multiply_adjoint operator instead of sqrt_inverse_multiply_adjoint
          ! Zilong is adding them in and will test
        END IF
      END IF

      IF (iflag .NE. 1) cost = cost + weights(1) * (BHalfW.DOT.BHalfW)
    END IF

    ! Note iflag .EQ.3, there is no need to multiply 0.5D0 as it is a Hessian product approximation
    IF (iflag .EQ. 0 .OR. iflag .EQ. 2) cost = 0.5D0 * cost
  END SUBROUTINE costFuncGrad_s
END MODULE costFuncGrad_m
