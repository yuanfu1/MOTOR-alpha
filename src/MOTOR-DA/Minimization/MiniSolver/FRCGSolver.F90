!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Tu Sanshan, 2021/1/26, @GBA-MWF, Shenzhen as SolverFRCG.F90
!   Modified by Zilong Qin (zilong.qin@gmail.com) and Yuanfu Xie, 2022/3/13, @GBA-MWF, Shenzhen
!     Add constraints in the Minimization, change the search strategy of step length, from recursively searching to only one try (1/2).
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022-10-08, @GBA-MWF, Shenzhen
!     Add output of components of Jo for debugging purpose.
!   Modified and changed by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024-05-13, @GBA-MWF, Shenzhen
!     Change to FRCGSolver.F90 as introducing lower bounds for controls and calling
!     a new function and gradient module costFuncGrad.F90
!!--------------------------------------------------------------------------------------------------

!> @brief
!  This is a version of Fletcher Reeves Conjugate Gradient method adapted with lower bounds
!  constraints using MOTOR-DA data structure.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE FRCGSolver_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE variationalSolvers_m, ONLY: variationalSolvers_t

  TYPE, EXTENDS(variationalSolvers_t) :: FRCGSolver_t
  CONTAINS
    PROCEDURE, PRIVATE :: boundProjection_s
    PROCEDURE, PUBLIC, PASS :: solve_s
  END TYPE FRCGSolver_t

CONTAINS
  ! Fletcher and Reeves CG solver:
  ! This is a version of scaled FRCG solver: See the reference paper:
  ! "painless-conjugate-gradient" page 32, eqns: 45-49
  SUBROUTINE solve_s(this, ctlBkgd, weights, initGuess, cost, sg, iterations)
    USE State_m, ONLY: State_t
    USE costFuncGrad_m, ONLY: costFuncGrad_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    USE parameters_m, ONLY: machineEps

    IMPLICIT NONE

    CLASS(FRCGSolver_t) :: this
    TYPE(State_t), INTENT(IN) :: ctlBkgd       ! Assuming bkgd is in control variable space
    REAL(r_kind), INTENT(IN) :: weights(:)  ! Weights for hybrid bkgd covariances
    TYPE(State_t), INTENT(INOUT) :: initGuess
    TYPE(costFuncGrad_t), INTENT(INOUT) :: cost
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), OPTIONAL :: iterations

    ! Local variables:
    INTEGER(i_kind) :: iflag, i, it, numControls, numSearch, MaxIterations
    INTEGER(i_kind) :: maxSearch = 20   ! Maximum line search allowed
    REAL(r_kind) :: func(2)     ! cost function values at current and plus steps
    REAL(r_kind) :: gxg(2)      ! Norms of gradients at current and plus steps
    REAL(r_kind) :: alpha       ! Line search and conjugate parameters
    TYPE(state_t) :: grad(2)    ! gradients of the cost function at current and plus steps
    TYPE(state_t) :: CNew       ! New iteration point
    REAL(r_kind) :: t1, t2 ! for timing purpose

    ! Debug:
    REAL(r_kind) :: rng_b = 0.45994843834260685D0

    IF (.NOT. sg%isActiveProc()) RETURN ! No action needed

    WRITE (*, 1) sg%gLevel, sg%mpddGlob%myrank
1   FORMAT('Starting FRCG solver at G', I2, ' processor: ', I6)

    ! Numbers of control variables:
    numControls = 0
    DO i = LBOUND(initGuess%fields, 1), UBOUND(initGuess%fields, 1)
      numControls = numControls + initGuess%fields(i)%numTotalMask
    END DO

    ! Initial step:
    iflag = 2   ! Evaluate both function and gradient
    PRINT*,'Evaluating the cost...'
    CALL cost%costFuncGrad_s(ctlBkgd, weights, initGuess, func(1), grad(1), iflag)
    gxg(1) = grad(1) .DOT.grad(1)

    MaxIterations = this%MaxOptStep
    IF (PRESENT(iterations)) THEN
      IF (iterations > MaxIterations) MaxIterations = iterations
    END IF
    IF (sg%isBaseProc()) WRITE (*, 2) func(1), SQRT(gxg(1)), sg%gLevel
2   FORMAT('At minimization step 0, f = ', E14.8, ' |g| = ', E14.8, ' Glevel: ', I3)

    IF (gxg(1) .LT. machineEps) RETURN    ! The current point is an minimizer

    ! Conjugate iterations:
    ! See Fletcher and Reeves 1963 paper: under references/da/optimization
    DO it = 1, MaxIterations - 1 ! Match with old SolverFRCG.F90

      ! Calculate the scale factor g^T A g:
      iflag = 3
      CNew = ctlBkgd%zeroCopy()
      CALL cost%costFuncGrad_s(ctlBkgd, weights, grad(1), func(2), grad(2), iflag)

      CNew = initGuess - grad(1) * (gxg(1) / func(2))   ! Grad(1) holds the current conjugate direction

      ! Project onto the simple bound constraints:
      CALL this%boundProjection_s(ctlBkgd, CNew, sg)

      ! Evaluate the cost and gradient at the new point:
      iflag = 0   ! Evaluate the cost function value only
      CALL cost%costFuncGrad_s(ctlBkgd, weights, CNew, func(2), grad(2), iflag)

      PRINT *, 'rng_b: ', rng_b, ' f_rng_b:', func(2), func(1)

      ! Line search if necessary:
      IF (func(2) .GE. func(1)) THEN
        alpha = 1.0D0
        numSearch = 0
        DO WHILE (func(2) .GE. func(1))
          alpha = alpha * 0.618D0
          CNew = initGuess - grad(1) * alpha

          ! Project onto the simple bound constraints:
          CALL this%boundProjection_s(ctlBkgd, CNew, sg)

          CALL cost%costFuncGrad_s(ctlBkgd, weights, CNew, func(2), grad(2), iflag)
          numSearch = numSearch + 1
          IF (numSearch .GT. maxSearch) THEN
            WRITE (*, 3) numSearch, maxSearch
3           FORMAT('FRCG line search failed after ', I3, &
                   ' exceeding ', I3, ' Use current')
            RETURN  ! Take the current initGuess value as a solution
          END IF
        END DO
        IF (sg%isBaseProc()) WRITE (*, 4) sg%gLevel, numSearch, it, func(2), func(1), alpha
4       FORMAT('FRCG line search: G', I3, ' numSearch: ', I3, &
               ' iter: ', I4, ' funcs: ', 2E14.5, ' Step: ', E12.4)
      END IF

      ! Accepting this new point:
      initGuess = CNew
      func(1) = func(2)
      iflag = 1   ! Calculate gradient only
      CALL cost%costFuncGrad_s(ctlBkgd, weights, initGuess, func(2), grad(2), iflag)
      gxg(2) = grad(2) .DOT.grad(2)
      alpha = gxg(2) / gxg(1)
      grad(1) = grad(2) + grad(1) * alpha
      gxg(1) = gxg(2)

      IF (sg%isBaseProc()) &
        WRITE (*, 5) sg%gLevel, it, func(1), SQRT(gxg(1)), alpha
5     FORMAT('Minimizing G', I1, ' step ', I3, ' f', E14.6, &
             ' |g|', E14.6, ' beta', E14.6)

      ! Minimizer is found:
      IF (gxg(1) .LT. machineEps) RETURN
    END DO
  END SUBROUTINE solve_s

  SUBROUTINE boundProjection_s(this, bkgd, C, sg)
    USE State_m, ONLY: State_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    USE costFuncGrad_m, ONLY: costFuncGrad_t

    IMPLICIT NONE

    CLASS(FRCGSolver_t) :: this
    TYPE(State_t), INTENT(IN) :: bkgd       ! Assuming bkgd is in control variable space
    TYPE(State_t), INTENT(INOUT) :: C
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    ! Local variables:
    INTEGER(i_kind) :: i, tIdx, hIdx, vIdx

    ! Apply low bound onstraints
    DO i = LBOUND(C%fields, 1), UBOUND(C%fields, 1)
      IF (C%fields(i)%lowBound .GT. 0) THEN
        FORALL (tIdx=1:sg%tSlots, vIdx=1:sg%vLevel, hIdx=1:sg%num_cell, &
                C%Fields(i)%DATA(vIdx, hIdx, tIdx) .LT. -bkgd%fields(i)%DATA(vIdx, hIdx, tIdx))
          C%Fields(i)%DATA(vIdx, hIdx, tIdx) = -bkgd%Fields(i)%DATA(vIdx, hIdx, tIdx)
        END FORALL
      END IF
    END DO
  END SUBROUTINE boundProjection_s
END MODULE FRCGSolver_m
