!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/21, @GBA-MWF, Shenzhen
!     Changed its original name to variationalSolver.F90 for separate testings
!     purpose.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE variationalSolvers_m
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE variationalSolvers_t
    INTEGER(i_kind) :: MaxOptStep = 10

  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    FINAL :: destructor
    PROCEDURE, PRIVATE :: boundProjection_s
    PROCEDURE, PUBLIC :: solve1_s
  END TYPE variationalSolvers_t

CONTAINS

  SUBROUTINE initialize_s(this, configFile)
    USE YAMLRead_m
    IMPLICIT NONE
    CLASS(variationalSolvers_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind) :: ifile

    ifile = yaml_get_var(TRIM(configFile), 'Minimization', 'MaxOptStep', this%MaxOptStep)
  END SUBROUTINE initialize_s

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(variationalSolvers_t), INTENT(INOUT) :: this
  END SUBROUTINE destructor

  SUBROUTINE boundProjection_s(this, bkgd, C, sg)
    USE State_m, ONLY: State_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    USE costFuncGrad_m, ONLY: costFuncGrad_t

    IMPLICIT NONE

    CLASS(variationalSolvers_t) :: this
    TYPE(State_t), INTENT(IN) :: bkgd       ! Assuming bkgd is in control variable space
    TYPE(State_t), INTENT(INOUT) :: C
    TYPE(SingleGrid_t), INTENT(IN) :: sg

  END SUBROUTINE boundProjection_s

  SUBROUTINE solve1_s(this, ctlBkgd, weights, initGuess, cost, sg, iterations)
    USE State_m, ONLY: State_t
    USE SingleGrid_m, ONLY: SingleGrid_t
    USE costFuncGrad_m, ONLY: costFuncGrad_t

    IMPLICIT NONE

    CLASS(variationalSolvers_t) :: this
    TYPE(State_t), INTENT(IN) :: ctlBkgd    ! Assuming bkgd is in control variable space
    REAL(r_kind), INTENT(IN) :: weights(:)  ! Weights for hybrid bkgd covariances
    TYPE(State_t), INTENT(INOUT) :: initGuess
    TYPE(costFuncGrad_t), INTENT(INOUT) :: cost
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind), OPTIONAL :: iterations

    ! Local variables:
    INTEGER(i_kind) :: i,j,k,numIterations
    REAL(r_kind) :: costValue
    TYPE(State_t) :: gradValue

    PRINT*, 'Starting variational solver...',sg%gLevel

    IF (PRESENT(iterations)) THEN
      numIterations = iterations
    ELSE
      numIterations = this%MaxOptStep ! Default to MaxOptStep if not provided
    END IF

    ! Minimize the cost function using a simple gradient descent method
    gradValue = initGuess  ! Initialize gradient value
    DO i = 1, numIterations
      ! Calculate the cost function and its gradient
      CALL cost%costFuncGrad_s(ctlBkgd, weights, cost%C%analysis, costValue, gradValue, 2)
    END DO

  END SUBROUTINE solve1_s

END MODULE variationalSolvers_m
