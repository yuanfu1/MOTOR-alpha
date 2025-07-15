!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE MiniSolver_m
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE JFunc_m, ONLY: JFunc_t
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE MiniSolver_t
    INTEGER(i_kind) :: MaxOptStep = 10
    REAL(r_kind) :: J

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC, PASS :: run
  END TYPE MiniSolver_t

CONTAINS

  SUBROUTINE initialize(this, configFile)
    USE YAMLRead_m
    IMPLICIT NONE
    CLASS(MiniSolver_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind) :: ifile

    ifile = yaml_get_var(TRIM(configFile), 'Minimization', 'MaxOptStep', this%MaxOptStep)
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(MiniSolver_t), INTENT(INOUT) :: this
  END SUBROUTINE destructor

  SUBROUTINE run(this, X, JFunc, sg, iters)
    IMPLICIT NONE
    CLASS(MiniSolver_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(JFunc_t), INTENT(INOUT) :: JFunc
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind), OPTIONAL :: iters

    ! TYPE(State_t) :: Xb
  END SUBROUTINE run

END MODULE MiniSolver_m
