!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.possionSolver_Kp.psSolutn
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, in grids/gzGrid.F90
!     Created by Yuanfu Xie, 2019/04, @GBA-MWF, Shenzhen
!     Modified by Yuanfu Xie 2020-4 for adding boundary conditions
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/17, @GBA-MWF, Shenzhen, add parallelization
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE psSolutn_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  TYPE psSolutn_t
    REAL(r_kind), ALLOCATABLE :: rights(:, :)
    REAL(r_kind), ALLOCATABLE :: solutn(:, :)
  CONTAINS
    FINAL :: destructor
  END TYPE psSolutn_t

  INTERFACE psSolutn_t
    PROCEDURE :: constructor
  END INTERFACE psSolutn_t

CONTAINS

  FUNCTION constructor(vLevel, num_cell) RESULT(this)
    IMPLICIT NONE
    TYPE(psSolutn_t) :: this
    INTEGER(i_kind), INTENT(in) :: vlevel, num_cell

    ALLOCATE (this%rights(vlevel, num_cell), &
              this%solutn(vlevel, num_cell))

  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(psSolutn_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%rights)) DEALLOCATE (this%rights)
    IF (ALLOCATED(this%solutn)) DEALLOCATE (this%solutn)

  END SUBROUTINE destructor

END MODULE psSolutn_m
