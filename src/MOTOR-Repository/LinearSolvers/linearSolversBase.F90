!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.linearSolversBase.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2025/05/19, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines the base class for linear solvers in the MOTOR framework.
!! It provides an abstract interface for loading matrices and initializing the solver.
MODULE linearSolversBase_m
  USE , INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE State_m, ONLY: State_t

  TYPE, ABSTRACT :: linearSolversBase_t
    INTEGER(c_int) :: nonzero  ! Number of non-zero elements in the matrix
    INTEGER(c_int), ALLOCATABLE :: irow(:), jcol(:)  ! Row and column indices for the matrix
    REAL(c_double), ALLOCATABLE :: sln(:,:),rhs(:,:) ! Vertical and horizontal solution/rhs
    REAL(c_double), ALLOCATABLE :: matrix(:,:)       ! Matrix for the linear solver loaded in vertical and horizontal
    TYPE(State_t), POINTER :: X
  CONTAINS
    ! Abstract methods to be implemented by derived types
    PROCEDURE(init), DEFERRED :: initialize_s
    PROCEDURE(load), DEFERRED :: loadMatrix_s
    PROCEDURE, PUBLIC :: solver_s
  END TYPE linearSolversBase_t

  ABSTRACT INTERFACE
    SUBROUTINE init(this, X)
      IMPORT :: c_int, linearSolversBase_t, State_t
      CLASS(linearSolversBase_t) :: this
      TYPE(State_t), TARGET, INTENT(IN) :: X
    END SUBROUTINE init

    SUBROUTINE load(this)
      IMPORT :: c_int, linearSolversBase_t, State_t
      CLASS(linearSolversBase_t) :: this
    END SUBROUTINE load

    ! SUBROUTINE solv(this, ib, ix, X)
    !   IMPORT :: c_int, linearSolversBase_t, State_t
    !   CLASS(linearSolversBase_t) :: this
    !   TYPE(State_t), INTENT(IN) :: X
    !   INTEGER(c_int), INTENT(IN) :: ib, ix  ! ib: index of right hand side, ix: index of the solution stored
    !   ! The solver should return the solution in this%sln
    ! END SUBROUTINE solv

  END INTERFACE

  CONTAINS
    SUBROUTINE solver_s(this, ib, ix, X, it)
      CLASS(linearSolversBase_t) :: this
      TYPE(State_t), INTENT(IN) :: X
      INTEGER(c_int), INTENT(IN) :: ib, ix  ! ib: index of right hand side, ix: index of the solution stored
      INTEGER(c_int), INTENT(IN) :: it  ! 0 Ax=b; 1 ATx=b
      ! The solver should return the solution in this%sln
      ! Implement the solver logic here
    END SUBROUTINE solver_s

END MODULE linearSolversBase_m