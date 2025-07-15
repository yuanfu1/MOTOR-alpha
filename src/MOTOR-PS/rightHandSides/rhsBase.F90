!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.rightHandSide.rhsBase
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-11   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an abstract data structure of right hand side module. It allows users switch
!! one model to another without knowning the implementation details.
MODULE rhsBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  TYPE, ABSTRACT :: rhsBase_t
    INTEGER(i_kind) :: num_eqns
    INTEGER(i_kind), ALLOCATABLE :: variableIdx(:)
    TYPE(poissonSolver_t), POINTER :: poissonSolver
  CONTAINS
    ! Initialize the base class
    PROCEDURE, PUBLIC :: initialize_s
    ! Deferred procedure for the right-hand side calculation
    PROCEDURE(model), DEFERRED, PUBLIC :: rightHandSide
  END TYPE rhsBase_t

  ! Define the abstract interface for the right-hand side calculation (model)
  ABSTRACT INTERFACE
    SUBROUTINE model(this, it, current, forward)
      IMPORT :: rhsBase_t, State_t, i_kind, r_kind
      CLASS(rhsBase_t) :: this
      INTEGER(i_kind), INTENT(IN) :: it
      TYPE(State_t), INTENT(IN)  :: current
      TYPE(State_t), INTENT(INOUT) :: forward
    END SUBROUTINE model
  END INTERFACE

CONTAINS

  ! Initialize the base class with default values for num_eqns and variableIdx
  SUBROUTINE initialize_s(this, sg, psolver)
    IMPLICIT NONE
    CLASS(rhsBase_t) :: this   ! INTENT(INOUT) matches derived classes
    TYPE(poissonSolver_t), TARGET :: psolver  ! Optional Poisson solver
    TYPE(SingleGrid_t), INTENT(IN) :: sg  ! Single grid information for initializing a model

    ! Set the default number of equations and allocate variableIdx
    this%num_eqns = 5
    ALLOCATE (this%variableIdx(1))
    PRINT *, 'Exiting rhsBase initialization'

    this%poissonSolver => psolver  ! Assign the Poisson solver if provided
  END SUBROUTINE initialize_s

END MODULE rhsBase_m
