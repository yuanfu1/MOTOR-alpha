MODULE rhsBaseAdjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t

  TYPE, ABSTRACT :: rhsBaseAdjoint_t
    INTEGER(i_kind) :: num_eqns
    INTEGER(i_kind), ALLOCATABLE :: variableIdx(:)
  CONTAINS
    ! Initialize the base class
    PROCEDURE, PUBLIC :: initialize_s
    ! Deferred procedure for the adjoint right-hand side calculation
    PROCEDURE(adj_model), DEFERRED, PUBLIC :: rightHandSideAdjoint
  END TYPE rhsBaseAdjoint_t

  ! Define the abstract interface for the adjoint right-hand side calculation
  ABSTRACT INTERFACE
    SUBROUTINE adj_model(this, X, it, adj_current, adj_forward)
      IMPORT :: rhsBaseAdjoint_t, State_t, i_kind, r_kind
      CLASS(rhsBaseAdjoint_t), INTENT(INOUT) :: this
      TYPE(State_t), INTENT(IN) :: X
      INTEGER(i_kind), INTENT(IN) :: it
      REAL(r_kind), TARGET, INTENT(IN) :: adj_current(X%sg%vLevel, X%sg%num_cell, UBOUND(X%fields, 1))  ! Declare adj_current as TARGET
      REAL(r_kind), INTENT(OUT) :: adj_forward(X%sg%vLevel, X%sg%num_cell, UBOUND(X%fields, 1))
    END SUBROUTINE adj_model
  END INTERFACE

CONTAINS

  ! Initialize the adjoint base class with default values for num_eqns and variableIdx
  SUBROUTINE initialize_s(this)
    IMPLICIT NONE
    CLASS(rhsBaseAdjoint_t), INTENT(INOUT) :: this

    ! Set the default number of equations and allocate variableIdx
    this%num_eqns = 5
    ALLOCATE (this%variableIdx(1))
  END SUBROUTINE initialize_s

END MODULE rhsBaseAdjoint_m
