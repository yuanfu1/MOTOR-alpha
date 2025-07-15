MODULE rhsAdjointBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t

  TYPE, ABSTRACT :: rhsAdjointBase_t
    INTEGER(i_kind) :: num_eqns
    INTEGER(i_kind), ALLOCATABLE :: variableIdx(:)

  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE(adj_model), DEFERRED :: rightHandSideAdjoint
  END TYPE rhsAdjointBase_t

  ABSTRACT INTERFACE
    SUBROUTINE adj_model(this, adj_X, it, adj_current, adj_forward)
      IMPORT :: rhsAdjointBase_t, State_t, i_kind, r_kind
      CLASS(rhsAdjointBase_t) :: this
      TYPE(State_t), INTENT(IN) :: adj_X
      INTEGER(i_kind), INTENT(IN) :: it
      REAL(r_kind), INTENT(IN) :: adj_current(adj_X%sg%vLevel, adj_X%sg%num_cell, UBOUND(adj_X%fields, 1))
      REAL(r_kind), INTENT(OUT) :: adj_forward(adj_X%sg%vLevel, adj_X%sg%num_cell, UBOUND(adj_X%fields, 1))
    END SUBROUTINE adj_model
  END INTERFACE

CONTAINS
  SUBROUTINE initialize_s(this)
    CLASS(rhsAdjointBase_t) :: this
    this%num_eqns = 1
    ALLOCATE (this%variableIdx(1))
  END SUBROUTINE initialize_s
END MODULE rhsAdjointBase_m
