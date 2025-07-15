MODULE TimeIntegrationBaseAdjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE rhsBaseAdjoint_m, ONLY: rhsBaseAdjoint_t
  IMPLICIT NONE

  TYPE, ABSTRACT :: TimeIntegrationBaseAdjoint_t
    CLASS(rhsBaseAdjoint_t), POINTER :: rhsAdjoint
  CONTAINS
    PROCEDURE(initl_adj), DEFERRED :: initializationAdjoint_s
    PROCEDURE(match_adj), DEFERRED :: TimeIntegrationBaseAdjoint_s
  END TYPE TimeIntegrationBaseAdjoint_t

  ABSTRACT INTERFACE
    SUBROUTINE initl_adj(this, rhsAdjoint, X)
      IMPORT :: TimeIntegrationBaseAdjoint_t, State_t, rhsBaseAdjoint_t
      CLASS(TimeIntegrationBaseAdjoint_t), INTENT(INOUT) :: this
      CLASS(rhsBaseAdjoint_t), TARGET, INTENT(INOUT) :: rhsAdjoint
      TYPE(State_t), INTENT(IN) :: X
    END SUBROUTINE initl_adj

    SUBROUTINE match_adj(this, it, dt, adj_X)
      IMPORT :: TimeIntegrationBaseAdjoint_t, State_t, i_kind, r_kind
      CLASS(TimeIntegrationBaseAdjoint_t), INTENT(INOUT) :: this
      INTEGER(i_kind), INTENT(IN) :: it
      REAL(r_kind), INTENT(IN) :: dt
      TYPE(State_t), INTENT(INOUT) :: adj_X
    END SUBROUTINE match_adj
  END INTERFACE
END MODULE TimeIntegrationBaseAdjoint_m
