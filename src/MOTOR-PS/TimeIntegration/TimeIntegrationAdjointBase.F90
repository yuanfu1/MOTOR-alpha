MODULE TimeIntegrationAdjointBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE rhsAdjointBase_m, ONLY: rhsAdjointBase_t

  TYPE, ABSTRACT :: TimeIntegrationAdjointBase_t
    CLASS(rhsAdjointBase_t), POINTER :: rhsAdj
  CONTAINS
    PROCEDURE(initlAdj), DEFERRED :: initialization_s
    PROCEDURE(matchAdj), DEFERRED :: TimeIntegrationAdjointBase_s
  END TYPE TimeIntegrationAdjointBase_t

  ABSTRACT INTERFACE
    SUBROUTINE initlAdj(this, rhsAdj, adj_X)
      IMPORT :: TimeIntegrationAdjointBase_t, State_t, rhsAdjointBase_t
      CLASS(TimeIntegrationAdjointBase_t) :: this
      CLASS(rhsAdjointBase_t), TARGET :: rhsAdj
      TYPE(State_t), INTENT(IN) :: adj_X
    END SUBROUTINE initlAdj

    SUBROUTINE matchAdj(this, it, dt, adj_X)
      IMPORT :: TimeIntegrationAdjointBase_t, State_t, i_kind, r_kind
      CLASS(TimeIntegrationAdjointBase_t) :: this
      INTEGER(i_kind), INTENT(IN) :: it
      REAL(r_kind), INTENT(IN) :: dt
      TYPE(State_t), INTENT(INOUT) :: adj_X
    END SUBROUTINE matchAdj
  END INTERFACE

END MODULE TimeIntegrationAdjointBase_m
