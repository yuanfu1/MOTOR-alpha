!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.TimeIntegration.TimeIntegrationBase
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-11   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an abstract data structure of time integration schemes. It allows users switch
!! time integration schemes without knowning the implementation details.
MODULE TimeIntegrationBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  ! USE C2MBase_m, ONLY: C2MBase_t
  USE obsMG_m, ONLY: obsMG_t
  ! USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE rhsBase_m, ONLY: rhsBase_t
  IMPLICIT NONE

  TYPE, ABSTRACT :: TimeIntegrationBase_t
    CLASS(rhsBase_t), POINTER :: rhs
    CHARACTER(LEN=1024) :: yamlFile
    TYPE(State_t) :: ab3_adj
  CONTAINS
    PROCEDURE(initl), DEFERRED :: initialization_s
    PROCEDURE(match), DEFERRED :: TimeIntegrationBase_fwd
    PROCEDURE(adjnt), DEFERRED :: TimeIntegrationBase_adj
  END TYPE TimeIntegrationBase_t

  ABSTRACT INTERFACE
    SUBROUTINE initl(this, rhs, X, yamlFile)
      IMPORT :: TimeIntegrationBase_t, State_t, rhsBase_t
      CLASS(TimeIntegrationBase_t) :: this
      CLASS(rhsBase_t), TARGET :: rhs
      TYPE(State_t), INTENT(IN) :: X
      CHARACTER(LEN=1024), INTENT(IN) :: yamlFile
    END SUBROUTINE initl
    SUBROUTINE match(this, dt, int, it, iC, X, yk)
      IMPORT :: TimeIntegrationBase_t, State_t, i_kind, r_kind
      CLASS(TimeIntegrationBase_t) :: this
      REAL(r_kind), INTENT(IN) :: dt            !< Time step length
      INTEGER(i_kind), INTENT(IN) :: int, it    !< intermedimum step and Current model time frame
      !< IC holds the model states of a given time integration scheme, e.g.,
      !< a AB3, its first 3 time frames of IC holds the 3 model states needed
      !< for AB3 in it-3*dt, it-2*dt, it-dt to update it. Noted by Yuanfu Xie
      !< It may be initialized as the number of time stepping scheme, AB3 for 3
      !< RK4 for 1 to save memory. Currently, we use the same time frame as X.
      TYPE(State_t), INTENT(INOUT) :: iC        !< Model initial condition of the current time frame
      TYPE(State_t), INTENT(INOUT) :: X         !< Current model states
      TYPE(State_t), OPTIONAL, INTENT(INOUT) :: yk     !< For multiple time stepping schemes, save the yk: F(Xk) for dX/dt = F(X)
    END SUBROUTINE match

    ! The adjoint only marches the adjoint variables between time frame of a state
    ! since there is no observation forcing term. The 4DVar in C2M will add the forcing
    SUBROUTINE adjnt(this, dt, it, X, Xadj)
      IMPORT :: TimeIntegrationBase_t, State_t, obsMG_t, &
        BMatrix_t, i_kind, r_kind
      CLASS(TimeIntegrationBase_t) :: this
      REAL(r_kind), INTENT(IN) :: dt            !< Time step length
      INTEGER(i_kind), INTENT(IN) :: it         !< the current time frame of a state
      TYPE(State_t), INTENT(IN) :: X            !< Current model states
      TYPE(State_t), INTENT(INOUT) :: Xadj      !< the adjoint model states
    END SUBROUTINE adjnt
  END INTERFACE
END MODULE TimeIntegrationBase_m
