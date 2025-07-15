!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/10/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines abstract dynamic core
!!
MODULE dyCoresBase_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  TYPE, ABSTRACT :: dyCoresBase_t
    INTEGER(i_kind) :: num_eqns
    INTEGER(i_kind), ALLOCATABLE :: variableIdx(:)
    CLASS(TimeIntegrationBase_t), POINTER :: tim    !< Time integration scheme
    CLASS(rhsBase_t), POINTER :: rhs                !< DU/Dt = F, rhs = F
    CLASS(SingleGrid_t), POINTER :: sg               !< Geometry of the model domain
    CHARACTER(LEN=1024) :: yamlFile
    INTEGER(i_kind) :: intval     !< Number of time steps between two adjacent time frame of State
    REAL(r_kind) :: dt  !< Time integration steplength
  CONTAINS
    ! Initialize the base class
    PROCEDURE(init), DEFERRED, PUBLIC :: initialize_s
    ! Deferred procedure for the right-hand side calculation
    PROCEDURE(model), DEFERRED, PUBLIC :: dycore_s
    PROCEDURE(adjoint), DEFERRED, PUBLIC :: dyAdjoint_s
  END TYPE dyCoresBase_t

  ! Define the abstract interface for the right-hand side calculation (model)
  ABSTRACT INTERFACE
    SUBROUTINE init(this, yamlFile, tim, rhs, sg, dt)
      IMPORT :: dyCoresBase_t, TimeIntegrationBase_t, rhsBase_t, i_kind, r_kind, SingleGrid_t
      CLASS(dyCoresBase_t) :: this
      CHARACTER(LEN=1024), INTENT(IN) :: yamlFile
      CLASS(TimeIntegrationBase_t), TARGET :: tim
      CLASS(rhsBase_t), TARGET :: rhs
      CLASS(SingleGrid_t), TARGET :: sg
      REAL(r_kind), INTENT(IN) :: dt
    END SUBROUTINE init

    SUBROUTINE model(this, iC, X)
      IMPORT :: dyCoresBase_t, State_t, i_kind, r_kind
      CLASS(dyCoresBase_t) :: this
      TYPE(State_t), INTENT(INOUT) :: iC
      TYPE(State_t), INTENT(INOUT) :: X
    END SUBROUTINE model

    SUBROUTINE adjoint(this, it, iC, Xadj)
      IMPORT :: dyCoresBase_t, State_t, i_kind, r_kind
      CLASS(dyCoresBase_t) :: this
      INTEGER(i_kind), INTENT(IN) :: it
      TYPE(State_t), INTENT(INOUT) :: iC
      TYPE(State_t), INTENT(INOUT) :: Xadj
    END SUBROUTINE adjoint
  END INTERFACE
END MODULE dyCoresBase_m
