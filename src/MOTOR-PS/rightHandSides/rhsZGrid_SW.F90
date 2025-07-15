!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.rightHandSide.rhsZGrid_SW.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-12-02   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a shallow water equation system using a Z-grid model.
MODULE rhsZGrid_SW_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE rhsBase_m, ONLY: rhsBase_t
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE geometry_m, ONLY: geometry_t
  USE parameters_m, ONLY: Omega,g

  IMPLICIT NONE

  TYPE, EXTENDS(rhsBase_t) :: rhsZGrid_SW_t
    TYPE(gzm_t) :: gzm
    CONTAINS
      PROCEDURE, PUBLIC :: initialize_s => initializeZGrid_SW
      PROCEDURE, PUBLIC :: rightHandSide => rhsZGrid_SW_s
  END TYPE rhsZGrid_SW_t
  
  CONTAINS
  
    SUBROUTINE initializeZGrid_SW(this, sg, psolver)
      CLASS(rhsZGrid_SW_t) :: this
      TYPE(poissonSolver_t), TARGET :: psolver
      TYPE(SingleGrid_t), INTENT(IN) :: sg

      this%gzm = gzm_t(sg)
      this%poissonSolver => psolver
    END SUBROUTINE initializeZGrid_SW

    ! Assume the current states of absolute vorticity, divergence and height are given
    ! streamfunction and velocity potential are computed based on the absolute vorticity
    ! and divergence. Note relative vorticity is also calculated from the absolute one.
    SUBROUTINE rhsZGrid_SW_s(this, it, current, forward)
      CLASS(rhsZGrid_SW_t) :: this
      INTEGER(i_kind), INTENT(IN) :: it
      TYPE(State_t), INTENT(IN) :: current
      TYPE(State_t), INTENT(INOUT) :: forward

      ! Local variables:
      INTEGER(i_kind) :: i
      REAL(r_kind), ALLOCATABLE :: relativeVorticity(:,:),div(:,:),jac(:,:)

      ! 1. Relative voriticity:
      ALLOCATE(relativeVorticity(current%sg%vLevel,current%sg%num_cell))
      DO i=1,current%sg%vLevel
        relativeVorticity(i,:) = current%fields(1)%DATA(i,:,it)- &
          2.0D0*Omega*DSIN(current%sg%cell_cntr(1,:))
      END DO

      ! 2. Streamfunction and velocity potential:
      CALL this%poissonSolver%PoissonSol_sphere(current%sg%gLevel, &
        current%sg%vLevel, relativeVorticity, forward%fields(4)%DATA(:,:,it))
      CALL this%poissonSolver%PoissonSol_sphere(current%sg%gLevel, &
        current%sg%vLevel, current%fields(2)%DATA(:,:,it), forward%fields(5)%DATA(:,:,it))

      ! 3. RHS absolute vorticity:
      ALLOCATE(div(current%sg%vLevel,current%sg%num_cell),jac(current%sg%vLevel,current%sg%num_cell))
      CALL this%gzm%Divergen(current%fields(1)%DATA(:,:,it),forward%fields(5)%DATA(:,:,it),div)
      CALL this%gzm%Jacobian(current%fields(1)%DATA(:,:,it),forward%fields(4)%DATA(:,:,it),jac)
      forward%fields(1)%DATA(:,:,it) = jac-div

      ! 4. RHS divergence:
      ! 4.1 jac+div
      CALL this%gzm%Divergen(current%fields(1)%DATA(:,:,it),forward%fields(4)%DATA(:,:,it),div)
      CALL this%gzm%Jacobian(current%fields(1)%DATA(:,:,it),forward%fields(5)%DATA(:,:,it),jac)
      forward%fields(2)%DATA(:,:,it) = jac+div
      ! 4.2 Kinetic energy: temporarily saved in div:
      div = forward%fields(4)%DATA(:,:,it)*relativeVorticity + &
            forward%fields(5)%DATA(:,:,it)*current%fields(2)%DATA(:,:,it)
      ! Use jac as a temporary memory to save intermediate quantities:
      CALL this%gzm%Divergen(forward%fields(4)%DATA(:,:,it),forward%fields(4)%DATA(:,:,it),jac)
      div = jac - div
      CALL this%gzm%Divergen(forward%fields(5)%DATA(:,:,it),forward%fields(5)%DATA(:,:,it),jac)
      div = 0.5D0*(div + jac)
      CALL this%gzm%Jacobian(forward%fields(4)%DATA(:,:,it),forward%fields(5)%DATA(:,:,it),jac)
      div = div + jac
      ! 4.3 \nable^2 (K + gh)
      div = div + g*current%fields(3)%DATA(:,:,it)
      CALL this%gzm%Laplacia(div,jac)
      forward%fields(2)%DATA(:,:,it) = forward%fields(2)%DATA(:,:,it) - div

      ! 5. Height:
      CALL this%gzm%Divergen(current%fields(3)%DATA(:,:,it),forward%fields(5)%DATA(:,:,it),div)
      CALL this%gzm%Jacobian(current%fields(3)%DATA(:,:,it),forward%fields(4)%DATA(:,:,it),jac)
      forward%fields(3)%DATA(:,:,it) = jac - div

      ! Deallocate memory:
      DEALLOCATE(relativeVorticity,div,jac)

    END SUBROUTINE rhsZGrid_SW_s

END MODULE rhsZGrid_SW_m