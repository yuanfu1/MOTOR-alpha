!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-DA.possionSolver_mixPsiChiSolver_Dirichlet
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, on 2025-04-17
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module contains a solver of streamfunction and velocity potential from vorticity and divergence
!! with u and v as boundary conditions. It calls PoissonSolver_t to solve the Poisson equation with
!! Dirichlet boundary conditions.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
MODULE mixedPsiChiSolver_Dirichlet_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE geometry_m, ONLY: geometry_t
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: singleGrid_t
  USE State_m, ONLY: State_t
  USE poissonSolver_m, ONLY: poissonSolver_t

  IMPLICIT NONE

  TYPE :: mixedPsiChiSolver_Dirichlet_t
    TYPE(poissonSolver_t) :: ps
    TYPE(gzm_t) :: gzm
    TYPE(State_t) :: X
    INTEGER(c_int) :: mgStart, mgEnd

    CONTAINS
      PROCEDURE :: initialize_s
      PROCEDURE :: solve_s
  END TYPE mixedPsiChiSolver_Dirichlet_t

  CONTAINS

  SUBROUTINE initialize_s(this, mgStart, mgEnd, geo)
    CLASS(mixedPsiChiSolver_Dirichlet_t) :: this
    INTEGER(c_int), INTENT(IN) :: mgStart, mgEnd
    TYPE(geometry_t), INTENT(IN) :: geo

    ! Local variables:
    CHARACTER(LEN=1024) :: configFile

    ! Initialize the poissonSolver
    CALL getarg(1, configFile)
    this%ps = poissonSolver_t(configFile, geo)
    PRINT*,'Poisson solver is initiated: ',TRIM(configFile)

    ! Allocate memory and gzms:
    IF (mgStart .LT. geo%mg%mg_coarsest .OR. mgEnd .GT. geo%mg%mg_finest) THEN
      WRITE(*,1) mgStart,mgEnd,geo%mg%mg_coarsest,geo%mg%mg_finest
1       FORMAT('Multigrid levels: ',2I3,' requested are out of the geometry setting: ',2I3)
      STOP
    END IF

    this%mgStart = mgStart
    this%mgEnd = mgEnd
  END SUBROUTINE initialize_s

  SUBROUTINE solve_s(this, vor, div, u, v, psi, chi, sg)
    CLASS(mixedPsiChiSolver_Dirichlet_t) :: this
    REAL(c_double), INTENT(IN) :: vor(:,:), div(:,:)
    REAL(c_double), INTENT(IN) :: u(:,:), v(:,:)
    REAL(c_double), INTENT(OUT) :: psi(:,:), chi(:,:)
    TYPE(singleGrid_t), INTENT(IN) :: sg

    ! Local variables:
    INTEGER(c_int) :: i, j
    REAL(c_double), ALLOCATABLE :: rhs(:,:),uv(:,:,:)

    ALLOCATE(rhs(sg%vLevel,sg%num_cell),uv(sg%vLevel,sg%num_cell,2))
    rhs = vor
    ! Assign the BC value to zeros:
    DO i=1,sg%num_cell
      IF (sg%cell_type(i) .EQ. 1 .OR. sg%cell_type(i) .EQ. 2) THEN
        rhs(:,i) = 0.0D0
      END IF
    END DO

    ! Solve the mixed Psi problem
    CALL this%ps%PoissonSol(sg%gLevel, sg%vLevel, rhs, psi, 'forward')

    rhs = div
    ! Assign the BC value to zeros:
    DO i=1,sg%num_cell
      IF (sg%cell_type(i) .EQ. 1 .OR. sg%cell_type(i) .EQ. 2) THEN
        rhs(:,i) = 0.0D0
      END IF
    END DO

    ! Solve the mixed Psi and Chi problem
    CALL this%ps%PoissonSol(sg%gLevel, sg%vLevel, rhs, chi, 'forward')

    ! Determining the general solution:
    CALL addLinearSolutions(psi, chi, u, v, sg)

    ! Recalculate the boundary values after adding the general solution:
    CALL uvVelocityOnInterior(psi, chi, uv(:,:,1), uv(:,:,2), sg)
  END SUBROUTINE solve_s
END MODULE mixedPsiChiSolver_Dirichlet_m