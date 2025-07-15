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
MODULE PoissonDirichlet_m
  USE , INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE linearSolversBase_m, ONLY: linearSolversBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t

  IMPLICIT NONE

  TYPE, EXTENDS(linearSolversBase_t) :: PoissonDirichlet_t
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s => initializePD_s
    PROCEDURE, PUBLIC :: loadMatrix_s => loadMatrixPD_s
  END TYPE PoissonDirichlet_t
  CONTAINS
    SUBROUTINE initializePD_s(this, X)
      CLASS(PoissonDirichlet_t) :: this
      TYPE(State_t), TARGET, INTENT(IN) :: X

      ! Allocate the sln and rhs based on the size of X
      this%X => X
      this%nonzero = X%sg%num_cell*X%sg%num_stcl
      ALLOCATE(this%sln(X%sg%vLevel,X%sg%num_cell),this%rhs(X%sg%vLevel,X%sg%num_cell))
      ALLOCATE(this%matrix(X%sg%vLevel, this%nonzero))  ! num_cell*num_stcl
      ALLOCATE(this%irow(this%nonzero), this%jcol(this%nonzero))
      ! Initialize the matrix and indices
      this%sln = 0.0_c_double
      this%rhs = 0.0_c_double
      this%matrix = 0.0_c_double
      this%irow = 0
      this%jcol = 0
    END SUBROUTINE initializePD_s

    SUBROUTINE loadMatrixPD_s(this)
      CLASS(PoissonDirichlet_t) :: this

      ! Local variables:
      INTEGER(c_int) :: i, j, ig, is, nz

      ASSOCIATE (sg => this%X%sg)
        nz = 0
        DO i=1,sg%num_cell
          ! Load the matrix for the Poisson Dirichlet solver: interior points
          IF (sg%cell_type(i) .EQ. 0) THEN
            DO j=1,sg%num_edge
              ! Fill the matrix with appropriate values
              DO ig = 1, sg%numQuadPerEdge
                DO is=1, sg%num_stcl
                  nz = nz + 1
                  this%irow(nz) = i
                  this%jcol(nz) = sg%edge_stcl(i, j, i)
                  this%matrix(i,nz) = sg%coef_norm(is, ig, j, i)
                END DO
              END DO
            END DO
          ELSE IF (sg%cell_type(i) .EQ. 1 .OR. sg%cell_type(i) .EQ. 2) THEN
            ! Dirichlet boundary condition
            nz = nz + 1
            this%irow(nz) = i
            this%jcol(nz) = i
            this%matrix(i,nz) = 1.0_c_double
          END IF
        END DO

        ! Load the matrix for the Poisson Dirichlet solver: boundary points

      END ASSOCIATE
    END SUBROUTINE loadMatrixPD_s

  END MODULE PoissonDirichlet_m