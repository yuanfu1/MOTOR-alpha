!!--------------------------------------------------------------------------------------------------
! PROJECT           : gzm - generalized Z-grid model
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2023-03-01 transferred by Yuanfu Xie from his own Z-grid model
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2023/03/01, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides the set of operators needed for a Z-grid model

MODULE gzm_Xie_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Multigrid_m, ONLY: MultiGrid_t
  USE state_m, ONLY: state_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  ! USE vertexSpace_m, ONLY : vertexSpace_t

  IMPLICIT NONE

  TYPE(geometry_t), TARGET :: geometry

  TYPE :: gzm_Xie_t
    REAL(r_kind) :: coefJacobian(4, 4)
    REAL(r_kind), ALLOCATABLE :: fv(:, :, :) ! Function values: vlvl x 7 x num_vrtx
    !TYPE(SingleGrid_t), POINTER :: grid
    !TYPE(Multigrid_t), POINTER :: mg
    TYPE(state_t) :: states
    ! TYPE(mpddGlob_t) :: mpddGlob
    TYPE(geometry_t) :: geometry
  CONTAINS
    PROCEDURE :: initial
    PROCEDURE :: destroy
    PROCEDURE :: Divergen
    PROCEDURE :: DivergenAdjoint  ! Add this line
    PROCEDURE :: Jacobian
    PROCEDURE :: JacobianAdjoint1  ! Adjoint with respect to first variable
    PROCEDURE :: JacobianAdjoint2  ! Adjoint with respect to second variable
    PROCEDURE :: Laplacia
    PROCEDURE :: LaplaciaAdjoint
    PROCEDURE :: LaplaciaAdjoint_copilot
  END TYPE gzm_Xie_t

CONTAINS

  SUBROUTINE initial(this, configFile, mpddGlob)
    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(mpddGlob_t) :: mpddGlob

    ! Local variables:
    INTEGER(i_kind) :: i, j

    ! Pre-calculate the coefficient for a line integral
    ! of Jacobian operator:
!       CALL JacobianCoef(this%coefJacobian)
!       DO i=1,4
!         WRITE(*,4) i,(this%coefJacobian(i,j),j=1,4)
! 4       FORMAT('CoeffJacobian row: ',I2,4D12.4)
!       END DO
!       DO i=1,4
!         WRITE(*,2) i,SUM(this%coefJacobian(i,:)),SUM(this%coefJacobian(:,i))
! 2       FORMAT('Row sum: ',I1,' is: ',D12.4,' column sum: ',D12.4)
!       END DO

    ! mpdd initialization:
    ! CALL this%mpddGlob%initialize()

    ! Initialize a grid geometry:
    CALL this%geometry%initialize(configFile, mpddGlob)
    PRINT *, 'Finish initialization of mpddGlob: ', mpddGlob%myrank

    ! ! The model grid is at the finest multigrid level:
    !this%grid => this%geometry%mg%sg(this%geometry%mg%mg_finest)
    !this%mg => geometry%mg
    PRINT *, 'Debugging here ' !,this%grid%num_vrtx
  END SUBROUTINE initial

  SUBROUTINE destroy(this)

    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this

    PRINT *, 'gzm destructor works'
    IF (ALLOCATED(this%fv)) DEALLOCATE (this%fv)

    CALL this%geometry%destroy
    ! CALL this%mpddGlob%barrier
    ! CALL this%mpddGlob%finalize
  END SUBROUTINE destroy

  SUBROUTINE Divergen(this, i_glvl, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    ! These operands are in cell space:
    REAL(r_kind), INTENT(IN) :: &
      opr_left(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))

      ALLOCATE (edge_value(sg%vLevel), left_value(sg%vLevel), rightvalue(sg%vLevel))

      ! Loop through all cells:
      DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (sg%cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only
        VALUE(:, ic) = 0.0D0
        DO ie = 1, sg%num_edge
          edge_value = 0.0D0
          left_value = 0.0D0
          rightvalue = 0.0D0
          DO ig = 1, sg%numQuadPerEdge
            DO i = 1, UBOUND(sg%edge_stcl, 1)
              left_value = left_value + &
                           opr_left(:, sg%edge_stcl(i, ie, ic)) * sg%coef_func(i, ig, ie, ic)
              rightvalue = rightvalue + &
                           opr_right(:, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
            END DO
            edge_value = edge_value + sg%coef_gl(ig) * left_value * rightvalue
          END DO
          VALUE(:, ic) = VALUE(:, ic) + edge_value * sg%edge_leng(ie, ic)
        END DO
        VALUE(:, ic) = VALUE(:, ic) / sg%cell_area(ic)
      END DO

      ! Deallocate memory:
      DEALLOCATE (edge_value, left_value, rightvalue)

    END ASSOCIATE
  END SUBROUTINE Divergen

  SUBROUTINE Jacobian(this, i_glvl, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    ! These operands are in cell space:
    REAL(r_kind), INTENT(IN) :: &
      opr_left(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))

      ALLOCATE (edge_value(sg%vLevel), left_value(sg%vLevel), rightvalue(sg%vLevel))

      ! Loop through all cells:
      DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (sg%cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only
        VALUE(:, ic) = 0.0D0
        DO ie = 1, sg%num_edge
          edge_value = 0.0D0
          left_value = 0.0D0
          rightvalue = 0.0D0
          DO ig = 1, sg%numQuadPerEdge
            DO i = 1, UBOUND(sg%edge_stcl, 1)
              left_value = left_value + &
                           opr_left(:, sg%edge_stcl(i, ie, ic)) * sg%coef_func(i, ig, ie, ic)
              rightvalue = rightvalue + &
                           opr_right(:, sg%edge_stcl(i, ie, ic)) * sg%coef_tang(i, ig, ie, ic)
            END DO
            edge_value = edge_value + sg%coef_gl(ig) * left_value * rightvalue
          END DO
          VALUE(:, ic) = VALUE(:, ic) + edge_value * sg%edge_leng(ie, ic)
        END DO
        VALUE(:, ic) = VALUE(:, ic) / sg%cell_area(ic)
      END DO

      ! Deallocate memory:
      DEALLOCATE (edge_value, left_value, rightvalue)

    END ASSOCIATE
  END SUBROUTINE Jacobian

  SUBROUTINE JacobianAdjoint1(this, i_glvl, opr_left, opr_right, VALUE)
    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    REAL(r_kind), INTENT(IN) :: &
      opr_left(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, icell
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:), left_value(:), right_value(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))
      ALLOCATE(edge_value(sg%vLevel), quad_value(sg%vLevel), &
               left_value(sg%vLevel), right_value(sg%vLevel))
      
      VALUE = 0.0D0

      ! Loop through cells for adjoint contributions
      DO ic = 1, sg%num_icell
        IF (sg%cell_type(ic) .NE. 0) CYCLE

        ! Scaled input for adjoint
        quad_value = opr_right(:,ic) / sg%cell_area(ic)

        DO ie = 1, sg%num_edge
          DO ig = 1, sg%numQuadPerEdge
            ! Get function value at quadrature point
            left_value = 0.0D0
            DO i = 1, UBOUND(sg%edge_stcl,1)
              icell = sg%edge_stcl(i,ie,ic)
              IF (icell > 0) THEN
                left_value = left_value + opr_left(:,icell) * sg%coef_func(i,ig,ie,ic)
              END IF
            END DO

            ! Accumulate adjoint contributions
            edge_value = -left_value * quad_value * sg%coef_gl(ig) * sg%edge_leng(ie,ic)
            
            ! Distribute using tangential coefficients
            DO i = 1, UBOUND(sg%edge_stcl,1)
              icell = sg%edge_stcl(i,ie,ic)
              IF (icell > 0) THEN
                VALUE(:,icell) = VALUE(:,icell) + edge_value * sg%coef_tang(i,ig,ie,ic)
              END IF
            END DO
          END DO
        END DO
      END DO

      DEALLOCATE(edge_value, quad_value, left_value, right_value)
    END ASSOCIATE
  END SUBROUTINE JacobianAdjoint1

  SUBROUTINE JacobianAdjoint2(this, i_glvl, opr_left, opr_right, VALUE)
    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    REAL(r_kind), INTENT(IN) :: &
      opr_left(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, icell
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:), left_value(:), right_value(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))
      ALLOCATE(edge_value(sg%vLevel), quad_value(sg%vLevel), &
               left_value(sg%vLevel), right_value(sg%vLevel))
      
      VALUE = 0.0D0

      ! Loop through cells for adjoint contributions
      DO ic = 1, sg%num_icell
        IF (sg%cell_type(ic) .NE. 0) CYCLE

        ! Scaled input for adjoint
        quad_value = opr_left(:,ic) / sg%cell_area(ic)

        DO ie = 1, sg%num_edge
          DO ig = 1, sg%numQuadPerEdge
            ! Get tangential derivative at quadrature point
            right_value = 0.0D0
            DO i = 1, UBOUND(sg%edge_stcl,1)
              icell = sg%edge_stcl(i,ie,ic)
              IF (icell > 0) THEN
                right_value = right_value + opr_right(:,icell) * sg%coef_tang(i,ig,ie,ic)
              END IF
            END DO

            ! Accumulate adjoint contributions
            edge_value = right_value * quad_value * sg%coef_gl(ig) * sg%edge_leng(ie,ic)
            
            ! Distribute using function coefficients
            DO i = 1, UBOUND(sg%edge_stcl,1)
              icell = sg%edge_stcl(i,ie,ic)
              IF (icell > 0) THEN
                VALUE(:,icell) = VALUE(:,icell) + edge_value * sg%coef_func(i,ig,ie,ic)
              END IF
            END DO
          END DO
        END DO
      END DO

      DEALLOCATE(edge_value, quad_value, left_value, right_value)
    END ASSOCIATE
  END SUBROUTINE JacobianAdjoint2

  SUBROUTINE Laplacia(this, i_glvl, oprand, VALUE)

    IMPLICIT NONE

    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    ! These operands are in cell space:
    REAL(r_kind), INTENT(IN) :: &
      oprand(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))

      ALLOCATE (edge_value(sg%vLevel), quad_value(sg%vLevel))

      ! Loop through all cells:
      DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (sg%cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only
        VALUE(:, ic) = 0.0D0
        DO ie = 1, sg%num_edge
          edge_value = 0.0D0
          DO ig = 1, sg%numQuadPerEdge
            quad_value = 0.0D0
            DO i = 1, UBOUND(sg%edge_stcl, 1)
              quad_value = quad_value + &
                oprand(:, sg%edge_stcl(i, ie, ic)) * &
                sg%coef_norm(i, ig, ie, ic)
            END DO
            edge_value = edge_value + sg%coef_gl(ig) * quad_value
          END DO
          VALUE(:, ic) = VALUE(:, ic) + edge_value * sg%edge_leng(ie, ic)
        END DO
        VALUE(:, ic) = VALUE(:, ic) / sg%cell_area(ic)
      END DO

      ! Deallocate memory:
      DEALLOCATE (edge_value, quad_value)
    END ASSOCIATE
  END SUBROUTINE Laplacia

  SUBROUTINE LaplaciaAdjoint(this, i_glvl, oprand, VALUE)
    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    REAL(r_kind), INTENT(IN) :: &
      oprand(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, icell
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))
      ALLOCATE(edge_value(sg%vLevel), quad_value(sg%vLevel))
      
      ! Initialize VALUE to zero
      VALUE = 0.0D0

      ! Loop through cells that give contributions
      DO ic = 1, sg%num_icell
        IF (sg%cell_type(ic) .NE. 0) CYCLE

        ! Apply inverse area scaling to input
        quad_value = oprand(:,ic) / sg%cell_area(ic)

        ! Loop through edges
        DO ie = 1, sg%num_edge
          ! Scale by edge length
          edge_value = quad_value * sg%edge_leng(ie,ic)
          
          ! Add contributions to stencil points
          DO i = 1, UBOUND(sg%edge_stcl,1)
            icell = sg%edge_stcl(i,ie,ic)
            IF (icell > 0) THEN
              DO ig = 1, sg%numQuadPerEdge
                VALUE(:,icell) = VALUE(:,icell) + &
                  edge_value * sg%coef_gl(ig) * sg%coef_norm(i,ig,ie,ic)
              END DO
            END IF
          END DO
        END DO
      END DO

      DEALLOCATE(edge_value, quad_value)
    END ASSOCIATE
  END SUBROUTINE LaplaciaAdjoint

  SUBROUTINE LaplaciaAdjoint_copilot(this, i_glvl, oprand, VALUE)
    IMPLICIT NONE
    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    REAL(r_kind), INTENT(IN) :: &
      oprand(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, icell
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))
      ALLOCATE(edge_value(sg%vLevel), quad_value(sg%vLevel))
      
      ! Initialize VALUE to zero
      VALUE = 0.0D0

      ! Loop through all cells that contribute to the adjoint
      DO ic = 1, sg%num_icell  ! Changed back to num_icell
        IF (sg%cell_type(ic) .NE. 0) CYCLE ! Interior points only
        
        DO i = 1, UBOUND(sg%edge_stcl,1)
          DO ie = 1, sg%num_edge
            edge_value = 0.0D0
            
            DO ig = 1, sg%numQuadPerEdge
              ! The adjoint of the normal derivative
              quad_value = oprand(:,ic) * sg%coef_norm(i,ig,ie,ic) * &
                        sg%edge_leng(ie,ic) / sg%cell_area(ic)
              
              edge_value = edge_value + sg%coef_gl(ig) * quad_value
            END DO
            
            ! Accumulate contributions to stencil points
            icell = sg%edge_stcl(i,ie,ic)
            IF (icell > 0) THEN
              VALUE(:,icell) = VALUE(:,icell) + edge_value
            END IF
          END DO
        END DO
      END DO

      DEALLOCATE(edge_value, quad_value)
    END ASSOCIATE
  END SUBROUTINE LaplaciaAdjoint_copilot

  SUBROUTINE DivergenAdjoint(this, i_glvl, opr_left, opr_right, VALUE, opr_left_adj, opr_right_adj, VALUE_adj)
    ! filepath: /Users/xiey/developments/da/motor_training/test/MOTOR/src/MOTOR-PS/ZGrid/gzm_Xie.F90
    IMPLICIT NONE

    CLASS(gzm_Xie_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i_glvl
    REAL(r_kind), INTENT(IN) :: &
      opr_left(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      VALUE(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(IN) :: &
      VALUE_adj(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)
    REAL(r_kind), INTENT(OUT) :: &
      opr_left_adj(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell), &
      opr_right_adj(this%geometry%mg%sg(i_glvl)%vLevel, this%geometry%mg%sg(i_glvl)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ASSOCIATE (sg => this%geometry%mg%sg(i_glvl))

      ALLOCATE (edge_value(sg%vLevel), left_value(sg%vLevel), rightvalue(sg%vLevel))

      ! Initialize adjoint variables
      opr_left_adj = 0.0D0
      opr_right_adj = 0.0D0

      ! Loop through all cells:
      DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (sg%cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only

        ! Adjoint of the area division
        edge_value = VALUE_adj(:, ic) / sg%cell_area(ic)

        DO ie = 1, sg%num_edge
          ! Adjoint of edge length multiplication
          rightvalue = edge_value * sg%edge_leng(ie, ic)

          DO ig = 1, sg%numQuadPerEdge
            ! Adjoint of the Gauss-Legendre coefficient multiplication
            left_value = rightvalue * sg%coef_gl(ig)

            DO i = 1, UBOUND(sg%edge_stcl, 1)
              ! Distribute to the adjoint of the input arrays
              opr_left_adj(:, sg%edge_stcl(i, ie, ic)) = opr_left_adj(:, sg%edge_stcl(i, ie, ic)) + &
                                                          left_value * opr_right(:, sg%edge_stcl(i, ie, ic)) * sg%coef_func(i, ig, ie, ic)
              opr_right_adj(:, sg%edge_stcl(i, ie, ic)) = opr_right_adj(:, sg%edge_stcl(i, ie, ic)) + &
                                                           left_value * opr_left(:, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
            END DO
          END DO
        END DO
      END DO

      ! Deallocate memory:
      DEALLOCATE (edge_value, left_value, rightvalue)

    END ASSOCIATE

  END SUBROUTINE DivergenAdjoint
END MODULE gzm_Xie_m
