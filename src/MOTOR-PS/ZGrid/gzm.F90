!!--------------------------------------------------------------------------------------------------
! PROJECT           : gzm - generalized Z-grid model
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie, Jilong Chen
! VERSION           : V 0.0
! HISTORY           : 2023-03-01 transferred by Yuanfu Xie from his own Z-grid model
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2023/03/01, @GBA-MWF, Shenzhen
!   Modified by Jilong Chen followed by Yuanfu Xie, 2024/10/21, for its 3D extension with vector
!   calculation and other efficiency considerations.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides the set of operators needed for a Z-grid model

MODULE gzm_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE
  TYPE :: gzm_t
    REAL(r_kind), ALLOCATABLE :: areaRatio(:)   !< changes with sigma heights
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE, PASS :: Divergen
    PROCEDURE, PASS :: Jacobian
    PROCEDURE, PASS :: Laplacia
    PROCEDURE, PASS :: Gradient  ! New subroutine for calculating gradients
    FINAL :: destructor
  END TYPE gzm_t

  INTERFACE gzm_t
    PROCEDURE :: constructor
  END INTERFACE gzm_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    IMPLICIT NONE
    TYPE(gzm_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg

    this%sg => sg

    ! Area ratio changing with sigma heights:
    ALLOCATE(this%areaRatio(this%sg%vLevel))
    this%areaRatio = (EarthRadius / (EarthRadius + this%sg%sigma))**2
  END FUNCTION constructor

  SUBROUTINE Divergen(this, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(gzm_t) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))
    ! Loop through all cells:
    VALUE = 0.0D0
    DO ic = 1, this%sg%num_icell  ! Apply to inner points only, icell has halo points etc.
      IF (this%sg%cell_type(ic) .GE. 1) CYCLE ! Apply to interior points only
      DO ie = 1, this%sg%num_edge
        edge_value = 0.0D0
        DO ig = 1, this%sg%numQuadPerEdge
          left_value = 0.0D0
          rightvalue = 0.0D0
          DO i = 1, UBOUND(this%sg%edge_stcl, 1)
            left_value = left_value + &
              opr_left(:, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
            rightvalue = rightvalue + &
              opr_right(:, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
          END DO
          edge_value = edge_value + this%sg%coef_gl(ig) * left_value * rightvalue
        END DO
        VALUE(:, ic) = VALUE(:, ic) + edge_value * this%sg%edge_leng(ie, ic)
      END DO
      VALUE(:, ic) = VALUE(:, ic) / this%sg%cell_area(ic) * this%areaRatio
    END DO

    ! Deallocate memory:
    DEALLOCATE (edge_value, left_value, rightvalue)

  END SUBROUTINE Divergen

  SUBROUTINE Jacobian(this, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(gzm_t) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    ! Loop through all cells:
    VALUE = 0.0D0
    DO ic = 1, this%sg%num_icell  ! Apply to inner points only, icell has halo points etc.
      IF (this%sg%cell_type(ic) .GE. 1) CYCLE ! Apply to interior points only
      DO ie = 1, this%sg%num_edge
        edge_value = 0.0D0
        DO ig = 1, this%sg%numQuadPerEdge
          left_value = 0.0D0
          rightvalue = 0.0D0
          DO i = 1, UBOUND(this%sg%edge_stcl, 1)
            left_value = left_value + &
              opr_left(:, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
            rightvalue = rightvalue + &
              opr_right(:, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_tang(i, ig, ie, ic)
          END DO
          edge_value = edge_value + this%sg%coef_gl(ig) * left_value * rightvalue
        END DO
        VALUE(:, ic) = VALUE(:, ic) + edge_value * this%sg%edge_leng(ie, ic)
      END DO
      VALUE(:, ic) = VALUE(:, ic) / this%sg%cell_area(ic) * this%areaRatio
    END DO

    ! Deallocate memory:
    DEALLOCATE (edge_value, left_value, rightvalue)

  END SUBROUTINE Jacobian

  SUBROUTINE Laplacia(this, oprand, VALUE)

    IMPLICIT NONE

    CLASS(gzm_t) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%sg)

      ALLOCATE (edge_value(sg%vLevel), quad_value(sg%vLevel))
      ! Loop through all cells:
      VALUE = 0.0D0
      DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (sg%cell_type(ic) .EQ. 0) THEN ! Apply to interior
          DO ie = 1, sg%num_edge
            edge_value = 0.0D0
            DO ig = 1, sg%numQuadPerEdge
              quad_value = 0.0D0
              DO i = 1, UBOUND(sg%edge_stcl, 1)
                quad_value = quad_value + &
                  oprand(:, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
              END DO
              edge_value = edge_value + sg%coef_gl(ig) * quad_value
            END DO
            VALUE(:, ic) = VALUE(:, ic) + edge_value * sg%edge_leng(ie, ic)
          END DO
        END IF
        VALUE(:, ic) = VALUE(:, ic) / sg%cell_area(ic) * this%areaRatio
      END DO

    END ASSOCIATE

    ! Deallocate memory:
    DEALLOCATE (edge_value, quad_value)

  END SUBROUTINE Laplacia

  SUBROUTINE Gradient(this, oprand, grad_lat, grad_lon)
    !-----------------------------------------------------------------
    ! This subroutine calculates the gradient of a field (stream function or velocity potential)
    ! in both latitude (φ) and longitude (λ) directions using predefined SingleGrid coefficients.
    !
    ! INPUT:
    !   - oprand: The field for which the gradient is being calculated (e.g., stream function or velocity potential)
    ! OUTPUT:
    !   - grad_lat: The gradient in the latitude direction (∂/∂φ)
    !   - grad_lon: The gradient in the longitude direction (∂/∂λ)
    !-----------------------------------------------------------------
    CLASS(gzm_t) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)          ! Input scalar field (e.g., stream function or velocity potential)
    REAL(r_kind), INTENT(OUT) :: grad_lat(:, :), grad_lon(:, :)  ! Gradients in the latitude and longitude directions

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio    ! Ratio for Earth radius scaling
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    PRINT*,'Entered gradient... ',this%sg%vLevel

    ! Allocate memory for temporary arrays
    ALLOCATE (edge_value(this%sg%vLevel), quad_value(this%sg%vLevel))

  ! Loop over through all cells:
    grad_lat = 0.0D0
    grad_lon = 0.0D0
    PRINT*,'Gradient icell: ',this%sg%num_icell,this%sg%num_edge
    DO ic = 1, this%sg%num_icell  ! Skip boundary points
      ! Loop over edges of each cell
      DO ie = 1, this%sg%num_edge
        edge_value = 0.0D0
        ! DO ig = 1, this%sg%numQuadPerEdge
        !   quad_value = 0.0D0

        !   ! Loop through stencils of the current edge
        !     quad_value = quad_value + oprand(:, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
        !   END DO

        !   ! Accumulate edge contribution
        !   edge_value = edge_value + this%sg%coef_gl(ig) * quad_value
        ! END DO
        ! DO i = 1, UBOUND(this%sg%edge_stcl, 1)
        !   grad_lat(:,ic) = grad_lat(:,ic) + oprand(:,this%sg%edge_stcl(i, ie, ic))*this%sg%coef_norm(i,ie,ic)
        ! END DO

        ! Calculate gradients using precomputed normal vectors, edge lengths, and considering spherical geometry
        grad_lat(:, ic) = grad_lat(:, ic) + edge_value * this%sg%edgeNorm2(1, ie, ic) * this%sg%edge_leng(ie, ic) * this%areaRatio
        grad_lon(:, ic) = grad_lon(:, ic) + edge_value * this%sg%edgeNorm2(2, ie, ic) * this%sg%edge_leng(ie, ic) * this%areaRatio
      END DO

      ! Normalize gradients by cell area: WHY asked by Yuanfu Xie on 2024-10-21???
      grad_lat(:, ic) = grad_lat(:, ic) / this%sg%cell_area(ic)
      grad_lon(:, ic) = grad_lon(:, ic) / this%sg%cell_area(ic)
    END DO

    ! Deallocate temporary arrays
    DEALLOCATE (edge_value, quad_value)

  END SUBROUTINE Gradient

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(gzm_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

    IF (ALLOCATED(this%areaRatio)) DEALLOCATE(this%areaRatio)

  END SUBROUTINE destructor
END MODULE gzm_m
