MODULE gzm_tlm_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE

  TYPE, EXTENDS(gzm_t) :: gzm_tlm_t
  CONTAINS
    PROCEDURE :: Divergen_tlm
    PROCEDURE :: Jacobian_tlm
    PROCEDURE :: Laplacia_tlm
  END TYPE gzm_tlm_t

CONTAINS

  SUBROUTINE Divergen_tlm(this, opr_left, opr_right, opr_left_tlm, opr_right_tlm, value_tlm)
    IMPLICIT NONE
    CLASS(gzm_tlm_t) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(IN) :: opr_left_tlm(:, :), opr_right_tlm(:, :)
    REAL(r_kind), INTENT(OUT) :: value_tlm(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value_tlm(:), left_value_tlm(:), rightvalue_tlm(:)
    REAL(r_kind), ALLOCATABLE :: left_value(:), rightvalue(:)

    ! Allocate memory for local variables
    ALLOCATE (edge_value_tlm(this%sg%vLevel), left_value_tlm(this%sg%vLevel), rightvalue_tlm(this%sg%vLevel))
    ALLOCATE (left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    ! Loop over levels and cells
    DO k = 1, this%sg%vLevel
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius
      DO ic = 1, this%sg%num_icell
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE
        value_tlm(k, ic) = 0.0D0
        DO ie = 1, this%sg%num_edge
          edge_value_tlm(k) = 0.0D0
          DO ig = 1, this%sg%numQuadPerEdge
            left_value(k) = 0.0D0
            rightvalue(k) = 0.0D0
            left_value_tlm(k) = 0.0D0
            rightvalue_tlm(k) = 0.0D0

            ! Compute left_value, rightvalue and their tangent linear counterparts
            DO i = 1, UBOUND(this%sg%edge_stcl, 1)
              left_value(k) = left_value(k) + opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
              left_value_tlm(k) = left_value_tlm(k) + opr_left_tlm(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue_tlm(k) = rightvalue_tlm(k) + opr_right_tlm(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
            END DO

            ! Apply the chain rule
            edge_value_tlm(k) = edge_value_tlm(k) + this%sg%coef_gl(ig) * (left_value(k) * rightvalue_tlm(k) + left_value_tlm(k) * rightvalue(k))
          END DO
          value_tlm(k, ic) = value_tlm(k, ic) + edge_value_tlm(k) * this%sg%edge_leng(ie, ic)
        END DO
        value_tlm(k, ic) = value_tlm(k, ic) / this%sg%cell_area(ic) / ratio**2.0D0
      END DO
    END DO

    ! Deallocate memory
    DEALLOCATE (edge_value_tlm, left_value_tlm, rightvalue_tlm, left_value, rightvalue)
  END SUBROUTINE Divergen_tlm

  SUBROUTINE Jacobian_tlm(this, opr_left, opr_right, opr_left_tlm, opr_right_tlm, value_tlm)
    IMPLICIT NONE
    CLASS(gzm_tlm_t) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(IN) :: opr_left_tlm(:, :), opr_right_tlm(:, :)
    REAL(r_kind), INTENT(OUT) :: value_tlm(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value_tlm(:), left_value_tlm(:), rightvalue_tlm(:)
    REAL(r_kind), ALLOCATABLE :: left_value(:), rightvalue(:)

    ! Allocate memory for local variables
    ALLOCATE (edge_value_tlm(this%sg%vLevel), left_value_tlm(this%sg%vLevel), rightvalue_tlm(this%sg%vLevel))
    ALLOCATE (left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    ! Loop over levels and cells
    DO k = 1, this%sg%vLevel
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius
      DO ic = 1, this%sg%num_icell
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE
        value_tlm(k, ic) = 0.0D0
        DO ie = 1, this%sg%num_edge
          edge_value_tlm(k) = 0.0D0
          DO ig = 1, this%sg%numQuadPerEdge
            left_value(k) = 0.0D0
            rightvalue(k) = 0.0D0
            left_value_tlm(k) = 0.0D0
            rightvalue_tlm(k) = 0.0D0

            ! Compute left_value, rightvalue and their tangent linear counterparts
            DO i = 1, UBOUND(this%sg%edge_stcl, 1)
              left_value(k) = left_value(k) + opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_tang(i, ig, ie, ic)
              left_value_tlm(k) = left_value_tlm(k) + opr_left_tlm(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue_tlm(k) = rightvalue_tlm(k) + opr_right_tlm(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_tang(i, ig, ie, ic)
            END DO

            ! Apply the chain rule
            edge_value_tlm(k) = edge_value_tlm(k) + this%sg%coef_gl(ig) * (left_value(k) * rightvalue_tlm(k) + left_value_tlm(k) * rightvalue(k))
          END DO
          value_tlm(k, ic) = value_tlm(k, ic) + edge_value_tlm(k) * this%sg%edge_leng(ie, ic)
        END DO
        value_tlm(k, ic) = value_tlm(k, ic) / this%sg%cell_area(ic) / ratio**2.0D0
      END DO
    END DO

    ! Deallocate memory
    DEALLOCATE (edge_value_tlm, left_value_tlm, rightvalue_tlm, left_value, rightvalue)
  END SUBROUTINE Jacobian_tlm

  SUBROUTINE Laplacia_tlm(this, oprand, oprand_tlm, value_tlm)
    IMPLICIT NONE
    CLASS(gzm_tlm_t) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)
    REAL(r_kind), INTENT(IN) :: oprand_tlm(:, :)
    REAL(r_kind), INTENT(OUT) :: value_tlm(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value_tlm(:), quad_value_tlm(:)
    REAL(r_kind), ALLOCATABLE :: quad_value(:)

    ! Associate the grid structure and allocate memory
    ASSOCIATE (sg => this%sg)
      ALLOCATE (edge_value_tlm(sg%vLevel), quad_value_tlm(sg%vLevel))
      ALLOCATE (quad_value(sg%vLevel))

      ! Loop over levels and cells
      DO k = 1, sg%vLevel
        ratio = (EarthRadius + sg%sigma(k)) / EarthRadius
        DO ic = 1, sg%num_icell
          IF (sg%bdy_type(ic) .GE. 1) CYCLE
          value_tlm(k, ic) = 0.0D0
          DO ie = 1, sg%num_edge
            edge_value_tlm(k) = 0.0D0
            DO ig = 1, sg%numQuadPerEdge
              quad_value(k) = 0.0D0
              quad_value_tlm(k) = 0.0D0

              ! Compute quad_value and its tangent linear counterpart
              DO i = 1, UBOUND(sg%edge_stcl, 1)
                quad_value(k) = quad_value(k) + oprand(k, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
                quad_value_tlm(k) = quad_value_tlm(k) + oprand_tlm(k, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
              END DO

              ! Apply the chain rule
              edge_value_tlm(k) = edge_value_tlm(k) + sg%coef_gl(ig) * (quad_value(k) * quad_value_tlm(k))
            END DO
            value_tlm(k, ic) = value_tlm(k, ic) + edge_value_tlm(k) * sg%edge_leng(ie, ic)
          END DO
          value_tlm(k, ic) = value_tlm(k, ic) / sg%cell_area(ic) / ratio**2.0D0
        END DO
      END DO
    END ASSOCIATE

    ! Deallocate memory
    DEALLOCATE (edge_value_tlm, quad_value_tlm, quad_value)
  END SUBROUTINE Laplacia_tlm

END MODULE gzm_tlm_m
