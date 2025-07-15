MODULE gzm_adj_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: EarthRadius
  IMPLICIT NONE

  TYPE :: gzm_adj_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE :: constructor
    PROCEDURE :: Divergen_AD
    PROCEDURE :: Jacobian_AD
    PROCEDURE :: Laplacia_AD
    PROCEDURE :: Grad_Lat_AD
    PROCEDURE :: Grad_Lon_AD
    FINAL :: destructor
  END TYPE gzm_adj_t

CONTAINS

  ! Constructor for gzm_adj_t type
  SUBROUTINE constructor(this, sg_input)
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg_input  ! sg_input must be a TARGET

    ! Initialize the SingleGrid_t pointer in gzm_adj_t
    this%sg => sg_input

  END SUBROUTINE constructor

  SUBROUTINE Divergen_AD(this, opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    IMPLICIT NONE
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    DO k = this%sg%vLevel, 1, -1
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius

      DO ic = this%sg%num_icell, 1, -1
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE

        DO ie = this%sg%num_edge, 1, -1
          edge_value(k) = adj_value(k, ic) / (this%sg%cell_area(ic) * ratio**2.0) * this%sg%edge_leng(ie, ic)

          DO ig = this%sg%numQuadPerEdge, 1, -1
            left_value(k) = 0.0_R_KIND
            rightvalue(k) = 0.0_R_KIND

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              left_value(k) = left_value(k) + opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
            END DO

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              adj_opr_left(k, this%sg%edge_stcl(i, ie, ic)) = adj_opr_left(k, this%sg%edge_stcl(i, ie, ic)) + &
                                                              edge_value(k) * this%sg%coef_gl(ig) * rightvalue(k) * this%sg%coef_func(i, ig, ie, ic)
              adj_opr_right(k, this%sg%edge_stcl(i, ie, ic)) = adj_opr_right(k, this%sg%edge_stcl(i, ie, ic)) + &
                                                               edge_value(k) * this%sg%coef_gl(ig) * left_value(k) * this%sg%coef_norm(i, ig, ie, ic)
            END DO
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (edge_value, left_value, rightvalue)
  END SUBROUTINE Divergen_AD

  SUBROUTINE Jacobian_AD(this, opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    IMPLICIT NONE
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    DO k = this%sg%vLevel, 1, -1
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius

      DO ic = this%sg%num_icell, 1, -1
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE

        DO ie = this%sg%num_edge, 1, -1
          edge_value(k) = adj_value(k, ic) / (this%sg%cell_area(ic) * ratio**2.0) * this%sg%edge_leng(ie, ic)

          DO ig = this%sg%numQuadPerEdge, 1, -1
            left_value(k) = 0.0_R_KIND
            rightvalue(k) = 0.0_R_KIND

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              left_value(k) = left_value(k) + opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_tang(i, ig, ie, ic)
            END DO

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              adj_opr_left(k, this%sg%edge_stcl(i, ie, ic)) = adj_opr_left(k, this%sg%edge_stcl(i, ie, ic)) + &
                                                              edge_value(k) * this%sg%coef_gl(ig) * rightvalue(k) * this%sg%coef_func(i, ig, ie, ic)
              adj_opr_right(k, this%sg%edge_stcl(i, ie, ic)) = adj_opr_right(k, this%sg%edge_stcl(i, ie, ic)) + &
                                                               edge_value(k) * this%sg%coef_gl(ig) * left_value(k) * this%sg%coef_tang(i, ig, ie, ic)
            END DO
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (edge_value, left_value, rightvalue)
  END SUBROUTINE Jacobian_AD

  SUBROUTINE Laplacia_AD(this, oprand, adj_oprand, adj_value)
    IMPLICIT NONE
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_oprand(:, :)
    REAL(r_kind), INTENT(IN) :: adj_value(:, :)

    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ALLOCATE (edge_value(this%sg%vLevel), quad_value(this%sg%vLevel))

    DO k = this%sg%vLevel, 1, -1
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius

      DO ic = this%sg%num_icell, 1, -1
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE

        DO ie = this%sg%num_edge, 1, -1
          edge_value(k) = adj_value(k, ic) / (this%sg%cell_area(ic) * ratio**2.0) * this%sg%edge_leng(ie, ic)

          DO ig = this%sg%numQuadPerEdge, 1, -1
            quad_value(k) = 0.0_R_KIND

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              quad_value(k) = quad_value(k) + oprand(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
            END DO

            DO i = UBOUND(this%sg%edge_stcl, 1), 1, -1
              adj_oprand(k, this%sg%edge_stcl(i, ie, ic)) = adj_oprand(k, this%sg%edge_stcl(i, ie, ic)) + &
                                                            edge_value(k) * this%sg%coef_gl(ig) * this%sg%coef_norm(i, ig, ie, ic)
            END DO
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (edge_value, quad_value)
  END SUBROUTINE Laplacia_AD

  ! Adjoint of Latitude gradient operator (∇_lat)
  SUBROUTINE Grad_Lat_AD(this, oprand, adj_oprand, adj_grad_lat)
    IMPLICIT NONE
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_oprand(:, :)
    REAL(r_kind), INTENT(IN) :: adj_grad_lat(:, :)

    INTEGER(i_kind) :: i, ic, ie, k
    REAL(r_kind) :: ratio

    DO k = this%sg%vLevel, 1, -1
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius

      DO ic = this%sg%num_icell, 1, -1
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE

        DO ie = this%sg%num_edge, 1, -1
          adj_oprand(k, this%sg%edge_stcl(1, ie, ic)) = adj_oprand(k, this%sg%edge_stcl(1, ie, ic)) + &
                                                        adj_grad_lat(k, ic) * this%sg%edgeNorm2(1, ie, ic) * this%sg%edge_leng(ie, ic) / ratio**2.0
        END DO
      END DO
    END DO
  END SUBROUTINE Grad_Lat_AD

  ! Adjoint of Longitude gradient operator (∇_lon)
  SUBROUTINE Grad_Lon_AD(this, oprand, adj_oprand, adj_grad_lon)
    IMPLICIT NONE
    CLASS(gzm_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_oprand(:, :)
    REAL(r_kind), INTENT(IN) :: adj_grad_lon(:, :)

    INTEGER(i_kind) :: i, ic, ie, k
    REAL(r_kind) :: ratio, cos_phi

    DO k = this%sg%vLevel, 1, -1
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius

      DO ic = this%sg%num_icell, 1, -1
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE

        cos_phi = COS(this%sg%cell_cntr(1, ic))

        DO ie = this%sg%num_edge, 1, -1
          adj_oprand(k, this%sg%edge_stcl(2, ie, ic)) = adj_oprand(k, this%sg%edge_stcl(2, ie, ic)) + &
                                                        adj_grad_lon(k, ic) * this%sg%edgeNorm2(2, ie, ic) * this%sg%edge_leng(ie, ic) / (ratio**2.0 * cos_phi)
        END DO
      END DO
    END DO
  END SUBROUTINE Grad_Lon_AD

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(gzm_adj_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor

END MODULE gzm_adj_m
