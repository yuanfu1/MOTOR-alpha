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
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides the set of operators needed for a Z-grid model

MODULE gzm_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE
  TYPE :: gzm_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE :: Divergen
    PROCEDURE :: Jacobian
    PROCEDURE :: Laplacia
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

    ! PRINT *, 'this%sg%bdy_type is ', this%sg%vLevel, this%sg%bdy_type(1)
  END FUNCTION constructor

  SUBROUTINE Divergen(this, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(gzm_t) :: this
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, ic, ie, ig, k
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))
    DO k = 1, this%sg%vLevel
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius
      ! Loop through all cells:
      DO ic = 1, this%sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE ! Apply to interior points only
        VALUE(k, ic) = 0.0D0
        DO ie = 1, this%sg%num_edge
          edge_value(k) = 0.0D0
          DO ig = 1, this%sg%numQuadPerEdge
            left_value(k) = 0.0D0
            rightvalue(k) = 0.0D0
            DO i = 1, UBOUND(this%sg%edge_stcl, 1)
              left_value(k) = left_value(k) + &
                              opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + &
                              opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_norm(i, ig, ie, ic)
            END DO
            edge_value(k) = edge_value(k) + this%sg%coef_gl(ig) * left_value(k) * rightvalue(k)
          END DO
          VALUE(k, ic) = VALUE(k, ic) + edge_value(k) * this%sg%edge_leng(ie, ic)
        END DO
        VALUE(k, ic) = VALUE(k, ic) / this%sg%cell_area(ic) / ratio**2.0D0
      END DO
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
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ALLOCATE (edge_value(this%sg%vLevel), left_value(this%sg%vLevel), rightvalue(this%sg%vLevel))

    ! Loop through all cells:
    DO k = 1, this%sg%vLevel
      ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius
      DO ic = 1, this%sg%num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (this%sg%bdy_type(ic) .GE. 1) CYCLE ! Apply to interior points only
        VALUE(k, ic) = 0.0D0
        DO ie = 1, this%sg%num_edge
          edge_value(k) = 0.0D0
          DO ig = 1, this%sg%numQuadPerEdge
            left_value(k) = 0.0D0
            rightvalue(k) = 0.0D0
            DO i = 1, UBOUND(this%sg%edge_stcl, 1)
              left_value(k) = left_value(k) + &
                              opr_left(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_func(i, ig, ie, ic)
              rightvalue(k) = rightvalue(k) + &
                              opr_right(k, this%sg%edge_stcl(i, ie, ic)) * this%sg%coef_tang(i, ig, ie, ic)
            END DO
            IF (ic .EQ. 501) THEN
              PRINT *, 'ig and ie are: ', ig, ie
              PRINT *, 'coef_tang is', this%sg%coef_tang(:, ig, ie, ic)
              PRINT *, 'cell stacl lat are: ', this%sg%cell_cntr(1, this%sg%edge_stcl(:, ie, ic))
              PRINT *, 'cell stacl lon are: ', this%sg%cell_cntr(2, this%sg%edge_stcl(:, ie, ic))
              PRINT *, 'EarthRadius: ', EarthRadius
            END IF
            edge_value(k) = edge_value(k) + this%sg%coef_gl(ig) * left_value(k) * rightvalue(k)
          END DO
          VALUE(k, ic) = VALUE(k, ic) + edge_value(k) * this%sg%edge_leng(ie, ic)
          IF (ic .EQ. 501) THEN
            PRINT *, 'ie is: ', ie
            PRINT *, 'this%sg%coef_gl(ig) is', this%sg%coef_gl(:)
            PRINT *, 'sg%edge_leng(ie, ic) is', this%sg%edge_leng(ie, ic)
          END IF
        END DO
        VALUE(k, ic) = VALUE(k, ic) / this%sg%cell_area(ic) / ratio**2.0D0
        IF (ic .EQ. 501) THEN
          PRINT *, 'ic is: ', ic
          PRINT *, 'this%sg%cell_area is', this%sg%cell_area(ic)
        END IF
      END DO
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
    REAL(r_kind) :: ratio
    REAL(r_kind), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%sg)

      ALLOCATE (edge_value(sg%vLevel), quad_value(sg%vLevel))
      DO k = 1, this%sg%vLevel
        ratio = (EarthRadius + this%sg%sigma(k)) / EarthRadius
        ! Loop through all cells:
        DO ic = 1, sg%num_icell  ! Apply to inner points only, icell has halo points etc.
          IF (sg%bdy_type(ic) .GE. 1) CYCLE ! Apply to interior points only
          VALUE(k, ic) = 0.0D0
          DO ie = 1, sg%num_edge
            edge_value(k) = 0.0D0
            DO ig = 1, sg%numQuadPerEdge
              quad_value(k) = 0.0D0
              DO i = 1, UBOUND(sg%edge_stcl, 1)
                quad_value(k) = quad_value(k) + &
                                oprand(k, sg%edge_stcl(i, ie, ic)) * sg%coef_norm(i, ig, ie, ic)
              END DO
              edge_value(k) = edge_value(k) + sg%coef_gl(ig) * quad_value(k)
              IF (ic .EQ. 501) THEN
                PRINT *, 'ig and ie are: ', ig, ie
                PRINT *, 'coef_norm is', this%sg%coef_norm(:, ig, ie, ic)
                PRINT *, 'cell stacl lat are: ', this%sg%cell_cntr(1, this%sg%edge_stcl(:, ie, ic))
                PRINT *, 'cell stacl lon are: ', this%sg%cell_cntr(2, this%sg%edge_stcl(:, ie, ic))
              END IF
            END DO
            VALUE(k, ic) = VALUE(k, ic) + edge_value(k) * sg%edge_leng(ie, ic)
            IF (ic .EQ. 501) THEN
              PRINT *, 'ie is: ', ie
              PRINT *, 'this%sg%coef_gl(ig) is', this%sg%coef_gl(:)
              PRINT *, 'sg%edge_leng(ie, ic) is', this%sg%edge_leng(ie, ic)
            END IF
          END DO
          VALUE(k, ic) = VALUE(k, ic) / sg%cell_area(ic) / ratio**2.0D0

          IF (ic .EQ. 501) THEN
            PRINT *, 'ic is: ', ic
            PRINT *, 'this%sg%cell_area is', this%sg%cell_area(ic)
          END IF

        END DO
      END DO

    END ASSOCIATE

    ! Deallocate memory:
    DEALLOCATE (edge_value, quad_value)

  END SUBROUTINE Laplacia

  ! SUBROUTINE CalNormRatio(this, cellIdx, Height, ratio)
  !    IMPLICIT NONE

  !    CLASS(gzm_t) :: this
  !    INTEGER(i_kind), INTENT(IN) :: cellIdx
  !    REAL(r_kind), INTENT(IN) :: Height
  !    REAL(r_kind), INTENT(OUT) :: ratio(4)
  !    REAL(r_kind) :: lat(6), lon(6), clat, dlat, dlon
  !    INTEGER(i_kind) :: i

  !    DO i = 1, 4
  !       lat = this%sg%cell_cntr(1, this%sg%edge_stcl(:, i, cellIdx))
  !       lon = this%sg%cell_cntr(2, this%sg%edge_stcl(:, i, cellIdx))
  !       clat = lat(3) + lat(4)
  !       IF ((i .EQ. 1) .OR. (i .EQ. 3)) THEN
  !          ! dlat=dabs(lat(4)-lat(3))
  !          ratio(i)=EarthRadius/(EarthRadius+Height)
  !       ELSE
  !          clat = (lat(3) + lat(4))/2.0D0
  !          ! dlon=dabs(lon(4)-lon(3))
  !          ratio(i)=EarthRadius/(EarthRadius+Height)*

  !       END IF

  !    END DO

  ! END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(gzm_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor
END MODULE gzm_m
