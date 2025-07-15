MODULE CalVerDer_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE
  TYPE :: CalVerDer_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE :: FirstOrder
    PROCEDURE :: SecondOrder
    PROCEDURE :: FirstOrder_half
    FINAL :: destructor
  END TYPE

  INTERFACE CalVerDer_t
    PROCEDURE :: constructor
  END INTERFACE CalVerDer_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(CalVerDer_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    this%sg => sg
  END FUNCTION constructor

  SUBROUTINE FirstOrder(this, A, parA_sigma)
    IMPLICIT NONE
    CLASS(CalVerDer_t) :: this
    REAL(r_kind), INTENT(IN) :: A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :)
    INTEGER(i_kind) :: i, k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               coef => this%sg%coef_fstdif)
      parA_sigma = 0.0D0
      DO i = 1, num_cell
        ! IF (cell_type(i) .EQ. 0) THEN
        DO k = 1, vLevel
          IF (k .EQ. 1) THEN
            parA_sigma(k, i) = SUM(coef(k, i, :) * A(k:k + 2, i))
          ELSEIF (k .EQ. vLevel) THEN
            parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 2:k, i))
          ELSE
            parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 1, i))
          END IF
        END DO
        ! END IF
      END DO
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, parA_sigma)
    END ASSOCIATE
  END SUBROUTINE FirstOrder

  SUBROUTINE SecondOrder(this, A, parA_sigma_2nd)
    IMPLICIT NONE

    CLASS(CalVerDer_t) :: this
    REAL(r_kind), INTENT(IN) :: A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma_2nd(:, :)
    INTEGER(i_kind) :: i, k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               coef => this%sg%coef_scddif)

      parA_sigma_2nd = 0.0D0
      DO i = 1, num_cell
        ! IF (cell_type(i) .EQ. 0) THEN
        DO k = 1, vLevel
          IF (k .EQ. 1) THEN
            parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k:k + 4, i))
          ELSEIF (k .EQ. 2) THEN
            parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 3, i))
          ELSEIF (k .EQ. vLevel) THEN
            parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 4:k, i))
          ELSEIF (k .EQ. vLevel - 1) THEN
            parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 3:k + 1, i))
          ELSE
            parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 2:k + 2, i))
          END IF
        END DO
        ! END IF
      END DO
    END ASSOCIATE
    ! CALL sg%ExchangeMatOnHalo2D(vLevel, parA_sigma_2nd)
  END SUBROUTINE

  SUBROUTINE FirstOrder_half(this, A, parA_sigma_half)
    IMPLICIT NONE
    CLASS(CalVerDer_t) :: this
    REAL(r_kind), INTENT(IN) :: A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma_half(:, :)
    INTEGER(i_kind) :: i, k, sz

    ASSOCIATE (num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               coef => this%sg%coef_fstdif_half)
      parA_sigma_half = 0.0D0
      sz = SIZE(A, 1)
      DO i = 1, num_cell
        ! IF (cell_type(i) .EQ. 0) THEN
        IF (vLevel .GT. sz) THEN
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k:k + 3, i))
            ELSEIF (k .EQ. 2) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 2, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 4:k - 1, i))
            ELSEIF (k .EQ. vLevel - 1) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 3:k, i))
            ELSE
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 2:k + 1, i))
            END IF
          END DO
        ELSE
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k:k + 3, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 2:k + 1, i))
            ELSE
              parA_sigma_half(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 2, i))
            END IF
          END DO

        END IF
        ! END IF
      END DO
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, parA_sigma)
    END ASSOCIATE
  END SUBROUTINE FirstOrder_half

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(CalVerDer_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor
END MODULE CalVerDer_m
