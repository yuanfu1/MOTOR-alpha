MODULE CalVerDer_adj_m
  USE kinds_m, ONLY: r_kind, i_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  IMPLICIT NONE

  TYPE :: CalVerDer_adj_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE :: FirstOrder_AD
    PROCEDURE :: SecondOrder_AD
    PROCEDURE :: FirstOrder_half_AD
    FINAL :: destructor_adj
  END TYPE CalVerDer_adj_t

  INTERFACE CalVerDer_adj_t
    PROCEDURE :: constructor_adj
  END INTERFACE CalVerDer_adj_t

CONTAINS

  FUNCTION constructor_adj(sg) RESULT(this)
    TYPE(CalVerDer_adj_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    this%sg => sg
  END FUNCTION constructor_adj

  ! First-order adjoint computation
  SUBROUTINE FirstOrder_AD(this, A_adj, parA_sigma_adj, rhs_adj)
    CLASS(CalVerDer_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), POINTER :: A_adj(:, :), parA_sigma_adj(:, :), rhs_adj(:, :)
    INTEGER(i_kind) :: i, k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               vLevel => this%sg%vLevel, &
               coef => this%sg%coef_fstdif)

      ! Initialize A_adj to zero
      A_adj = 0.0_R_KIND

      ! Adjoint computation: propagate the gradient from parA_sigma_adj and rhs_adj to A_adj
      DO i = 1, num_cell
        DO k = 1, vLevel
          IF (k == 1) THEN
            A_adj(k:k + 2, i) = A_adj(k:k + 2, i) + coef(k, i, :) * parA_sigma_adj(k, i) + rhs_adj(k, i)
          ELSEIF (k == vLevel) THEN
            A_adj(k - 2:k, i) = A_adj(k - 2:k, i) + coef(k, i, :) * parA_sigma_adj(k, i) + rhs_adj(k, i)
          ELSE
            A_adj(k - 1:k + 1, i) = A_adj(k - 1:k + 1, i) + coef(k, i, :) * parA_sigma_adj(k, i) + rhs_adj(k, i)
          END IF
        END DO
      END DO
    END ASSOCIATE
  END SUBROUTINE FirstOrder_AD

  ! Second-order adjoint computation
  SUBROUTINE SecondOrder_AD(this, A_adj, para_sigma_2nd_adj, rhs_adj)
    CLASS(CalVerDer_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), POINTER :: A_adj(:, :), para_sigma_2nd_adj(:, :), rhs_adj(:, :)
    INTEGER(i_kind) :: i, k

    ! Initialize A_adj to zero
    A_adj = 0.0_R_KIND

    ! Adjoint computation: propagate the gradient from para_sigma_2nd_adj and rhs_adj to A_adj
    DO i = 1, this%sg%num_cell
      DO k = 1, this%sg%vLevel
        IF (k == 1) THEN
          A_adj(k:k + 2, i) = A_adj(k:k + 2, i) + this%sg%coef_scddif(k, i, :) * para_sigma_2nd_adj(k, i) + rhs_adj(k, i)
        ELSEIF (k == this%sg%vLevel) THEN
          A_adj(k - 2:k, i) = A_adj(k - 2:k, i) + this%sg%coef_scddif(k, i, :) * para_sigma_2nd_adj(k, i) + rhs_adj(k, i)
        ELSE
          A_adj(k - 1:k + 1, i) = A_adj(k - 1:k + 1, i) + this%sg%coef_scddif(k, i, :) * para_sigma_2nd_adj(k, i) + rhs_adj(k, i)
        END IF
      END DO
    END DO
  END SUBROUTINE SecondOrder_AD

  ! First-order adjoint computation with half-order accuracy
  SUBROUTINE FirstOrder_half_AD(this, A_adj, parA_sigma_half_adj, rhs_adj)
    IMPLICIT NONE
    CLASS(CalVerDer_adj_t), INTENT(INOUT) :: this
    REAL(r_kind), POINTER :: A_adj(:, :), parA_sigma_half_adj(:, :), rhs_adj(:, :)
    INTEGER(i_kind) :: i, k
    INTEGER :: sz

    ! Get the size of the first dimension of A_adj
    sz = SIZE(A_adj, 1)

    ASSOCIATE (num_cell => this%sg%num_cell, &
               vLevel => this%sg%vLevel, &
               coef => this%sg%coef_fstdif_half)

      ! Initialize parA_sigma_half_adj to zero
      parA_sigma_half_adj = 0.0_R_KIND

      ! Adjoint computation: propagate the gradient from A_adj and rhs_adj to parA_sigma_half_adj
      DO i = 1, num_cell
        DO k = 1, vLevel
          IF (sz >= 4) THEN
            IF (k == 1) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, :) * A_adj(k:k + 3, i)) + rhs_adj(k, i)
            ELSEIF (k == 2) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, :) * A_adj(k - 1:k + 2, i)) + rhs_adj(k, i)
            ELSEIF (k == vLevel) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, :) * A_adj(k - 4:k - 1, i)) + rhs_adj(k, i)
            ELSEIF (k == vLevel - 1) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, :) * A_adj(k - 3:k, i)) + rhs_adj(k, i)
            ELSE
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, :) * A_adj(k - 2:k + 1, i)) + rhs_adj(k, i)
            END IF
          ELSEIF (sz == 3) THEN
            IF (k == 1) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, 1:3) * A_adj(k:k + 2, i)) + rhs_adj(k, i)
            ELSEIF (k == vLevel) THEN
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, 1:3) * A_adj(k - 2:k, i)) + rhs_adj(k, i)
            ELSE
              parA_sigma_half_adj(k, i) = SUM(coef(k, i, 1:3) * A_adj(k - 1:k + 1, i)) + rhs_adj(k, i)
            END IF
          ELSE
            WRITE (*, *) "Error: Array size not supported in FirstOrder_half_AD subroutine"
            STOP
          END IF
        END DO
      END DO
    END ASSOCIATE
  END SUBROUTINE FirstOrder_half_AD

  ! Destructor for CalVerDer_adj_t type
  IMPURE ELEMENTAL SUBROUTINE destructor_adj(this)
    TYPE(CalVerDer_adj_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
  END SUBROUTINE destructor_adj

END MODULE CalVerDer_adj_m
