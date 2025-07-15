MODULE CalVerDer_TL_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  IMPLICIT NONE

  TYPE :: CalVerDer_TL_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE :: FirstOrder_TL
    PROCEDURE :: SecondOrder_TL
    PROCEDURE :: FirstOrder_half_TL
    PROCEDURE :: initialize
    FINAL :: destructor
  END TYPE CalVerDer_TL_t

CONTAINS

  SUBROUTINE initialize(this)
    CLASS(CalVerDer_TL_t), INTENT(INOUT) :: this
    IF (.NOT. ALLOCATED(this%sg%coef_fstdif)) THEN
      ALLOCATE (this%sg%coef_fstdif(3, this%sg%vLevel, this%sg%num_cell))
      this%sg%coef_fstdif = 1.0_R_KIND
    END IF
    IF (.NOT. ALLOCATED(this%sg%coef_scddif)) THEN
      ALLOCATE (this%sg%coef_scddif(5, this%sg%vLevel, this%sg%num_cell))
      this%sg%coef_scddif = 1.0_R_KIND
    END IF
    IF (.NOT. ALLOCATED(this%sg%coef_fstdif_half)) THEN
      ALLOCATE (this%sg%coef_fstdif_half(4, this%sg%vLevel, this%sg%num_cell))
      this%sg%coef_fstdif_half = 1.0_R_KIND
    END IF
  END SUBROUTINE initialize

  SUBROUTINE FirstOrder_TL(this, x, perturbation, RESULT)
    CLASS(CalVerDer_TL_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: x(:, :)
    REAL(r_kind), INTENT(IN) :: perturbation(:, :)
    REAL(r_kind), INTENT(OUT) :: RESULT(SIZE(x, 1), SIZE(x, 2))

    INTEGER(i_kind) :: i, k

    CALL this%initialize()

    RESULT = 0.0_R_KIND
    DO i = 1, this%sg%num_cell
      DO k = 1, this%sg%vLevel
        IF (k .EQ. 1) THEN
          RESULT(k, i) = SUM(this%sg%coef_fstdif(:, k, i) * perturbation(k:k + 2, i))
        ELSEIF (k .EQ. this%sg%vLevel) THEN
          RESULT(k, i) = SUM(this%sg%coef_fstdif(:, k, i) * perturbation(k - 2:k, i))
        ELSE
          RESULT(k, i) = SUM(this%sg%coef_fstdif(:, k, i) * perturbation(k - 1:k + 1, i))
        END IF
      END DO
    END DO
  END SUBROUTINE FirstOrder_TL

  SUBROUTINE SecondOrder_TL(this, x, perturbation, RESULT)
    CLASS(CalVerDer_TL_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: x(:, :)
    REAL(r_kind), INTENT(IN) :: perturbation(:, :)
    REAL(r_kind), INTENT(OUT) :: RESULT(SIZE(x, 1), SIZE(x, 2))

    INTEGER(i_kind) :: i, k

    CALL this%initialize()

    RESULT = 0.0_R_KIND
    DO i = 1, this%sg%num_cell
      DO k = 1, this%sg%vLevel
        IF (k .EQ. 1) THEN
          RESULT(k, i) = SUM(this%sg%coef_scddif(:, k, i) * perturbation(k:k + 4, i))
        ELSEIF (k .EQ. 2) THEN
          RESULT(k, i) = SUM(this%sg%coef_scddif(:, k, i) * perturbation(k - 1:k + 3, i))
        ELSEIF (k .EQ. this%sg%vLevel) THEN
          RESULT(k, i) = SUM(this%sg%coef_scddif(:, k, i) * perturbation(k - 4:k, i))
        ELSEIF (k .EQ. this%sg%vLevel - 1) THEN
          RESULT(k, i) = SUM(this%sg%coef_scddif(:, k, i) * perturbation(k - 3:k + 1, i))
        ELSE
          RESULT(k, i) = SUM(this%sg%coef_scddif(:, k, i) * perturbation(k - 2:k + 2, i))
        END IF
      END DO
    END DO
  END SUBROUTINE SecondOrder_TL

  SUBROUTINE FirstOrder_half_TL(this, x, perturbation, RESULT)
    CLASS(CalVerDer_TL_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN) :: x(:, :)
    REAL(r_kind), INTENT(IN) :: perturbation(:, :)
    REAL(r_kind), INTENT(OUT) :: RESULT(SIZE(x, 1), SIZE(x, 2))

    INTEGER(i_kind) :: i, k

    CALL this%initialize()

    RESULT = 0.0_R_KIND
    DO i = 1, this%sg%num_cell
      DO k = 1, this%sg%vLevel
        IF (k .EQ. 1) THEN
          RESULT(k, i) = SUM(this%sg%coef_fstdif_half(:, k, i) * perturbation(k:k + 3, i))
        ELSEIF (k .EQ. this%sg%vLevel) THEN
          RESULT(k, i) = SUM(this%sg%coef_fstdif_half(:, k, i) * perturbation(k - 2:k + 1, i))
        ELSE
          RESULT(k, i) = SUM(this%sg%coef_fstdif_half(:, k, i) * perturbation(k - 1:k + 2, i))
        END IF
      END DO
    END DO
  END SUBROUTINE FirstOrder_half_TL

  SUBROUTINE destructor(this)
    TYPE(CalVerDer_TL_t), INTENT(INOUT) :: this

    ! Deallocate and nullify pointers
    IF (ASSOCIATED(this%sg)) THEN
      IF (ALLOCATED(this%sg%coef_fstdif)) DEALLOCATE (this%sg%coef_fstdif)
      IF (ALLOCATED(this%sg%coef_scddif)) DEALLOCATE (this%sg%coef_scddif)
      IF (ALLOCATED(this%sg%coef_fstdif_half)) DEALLOCATE (this%sg%coef_fstdif_half)
      NULLIFY (this%sg)
    END IF
  END SUBROUTINE destructor

END MODULE CalVerDer_TL_m
