MODULE InterpValue_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  TYPE :: InterpValue_t
    TYPE(SingleGrid_t), POINTER :: sg
  CONTAINS
    PROCEDURE, NOPASS :: UpdateSgInterpCoef
    PROCEDURE :: InterpValue
    FINAL :: destructor
  END TYPE InterpValue_t

  INTERFACE InterpValue_t
    PROCEDURE :: constructor
  END INTERFACE InterpValue_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    IMPLICIT NONE
    TYPE(InterpValue_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(INOUT) :: sg

    this%sg => sg

  END FUNCTION constructor

  SUBROUTINE UpdateSgInterpCoef(sg_sigma_in, sg_to_update)
    USE GenInterpCoef_m

    IMPLICIT NONE

    CLASS(SingleGrid_t), INTENT(IN) :: sg_sigma_in
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_to_update

    ALLOCATE (sg_to_update%InterpCoef(sg_to_update%vLevel, sg_to_update%num_cell, 4), &
              sg_to_update%ks_interp(sg_to_update%vLevel, sg_to_update%num_cell))
    sg_to_update%InterpCoef = 0.0D0
    sg_to_update%ks_interp = 0
    CALL GenInterpCoef(sg_sigma_in%sigma, sg_to_update%sigma, sg_to_update%InterpCoef, sg_to_update%ks_interp)
    ! Yuanfu Xie temporarily turned these output as they serve as debugging on 2025-02-13
    ! PRINT *, 'sg%sigma: ', sg_sigma_in%sigma, '===', sg_to_update%sigma
    ! PRINT *, 'ks is: ', sg_to_update%ks_interp(:, 10), sg_to_update%InterpCoef(2, 10, :)

  END SUBROUTINE

  SUBROUTINE InterpValue(this, valuein, valueout)

    IMPLICIT NONE
    CLASS(InterpValue_t) :: this
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: valuein
    REAL(r_kind), INTENT(OUT) :: valueout(:, :)
    INTEGER(i_kind) :: i, k, sz

    ASSOCIATE (vLevel => this%sg%vLevel, &
               ks => this%sg%ks_interp, &
               num_icell => this%sg%num_icell)

      sz = SIZE(valuein(:, 1))

      ! PRINT *, 'sz is', sz

      ! IF (sz .LT. vLevel) THEN
      !    DO k = 1, vLevel
      !       IF (k .EQ. 1) THEN
      !          DO i = 1, num_icell
      !             valueout(k, i) = sum(valuein(k:k + 2, i)*this%sg%InterpCoef(k, i, :))
      !          END DO
      !       ELSEIF (k .EQ. vLevel) THEN
      !          DO i = 1, num_icell
      !             valueout(k, i) = sum(valuein(k - 3:k - 1, i)*this%sg%InterpCoef(k, i, :))
      !          END DO
      !       ELSE
      !          DO i = 1, num_icell
      !             valueout(k, i) = sum(valuein(k - 1:k, i)*this%sg%InterpCoef(k, i, 1:2))
      !          END DO
      !       END IF
      !    END DO
      ! ELSE
      !    DO k = 1, vLevel
      !       DO i = 1, num_icell
      !          valueout(k, i) = sum(valuein(k:k + 1, i)*this%sg%InterpCoef(k, i, 1:2))
      !       END DO
      !    END DO
      ! END IF
      PRINT *, 'sg%ks_interp is: ', SIZE(ks, 1), SIZE(ks, 2), MAXVAL(ks), MINVAL(ks)
      DO k = 1, vLevel
        DO i = 1, num_icell
          IF (ks(k, i) .EQ. 0) THEN
            PRINT *, 'ks has 0 in ', k, i
            PRINT *, ks(:, i - 1)
            STOP
          END IF
          valueout(k, i) = SUM(valuein(ks(k, i):ks(k, i) + 3, i) * this%sg%InterpCoef(k, i, :))
        END DO
      END DO
    END ASSOCIATE

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(InterpValue_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor

END MODULE InterpValue_m
