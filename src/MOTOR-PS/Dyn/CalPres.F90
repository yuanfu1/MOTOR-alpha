MODULE CalPres_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: k_d, r_d
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE InterpValue_m, ONLY: InterpValue_t

  IMPLICIT NONE
  TYPE :: CalPres_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
  CONTAINS
    PROCEDURE :: CalPres
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE

  INTERFACE CalPres_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS
  FUNCTION constructor(sg) RESULT(this)
    TYPE(CalPres_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg

    this%sg => sg

    this%InterpValue = InterpValue_t(this%sg)

  END FUNCTION constructor
  SUBROUTINE CalPres(this, rho, theta, pres)
    IMPLICIT NONE
    CLASS(CalPres_t) :: this
    REAL(r_kind), INTENT(IN) :: rho(:, :), theta(:, :)
    REAL(r_kind), INTENT(inout) :: pres(:, :)
    REAL(r_kind), ALLOCATABLE :: theta2(:, :)

    ALLOCATE (theta2(this%sg%vLevel, this%sg%num_cell))

    theta2 = 2.0D0

    CALL this%InterpValue%InterpValue(theta, theta2)

    pres(:, 1:this%sg%num_iCell) = 10.0D0**(-5.0D0 * k_d / (1.0D0 - k_d)) &
                                   * (rho(:, 1:this%sg%num_iCell) &
                                      * r_d * theta2(:, 1:this%sg%num_iCell))**(1 / (1.0D0 - k_d))

    DEALLOCATE (theta2)

  END SUBROUTINE

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(CalPres_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destroy

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(CalPres_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor
END MODULE
