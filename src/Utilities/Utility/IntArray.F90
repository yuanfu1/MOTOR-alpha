
MODULE intArray_m
  USE kinds_m, ONLY: i_kind
  TYPE :: intArray_t
    INTEGER(i_kind) :: num_elements
    INTEGER(i_kind), ALLOCATABLE :: iarray(:)
  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: initialize_s
  END TYPE intArray_t

  INTERFACE intArray_t
    PROCEDURE :: constructor
  END INTERFACE intArray_t

CONTAINS
  SUBROUTINE initialize_s(this, num)
    CLASS(intArray_t) :: this
    INTEGER(i_kind), INTENT(IN) :: num

    this%num_elements = num
    ALLOCATE (this%iarray(num))
  END SUBROUTINE initialize_s

  FUNCTION constructor(num) RESULT(this)
    IMPLICIT NONE
    TYPE(intArray_t) :: this
    INTEGER(i_kind), INTENT(IN) :: num

    this%num_elements = num
    ALLOCATE (this%iarray(num))
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(intArray_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%iarray)) DEALLOCATE (this%iarray)
  END SUBROUTINE destructor
END MODULE intArray_m
