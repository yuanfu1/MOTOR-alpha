MODULE Matrix_Utils
  USE kinds_m, ONLY: r_kind
  IMPLICIT NONE

CONTAINS

  FUNCTION IDENTITY_MATRIX(size) RESULT(matrix)
    INTEGER, INTENT(IN) :: size
    REAL(r_kind) :: matrix(size, size)
    INTEGER :: i

    matrix = 0.0_R_KIND
    DO i = 1, size
      matrix(i, i) = 1.0_R_KIND
    END DO
  END FUNCTION IDENTITY_MATRIX

END MODULE Matrix_Utils
