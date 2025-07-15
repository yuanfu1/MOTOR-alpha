MODULE ilu_precond_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ilu_t, create_ilu, apply_ilu

  TYPE :: ilu_t
    INTEGER(c_int) :: n           ! Matrix size
    INTEGER(c_int), ALLOCATABLE :: ia(:)  ! Row pointers
    INTEGER(c_int), ALLOCATABLE :: ja(:)  ! Column indices
    REAL(c_double), ALLOCATABLE :: l(:)   ! Lower triangular part
    REAL(c_double), ALLOCATABLE :: u(:)   ! Upper triangular part
    REAL(c_double), ALLOCATABLE :: d(:)   ! Diagonal
  END TYPE ilu_t

CONTAINS
  SUBROUTINE create_ilu(ilu, ia, ja, a, n)
    TYPE(ilu_t), INTENT(OUT) :: ilu
    INTEGER(c_int), INTENT(IN) :: n, ia(:), ja(:)
    REAL(c_double), INTENT(IN) :: a(:)
    INTEGER(c_int) :: i, j, k
    REAL(c_double) :: sum

    ! Initialize ILU structure
    ilu%n = n
    ALLOCATE(ilu%ia(n+1), ilu%ja(SIZE(ja)))
    ALLOCATE(ilu%l(SIZE(ja)), ilu%u(SIZE(ja)), ilu%d(n))
    
    ! Copy matrix structure
    ilu%ia = ia
    ilu%ja = ja

    ! Perform ILU(0) factorization
    DO i = 1, n
      ! Copy row i
      DO k = ia(i), ia(i+1)-1
        j = ja(k)
        IF (j < i) THEN
          ilu%l(k) = a(k)  ! Lower part
        ELSE IF (j == i) THEN
          ilu%d(i) = a(k)  ! Diagonal
        ELSE
          ilu%u(k) = a(k)  ! Upper part
        END IF
      END DO

      ! Eliminate entries in row i
      DO k = ia(i), ia(i+1)-1
        IF (ja(k) < i) THEN
          sum = ilu%l(k) / ilu%d(ja(k))
          ilu%l(k) = sum
          ! Update remaining entries in row i
          DO j = k+1, ia(i+1)-1
            ilu%u(j) = ilu%u(j) - sum * ilu%u(ia(ja(k))+j-k)
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE create_ilu

  ! Apply ILU preconditioner: M^(-1)v where M = LU
  SUBROUTINE apply_ilu(ilu, v, result)
    TYPE(ilu_t), INTENT(IN) :: ilu
    REAL(c_double), INTENT(IN) :: v(:)
    REAL(c_double), INTENT(OUT) :: result(:)
    REAL(c_double), ALLOCATABLE :: y(:)
    INTEGER(c_int) :: i, j, k

    ALLOCATE(y(ilu%n))
    
    ! Forward solve Ly = v
    DO i = 1, ilu%n
      y(i) = v(i)
      DO k = ilu%ia(i), ilu%ia(i+1)-1
        j = ilu%ja(k)
        IF (j < i) THEN
          y(i) = y(i) - ilu%l(k) * y(j)
        END IF
      END DO
    END DO

    ! Backward solve Ux = y
    DO i = ilu%n, 1, -1
      result(i) = y(i) / ilu%d(i)
      DO k = ilu%ia(i), ilu%ia(i+1)-1
        j = ilu%ja(k)
        IF (j > i) THEN
          result(i) = result(i) - ilu%u(k) * result(j) / ilu%d(i)
        END IF
      END DO
    END DO

    DEALLOCATE(y)
  END SUBROUTINE apply_ilu
END MODULE ilu_precond_m
