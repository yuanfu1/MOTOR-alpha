MODULE Derivative_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  TYPE :: Derivative_t
  CONTAINS
    PROCEDURE :: merge_sort_index
  END TYPE
  TYPE, EXTENDS(Derivative_t) :: HalfLayerShortestDist_t
  CONTAINS
    PROCEDURE :: FirstOrderHalfLayerShortestDist
    PROCEDURE :: FirstOrder => FirstOrderHalfLayerShortestDist
  END TYPE

  TYPE, EXTENDS(Derivative_t) :: HalfLayerMinNorm_t
  CONTAINS
    PROCEDURE, NOPASS :: FirstOrderHalfLayerMinNorm
    PROCEDURE, NOPASS :: FirstOrder => FirstOrderHalfLayerMinNorm
  END TYPE

  TYPE, EXTENDS(Derivative_t) :: FirstOrder_t
  CONTAINS
    PROCEDURE, NOPASS :: FirstOrderFullLayer
    PROCEDURE, NOPASS :: FirstOrder => FirstOrderFullLayer
  END TYPE

  TYPE, EXTENDS(Derivative_t) :: SecOrderShortestDist_t
  CONTAINS
    PROCEDURE :: SecondOrderShortestDist
    PROCEDURE :: SecondOrder => SecondOrderShortestDist
  END TYPE
  TYPE, EXTENDS(Derivative_t) :: SecOrderMinNorm_t
  CONTAINS
    PROCEDURE, NOPASS :: SecondOrderMinNorm
    PROCEDURE, NOPASS :: SecondOrder => SecondOrderMinNorm
  END TYPE

CONTAINS

  SUBROUTINE FirstOrderHalfLayerShortestDist(this, f, z, ks)

    IMPLICIT NONE

    CLASS(HalfLayerShortestDist_t) :: this
    REAL(r_kind), INTENT(IN) :: z(:)
    REAL(r_kind), INTENT(INOUT) :: f(:)
    INTEGER(i_kind), INTENT(IN) :: ks

    INTEGER(i_kind) :: i, j, sz, n, info
    INTEGER(i_kind), ALLOCATABLE :: id(:), ipiv(:)

    REAL(r_kind) :: zm
    REAL(r_kind), ALLOCATABLE :: Y(:, :), dist(:), dist2(:), A(:, :)

    n = SIZE(z)
    sz = n - 1
    f = 0.0D0

    ALLOCATE (id(n), ipiv(sz), Y(sz, 1), dist(n), dist2(n), A(sz, sz))
    IF (ks .EQ. -1) THEN
      zm = z(1) - (z(2) - z(1)) / 2.0D0
    ELSEIF (ks .EQ. -2) THEN
      zm = z(n) + (z(n) - z(n - 1)) / 2.0D0
    ELSE
      zm = (z(ks) + z(ks + 1)) / 2.0D0
    END IF

    dist = z - zm
    ! choose the grid points with the nearest distance
    dist2 = DABS(dist)

    ! print *, dist2
    CALL this%merge_sort_index(dist2, id)
    ! print *, 'id:', id

    dist2 = dist(id)
    ! PRINT *, 'dist2:', dist2
    ! creat A matrix
    DO i = 1, sz
      A(i, :) = dist2(1:sz)**(i - 1)
    END DO
    Y(:, 1) = [0.0D0, 1.0D0, 0.0D0]
    CALL DGESV(sz, 1, A, sz, ipiv, Y, sz, info)
    ! PRINT *, 'Y:', Y
    DO i = 1, sz
      f(id(i)) = Y(i, 1)
    END DO
    ! id2 = id(1:sz) - n/2

    DEALLOCATE (id, ipiv, Y, dist, dist2, A)

  END SUBROUTINE

  SUBROUTINE FirstOrderHalfLayerMinNorm(f, z, ks)

    IMPLICIT NONE

    ! CLASS(HalfLayerMinNorm) :: this
    REAL(r_kind), INTENT(IN) :: z(:)
    REAL(r_kind), INTENT(INOUT) :: f(:)
    INTEGER(i_kind), INTENT(IN) :: ks
    INTEGER(i_kind) :: i, j, sz, info
    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)
    REAL(r_kind) :: zm

    REAL(r_kind), ALLOCATABLE :: Yt(:), Y(:, :), dist(:), A(:, :), B(:, :), &
                                 sigma(:, :)

    sz = SIZE(z)
    f = 0.0D0

    ALLOCATE (ipiv(2 * sz - 1), Yt(sz - 1), Y(2 * sz - 1, 1), dist(sz), A(sz - 1, sz), &
              sigma(sz, sz), B(2 * sz - 1, 2 * sz - 1))

    sigma = 0.0D0
    A = 0.0D0
    B = 0.0D0
    Y = 0.0D0

    IF (ks .EQ. -1) THEN
      zm = z(1) - (z(2) - z(1)) / 2.0D0
    ELSEIF (ks .EQ. -2) THEN
      zm = z(sz) + (z(sz) - z(sz - 1)) / 2.0D0
    ELSE
      zm = (z(ks) + z(ks + 1)) / 2.0D0
    END IF

    dist = z - zm
    ! PRINT *, 'dist: ', dist
    ! creat A matrix
    DO i = 1, sz - 1
      A(i, :) = dist(1:sz)**(i - 1)
    END DO
    DO i = 1, sz
      sigma(i, i) = (dist(i)**3)
    END DO
    B(1:sz, 1:sz) = 2.0D0 * sigma
    B(sz + 1:2 * sz - 1, 1:sz) = A
    B(1:sz, sz + 1:2 * sz - 1) = TRANSPOSE(A)
    Yt = [0.0D0, 1.0D0, 0.0D0]

    Y(sz + 1:2 * sz - 1, 1) = Yt
    ! DO i = 1, 2*sz - 1
    !    PRINT *, 'B(i):', B(i, :)
    ! END DO
    CALL DGESV(2 * sz - 1, 1, B, 2 * sz - 1, ipiv, Y, 2 * sz - 1, info)
    ! PRINT *, 'ipiv:', ipiv
    ! PRINT *, 'Y:', Y
    f = Y(1:sz, 1)
    ! print *, 'f:', f
    ! id2 = id(1:sz) - n/2

    ! PRINT *, 'sizes of the parameters: ', size(ipiv), size(Yt), size(Y), size(dist), size(A), size(sigma), size(B)

    DEALLOCATE (ipiv, Yt, Y, dist, A, sigma, B)

  END SUBROUTINE

  SUBROUTINE FirstOrderFullLayer(f, z, ks)

    IMPLICIT NONE

    ! CLASS(FirstOrder) :: this
    REAL(r_kind), INTENT(IN) :: z(:)
    REAL(r_kind), INTENT(INOUT) :: f(:)
    INTEGER(i_kind), INTENT(IN) :: ks

    INTEGER(i_kind) :: i, j, sz, info
    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)

    REAL(r_kind), ALLOCATABLE :: Y(:, :), dist(:), A(:, :)

    sz = SIZE(z)
    f = 0.0D0

    ALLOCATE (ipiv(sz), Y(sz, 1), dist(sz), A(sz, sz))
    dist = z - z(ks)
    ! creat A matrix
    DO i = 1, sz
      A(i, :) = dist(1:sz)**(i - 1)
    END DO
    Y(:, 1) = [0.0D0, 1.0D0, 0.0D0]
    CALL DGESV(sz, 1, A, sz, ipiv, Y, sz, info)
    ! PRINT *, 'Y:', Y
    f = Y(:, 1)
    ! id2 = id(1:sz) - n/2

    DEALLOCATE (ipiv, Y, dist, A)

  END SUBROUTINE

  SUBROUTINE SecondOrderShortestDist(this, f, z, ks)

    IMPLICIT NONE

    CLASS(SecOrderShortestDist_t) :: this
    REAL(r_kind), INTENT(IN) :: z(:)
    REAL(r_kind), INTENT(INOUT) :: f(:)
    INTEGER(i_kind), INTENT(IN) :: ks

    INTEGER(i_kind) :: i, j, sz, n, info
    INTEGER(i_kind), ALLOCATABLE :: id(:), ipiv(:)

    ! REAL(r_kind) :: zm
    REAL(r_kind), ALLOCATABLE :: Y(:, :), dist(:), dist2(:), A(:, :)

    n = SIZE(z)
    sz = n - 1
    f = 0.0D0

    ALLOCATE (id(n), ipiv(sz), Y(sz, 1), dist(n), dist2(n), A(sz, sz))
    dist = z - z(ks)
    ! choose the grid points with the nearest distance
    dist2 = DABS(dist)

    ! print *, dist2
    CALL this%merge_sort_index(dist2, id)
    ! print *, 'id:', id

    dist2 = dist(id)
    ! PRINT *, 'dist2:', dist2
    ! creat A matrix
    DO i = 1, sz
      A(i, :) = dist2(1:sz)**(i - 1)
    END DO
    Y(:, 1) = [0.0D0, 0.0D0, 2.0D0, 0.0D0]
    CALL DGESV(sz, 1, A, sz, ipiv, Y, sz, info)
    ! PRINT *, 'Y:', Y
    DO i = 1, sz
      f(id(i)) = Y(i, 1)
    END DO
    ! id2 = id(1:sz) - n/2

    DEALLOCATE (id, ipiv, Y, dist, dist2, A)

  END SUBROUTINE

  SUBROUTINE SecondOrderMinNorm(f, z, ks)

    IMPLICIT NONE

    ! CLASS(SecOrderMinNorm) :: this
    REAL(r_kind), INTENT(IN) :: z(:)
    REAL(r_kind), INTENT(INOUT) :: f(:)
    INTEGER(i_kind), INTENT(IN) :: ks
    INTEGER(i_kind) :: i, j, sz, info
    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)

    REAL(r_kind), ALLOCATABLE :: Yt(:), Y(:, :), dist(:), A(:, :), B(:, :), &
                                 sigma(:, :)

    sz = SIZE(z)
    f = 0.0D0

    ALLOCATE (ipiv(2 * sz - 1), Yt(sz - 1), Y(2 * sz - 1, 1), dist(sz), A(sz - 1, sz), &
              sigma(sz, sz), B(2 * sz - 1, 2 * sz - 1))

    ! PRINT *, 'sizes of the parameters: ', size(ipiv), size(Yt), size(Y), size(dist), size(A), size(sigma), size(B)

    sigma = 0.0D0
    A = 0.0D0
    B = 0.0D0
    Y = 0.0D0
    dist = z - z(ks)
    ! PRINT *, 'dist: ', dist
    ! creat A matrix
    DO i = 1, sz - 1
      A(i, :) = dist(1:sz)**(i - 1)
    END DO
    DO i = 1, sz
      sigma(i, i) = dist(i)**4
    END DO
    B(1:sz, 1:sz) = 2.0D0 * sigma
    B(sz + 1:2 * sz - 1, 1:sz) = A
    B(1:sz, sz + 1:2 * sz - 1) = TRANSPOSE(A)
    Yt = [0.0D0, 0.0D0, 2.0D0, 0.0D0]

    Y(sz + 1:2 * sz - 1, 1) = Yt
    ! DO i = 1, 2*sz - 1
    !    PRINT *, 'B(i):', B(i, :)
    ! END DO
    CALL DGESV(2 * sz - 1, 1, B, 2 * sz - 1, ipiv, Y, 2 * sz - 1, info)
    ! PRINT *, 'ipiv:', ipiv
    ! PRINT *, 'Y:', Y
    f = Y(1:sz, 1)
    ! print *, 'f:', f
    ! id2 = id(1:sz) - n/2

    ! PRINT *, 'sizes of the parameters: ', size(ipiv), size(Yt), size(Y), size(dist), size(A), size(sigma), size(B)

    DEALLOCATE (ipiv, Yt, Y, dist, A, sigma, B)

  END SUBROUTINE

  SUBROUTINE merge_sort_index(this, a, idx)
    IMPLICIT NONE
    CLASS(Derivative_t) :: this
    REAL(r_kind), INTENT(inout) :: a(:)
    INTEGER(i_kind), INTENT(out) :: idx(:)
    INTEGER(i_kind) :: n, i

    n = SIZE(a, 1)
    IF (SIZE(a) /= SIZE(idx)) THEN
      PRINT *, 'ERROR in sorting, size of a is not equal to size of idx'
      STOP
    END IF

    DO i = 1, n
      idx(i) = i
    END DO

    CALL merge_sort(a, idx, n)
  END SUBROUTINE

  RECURSIVE SUBROUTINE merge_sort(a, idx, n)
    REAL(r_kind), INTENT(inout) :: a(n)
    INTEGER(i_kind), INTENT(inout) :: idx(n)
    INTEGER(i_kind), INTENT(in) :: n
    REAL(r_kind), DIMENSION(n) :: temp
    INTEGER(i_kind), DIMENSION(n) :: temp_idx
    INTEGER(i_kind) :: i, j, k, mid

    IF (n <= 1) THEN
      RETURN
    END IF

    mid = n / 2

    CALL merge_sort(a, idx, mid)
    CALL merge_sort(a(mid + 1), idx(mid + 1), n - mid)

    i = 1
    j = mid + 1
    DO k = 1, n
      ! print *, k
      IF (j > n) THEN
        temp(k) = a(i)
        temp_idx(k) = idx(i)
        i = i + 1
      ELSE IF (i <= mid .AND. a(i) <= a(j)) THEN
        temp(k) = a(i)
        temp_idx(k) = idx(i)
        i = i + 1
      ELSE
        temp(k) = a(j)
        temp_idx(k) = idx(j)
        j = j + 1
      END IF
    END DO

    a(1:n) = temp(1:n)
    idx(1:n) = temp_idx(1:n)
  END SUBROUTINE merge_sort

END MODULE Derivative_m
