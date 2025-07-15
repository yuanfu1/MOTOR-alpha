MODULE InterpCoef_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  ! TYPE :: InterpCoef_t
  ! CONTAINS
  !    PROCEDURE, NOPASS :: InterpCoef_3points
  !    PROCEDURE, NOPASS :: InterpCoef => InterpCoef_3points
  ! END TYPE

  ! TYPE :: InterpCoefMinNorm_t
  ! CONTAINS
  !    PROCEDURE, NOPASS :: InterpCoefMinNorm
  !    PROCEDURE, NOPASS :: InterpCoef => InterpCoefMinNorm
  ! END TYPE

CONTAINS

  ! SUBROUTINE InterpCoef_3points(z, zinterp, f)

  !    IMPLICIT NONE

  !    REAL(r_kind), INTENT(IN) :: z(:), zinterp
  !    REAL(r_kind), INTENT(INOUT) :: f(:)

  !    INTEGER(i_kind) :: i, j, sz, info
  !    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)

  !    REAL(r_kind), ALLOCATABLE :: Y(:, :), dist(:), A(:, :)

  !    sz = SIZE(z)
  !    f = 0.0D0

  !    ALLOCATE (ipiv(sz), Y(sz, 1), dist(sz), A(sz, sz))
  !    dist = z - zinterp
  !    ! creat A matrix
  !    DO i = 1, sz
  !       A(i, :) = dist(1:sz)**(i - 1)
  !    END DO
  !    Y(:, 1) = [1.0D0, 0.0D0, 0.0D0]
  !    CALL DGESV(sz, 1, A, sz, ipiv, Y, sz, info)
  !    ! PRINT *, 'Y:', Y
  !    f = Y(:, 1)
  !    ! id2 = id(1:sz) - n/2

  !    DEALLOCATE (ipiv, Y, dist, A)

  ! END SUBROUTINE

  SUBROUTINE InterpCoef(z, zinterp, f)

    IMPLICIT NONE

    REAL(r_kind), INTENT(IN) :: z(:), zinterp
    REAL(r_kind), INTENT(INOUT) :: f(:)

    INTEGER(i_kind) :: i, j, sz, info
    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)

    REAL(r_kind), ALLOCATABLE :: Y(:, :), dist(:), A(:, :)

    sz = SIZE(z)
    f = 0.0D0

    ALLOCATE (ipiv(sz), Y(sz, 1), dist(sz), A(sz, sz))
    dist = z - zinterp
    ! creat A matrix
    DO i = 1, sz
      A(i, :) = dist(1:sz)**(i - 1)
    END DO
    Y(:, 1) = [1.0D0, 0.0D0, 0.0D0, 0.0D0]
    CALL DGESV(sz, 1, A, sz, ipiv, Y, sz, info)
    ! PRINT *, 'Y:', Y
    f = Y(:, 1)
    ! id2 = id(1:sz) - n/2

    DEALLOCATE (ipiv, Y, dist, A)

  END SUBROUTINE

  ! SUBROUTINE InterpCoef(z, zinterp, f)

  !    IMPLICIT NONE

  !    ! CLASS(HalfLayerMinNorm) :: this
  !    REAL(r_kind), INTENT(IN) :: z(:), zinterp
  !    REAL(r_kind), INTENT(INOUT) :: f(:)
  !    INTEGER(i_kind) :: i, j, sz, info
  !    INTEGER(i_kind), ALLOCATABLE :: ipiv(:)

  !    REAL(r_kind), ALLOCATABLE :: Yt(:), Y(:, :), dist(:), A(:, :), B(:, :), &
  !                                 sigma(:, :)

  !    sz = SIZE(z)
  !    f = 0.0D0

  !    ALLOCATE (ipiv(2*sz - 1), Yt(sz - 1), Y(2*sz - 1, 1), dist(sz), A(sz - 1, sz), &
  !              sigma(sz, sz), B(2*sz - 1, 2*sz - 1))

  !    sigma = 0.0D0
  !    A = 0.0D0
  !    B = 0.0D0
  !    Y = 0.0D0

  !    dist = z - zinterp
  !    ! PRINT *, 'dist: ', dist
  !    ! creat A matrix
  !    DO i = 1, sz - 1
  !       A(i, :) = dist(1:sz)**(i - 1)
  !    END DO
  !    DO i = 1, sz
  !       sigma(i, i) = (dist(i)**3)
  !    END DO
  !    B(1:sz, 1:sz) = 2.0D0*sigma
  !    B(sz + 1:2*sz - 1, 1:sz) = A
  !    B(1:sz, sz + 1:2*sz - 1) = TRANSPOSE(A)
  !    Yt = [1.0D0, 0.0D0, 0.0D0]

  !    Y(sz + 1:2*sz - 1, 1) = Yt
  !    ! DO i = 1, 2*sz - 1
  !    !    PRINT *, 'B(i):', B(i, :)
  !    ! END DO
  !    CALL DGESV(2*sz - 1, 1, B, 2*sz - 1, ipiv, Y, 2*sz - 1, info)
  !    ! PRINT *, 'ipiv:', ipiv
  !    ! PRINT *, 'Y:', Y
  !    f = Y(1:sz, 1)
  !    ! print *, 'f:', f
  !    ! id2 = id(1:sz) - n/2

  !    ! PRINT *, 'sizes of the parameters: ', size(ipiv), size(Yt), size(Y), size(dist), size(A), size(sigma), size(B)

  !    DEALLOCATE (ipiv, Yt, Y, dist, A, sigma, B)

  ! END SUBROUTINE
END MODULE
