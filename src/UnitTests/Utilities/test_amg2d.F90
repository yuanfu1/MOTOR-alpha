PROGRAM test_amg2d
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE amg2d_m, ONLY: amg2d_t, amg_level_t  ! Added amg_level_t
  IMPLICIT NONE

  ! Add pi as a parameter
  REAL(c_double), PARAMETER :: pi = 3.14159265358979323846_c_double

  TYPE(amg2d_t) :: amg
  INTEGER(c_int) :: nx, ny, n, i, j, k
  REAL(c_double), ALLOCATABLE :: b(:), x(:), exact(:)
  REAL(c_double) :: h, tol, error

  ! Set up test problem parameters
  nx = 32  ! Grid points in x
  ny = 32  ! Grid points in y
  n = nx * ny
  h = 1.0_c_double / (nx + 1)
  tol = 1.0E-10_c_double

  ! Allocate arrays
  ALLOCATE(b(n), x(n), exact(n))

  ! Initialize AMG solver
  CALL amg%initialize(n, tol, 100)

  ! Set up system matrix (2D Poisson equation)
  CALL setup_poisson_matrix(amg%levels(1), nx, ny, h)

  ! Set up right-hand side and exact solution
  DO j = 1, ny
    DO i = 1, nx
      k = (j-1)*nx + i
      ! Example: u(x,y) = sin(πx)sin(πy)
      exact(k) = SIN(pi*i*h) * SIN(pi*j*h)
      ! f = -Δu
      b(k) = 2.0_c_double * (pi**2) * exact(k)
    END DO
  END DO

  ! Setup AMG hierarchy
  CALL amg%setup()

  ! Solve system
  CALL amg%solve(b, x)

  ! Compute error
  error = MAXVAL(ABS(x - exact))
  PRINT *, 'Maximum error:', error

  ! Clean up
  CALL amg%destroy()
  DEALLOCATE(b, x, exact)

CONTAINS

  SUBROUTINE setup_poisson_matrix(level, nx, ny, h)
    TYPE(amg_level_t), INTENT(INOUT) :: level
    INTEGER(c_int), INTENT(IN) :: nx, ny
    REAL(c_double), INTENT(IN) :: h
    INTEGER(c_int) :: i, j, k, row

    ! 5-point stencil for 2D Poisson
    level%ia(1) = 1
    k = 1
    DO j = 1, ny
      DO i = 1, nx
        row = (j-1)*nx + i
        level%ia(row+1) = level%ia(row)

        ! Center point
        k = level%ia(row+1)
        level%ja(k) = row
        level%a(k) = 4.0_c_double/(h*h)
        level%ia(row+1) = level%ia(row+1) + 1

        ! Left neighbor
        IF (i > 1) THEN
          k = level%ia(row+1)
          level%ja(k) = row - 1
          level%a(k) = -1.0_c_double/(h*h)
          level%ia(row+1) = level%ia(row+1) + 1
        END IF

        ! Right neighbor
        IF (i < nx) THEN
          k = level%ia(row+1)
          level%ja(k) = row + 1
          level%a(k) = -1.0_c_double/(h*h)
          level%ia(row+1) = level%ia(row+1) + 1
        END IF

        ! Bottom neighbor
        IF (j > 1) THEN
          k = level%ia(row+1)
          level%ja(k) = row - nx
          level%a(k) = -1.0_c_double/(h*h)
          level%ia(row+1) = level%ia(row+1) + 1
        END IF

        ! Top neighbor
        IF (j < ny) THEN
          k = level%ia(row+1)
          level%ja(k) = row + nx
          level%a(k) = -1.0_c_double/(h*h)
          level%ia(row+1) = level%ia(row+1) + 1
        END IF
      END DO
    END DO
  END SUBROUTINE setup_poisson_matrix

END PROGRAM test_amg2d
