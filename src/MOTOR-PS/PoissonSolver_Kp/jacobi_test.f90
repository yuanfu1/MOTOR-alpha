
PROGRAM jacobi_test

  !>
  !!=================================================================
  !!  This is a unit test program of a Jacobi iteration method for a
  !!  linear system of equations Ax = b derived from a Poisson system
  !!
  !!
  !!  \author Yuanfu Xie
  !!  \b history
  !!    Created by Yuanfu Xie Sep. 2018
  !!    Modified by Yuanfu Xie Jun. 2019 for a new sparse matrix data
  !!    structure.
  !!=================================================================
  !

  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE unitPoisson_m

  IMPLICIT NONE

  INTEGER(i_kind) :: ndim, mitr, i, ngrid(3), io
  LOGICAL :: passed
  REAL(r_kind), ALLOCATABLE :: truth(:), omega(:), guess(:)

  ! These errors are set from an independent jacobi and Gauss-Seidel
  ! implementation under /Users/xieyuanfu/developments/models/gzm/src/util
  ! test_jacobi.f90
  ! If the test code matches these errors, we pass.
  REAL(r_kind) :: jacobi_err(2), gaussSeidel_err(2), err(2)

  TYPE(unitPoisson_t) :: unit_test

  passed = .TRUE.

  ! Note: the unit test is set for the following parameters according
  ! to the test code of test_jacobi.f90 mentioned above.
  ngrid(1) = 65
  ngrid(2) = 65
  ngrid(3) = 1
  jacobi_err(1) = 0.9781    ! Relative error
  jacobi_err(2) = 1.5193    ! maximum error
  gaussSeidel_err(1) = 0.9455   ! Relative error
  gaussSeidel_err(2) = 1.4698   ! maximum error

  io = 1

  ! Number of iterations requested:
  mitr = 2

  CALL unit_test%construct(ngrid)

  ALLOCATE (truth(unit_test%num_eqns), guess(unit_test%num_eqns), omega(mitr))

  CALL unit_test%loadMatrx()

  CALL unit_test%Analytics(truth)

  omega = 0.85
  omega(1) = 1.3895D0
  omega(2) = 0.5617D0

  ! Testing Jacobi method:
  guess = 0.0D0
  CALL jacobi(unit_test%num_eqns, unit_test%num_vlvl, unit_test%max_widt, &
              unit_test%num_elem, unit_test%idx_colu, io, &
              unit_test%row_elem, mitr, omega, &
              unit_test%rgt_hand, guess, &
              unit_test%sol_curr)

  CALL unit_test%CalResidu()

  ! Check the solution accuracy:
  CALL unit_test%CalErrors(truth, err)
  IF (ABS(err(1) - jacobi_err(1)) .GE. 1.0D-4 .OR. &
      ABS(err(2) - jacobi_err(2)) .GE. 1.0D-4) passed = .FALSE.

  DO i = 1, unit_test%num_eqns
    WRITE (11, *) i, truth(i), unit_test%sol_curr(1, i)
  END DO

  ! Testing Gauss-Seidel method:
  guess = 0.0D0
  CALL GaussSeidel(unit_test%num_eqns, unit_test%num_vlvl, unit_test%max_widt, &
                   unit_test%num_elem, unit_test%idx_colu, io, &
                   unit_test%row_elem, mitr, omega, &
                   unit_test%rgt_hand, guess, &
                   unit_test%sol_curr)

  CALL unit_test%CalResidu()

  ! Check the solution accuracy:
  CALL unit_test%CalErrors(truth, err)
  IF (ABS(err(1) - gaussSeidel_err(1)) .GE. 1.0D-4 .OR. &
      ABS(err(2) - gaussSeidel_err(2)) .GE. 1.0D-4) passed = .FALSE.

  DO i = 1, unit_test%num_eqns
    WRITE (12, *) i, truth(i), unit_test%sol_curr(1, i)
  END DO

  IF (passed) THEN
    PRINT *, 'Test passed'
  ELSE
    PRINT *, 'Test failed'
  END IF

  ! Deallocate memory:
  CALL unit_test%deconstru()
  DEALLOCATE (truth)

END PROGRAM jacobi_test
