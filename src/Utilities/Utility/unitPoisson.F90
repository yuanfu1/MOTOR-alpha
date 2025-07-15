MODULE unitPoisson_m

  !>
  !!=================================================================
  !!  This module defines various test cases for Relaxed Jacobi (RJ)
  !!  methods. It will serve as a unit test for these methods.
  !!
  !!  \author Yuanfu Xie
  !!  \b History
  !!    Created by Yuanfu Xie May 2019
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double

  IMPLICIT NONE

  PRIVATE

  PUBLIC unitPoisson_t

  TYPE :: unitPoisson_t
    INTEGER(i_kind) :: num_eqns, num_vlvl, max_widt
    INTEGER(i_kind) :: num_grid(2)
    INTEGER(i_kind), ALLOCATABLE :: num_elem(:), idx_colu(:, :)
    REAL(r_kind) :: del_dxdy(2)
    REAL(r_kind), ALLOCATABLE :: row_elem(:, :), sol_curr(:, :), &
                                 rgt_hand(:, :), res_dual(:, :)

  CONTAINS
    PROCEDURE, PUBLIC :: Analytics
    PROCEDURE, PUBLIC :: Construct
    PROCEDURE, PUBLIC :: Deconstru
    PROCEDURE, PUBLIC :: LoadMatrx
    PROCEDURE, PUBLIC :: CalResidu
    PROCEDURE, PUBLIC :: CalErrors
  END TYPE unitPoisson_t

CONTAINS

  !>
    !!===============================================================
    !! Constructor: explicit as I felt it has some advantages
    !!===============================================================
  !
  SUBROUTINE Construct(this, ngrd)
    CLASS(unitPoisson_t) :: this
    INTEGER(i_kind), INTENT(IN) :: ngrd(3)

    this%num_eqns = ngrd(1) * ngrd(2)
    this%num_vlvl = ngrd(3)
    this%max_widt = 5
    this%num_grid = ngrd(1:2)
    this%del_dxdy = 1.0D0 / DBLE(this%num_grid(1:2) + 1)

    ALLOCATE (this%num_elem(this%num_eqns), &
              this%idx_colu(this%max_widt, this%num_eqns), &
              this%row_elem(this%max_widt, this%num_eqns), &
              this%sol_curr(this%num_vlvl, this%num_eqns), &
              this%rgt_hand(this%num_vlvl, this%num_eqns), &
              this%res_dual(this%num_vlvl, this%num_eqns))
  END SUBROUTINE Construct

  !>
    !!===============================================================
    !! Deconstructor: explicit as I felt it has some advantages
    !!===============================================================
  !
  SUBROUTINE Deconstru(this)
    CLASS(unitPoisson_t) :: this

    DEALLOCATE (this%num_elem, this%idx_colu, this%row_elem, &
                this%sol_curr, this%rgt_hand, this%res_dual)
  END SUBROUTINE Deconstru

  !>
    !!===============================================================
    !! Analytic sol_curr: use a statement function for unit_test
    !!===============================================================
  !
  SUBROUTINE Analytics(this, tru)
    CLASS(unitPoisson_t) :: this
    REAL(r_kind), INTENT(OUT) :: tru(this%num_eqns)

    ! Local variables:
    INTEGER(i_kind) :: i, j, iequation
    REAL(r_kind) :: x, y, pi, truth, lapla

    truth(x, y) = SIN(2.0D0 * pi * x) * SIN(2.0D0 * pi * y) + SIN(4.0D0 * pi * x) * SIN(4.0D0 * pi * y)
    lapla(x, y) = -8.0D0 * pi * pi * SIN(2.0D0 * pi * x) * SIN(2.0D0 * pi * y) &
                  - 32.0D0 * pi * pi * SIN(4.0D0 * pi * x) * SIN(4.0D0 * pi * y)

    pi = 4.0D0 * ATAN(1.0D0)

    ! True sol_curr and right hand side:
    ! If the second dimension is 1, please select the analytic
    ! sol_curr such that f(x,y) = f(x), notice y=1/2 if num_grid(2)=1.
    DO j = 1, this%num_grid(2)
      y = DBLE(j) / DBLE(this%num_grid(2) + 1)
      DO i = 1, this%num_grid(1)
        x = DBLE(i) / DBLE(this%num_grid(1) + 1)

        iequation = i + this%num_grid(1) * (j - 1)
        tru(iequation) = truth(x, y)

        this%rgt_hand(:, iequation) = lapla(x, y)
      END DO
    END DO

    ! Boundary conditions:
    DO j = 1, this%num_grid(2)
      y = DBLE(j) / DBLE(this%num_grid(2) + 1)
      i = 0
      x = DBLE(i) / DBLE(this%num_grid(1) + 1)
      iequation = 1 + this%num_grid(1) * (j - 1)
      this%rgt_hand(:, iequation) = this%rgt_hand(:, iequation) - &
                                    truth(x, y) / (this%del_dxdy(1))**2

      i = this%num_grid(1) + 1
      x = DBLE(i) / DBLE(this%num_grid(1) + 1)
      iequation = i - 1 + this%num_grid(1) * (j - 1)
      this%rgt_hand(:, iequation) = this%rgt_hand(:, iequation) - &
                                    truth(x, y) / (this%del_dxdy(1))**2
    END DO
    DO i = 1, this%num_grid(1)
      x = DBLE(i) / DBLE(this%num_grid(1) + 1)
      j = 0
      y = DBLE(j) / DBLE(this%num_grid(2) + 1)
      iequation = i + this%num_grid(1) * (j + 1 - 1)
      this%rgt_hand(:, iequation) = this%rgt_hand(:, iequation) - &
                                    truth(x, y) / (this%del_dxdy(2))**2

      j = this%num_grid(2) + 1
      y = DBLE(j) / DBLE(this%num_grid(2) + 1)
      iequation = i + this%num_grid(1) * (j - 1 - 1)
      this%rgt_hand(:, iequation) = this%rgt_hand(:, iequation) - &
                                    truth(x, y) / (this%del_dxdy(2))**2
    END DO

  END SUBROUTINE Analytics

  SUBROUTINE LoadMatrx(this)
    CLASS(unitPoisson_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i, j, iequation

    ! Default band width: at boundaries, this needs to reassign
    this%num_elem = 0
    DO j = 1, this%num_grid(2)
      DO i = 1, this%num_grid(1)
        iequation = i + this%num_grid(1) * (j - 1)
        this%num_elem(iequation) = this%num_elem(iequation) + 1
        this%idx_colu(this%num_elem(iequation), iequation) = iequation
        this%row_elem(1, iequation) = -2.0D0 / (this%del_dxdy(1))**2
        IF (this%num_grid(2) .GT. 1) &
          this%row_elem(1, iequation) = &
          this%row_elem(1, iequation) - 2.0D0 / (this%del_dxdy(2))**2

        IF (i .GT. 1) THEN
          this%num_elem(iequation) = this%num_elem(iequation) + 1
          this%idx_colu(this%num_elem(iequation), iequation) = iequation - 1
          this%row_elem(this%num_elem(iequation), iequation) = 1.0D0 / (this%del_dxdy(1))**2
        END IF
        IF (i .LT. this%num_grid(1)) THEN
          this%num_elem(iequation) = this%num_elem(iequation) + 1
          this%idx_colu(this%num_elem(iequation), iequation) = iequation + 1
          this%row_elem(this%num_elem(iequation), iequation) = 1.0D0 / (this%del_dxdy(1))**2
        END IF

        IF (j .GT. 1) THEN
          this%num_elem(iequation) = this%num_elem(iequation) + 1
          this%idx_colu(this%num_elem(iequation), iequation) = i + this%num_grid(1) * (j - 2)
          this%row_elem(this%num_elem(iequation), iequation) = 1.0D0 / (this%del_dxdy(2))**2
        END IF
        IF (j .LT. this%num_grid(2)) THEN
          this%num_elem(iequation) = this%num_elem(iequation) + 1
          this%idx_colu(this%num_elem(iequation), iequation) = i + this%num_grid(1) * (j)
          this%row_elem(this%num_elem(iequation), iequation) = 1.0D0 / (this%del_dxdy(2))**2
        END IF
      END DO
    END DO

  END SUBROUTINE LoadMatrx

  SUBROUTINE CalResidu(this)
    CLASS(unitPoisson_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i, ie
    REAL(r_kind) :: rsd_max

    rsd_max = 0.0D0
    ! rhs-diag*sol:
    DO i = 1, this%num_vlvl
      this%res_dual(i, :) = this%rgt_hand(i, :) - &
                            this%row_elem(1, :) * this%sol_curr(i, :)
    END DO
    ! Off diagonal elements:
    DO i = 1, this%num_eqns
      DO ie = 2, this%num_elem(i) ! Off diagonal: starts from 2
        this%res_dual(:, i) = this%res_dual(:, i) - &
                              this%row_elem(ie, i) * this%sol_curr(:, this%idx_colu(ie, i))
      END DO
    END DO
    WRITE (*, 1)
1   FORMAT('--------------------------------------------------')
    WRITE (*, 2) MAXVAL(this%res_dual), MINVAL(this%res_dual)
2   FORMAT('Max/min residuals: ', 2E12.4)
    WRITE (*, 1)
  END SUBROUTINE CalResidu

  SUBROUTINE CalErrors(this, tru, err)
    CLASS(unitPoisson_t) :: this
    REAL(r_kind), INTENT(IN) :: tru(this%num_eqns)
    REAL(r_kind), INTENT(OUT) :: err(2) ! err(1) relative; err(2) maximum

    ! Local variables:
    INTEGER(i_kind) :: i, j, i_err, j_err, iequation
    REAL(r_kind) :: err_l2(this%num_vlvl), &
                    err_lm(this%num_vlvl), fun(this%num_vlvl)

    err_l2 = 0.0D0
    err_lm = 0.0D0
    fun = 0.0D0
    DO j = 1, this%num_grid(2)
      DO i = 1, this%num_grid(1)
        iequation = i + this%num_grid(1) * (j - 1)
        err_l2 = err_l2 + (this%sol_curr(:, iequation) - tru(iequation))**2
        fun = fun + tru(iequation)**2

        IF (err_lm(1) .LT. ABS(this%sol_curr(1, iequation) - tru(iequation))) THEN
          err_lm(1) = ABS(this%sol_curr(1, iequation) - tru(iequation))
          i_err = i
          j_err = j
        END IF
      END DO
    END DO

    err(1) = SQRT(err_l2(1) / fun(1))
    err(2) = err_lm(1)

    WRITE (*, 1) SQRT(err_l2(1) / fun(1))
    WRITE (*, 2) err_lm(1), i_err, j_err
1   FORMAT('Relative error: ', F12.4)
2   FORMAT(' Maximum error: 'F12.4, ' at: ', 2I4)
  END SUBROUTINE CalErrors

END MODULE unitPoisson_m
