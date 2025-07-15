!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by YUANFU XIE  (yuanfu_xie@yahoo.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
! A unit test of regular grid utility:
!
! An interesting note: Usually, when we applied a bilinear interpolation, we always expect
! the interpolation errors are reduced by a factor of 4 when the gridspace is halved. I had
! hard time to pass this unit test as the error reduction is less than a factor of 4 at some
! levels. The reason is that when I used an interpolation of grid points to a set of fixed
! observations, the Taylor expansion ends up with truncation errors multiplied by the
! interpolation coefficients. These coefficients make the error reductions higher than 4 at
! some levels but less than 4 at others. But THE IMPORTANT NOTE IS THAT the overall reductions
! are averaged by a factor of 4.
!
! Therefore, I used a total reduction comparing to 4**nr for this unit test!
PROGRAM Test_RegularGrid
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE parameters_m, ONLY: degree2radian
  USE RectangleGridUtility_m

  IMPLICIT NONE

  TYPE(RectangleGridUtility_t) :: regular

  ! Local variables:
  INTEGER(i_kind) :: numsrc, numtgt, i, j, ng(2), k, nr, iproc
  INTEGER(i_kind), ALLOCATABLE :: idxtgt(:, :)
  REAL(r_kind), PARAMETER :: domain(4) = (/22.0D0, 40.0D0, 105.0D0, 150D0/)
  REAL(r_kind), ALLOCATABLE :: srcLL(:, :), tgtLL(:, :), coetgt(:, :), reg(:, :, :)
  REAL(r_kind), ALLOCATABLE :: srcFun(:), tgtFun(:), err(:, :)

  REAL(r_kind) :: pi, dlat, dlon, anal(2)

  pi = 4.0D0 * ATAN(1.0D0)

  nr = 8
  iproc = 1 ! Single processor test

  anal(1) = 1.0D0
  anal(2) = 1.0D0

  ! Observations:
  !numtgt = 2
  ! Test halos:
  numtgt = 8

  ALLOCATE (err(numtgt, nr))
  err = 0.0D0

  ALLOCATE (tgtll(numtgt, 2), idxtgt(4, numtgt), coetgt(4, numtgt))
  ALLOCATE (tgtFun(numtgt))

  !tgtll(1,1) = 35.0D0*degree2radian
  !tgtll(1,2) = 110.0D0*degree2radian
  !tgtll(2,1) = 24.0D0*degree2radian
  !tgtll(2,2) = 136.2D0*degree2radian
  ! 1. bottom left
  tgtll(1, 1) = domain(1) * degree2radian + 2.0D-4
  tgtll(1, 2) = domain(3) * degree2radian + 0.001
  ! 2. bottom middle:
  tgtll(2, 1) = domain(1) * degree2radian + 2.0D-4
  tgtll(2, 2) = 0.5D0 * (domain(3) + domain(4)) * degree2radian
  ! 3. bottom right:
  tgtll(3, 1) = domain(1) * degree2radian + 2.0D-4
  tgtll(3, 2) = domain(4) * degree2radian - 0.001
  ! 4. left middle:
  tgtll(4, 1) = 0.5D0 * (domain(1) + domain(2)) * degree2radian
  tgtll(4, 2) = domain(3) * degree2radian + 0.001
  ! 5. right middle:
  tgtll(5, 1) = 0.5D0 * (domain(1) + domain(2)) * degree2radian
  tgtll(5, 2) = domain(4) * degree2radian - 0.001
  ! 6. top left:
  tgtll(6, 1) = domain(2) * degree2radian - 2.0D-4
  tgtll(6, 2) = domain(3) * degree2radian + 0.001
  ! 7. top middle:
  tgtll(7, 1) = domain(2) * degree2radian - 0.001D0
  tgtll(7, 2) = 0.5D0 * (domain(3) + domain(4)) * degree2radian
  ! 8. Top right:
  tgtll(8, 1) = domain(2) * degree2radian - 2.0D-4
  tgtll(8, 2) = domain(4) * degree2radian - 0.001

  DO k = 1, nr ! nr resolutions to test
    ! ALlocate memory:
    ng(1) = 5 * 2**(k - 1) + 1
    ng(2) = 5 * 2**(k - 1) + 1

    numsrc = ng(1) * ng(2)

    ! Allocate memory:
    ALLOCATE (srcLL(numsrc, 2), reg(ng(1), ng(2), 2))
    ALLOCATE (srcFun(numsrc))

    ! Lat-Lon grid:
    dlat = (domain(2) - domain(1)) / (DBLE(ng(2) - 1))
    dlon = (domain(4) - domain(3)) / (DBLE(ng(1) - 1))

    srcLL = 0.0D0
    IF (.FALSE.) THEN ! No halo points
      DO j = 1, ng(2)
      DO i = 1, ng(1)
        reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
        reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)

        srcll(i + ng(1) * (j - 1), 1) = reg(i, j, 1) * degree2radian
        srcll(i + ng(1) * (j - 1), 2) = reg(i, j, 2) * degree2radian

        ! Analytic src:
        srcFun(i + ng(1) * (j - 1)) = anal(1) * srcll(i + ng(1) * (j - 1), 1)**2 + &
                                      anal(2) * srcll(i + ng(1) * (j - 1), 2)**2
      END DO
      END DO
    ELSE    ! Halo around all edges
      ! Interior points:
      DO j = 2, ng(2) - 1
        DO i = 2, ng(1) - 1
          reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
          reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)

          srcll(i - 1 + (ng(1) - 2) * (j - 2), 1) = reg(i, j, 1) * degree2radian
          srcll(i - 1 + (ng(1) - 2) * (j - 2), 2) = reg(i, j, 2) * degree2radian

          ! Analytic src:
          srcFun(i - 1 + (ng(1) - 2) * (j - 2)) = anal(1) * srcll(i - 1 + (ng(1) - 2) * (j - 2), 1)**2 + &
                                                  anal(2) * srcll(i - 1 + (ng(1) - 2) * (j - 2), 2)**2
        END DO
      END DO
      ! Bottom side:
      j = 1
      DO i = 1, ng(1)
        reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
        reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)
        srcll((ng(1) - 2) * (ng(2) - 2) + i, 1) = reg(i, j, 1) * degree2radian
        srcll((ng(1) - 2) * (ng(2) - 2) + i, 2) = reg(i, j, 2) * degree2radian

        ! Analytic src:
        srcFun((ng(1) - 2) * (ng(2) - 2) + i) = anal(1) * srcll((ng(1) - 2) * (ng(2) - 2) + i, 1)**2 + &
                                                anal(2) * srcll((ng(1) - 2) * (ng(2) - 2) + i, 2)**2
      END DO
      ! Right side:
      i = ng(1)
      DO j = 2, ng(2) - 1
        reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
        reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)
        srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + j - 1, 1) = reg(i, j, 1) * degree2radian
        srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + j - 1, 2) = reg(i, j, 2) * degree2radian

        ! Analytic src:
        srcFun((ng(1) - 2) * (ng(2) - 2) + ng(1) + j - 1) = anal(1) * srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + j - 1, 1)**2 + &
                                                            anal(2) * srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + j - 1, 2)**2
      END DO
      ! Top side:
      j = ng(2)
      DO i = 1, ng(1)
        reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
        reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)
        srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + ng(2) - 2 + i, 1) = reg(i, j, 1) * degree2radian
        srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + ng(2) - 2 + i, 2) = reg(i, j, 2) * degree2radian

        ! Analytic src:
        srcFun((ng(1) - 2) * (ng(2) - 2) + ng(1) + ng(2) - 2 + i) = &
          anal(1) * srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + ng(2) - 2 + i, 1)**2 + &
          anal(2) * srcll((ng(1) - 2) * (ng(2) - 2) + ng(1) + ng(2) - 2 + i, 2)**2
      END DO
      ! Left side:
      i = 1
      DO j = 2, ng(2) - 1
        reg(i, j, 1) = domain(1) + dlat * DBLE(j - 1)
        reg(i, j, 2) = domain(3) + dlon * DBLE(i - 1)
        srcll((ng(1) - 2) * (ng(2) - 2) + 2 * ng(1) + ng(2) - 2 + j - 1, 1) = reg(i, j, 1) * degree2radian
        srcll((ng(1) - 2) * (ng(2) - 2) + 2 * ng(1) + ng(2) - 2 + j - 1, 2) = reg(i, j, 2) * degree2radian

        ! Analytic src:
        srcFun((ng(1) - 2) * (ng(2) - 2) + 2 * ng(1) + ng(2) - 2 + j - 1) = &
          anal(1) * srcll((ng(1) - 2) * (ng(2) - 2) + 2 * ng(1) + ng(2) - 2 + j - 1, 1)**2 + &
          anal(2) * srcll((ng(1) - 2) * (ng(2) - 2) + 2 * ng(1) + ng(2) - 2 + j - 1, 2)**2
      END DO
    END IF

    ! Interpolation:
    CALL regular%RectangleHorizontalIntp(srcll, numsrc, tgtll, numtgt, idxtgt, coetgt, iproc)

    ! Calculate errors:
    DO i = 1, numtgt
      ! Analytic at obs locations:
      tgtFun(i) = anal(1) * tgtll(i, 1)**2 + anal(2) * tgtll(i, 2)**2
      err(i, k) = 0.0D0
      IF (MINVAL(idxtgt(:, i)) .GT. 0) THEN
        DO j = 1, 4
          err(i, k) = err(i, k) + srcFun(idxtgt(j, i)) * coetgt(j, i)
        END DO
        err(i, k) = err(i, k) - tgtFun(i)
      END IF
    END DO

    IF (k .GT. 1) WRITE (*, 2) k, MAXVAL(ABS(err(:, k))), MAXVAL(err(:, k - 1)) / MAXVAL(err(:, k))
2   FORMAT('Maximum error at level: ', I1, ' is: ', D12.4, ' Ratio: ', D12.4)

    DEALLOCATE (srcll, reg)
    DEALLOCATE (srcFun)
  END DO

  ! Unit test result:
  IF (ABS(MAXVAL(err(:, 1)) / MAXVAL(err(:, nr))) .GE. 4.0D0 * (nr - 1) - 1.0D-3) THEN
    WRITE (*, 3) ABS(MAXVAL(err(:, 1)) / MAXVAL(err(:, nr))), 4.0D0**(nr - 1)
3   FORMAT('Test passed ', 2D14.7)
  ELSE
    WRITE (*, 4) ABS(MAXVAL(err(:, 1)) / MAXVAL(err(:, nr))), 4.0D0**(nr - 1) + 1.0D-3
4   FORMAT('Test failed', 2D14.7)
  END IF

  DEALLOCATE (err, tgtLL, idxtgt, coetgt, tgtFun)

END PROGRAM Test_RegularGrid
