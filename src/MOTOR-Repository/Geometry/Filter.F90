!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE Filter_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double, r_single

  IMPLICIT NONE

CONTAINS
  ! Yuanfu Xie added smootherAvg: using distance weighted average to smooth:
  SUBROUTINE smoothField(sg, vlvl, tslot, numSmooth, toSmooth, smoothed)
    IMPLICIT NONE
    TYPE(SingleGrid_t) sg
    INTEGER(i_kind), INTENT(IN) :: vlvl, tslot, numSmooth ! Number of smooths
    REAL(r_kind), INTENT(IN) :: toSmooth(vlvl, sg%num_cell, tslot)
    REAL(r_kind), INTENT(OUT) :: smoothed(vlvl, sg%num_cell, tslot)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, nstcl
    REAL(r_kind) :: totalDist, relaxation, sigma
    REAL(r_kind), ALLOCATABLE :: copy(:, :, :), weight(:), SAVE(:, :, :)

    ! Gaussian convolution relaxation parameter:
    relaxation = 0.05D0

    nstcl = UBOUND(sg%cell_stcl, 1)
    ALLOCATE (copy(vlvl, sg%num_cell, tslot), weight(nstcl))

    ALLOCATE (SAVE(vlvl, sg%num_cell, tslot))

    ! Repeat smooths:
    copy = toSmooth
    DO k = 1, numSmooth
      DO j = 1, sg%num_icell

        smoothed(:, j, :) = 0.0D0
        totalDist = 0.0D0
        weight = 0.0D0
        sigma = DSQRT(2.0D0) * MAXVAL(sg%cell_dist(:, j)) ! for 2*sigma**2
        DO i = 1, nstcl
          IF (sg%cell_stcl(i, j) .GT. 0 .AND. &
              sg%cell_stcl(i, j) .NE. j) THEN
            ! Use of interior or halo points only:
            IF (sg%cell_type(sg%cell_stcl(i, j)) .EQ. 0 .OR. &
                sg%cell_type(sg%cell_stcl(i, j)) .EQ. 3) THEN
              weight(i) = EXP(-(sg%cell_dist(i, j) / sigma)**2)
              smoothed(:, j, :) = smoothed(:, j, :) + copy(:, sg%cell_stcl(i, j), :) * weight(i)
              totalDist = totalDist + weight(i)

            END IF
          END IF
        END DO
        IF (sg%cell_type(j) .EQ. 0) THEN
          smoothed(:, j, :) = copy(:, j, :) * relaxation + (1.0D0 - relaxation) * smoothed(:, j, :) / totalDist
        ELSE
          smoothed(:, j, :) = smoothed(:, j, :) / totalDist ! Boundary points
        END IF
      END DO

      ! Exchange the filtered halo values:
      DO i = 1, tslot
        CALL sg%ExchangeMatOnHalo2D(vlvl, smoothed(:, :, i))
      END DO
      IF (k .LT. numSmooth) copy = smoothed
    END DO

    DEALLOCATE (copy, weight)
  END SUBROUTINE smoothField

  SUBROUTINE guidedfilter(I, p, r, eps, q)
    INTEGER :: hei, wid
    INTEGER(i_kind), INTENT(IN) :: r
    REAL(r_kind), INTENT(IN) :: I(:, :), p(:, :), eps
    REAL(r_kind), INTENT(OUT) :: q(SIZE(I, 1), SIZE(I, 2))

    REAL(r_kind) :: One(SIZE(I, 1), SIZE(I, 2)), &
                    NMat(SIZE(I, 1), SIZE(I, 2)), &
                    mean_I(SIZE(I, 1), SIZE(I, 2)), &
                    mean_p(SIZE(I, 1), SIZE(I, 2)), &
                    mean_Ip(SIZE(I, 1), SIZE(I, 2)), &
                    cov_Ip(SIZE(I, 1), SIZE(I, 2)), &
                    ! mean_II(size(I, 1), size(I, 2)), &
                    var_I(SIZE(I, 1), SIZE(I, 2)), &
                    a(SIZE(I, 1), SIZE(I, 2)), &
                    b(SIZE(I, 1), SIZE(I, 2)), &
                    mean_a(SIZE(I, 1), SIZE(I, 2)), &
                    mean_b(SIZE(I, 1), SIZE(I, 2))

    hei = SIZE(I, 1)
    wid = SIZE(I, 2)

    One = 1.0D0
    NMat = boxfilter(One, r)
    ! PRINT*, 'NMat:', NMat

    ! PRINT*, 'NMat:', MAXVAL(NMat), MINVAL(NMat)

    mean_I = boxfilter(I, r) / NMat
    PRINT *, 'mean_I:', MAXVAL(mean_I), MINVAL(mean_I)

    mean_p = boxfilter(p, r) / NMat
    mean_Ip = boxfilter(I * p, r) / NMat
    cov_Ip = mean_Ip - mean_I * mean_p
    PRINT *, 'mean_I:', MAXVAL(mean_I), MINVAL(mean_I)
    PRINT *, 'mean_p:', MAXVAL(mean_p), MINVAL(mean_p)
    PRINT *, 'mean_Ip:', MAXVAL(mean_Ip), MINVAL(mean_Ip)
    PRINT *, 'mean_I*mean_p:', MAXVAL(mean_I * mean_p), MINVAL(mean_I * mean_p)

    PRINT *, 'cov_Ip:', MAXVAL(cov_Ip), MINVAL(cov_Ip)

    ! mean_II = boxfilter(I*I, r)/N
    var_I = boxfilter(I * I, r) / NMat - mean_I * mean_I
    PRINT *, 'var_I:', MAXVAL(var_I), MINVAL(var_I)

    a = cov_Ip / (var_I + eps)
    b = mean_p - a * mean_I
    PRINT *, 'a:', MAXVAL(a), MINVAL(a)
    PRINT *, 'b:', MAXVAL(b), MINVAL(b)

    mean_a = boxfilter(a, r) / NMat
    mean_b = boxfilter(b, r) / NMat
    PRINT *, 'mean_a:', MAXVAL(boxfilter(a, r)), MINVAL(boxfilter(a, r))
    PRINT *, 'mean_b:', MAXVAL(boxfilter(b, r)), MINVAL(boxfilter(b, r))
    PRINT *, 'NMat:', MAXVAL(NMat), MINVAL(NMat)

    q = mean_a * I + mean_b
  END SUBROUTINE

  FUNCTION boxfilter(imSrc, r) RESULT(imDst)
    INTEGER(i_kind) :: hei, wid, i, j
    INTEGER(i_kind), INTENT(IN) :: r
    REAL(r_kind), INTENT(IN) :: imSrc(:, :)
    REAL(r_kind) :: imDst(SIZE(imSrc, 1), SIZE(imSrc, 2))
    REAL(r_kind) :: imCum(SIZE(imSrc, 1), SIZE(imSrc, 2))

    hei = SIZE(imSrc, 1)
    wid = SIZE(imSrc, 2)

    imDst = 0.0D0

    ! cumulative sum over X axis
    imCum = imSrc
    DO i = 2, wid
      imCum(:, i) = imCum(:, i) + imCum(:, i - 1)
    END DO

    ! difference over X axis
    imDst(:, 1:r + 1) = imCum(:, 1 + r:2 * r + 1); 
    imDst(:, r + 2:wid - r) = imCum(:, 2 * r + 2:wid) - imCum(:, 1:wid - 2 * r - 1); 
    DO j = wid - r + 1, wid
      imDst(:, j) = imCum(:, wid) - imCum(:, j - r - 1)
    END DO

    ! cumulative sum over Y axis
    imCum = imDst
    DO i = 2, hei
      imCum(i, :) = imCum(i, :) + imCum(i - 1, :)
    END DO

    ! difference over Y axis
    imDst(1:r + 1, :) = imCum(1 + r:2 * r + 1, :)
    imDst(r + 2:hei - r, :) = imCum(2 * r + 2:hei, :) - imCum(1:hei - 2 * r - 1, :); 
    DO i = hei - r + 1, hei
      imDst(i, :) = imCum(hei, :) - imCum(i - r - 1, :)
    END DO

  END FUNCTION

END MODULE
