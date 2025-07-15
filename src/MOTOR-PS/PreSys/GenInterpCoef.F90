!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.vertical interpolation coef calculation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2023/06/06, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
MODULE GenInterpCoef_m

  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpCoef_m

  IMPLICIT NONE

  INTERFACE GenInterpCoef
    MODULE PROCEDURE GenInterpCoef_1D
    MODULE PROCEDURE GenInterpCoef_2D
  END INTERFACE

CONTAINS

  SUBROUTINE GenInterpCoef_1D(sigmain, sigmaout, coefout, ks)
    IMPLICIT NONE
    ! TYPE(InterpCoef_t) :: InterpCoef
    REAL(r_kind), INTENT(IN) :: sigmain(:), sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: coefout(:, :)
    INTEGER(i_kind), INTENT(OUT) :: ks(:)
    INTEGER(i_kind) :: i, j, xin_size, xout_size
    REAL(r_kind) :: d1, d2, d3

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) .GT. sigmaout(i)) THEN
          IF (j .EQ. 1) THEN
            CALL InterpCoef(sigmain(1:3), sigmaout(i), coefout(i, :))
            ks(i) = 1
            EXIT
          ELSEIF (j .GT. 1 .AND. j .LT. xin_size) THEN
            CALL InterpCoef(sigmain(j - 1:j + 1), sigmaout(i), coefout(i, :))
            ks(i) = j - 1
            EXIT
          ELSE
            CALL InterpCoef(sigmain(j - 2:j), sigmaout(i), coefout(i, :))
            ks(i) = j - 2
          END IF
        ELSEIF (sigmain(j) < sigmaout(i) .AND. j == xin_size) THEN
          CALL InterpCoef(sigmain(j - 2:j), sigmaout(i), coefout(i, :))
          ks(i) = j - 2
          EXIT
        END IF
      END DO
    END DO

  END SUBROUTINE GenInterpCoef_1D

  SUBROUTINE GenInterpCoef_2D(sigmain, sigmaout, coefout, ks)
    IMPLICIT NONE
    ! TYPE(InterpCoefMinNorm_t) :: InterpCoef
    REAL(r_kind), INTENT(IN) :: sigmain(:), sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: coefout(:, :, :)
    INTEGER(i_kind), INTENT(OUT) :: ks(:, :)
    INTEGER(i_kind) :: i, j, hi, xin_size, xout_size, h_size

    xout_size = SIZE(sigmaout, 1)
    xin_size = SIZE(sigmain, 1)
    h_size = SIZE(ks, 2)

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) .GT. sigmaout(i)) THEN
          IF (j .LE. 2) THEN
            DO hi = 1, h_size
              CALL InterpCoef(sigmain(1:4), sigmaout(i), coefout(i, hi, :))
              ks(i, hi) = 1
              IF (ks(i, hi) .EQ. 0) THEN
                PRINT *, 'stop 1 and j=', j
                STOP
              END IF
            END DO
            EXIT
          ELSEIF (j .GT. 2 .AND. j .LT. xin_size) THEN
            DO hi = 1, h_size
              CALL InterpCoef(sigmain(j - 2:j + 1), sigmaout(i), coefout(i, hi, :))
              ks(i, hi) = j - 2
              IF (ks(i, hi) .EQ. 0) THEN
                PRINT *, 'stop 2 and j=', j
                STOP
              END IF
            END DO
            EXIT
          ELSE
            DO hi = 1, h_size
              CALL InterpCoef(sigmain(j - 3:j), sigmaout(i), coefout(i, hi, :))
              ks(i, hi) = j - 3
              IF (ks(i, hi) .EQ. 0) THEN
                PRINT *, 'stop 3 and j=', j
                STOP
              END IF
            END DO
          END IF
        ELSEIF (j == xin_size) THEN
          DO hi = 1, h_size
            CALL InterpCoef(sigmain(j - 3:j), sigmaout(i), coefout(i, hi, :))
            ks(i, hi) = j - 3
            IF (ks(i, hi) .EQ. 0) THEN
              PRINT *, 'stop 4 and j=', j
              STOP
            END IF
          END DO
          EXIT
        END IF
      END DO
    END DO

  END SUBROUTINE GenInterpCoef_2D

END MODULE GenInterpCoef_m
