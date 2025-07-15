!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE Interp1d2_m

CONTAINS

  SUBROUTINE interp1d(sigmain, xin, xin_size, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: xin_size
    INTEGER(4) :: i, j, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(8), INTENT(IN) :: sigmain(:), xin(:)
    REAL(8), INTENT(IN) :: sigmaout(:)
    REAL(8), INTENT(OUT) :: xout(:)
    REAL(8) :: xint(xin_size)

    xout_size = SIZE(sigmaout)

    IF (PRESENT(log_flag) .AND. log_flag == "LOG") THEN
      xint = dlog(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(i) = xint(j)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(i) = xint(j - 1) + &
                      (xint(j) - xint(j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                      (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(i) = xint(1) + &
                      (xint(2) - xint(1)) / (sigmain(2) - sigmain(1)) * &
                      (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size) THEN
          xout(i) = xint(j - 1) + &
                    (xint(j) - xint(j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                    (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag) .AND. log_flag == "LOG") xout = EXP(xout)
  END SUBROUTINE interp1d

END MODULE Interp1d2_m
