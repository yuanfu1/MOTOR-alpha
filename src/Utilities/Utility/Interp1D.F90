!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.utility.vertical interpolation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2021/10/26, @GBA-MWF, Shenzhen
!   Modified by Jilong Chen (jchen@link.cuhk.edu.hk) for linear extrapolation, 2022/1/18
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2022/3/7, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
MODULE Interp1D_m

  USE kinds_m, ONLY: i_kind, r_kind, r_double, r_single

  INTERFACE interp1d
    MODULE PROCEDURE interp1d_1D
    MODULE PROCEDURE interp1d_2D
    MODULE PROCEDURE interp1d_3D_idx3
    MODULE PROCEDURE interp1d_3D_idx3_single
    MODULE PROCEDURE interp1d_3D_idx2
    MODULE PROCEDURE interp1d_3D_idx2_single

  END INTERFACE

CONTAINS

  SUBROUTINE interp1d_1D(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_kind), INTENT(IN) :: sigmain(:), xin(:)
    REAL(r_kind), INTENT(IN) :: sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: xout(:)
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: xint

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    ALLOCATE (xint(1:xin_size))

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = dlog(xin)
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

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_1D

  SUBROUTINE interp1d_2D(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_kind), INTENT(IN) :: sigmain(:), xin(:, :)
    REAL(r_kind), INTENT(IN) :: sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: xout(:, :)
    REAL(r_kind), ALLOCATABLE :: xint(:, :)

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    xint = xin
    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = dlog(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(i, :) = xint(j, :)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(i, :) = xint(j - 1, :) + &
                         (xint(j, :) - xint(j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                         (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(i, :) = xint(1, :) + &
                         (xint(2, :) - xint(1, :)) / (sigmain(2) - sigmain(1)) * &
                         (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size) THEN
          xout(i, :) = xint(j - 1, :) + &
                       (xint(j, :) - xint(j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                       (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_2D

  SUBROUTINE interp1d_3D_idx3(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_kind), INTENT(IN) :: sigmain(:), xin(:, :, :)
    REAL(r_kind), INTENT(IN) :: sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: xout(:, :, :)
    REAL(r_kind), ALLOCATABLE :: xint(:, :, :)

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    xint = xin
    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = dlog(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(:, :, i) = xint(:, :, j)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(:, :, i) = xint(:, :, j - 1) + &
                            (xint(:, :, j) - xint(:, :, j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                            (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(:, :, i) = xint(:, :, 1) + &
                            (xint(:, :, 2) - xint(:, :, 1)) / (sigmain(2) - sigmain(1)) * &
                            (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size .AND. (SIZE(sigmain) .GT. 1)) THEN
          xout(:, :, i) = xint(:, :, j - 1) + &
                          (xint(:, :, j) - xint(:, :, j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                          (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_3D_idx3

  SUBROUTINE interp1d_3D_idx3_single(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_kind), INTENT(IN) :: sigmain(:)
    REAL(r_single), INTENT(IN) :: xin(:, :, :)
    REAL(r_kind), INTENT(IN) :: sigmaout(:)
    REAL(r_single), INTENT(OUT) :: xout(:, :, :)
    REAL(r_single), ALLOCATABLE :: xint(:, :, :)

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    xint = xin
    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = LOG(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(:, :, i) = xint(:, :, j)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(:, :, i) = xint(:, :, j - 1) + &
                            (xint(:, :, j) - xint(:, :, j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                            (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(:, :, i) = xint(:, :, 1) + &
                            (xint(:, :, 2) - xint(:, :, 1)) / (sigmain(2) - sigmain(1)) * &
                            (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size .AND. (SIZE(sigmain) .GT. 1)) THEN
          xout(:, :, i) = xint(:, :, j - 1) + &
                          (xint(:, :, j) - xint(:, :, j - 1)) / (sigmain(j) - sigmain(j - 1)) * &
                          (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_3D_idx3_single

  SUBROUTINE interp1d_3D_idx2_single(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_single), INTENT(IN) :: sigmain(:)
    REAL(r_single), INTENT(IN) :: xin(:, :, :)
    REAL(r_single), INTENT(IN) :: sigmaout(:)
    REAL(r_single), INTENT(OUT) :: xout(:, :, :)
    REAL(r_single), ALLOCATABLE :: xint(:, :, :)

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    xint = xin
    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = LOG(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(:, i, :) = xint(:, j, :)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(:, i, :) = xint(:, j - 1, :) + &
                            (xint(:, j, :) - xint(:, j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                            (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(:, i, :) = xint(:, 1, :) + &
                            (xint(:, 2, :) - xint(:, 1, :)) / (sigmain(2) - sigmain(1)) * &
                            (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size .AND. (SIZE(sigmain) .GT. 1)) THEN
          xout(:, i, :) = xint(:, j - 1, :) + &
                          (xint(:, j, :) - xint(:, j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                          (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_3D_idx2_single

  SUBROUTINE interp1d_3D_idx2(sigmain, xin, sigmaout, xout, log_flag)
    IMPLICIT NONE

    INTEGER(i_kind) :: i, j, xin_size, xout_size
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: log_flag
    REAL(r_single), INTENT(IN) :: sigmain(:)
    REAL(r_kind), INTENT(IN) :: xin(:, :, :)
    REAL(r_single), INTENT(IN) :: sigmaout(:)
    REAL(r_kind), INTENT(OUT) :: xout(:, :, :)
    REAL(r_kind), ALLOCATABLE :: xint(:, :, :)

    xout_size = SIZE(sigmaout)
    xin_size = SIZE(sigmain)

    xint = xin
    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xint = dlog(xin)
    ELSE
      xint = xin
    END IF

    DO i = 1, xout_size
      DO j = 1, xin_size
        IF (sigmain(j) == sigmaout(i)) THEN
          xout(:, i, :) = xint(:, j, :)
          EXIT
        ELSE IF (sigmain(j) > sigmaout(i) .AND. j <= xin_size) THEN
          IF (j > 1) THEN
            xout(:, i, :) = xint(:, j - 1, :) + &
                            (xint(:, j, :) - xint(:, j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                            (sigmaout(i) - sigmain(j - 1))
            EXIT
          ELSE
            xout(:, i, :) = xint(:, 1, :) + &
                            (xint(:, 2, :) - xint(:, 1, :)) / (sigmain(2) - sigmain(1)) * &
                            (sigmaout(i) - sigmain(1))
            EXIT
          END IF
        ELSE IF (sigmain(j) < sigmaout(i) .AND. j == xin_size .AND. (SIZE(sigmain) .GT. 1)) THEN
          xout(:, i, :) = xint(:, j - 1, :) + &
                          (xint(:, j, :) - xint(:, j - 1, :)) / (sigmain(j) - sigmain(j - 1)) * &
                          (sigmaout(i) - sigmain(j - 1))
        END IF
      END DO
    END DO

    IF (PRESENT(log_flag)) THEN
      IF (log_flag == "LOG") xout = EXP(xout)
    END IF
    IF (ALLOCATED(xint)) DEALLOCATE (xint)
  END SUBROUTINE interp1d_3D_idx2

END MODULE Interp1D_m
