!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/20, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
PROGRAM Test_QRSolver
  USE kinds_m, ONLY: r_kind, i_kind

  REAL(r_kind) :: sigma(4)
  REAL(r_kind) :: matrix(4, 4), RHS(4), TAU(4, 4), WORK(256), A(4, 4)
  INTEGER(i_kind) :: i, j, k, LWORK, INFO
  REAL(r_kind) :: h1, h2, h3, h4, a1, a2, a3, a4

  matrix = 0
  sigma = (/0, 100, 300, 1000/)

  DO i = 2, 3
    IF (i == 2) THEN
      h1 = sigma(i - 1) - sigma(i)
      h2 = 0
      h3 = sigma(i + 1) - sigma(i)
      h4 = sigma(i + 2) - sigma(i)
    END IF

    IF (i == 3) THEN
      h1 = sigma(i - 2) - sigma(i)
      h2 = sigma(i - 1) - sigma(i)
      h3 = 0
      h4 = sigma(i + 1) - sigma(i)
    END IF

    matrix(1, i) = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
    matrix(2, i) = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
    matrix(3, i) = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
    matrix(4, i) = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))
  END DO

  matrix(1, 1) = 1.0D0; 
  matrix(4, 4) = 1.0D0; 
  ! !! cautious: fortran is column-major order
  WRITE (*, '(4F20.15)') TRANSPOSE(matrix)
! 1001 FORMAT

  LWORK = 16
  A = matrix

  ! QR Decomposition
  CALL DGEQRF(4, 4, A, 4, TAU, WORK, LWORK, INFO)

  ! A * X=B
  RHS = 1.0D0
  CALL DORMQR('L', 'T', 4, 1, 4, A, 4, TAU, rhs, 4, WORK, LWORK, INFO)
  CALL DTRTRS('U', 'N', 'N', 4, 1, A, 4, rhs, 4, INFO)
  PRINT *, 'INFO: ', MATMUL(matrix, rhs)

  ! A**T * X=B
  RHS = 1.0D0
  CALL DTRTRS('U', 'T', 'N', 4, 1, A, 4, rhs, 4, INFO)
  CALL DORMQR('L', 'N', 4, 1, 4, A, 4, TAU, rhs, 4, WORK, LWORK, INFO)

  PRINT *, 'INFO: ', MATMUL(TRANSPOSE(matrix), rhs)

END PROGRAM Test_QRSolver
