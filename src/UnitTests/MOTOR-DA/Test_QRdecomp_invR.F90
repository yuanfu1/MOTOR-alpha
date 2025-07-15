!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.EnLoc.Test_EnLoc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0   Test of calculating the inverse of R matrix by using LAPACK
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/5/9, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_QRdecomp_invR

  USE kinds_m
  USE BKErr_m, ONLY: BKErr_t

  IMPLICIT NONE

  TYPE(BKErr_t)              ::  BKErr, BKErr_r
  REAL(r_kind)               ::  maxdiff
  REAL(r_kind), ALLOCATABLE  ::  A(:, :), A_recover(:, :)
  REAL(r_kind), ALLOCATABLE  ::  invR(:, :), invR_r(:, :)
  INTEGER(i_kind)            ::  m, n, status, i, info
  CHARACTER(LEN=1024)        ::  outputDir, datFileName

  ! Get the OUTPUT_DIR
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)
  datFileName = TRIM(outputDir)//'/QRdecomp_data.dat'

  ! Initialize A
  m = 10
  n = 5
  ALLOCATE (A(m, n), A_recover(m, n))

  A = RESHAPE((/92.0D-0, 99.0D-0, 1.0D-0, 8.0D-0, 15.0D-0, &
                98.0D-0, 80.0D-0, 7.0D-0, 14.0D-0, 16.0D-0, &
                4.0D-0, 81.0D-0, 88.0D-0, 20.0D-0, 22.0D-0, &
                85.0D-0, 87.0D-0, 19.0D-0, 21.0D-0, 3.0D-0, &
                86.0D-0, 93.0D-0, 25.0D-0, 2.0D-0, 9.0D-0, &
                17.0D-0, 24.0D-0, 76.0D-0, 83.0D-0, 90.0D-0, &
                23.0D-0, 5.0D-0, 82.0D-0, 89.0D-0, 91.0D-0, &
                79.0D-0, 6.0D-0, 13.0D-0, 95.0D-0, 97.0D-0, &
                10.0D-0, 12.0D-0, 94.0D-0, 96.0D-0, 78.0D-0, &
                11.0D-0, 18.0D-0, 100.0D-0, 77.0D-0, 84.0D-0 &
                /), shape=(/m, n/), order=(/2, 1/))

  ! A = Q * R, R^(-1) = inv(R), inv() is one function in MATLAB
  ! invR = inv(R), regarded as the true value, used for verification
  invR = RESHAPE((/-0.005000937763754D0, -0.008027107704547D0, -0.001634728280714D0, -0.026491504638836D0, 0.004628247869395D0, &
                   0.0D0, 0.009445671029026D0, 0.003753111462132D0, 0.026267495348697D0, -0.004291019971609D0, &
                   0.0D0, 0.0D0, -0.005580441398466D0, -0.027010954190419D0, 0.004559497333168D0, &
                   0.0D0, 0.0D0, 0.0D0, 0.027643828005699D0, 0.029717008716859D0, &
                   0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.034752719739647D0 &
                   /), shape=(/n, n/), order=(/2, 1/))

  ! QR Decomposition, return the whole Q and R matrixes
  BKErr%nc = m
  CALL BKErr%b_qrDecomposition(n, A, 0, status)
  PRINT *, "nc, nz:", BKErr%nc, BKErr%nz
  ! PRINT *,"done b_qrDecomposition"
  CALL BKErr%b_qrDecomp_save(datFileName, .TRUE., status)
  ! PRINT *,"done b_qrDecomp_save"
  CALL BKErr%free()

  CALL BKErr_r%b_qrDecomp_read(datFileName, status)
  CALL BKErr_r%b_qrReconstruction(status)

  ! rebuild A by using Q*R
  A_recover = MATMUL(BKErr_r%Q, BKErr_r%R)

  ! WRITE(*,*) "Q:"
  ! DO i=1,m
  !   PRINT "(100f10.4)", BKErr_r%Q(i,:)
  ! END DO

  ! WRITE(*,*) "R:"
  ! DO i=1,m
  !   PRINT "(100f10.4)", BKErr_r%R(i,:)
  ! END DO

  ! print out the max difference between A and A_recover
  PRINT *, "MAXVAL(ABS(A_recover - A_orig)):", MAXVAL(ABS(A_recover - A))

  invR_r = BKErr_r%R(1:n, 1:n)
  ! PRINT *, 'invR_r before DGETRI:'
  ! DO i=1,n
  !   PRINT "(100f10.4)", invR_r(i,:)
  ! END DO

  CALL DTRTRI('U', 'N', n, invR_r, n, info)
  PRINT *, "status of DTRTRI:", info
  ! PRINT *, 'invR_r after DTRTRI:'
  ! DO i=1,n
  !   PRINT "(100f10.4)", invR_r(i,:)
  ! END DO

  ! PRINT *, 'invR:'
  ! DO i=1,n
  !   PRINT "(100f10.4)", invR(i,:)
  ! END DO

  ! print out the max difference between invR and invR_r
  PRINT *, "MAXVAL(ABS(invR_cal - invR_true)):", MAXVAL(ABS(invR_r - invR))

END PROGRAM Test_QRdecomp_invR
