!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021/5/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!  This program is used for adjoint checking of EnLoc.
!!
PROGRAM Test_EnLoc_ADcheck
  USE kinds_m, ONLY: i_kind, r_kind
  USE BKErr_m, ONLY: BKErr_t

  TYPE(BKErr_t)           :: br
  CHARACTER(LEN=1024)     :: outputDir, datFileName
  CHARACTER(LEN=3)        :: gid
  INTEGER(i_kind)         :: status, nlvl, nhor, ntim, nz, nc, i, j, info

  REAL(r_kind), ALLOCATABLE :: qr_aux(:), qr_a(:, :)
  REAL(r_kind), ALLOCATABLE :: Rup(:, :), Rinv(:, :), Rinvtmp(:, :)
  REAL(r_kind) :: ADcheck_left, ADcheck_right, tb, te
  REAL(r_kind), ALLOCATABLE :: t1(:), t2(:)

  ! ! Here is the code of monk for simulating the fileinput of B mat files
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)
  gid = 'G03'
  datFileName = TRIM(outputDir)//'/BKErr_QRdecomp_temp_'//gid//'.dat'
  PRINT *, 'datFileName:', TRIM(datFileName)
  CALL CPU_TIME(tb)
  CALL br%b_qrDecomp_read(TRIM(datFileName), status)
  CALL CPU_TIME(te)
  PRINT *, "b_qrDecomp_read - Time cost: ", (te - tb), "s"

  nc = br%nc
  nz = br%nz

  ALLOCATE (qr_aux(nz), qr_a(nc, nz))
  qr_aux = br%qr_aux
  qr_a = br%qr_a

  ALLOCATE (Rup(nz, nz))
  Rup = 0.0D0
  DO i = 1, nz
    DO j = i, nz
      Rup(i, j) = qr_a(i, j)
    END DO
  END DO
  ! PRINT *, 'Rup:'
  ! DO i=1,nz
  !   WRITE(*, '(10f18.12)') Rup(i,:)
  ! END DO

  ! calculate inverse of R matrix
  ALLOCATE (Rinvtmp(nz, nz))
  Rinvtmp = Rup
  CALL DTRTRI('U', 'N', nz, Rinvtmp, nz, info)
  IF (info .EQ. 0) THEN
    ALLOCATE (Rinv(nz, nz))
    Rinv = Rinvtmp
  ELSE
    PRINT *, "ERROR! Failed in calculating the inverse of R matirx."
    STOP
  END IF

  ! PRINT *, 'Rinv:'
  ! DO i=1,nz
  !   WRITE(*, '(10f18.12)') Rinv(i,:)
  ! END DO
  DEALLOCATE (Rinvtmp)
  ALLOCATE (t1(nz), t2(nz))

  !! FOR sqrt_inverse_multiply
  BLOCK
    REAL(r_kind), ALLOCATABLE :: datatmp(:, :)
    REAL(r_kind), ALLOCATABLE :: bx(:, :), w(:, :), qtw(:, :), v(:, :)
    REAL(r_kind), ALLOCATABLE :: qtwtmp(:), x(:, :)
    INTEGER(i_kind) :: n, jj
    REAL(r_kind) :: tmp
    ! REAL(r_kind) :: tb, te, t2b, t2e

    ALLOCATE (w(nc, 1), qtw(nc, 1), v(nc, 1))
    ! ALLOCATE(qtwtmp(nc))
    ! w   = 1.0D0
    qtw = 0.0D0
    v = 0.0D0
    DO i = 1, nc
      w(i, 1) = i * 1.0D0
    END DO
    ! qtwtmp = 0.0D0

    ! CALL cpu_time(tb)
    DO n = 1, nz
      CALL CPU_TIME(tb)
      v(n, 1) = qr_aux(n)
      v(n + 1:nc, 1) = qr_a(n + 1:nc, n)
      IF (n .GT. 1) v(1:n - 1, 1) = 0.0D0
      ! PRINT *, 'sum(v):', SUM(v)
      ! PRINT *, 'qr_aux:', v(n,1)
      ! IF (n .EQ. nz) THEN
      !   CALL cpu_time(te)
      !   PRINT *, 'sqrt_inverse_multiply - Time cost in generating v vector:', (te-tb), 's'
      ! END IF
      CALL CPU_TIME(te)
      t1(n) = te - tb

      IF (n .EQ. 1) qtw = w

      ! x(i) = DOT_PRODUCT(tmp(i,:), bv)
      ! IF (n .EQ. 1) CALL cpu_time(t2b)
      CALL CPU_TIME(tb)
      ! qtwtmp = 0.0D0
      tmp = DOT_PRODUCT(v(:, 1), qtw(:, 1))
      DO i = 1, nc
        ! qtwtmp(i) = qtw(i,1) - v(i,1) * tmp / v(n,1)
        qtw(i, 1) = qtw(i, 1) - v(i, 1) * tmp / v(n, 1)
      END DO
      ! qtw(:,1) = qtwtmp
      ! PRINT *, 'sum(qtw):', SUM(qtw)
      ! IF (n .EQ. nz) THEN
      !   CALL cpu_time(t2e)
      !   PRINT *, 'sqrt_inverse_multiply - Time cost in calculating qtw:', (t2e-t2b), 's'
      ! END IF
      ! qtw = qtw - MATMUL(MATMUL(v, TRANSPOSE(v)), qtw) / v(i,1)
      CALL CPU_TIME(te)
      t2(n) = te - tb

    END DO
    ! CALL cpu_time(te)
    ! PRINT *, 'sqrt_inverse_multiply - Time cost:', (te-tb), 's'

    CALL CPU_TIME(tb)
    ALLOCATE (datatmp(nz, 1))
    datatmp = 0.0D0
    datatmp = MATMUL(Rinv, qtw(1:nz, :))
    CALL CPU_TIME(te)

    PRINT *, 'sqrt_inverse_multiply - Time cost in generating v vector:', SUM(t1), 's'
    PRINT *, 'sqrt_inverse_multiply - Time cost in calculating QTx:', SUM(t2), 's'
    PRINT *, 'sqrt_inverse_multiply - Time cost in calculating QnTx at n-th loop:'
    DO i = 1, nz
      PRINT *, t2(i), 's'
    END DO
    PRINT *, 'sqrt_inverse_multiply - Time cost in calculating RinvQTx:', (te - tb), 's'

    ALLOCATE (x(nz, 1))
    DO i = 1, nz
      x(i, 1) = i * 1.0D0
    END DO
    ! x = 1.0D0
    ADcheck_left = DOT_PRODUCT(datatmp(:, 1), x(:, 1))

    DEALLOCATE (datatmp, w, qtw, v, x)

    ! PRINT *, "PJM DEBUG --- Done sqrt_inverse_multiply"
  END BLOCK
  !! END sqrt_inverse_multiply

  !! FOR sqrt_inverse_multiply_adjoint
  BLOCK
    REAL(r_kind), ALLOCATABLE :: bv(:, :), v(:, :), y(:, :)
    REAL(r_kind), ALLOCATABLE :: RinvTx(:, :), QnRinvTx(:, :), QnRinvTxtmp(:)
    INTEGER(i_kind) :: n, mused
    REAL(r_kind) :: tmp
    ! REAL(r_kind) :: tb, te, t2b, t2e

    ! PRINT *, "PJM DEBUG --- Begin sqrt_inverse_multiply_adjoint"
    ALLOCATE (bv(nz, 1), RinvTx(nz, 1), v(nc, 1))
    ALLOCATE (QnRinvTx(nc, 1))
    ! ALLOCATE(QnRinvTxtmp(nc))
    DO i = 1, nz
      bv(i, 1) = i * 1.0D0
    END DO
    ! bv = 1.0D0

    RinvTx = 0.0D0
    RinvTx = MATMUL(TRANSPOSE(Rinv), bv)
    ! PRINT *, 'RinvTx:'
    ! PRINT *, RinvTx
    v = 0.0D0
    QnRinvTx = 0.0D0
    QnRinvTx(1:nz, 1) = RinvTx(:, 1)

    ! CALL cpu_time(tb)
    DO j = 1, nz
      ! generated v vector, v = [zeros(n-1,1), qraux(n), a(n+1:nc, n)]
      ! IF (j .EQ. 1) CALL cpu_time(tb)
      n = nz - j + 1
      CALL CPU_TIME(tb)
      v(n, 1) = qr_aux(n)
      v(n + 1:nc, 1) = qr_a(n + 1:nc, n)
      ! IF (j .EQ. nz) THEN
      !   CALL cpu_time(te)
      !   PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in generating v vector:', (te-tb), 's'
      ! END IF

      IF (n .GT. 1) v(1:n - 1, 1) = 0.0D0
      CALL CPU_TIME(te)
      t1(j) = te - tb

      ! w_0 = (inv(R))' * X, [n n]*[n 1] = [n 1]
      ! when i=1, w_1 = Q_n * w_0 = w_0 - v_n * v_n' * w_0 ./ tau,
      !           dimension is [m n]*[n n]*[n 1] = [m 1];
      ! when i>1, w_i = Q_n-i+1 * w_i = w_i - v_n-i+1 * v_n-i+1' * w_i ./ tau,
      !           dimension is [m m]*[m m]*[m 1] = [m 1].
      IF (j .EQ. 1) THEN
        mused = nz
      ELSE
        mused = nc
      END IF
      ! IF (j .EQ. 1) CALL cpu_time(t2b)
      CALL CPU_TIME(tb)
      tmp = DOT_PRODUCT(v(1:mused, 1), QnRinvTx(1:mused, 1))
      DO i = 1, nc
        ! QnRinvTxtmp(i) = QnRinvTx(i,1) - v(i,1) * tmp / v(n,1)
        QnRinvTx(i, 1) = QnRinvTx(i, 1) - v(i, 1) * tmp / v(n, 1)
      END DO
      ! QnRinvTx(:,1) = QnRinvTxtmp(:)
      ! IF (j .EQ. nz) THEN
      !   CALL cpu_time(t2e)
      !   PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in calculating QnRinvTx:', (t2e-t2b), 's'
      ! END IF
      CALL CPU_TIME(te)
      t2(j) = te - tb
    END DO
    ! CALL cpu_time(te)
    ! PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost:', (te-tb), 's'

    PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in generating v vector:', SUM(t1), 's'
    PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in calculating QRinvTx:', SUM(t2), 's'
    PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in calculating QnRinvTx at n-th loop:'
    DO i = 1, nz
      PRINT *, t2(i), 's'
    END DO
    ! PRINT *, 'sqrt_inverse_multiply_adjoint - Time cost in MATMUL at last step:', (te-tb), 's'

    ALLOCATE (y(nc, 1))
    ! y = 1.0D0
    DO i = 1, nc
      y(i, 1) = i * 1.0D0
    END DO
    ADcheck_right = DOT_PRODUCT(QnRinvTx(:, 1), y(:, 1))

    DEALLOCATE (bv, v, RinvTx, QnRinvTx, y)
    ! DEALLOCATE(QnRinvTxtmp)

    ! PRINT *, "PJM DEBUG --- DONE sqrt_inverse_multiply_adjoint"
  END BLOCK
  !! END sqrt_inverse_multiply_adjoint

  PRINT *, 'ADcheck_left:', ADcheck_left
  PRINT *, 'ADcheck_right:', ADcheck_right

END PROGRAM Test_EnLoc_ADcheck
