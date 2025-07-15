PROGRAM Simple_Test_GZM_TLM_ADJ
  USE kinds_m, ONLY: r_kind
  USE Verification_m
  IMPLICIT NONE

  REAL(r_kind), ALLOCATABLE :: x(:, :), dx(:, :), y(:, :), tlm_output(:, :), ad_output(:, :)
  REAL(r_kind) :: epsilon
  INTEGER :: n, m, i, j

  n = 10
  m = 10
  epsilon = 1.0E-3_R_KIND

  ALLOCATE (x(n, m), dx(n, m), y(n, m), tlm_output(n, m), ad_output(n, m))

  ! 初始化数据
  DO i = 1, n
    DO j = 1, m
      x(i, j) = REAL((i - 1) * m + j, r_kind)
      dx(i, j) = epsilon
    END DO
  END DO

  ! 简单的线性操作作为 TLM 和 AD
  y = x
  tlm_output = x + dx
  ad_output = x - dx

  ! 调试输出
  PRINT *, "Initial values of x: "
  DO i = 1, n
    PRINT *, x(i, :)
  END DO

  PRINT *, "Initial values of dx: "
  DO i = 1, n
    PRINT *, dx(i, :)
  END DO

  PRINT *, "Computed tlm_output: "
  DO i = 1, n
    PRINT *, tlm_output(i, :)
  END DO

  PRINT *, "Computed ad_output: "
  DO i = 1, n
    PRINT *, ad_output(i, :)
  END DO

  ! 验证 TLM 和 AD
  CALL Verify_TLM_AD(x, dx, y, tlm_output, ad_output, "Simple Test")

  ! 清理
  DEALLOCATE (x, dx, y, tlm_output, ad_output)

END PROGRAM Simple_Test_GZM_TLM_ADJ
