PROGRAM verify_jacobian_adjoint
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t
  USE gzm_adj_m, ONLY: gzm_adj_t
  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  REAL(r_kind), POINTER :: opr_left(:, :), opr_right(:, :), oprand(:, :), adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
  REAL(r_kind), POINTER :: VALUE(:, :), adj_value(:, :)
  REAL(r_kind) :: epsilon, max_diff, sum_error, sum_relative_error, sum_square_error, sum_square_relative_error
  REAL(r_kind) :: mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, diff, numeric_grad
  INTEGER(i_kind) :: num_cell, vLevel, i, k, num_levels, n, count_pass
  CHARACTER(LEN=100) :: file_name, error_file_name
  INTEGER :: unit_num, error_unit_num, total_points
  REAL(r_kind) :: tolerance, tolerance_rate
  LOGICAL :: global_pass
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Set parameters
  epsilon = 1.0_R_KIND * 1.0E-6
  tolerance = 1.0_R_KIND * 1.0E-12
  tolerance_rate = 0.9995_R_KIND  ! 99.95% pass rate

  file_name = "adjoint_jacobian_ad_results.txt"
  error_file_name = "jacobian_large_errors.txt"

  unit_num = 10
  error_unit_num = 20

  ! Open file for writing results
  OPEN (UNIT=unit_num, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')

  ! Open file for logging large errors
  OPEN (UNIT=error_unit_num, FILE=error_file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')

  ! Configure file paths
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testjacobian_tlad.yaml"
  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialization
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)

  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type
  gzm%sg => sg
  gzm_adj%sg => sg

  num_cell = sg%num_icell
  vLevel = sg%vLevel
  num_levels = vLevel  ! Number of vertical levels
  total_points = num_cell * vLevel

  ! Allocate arrays for operator data and adjoint verification
  ALLOCATE (opr_left(vLevel, num_cell), opr_right(vLevel, num_cell), oprand(vLevel, num_cell))
  ALLOCATE (adj_opr_left(vLevel, num_cell), adj_opr_right(vLevel, num_cell), adj_oprand(vLevel, num_cell))
  ALLOCATE (VALUE(vLevel, num_cell), adj_value(vLevel, num_cell))

  ! Generate test data
  CALL generate_test_data(opr_left, opr_right, oprand, vLevel, num_cell)

  ! Initialize metrics
  max_diff = 0.0_R_KIND
  sum_error = 0.0_R_KIND
  sum_square_error = 0.0_R_KIND
  sum_relative_error = 0.0_R_KIND
  sum_square_relative_error = 0.0_R_KIND
  count_pass = 0  ! Counter for passing grid points
  n = 0  ! Counter for grid points with significant error

  global_pass = .TRUE.

  adj_opr_left = 0.0_R_KIND
  adj_opr_right = 0.0_R_KIND
  adj_oprand = 0.0_R_KIND
  adj_value = 1.0_R_KIND  ! Set adjoint value to 1.0

  ! Forward calculation
  CALL gzm%Jacobian(opr_left, opr_right, VALUE)

  ! Adjoint calculation
  CALL gzm_adj%Jacobian_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)

  ! Loop over all cells and vertical levels
  DO i = 1, num_cell
    DO k = 1, vLevel

      ! Compute numeric gradient for each grid point
      CALL compute_numeric_gradient(gzm, opr_left, opr_right, oprand, VALUE, epsilon, k, i, numeric_grad, "Jacobian")

      ! Verify adjoint for each grid point
      CALL verify_adjoint(adj_opr_right(k, i), numeric_grad, "Jacobian", k, i, max_diff, sum_error, sum_square_error, &
                          sum_relative_error, sum_square_relative_error, n, count_pass, global_pass, tolerance, unit_num, error_unit_num)

    END DO
  END DO

  ! Apply tolerance rate check
  IF (REAL(count_pass) / REAL(total_points) < tolerance_rate) THEN
    global_pass = .FALSE.
  END IF

  ! Calculate final metrics
  IF (n > 0) THEN
    mean_error = sum_error / n
    std_dev_error = SQRT(sum_square_error / n - mean_error**2)
    mean_relative_error = sum_relative_error / n
    std_dev_relative_error = SQRT(sum_square_relative_error / n - mean_relative_error**2)
  ELSE
    mean_error = 0.0_R_KIND
    std_dev_error = 0.0_R_KIND
    mean_relative_error = 0.0_R_KIND
    std_dev_relative_error = 0.0_R_KIND
  END IF

  ! Print final metrics
  WRITE (unit_num, *) "Final Metrics for Jacobian:"
  CALL print_final_metrics("Jacobian", max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, global_pass, tolerance, unit_num)

  CLOSE (unit_num)
  CLOSE (error_unit_num)

  ! Deallocate arrays
  DEALLOCATE (opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value)

CONTAINS

  ! Subroutine to generate test data
  SUBROUTINE generate_test_data(opr_left, opr_right, oprand, vLevel, num_icell)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    REAL(r_kind), INTENT(OUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    INTEGER, INTENT(IN) :: vLevel, num_icell
    INTEGER :: i, j

    ! Initialize operator data
    DO i = 1, vLevel
      DO j = 1, num_icell
        opr_left(i, j) = SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                         EXP(-REAL(i) / REAL(vLevel)) + 0.1_R_KIND * REAL(j) / REAL(num_icell)

        opr_right(i, j) = COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                          EXP(-REAL(i) / REAL(vLevel)) + 0.1_R_KIND * REAL(i) / REAL(vLevel)

        oprand(i, j) = SIN(3.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel)) * &
                       COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell))
      END DO
    END DO
  END SUBROUTINE generate_test_data

  ! Subroutine to compute the numeric gradient
  SUBROUTINE compute_numeric_gradient(gzm, opr_left, opr_right, oprand, VALUE, epsilon, k, i, numeric_grad, order)
    USE kinds_m, ONLY: r_kind
    TYPE(gzm_t), INTENT(IN) :: gzm
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: VALUE(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    INTEGER, INTENT(IN) :: k, i
    REAL(r_kind), INTENT(OUT) :: numeric_grad
    CHARACTER(LEN=*), INTENT(IN) :: order

    REAL(r_kind) :: value_plus, value_minus

    IF (TRIM(order) == "Jacobian") THEN
      opr_right(k, i) = opr_right(k, i) + epsilon
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
      value_plus = SUM(VALUE)

      opr_right(k, i) = opr_right(k, i) - 2.0_R_KIND * epsilon
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
      value_minus = SUM(VALUE)

      opr_right(k, i) = opr_right(k, i) + epsilon

    END IF

    ! Compute numeric gradient using central difference
    numeric_grad = (value_plus - value_minus) / (2.0_R_KIND * epsilon)
  END SUBROUTINE compute_numeric_gradient

  ! Subroutine to verify adjoint gradients
  SUBROUTINE verify_adjoint(adjoint_grad, numeric_grad, order, k, i, max_diff, sum_error, sum_square_error, &
                            sum_relative_error, sum_square_relative_error, n, count_pass, global_pass, tolerance, unit_num, error_unit_num)
    USE kinds_m, ONLY: r_kind
    REAL(r_kind), INTENT(IN) :: adjoint_grad, numeric_grad, tolerance
    CHARACTER(LEN=*), INTENT(IN) :: order
    INTEGER, INTENT(IN) :: k, i, unit_num, error_unit_num
    REAL(r_kind), INTENT(INOUT) :: max_diff, sum_error, sum_square_error, sum_relative_error, sum_square_relative_error
    INTEGER, INTENT(INOUT) :: n, count_pass
    LOGICAL, INTENT(INOUT) :: global_pass
    REAL(r_kind) :: diff, rel_error

    diff = ABS(numeric_grad - adjoint_grad)
    IF (diff > max_diff) max_diff = diff

    rel_error = diff / MAX(ABS(numeric_grad), epsilon)

    ! Accumulate errors and increment counters for significant errors
    IF (diff > tolerance) THEN
      n = n + 1
      sum_error = sum_error + diff
      sum_square_error = sum_square_error + diff**2
      sum_relative_error = sum_relative_error + rel_error
      sum_square_relative_error = sum_square_relative_error + rel_error**2
      WRITE (error_unit_num, *) "Significant error -- ", order, " Cell: ", i, " Level: ", k, " Diff: ", diff, " Rel Error: ", rel_error
    ELSE
      count_pass = count_pass + 1
    END IF

    WRITE (unit_num, *) order, " Adjoint Verification -- Cell: ", i, " Level: ", k, &
      " Numeric Grad: ", numeric_grad, " Adjoint Grad: ", adjoint_grad, " Diff: ", diff, &
      " Relative Error: ", rel_error
  END SUBROUTINE verify_adjoint

  ! Subroutine to print final metrics and determine if a test passed
  SUBROUTINE print_final_metrics(order, max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, global_pass, tolerance, unit_num)
    CHARACTER(LEN=*), INTENT(IN) :: order
    REAL(r_kind), INTENT(IN) :: max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, tolerance
    LOGICAL, INTENT(IN) :: global_pass
    INTEGER, INTENT(IN) :: unit_num

    WRITE (unit_num, *) "Metrics for ", order, ":"
    WRITE (unit_num, *) "Maximum difference: ", max_diff
    WRITE (unit_num, *) "Mean error: ", mean_error
    WRITE (unit_num, *) "Standard deviation of error: ", std_dev_error
    !WRITE(unit_num, *) "Mean relative error: ", mean_relative_error
    !WRITE(unit_num, *) "Standard deviation of relative error: ", std_dev_relative_error

    IF (global_pass .OR. max_diff < tolerance .OR. mean_error < tolerance .OR. std_dev_error < tolerance) THEN
      WRITE (unit_num, *) order, " Adjoint Test: PASSED"
    ELSE
      WRITE (unit_num, *) order, " Adjoint Test: FAILED"
    END IF
  END SUBROUTINE print_final_metrics

END PROGRAM verify_jacobian_adjoint
