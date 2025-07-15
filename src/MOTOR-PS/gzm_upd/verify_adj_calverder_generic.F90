PROGRAM verify_adj_calverder_generic
  USE kinds_m, ONLY: i_kind, r_kind
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE CalVerDer_adj_m, ONLY: CalVerDer_adj_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE GenContainers_m, ONLY: GenContainers_t, updateSg, GenGeometry
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(CalVerDer_t) :: calverder
  TYPE(CalVerDer_adj_t) :: calverder_adj
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(GenContainers_t) :: containers
  TYPE(State_t) :: state

  REAL(r_kind), POINTER :: altitude(:), temperature(:), pressure(:), wind_speed(:), humidity(:)
  REAL(r_kind), POINTER :: A(:, :), parA_sigma(:, :), A_adj(:, :), parA_sigma_adj(:, :)
  REAL(r_kind), POINTER :: parA_sigma_half(:, :), parA_sigma_2nd(:, :)
  REAL(r_kind), POINTER :: parA_sigma_half_adj(:, :), parA_sigma_2nd_adj(:, :)
  REAL(r_kind) :: epsilon, max_diff_first, max_diff_second, max_diff_half, numeric_grad
  REAL(r_kind) :: mean_error_first, std_dev_error_first, mean_relative_error_first, std_dev_relative_error_first
  REAL(r_kind) :: mean_error_second, std_dev_error_second, mean_relative_error_second, std_dev_relative_error_second
  REAL(r_kind) :: mean_error_half, std_dev_error_half, mean_relative_error_half, std_dev_relative_error_half
  INTEGER(i_kind) :: num_cell, vLevel, i, k, blockSize, startCell, endCell, num_levels, total_points
  CHARACTER(LEN=100) :: file_name
  INTEGER :: unit_num
  REAL(r_kind) :: tolerance
  LOGICAL :: pass_first, pass_second, pass_half
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Set parameters
  epsilon = 1.0_R_KIND * 1.0E-6
  tolerance = 1.0_R_KIND * 1.0E-5
  blockSize = 10  ! Set block size for each processing block

  file_name = "adjoint_calverder_verification_results.txt"
  unit_num = 10

  ! Open file for writing results
  OPEN (UNIT=unit_num, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  containers = GenContainers_t(configFile)
  CALL containers%GenGeometry(geometry)

  sg => geometry%mg%sg(6)
  calverder%sg => sg

  CALL containers%updateSg(sg, sg)

  num_cell = sg%num_cell
  vLevel = sg%vLevel
  total_points = num_cell * vLevel
  num_levels = vLevel  ! Number of vertical levels

  ! Initialize pass flags
  pass_first = .TRUE.
  pass_second = .TRUE.
  pass_half = .TRUE.

  ! Allocate arrays for atmospheric profiles
  ALLOCATE (altitude(num_levels), temperature(num_levels), pressure(num_levels))
  ALLOCATE (wind_speed(num_levels), humidity(num_levels))

  ! Generate atmospheric profiles
  CALL generate_atmospheric_profiles(altitude, temperature, pressure, wind_speed, humidity)

  ! 第一阶验证
  ALLOCATE (A(vLevel, num_cell), parA_sigma(vLevel, num_cell), A_adj(vLevel, num_cell), parA_sigma_adj(vLevel, num_cell))
  max_diff_first = 0.0_R_KIND
  mean_error_first = 0.0_R_KIND
  std_dev_error_first = 0.0_R_KIND
  mean_relative_error_first = 0.0_R_KIND
  std_dev_relative_error_first = 0.0_R_KIND

  DO startCell = 1, num_cell, blockSize
    endCell = MIN(startCell + blockSize - 1, num_cell)
    DO i = startCell, endCell
      DO k = 1, vLevel
        CALL generate_test_data(A, temperature, pressure, wind_speed, humidity)
        CALL calverder%FirstOrder(A, parA_sigma)
        CALL calverder_adj%FirstOrder_AD(A_adj, parA_sigma_adj)
        CALL compute_numeric_gradient(calverder, A, parA_sigma, epsilon, k, i, numeric_grad, "FirstOrder")
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "FirstOrder", k, i, max_diff_first, mean_error_first, std_dev_error_first, &
                            mean_relative_error_first, std_dev_relative_error_first, pass_first, unit_num)
      END DO
    END DO
  END DO

  ! 打印并释放第一阶内存
  WRITE (unit_num, *) "Final Metrics for First Order:"
  CALL print_final_metrics("FirstOrder", max_diff_first, mean_error_first, std_dev_error_first, mean_relative_error_first, std_dev_relative_error_first, pass_first, tolerance, unit_num)
  DEALLOCATE (A, parA_sigma, A_adj, parA_sigma_adj)

  ! 第二阶验证
  ALLOCATE (A(vLevel, num_cell), parA_sigma_2nd(vLevel, num_cell), A_adj(vLevel, num_cell), parA_sigma_2nd_adj(vLevel, num_cell))
  max_diff_second = 0.0_R_KIND
  mean_error_second = 0.0_R_KIND
  std_dev_error_second = 0.0_R_KIND
  mean_relative_error_second = 0.0_R_KIND
  std_dev_relative_error_second = 0.0_R_KIND

  DO startCell = 1, num_cell, blockSize
    endCell = MIN(startCell + blockSize - 1, num_cell)
    DO i = startCell, endCell
      DO k = 1, vLevel
        CALL generate_test_data(A, temperature, pressure, wind_speed, humidity)
        CALL calverder%SecondOrder(A, parA_sigma_2nd)
        CALL calverder_adj%SecondOrder_AD(A_adj, parA_sigma_2nd_adj)
        CALL compute_numeric_gradient(calverder, A, parA_sigma_2nd, epsilon, k, i, numeric_grad, "SecondOrder")
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "SecondOrder", k, i, max_diff_second, mean_error_second, std_dev_error_second, &
                            mean_relative_error_second, std_dev_relative_error_second, pass_second, unit_num)
      END DO
    END DO
  END DO

  ! 打印并释放第二阶内存
  WRITE (unit_num, *) "Final Metrics for Second Order:"
  CALL print_final_metrics("SecondOrder", max_diff_second, mean_error_second, std_dev_error_second, mean_relative_error_second, std_dev_relative_error_second, pass_second, tolerance, unit_num)
  DEALLOCATE (A, parA_sigma_2nd, A_adj, parA_sigma_2nd_adj)

  ! 半阶验证
  ALLOCATE (A(vLevel, num_cell), parA_sigma_half(vLevel, num_cell), A_adj(vLevel, num_cell), parA_sigma_half_adj(vLevel, num_cell))
  max_diff_half = 0.0_R_KIND
  mean_error_half = 0.0_R_KIND
  std_dev_error_half = 0.0_R_KIND
  mean_relative_error_half = 0.0_R_KIND
  std_dev_relative_error_half = 0.0_R_KIND

  DO startCell = 1, num_cell, blockSize
    endCell = MIN(startCell + blockSize - 1, num_cell)
    DO i = startCell, endCell
      DO k = 1, vLevel
        CALL generate_test_data(A, temperature, pressure, wind_speed, humidity)
        CALL calverder%FirstOrder_half(A, parA_sigma_half)
        CALL calverder_adj%FirstOrder_half_AD(A_adj, parA_sigma_half_adj)
        CALL compute_numeric_gradient(calverder, A, parA_sigma_half, epsilon, k, i, numeric_grad, "FirstOrderHalf")
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "FirstOrderHalf", k, i, max_diff_half, mean_error_half, std_dev_error_half, &
                            mean_relative_error_half, std_dev_relative_error_half, pass_half, unit_num)
      END DO
    END DO
  END DO

  ! 打印并释放半阶内存
  WRITE (unit_num, *) "Final Metrics for First Order Half:"
  CALL print_final_metrics("FirstOrderHalf", max_diff_half, mean_error_half, std_dev_error_half, mean_relative_error_half, std_dev_relative_error_half, pass_half, tolerance, unit_num)
  DEALLOCATE (A, parA_sigma_half, A_adj, parA_sigma_half_adj)

  CLOSE (unit_num)

  ! Deallocate atmospheric profiles
  DEALLOCATE (altitude, temperature, pressure, wind_speed, humidity)

CONTAINS

  ! Subroutine to generate atmospheric profiles
  SUBROUTINE generate_atmospheric_profiles(altitude, temperature, pressure, wind_speed, humidity)
    REAL(r_kind), INTENT(OUT) :: altitude(:), temperature(:), pressure(:), wind_speed(:), humidity(:)
    INTEGER :: i, num_levels
    REAL(r_kind) :: lapse_rate, sea_level_temp, sea_level_pressure, temp, pres

    num_levels = SIZE(altitude, 1)
    lapse_rate = 6.5_R_KIND / 1000.0_R_KIND  ! Lapse rate in K/m
    sea_level_temp = 288.15_R_KIND  ! Sea level standard temperature in K
    sea_level_pressure = 101325.0_R_KIND  ! Sea level standard pressure in Pa

    DO i = 1, num_levels
      altitude(i) = REAL(i - 1, r_kind) * 200.0_R_KIND  ! Example: 0 to 20,000 meters
      temp = sea_level_temp - lapse_rate * altitude(i)
      pres = sea_level_pressure * (temp / sea_level_temp)**5.256_R_KIND

      temperature(i) = temp
      pressure(i) = pres
      wind_speed(i) = 5.0_R_KIND + 10.0_R_KIND * EXP(-((altitude(i) - 10000.0_R_KIND)**2) / (2000.0_R_KIND**2))  ! Jet stream example
      humidity(i) = 0.8_R_KIND * EXP(-altitude(i) / 2000.0_R_KIND)  ! Exponential decrease with altitude
    END DO
  END SUBROUTINE generate_atmospheric_profiles

  ! Subroutine to generate test data based on profiles
  SUBROUTINE generate_test_data(A, temperature, pressure, wind_speed, humidity)
    REAL(r_kind), INTENT(OUT) :: A(:, :)
    REAL(r_kind), INTENT(IN) :: temperature(:), pressure(:), wind_speed(:), humidity(:)
    INTEGER :: i, j, vLevel, num_cell

    vLevel = SIZE(A, 1)
    num_cell = SIZE(A, 2)

    DO i = 1, vLevel
      DO j = 1, num_cell
        A(i, j) = temperature(i) + pressure(i) * 1.0E-5_R_KIND + wind_speed(i) * 0.1_R_KIND + humidity(i) * 10.0_R_KIND
      END DO
    END DO
  END SUBROUTINE generate_test_data

  ! Subroutine to compute the numeric gradient
  SUBROUTINE compute_numeric_gradient(calverder, A, parA_sigma, epsilon, k, i, numeric_grad, order)
    TYPE(CalVerDer_t), INTENT(IN) :: calverder
    REAL(r_kind), INTENT(INOUT) :: A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    INTEGER, INTENT(IN) :: k, i
    REAL(r_kind), INTENT(OUT) :: numeric_grad
    CHARACTER(LEN=*), INTENT(IN) :: order

    REAL(r_kind) :: parA_sigma_plus, parA_sigma_minus

    IF (TRIM(order) == "FirstOrder") THEN
      A(k, i) = A(k, i) + epsilon
      CALL calverder%FirstOrder(A, parA_sigma)
      parA_sigma_plus = SUM(parA_sigma)

      A(k, i) = A(k, i) - 2.0_R_KIND * epsilon
      CALL calverder%FirstOrder(A, parA_sigma)
      parA_sigma_minus = SUM(parA_sigma)

    ELSE IF (TRIM(order) == "SecondOrder") THEN
      A(k, i) = A(k, i) + epsilon
      CALL calverder%SecondOrder(A, parA_sigma)
      parA_sigma_plus = SUM(parA_sigma)

      A(k, i) = A(k, i) - 2.0_R_KIND * epsilon
      CALL calverder%SecondOrder(A, parA_sigma)
      parA_sigma_minus = SUM(parA_sigma)

    ELSE IF (TRIM(order) == "FirstOrderHalf") THEN
      A(k, i) = A(k, i) + epsilon
      CALL calverder%FirstOrder_half(A, parA_sigma)
      parA_sigma_plus = SUM(parA_sigma)

      A(k, i) = A(k, i) - 2.0_R_KIND * epsilon
      CALL calverder%FirstOrder_half(A, parA_sigma)
      parA_sigma_minus = SUM(parA_sigma)
    END IF

    A(k, i) = A(k, i) + epsilon

    ! Compute numeric gradient using central difference
    numeric_grad = (parA_sigma_plus - parA_sigma_minus) / (2.0_R_KIND * epsilon)
  END SUBROUTINE compute_numeric_gradient

  ! Subroutine to verify adjoint gradients
  SUBROUTINE verify_adjoint(adjoint_grad, numeric_grad, order, k, i, max_diff, mean_error, std_dev_error, &
                            mean_relative_error, std_dev_relative_error, PASS, unit_num)
    REAL(r_kind), INTENT(IN) :: adjoint_grad, numeric_grad
    CHARACTER(LEN=*), INTENT(IN) :: order
    INTEGER, INTENT(IN) :: k, i, unit_num
    REAL(r_kind), INTENT(INOUT) :: max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error
    LOGICAL, INTENT(INOUT) :: PASS
    REAL(r_kind) :: diff, rel_error

    diff = ABS(numeric_grad - adjoint_grad)
    IF (diff > max_diff) max_diff = diff

    rel_error = diff / MAX(ABS(numeric_grad), 1.0_R_KIND * 1.0E-8)
    mean_error = mean_error + diff
    mean_relative_error = mean_relative_error + rel_error
    std_dev_error = std_dev_error + diff**2
    std_dev_relative_error = std_dev_relative_error + rel_error**2

    IF (diff > 1.0_R_KIND * 1.0E-5) THEN
      PASS = .FALSE.
    END IF

    WRITE (unit_num, *) order, " Adjoint Verification -- Cell: ", i, " Level: ", k, &
      " Numeric Grad: ", numeric_grad, " Adjoint Grad: ", adjoint_grad, " Diff: ", diff, &
      " Relative Error: ", rel_error
  END SUBROUTINE verify_adjoint

  ! Subroutine to print final metrics and determine if a test passed
  SUBROUTINE print_final_metrics(order, max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, PASS, tolerance, unit_num)
    CHARACTER(LEN=*), INTENT(IN) :: order
    REAL(r_kind), INTENT(IN) :: max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, tolerance
    LOGICAL, INTENT(IN) :: PASS
    INTEGER, INTENT(IN) :: unit_num

    WRITE (unit_num, *) "Metrics for ", order, ":"
    WRITE (unit_num, *) "Maximum difference: ", max_diff
    WRITE (unit_num, *) "Mean error: ", mean_error
    WRITE (unit_num, *) "Standard deviation of error: ", std_dev_error
    WRITE (unit_num, *) "Mean relative error: ", mean_relative_error
    WRITE (unit_num, *) "Standard deviation of relative error: ", std_dev_relative_error

    IF (max_diff > tolerance .OR. mean_error > tolerance .OR. std_dev_error > tolerance .OR. &
        mean_relative_error > tolerance .OR. std_dev_relative_error > tolerance) THEN
      WRITE (unit_num, *) order, " Adjoint Test: FAILED"
    ELSE
      WRITE (unit_num, *) order, " Adjoint Test: PASSED"
    END IF
  END SUBROUTINE print_final_metrics

END PROGRAM verify_adj_calverder_generic
