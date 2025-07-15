PROGRAM verify_adj_calverder_firstorder
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
  REAL(r_kind) :: epsilon, diff, max_diff, numeric_grad
  REAL(r_kind) :: mean_error, std_dev_error, mean_relative_error, std_dev_relative_error
  REAL(r_kind) :: sum_error, sum_relative_error, sum_square_error, sum_square_relative_error
  INTEGER(i_kind) :: num_cell, vLevel, i, k, num_levels, n, count_pass
  CHARACTER(LEN=100) :: file_name
  INTEGER :: unit_num, total_points
  REAL(r_kind) :: tolerance, tolerance_rate, diff_sum
  LOGICAL :: global_pass
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Set parameters
  epsilon = 1.0_R_KIND * 1.0E-6
  tolerance = 1.0_R_KIND * 1.0E-5
  tolerance_rate = 0.9995_R_KIND  ! 99.95% pass rate

  file_name = "adjoint_calverder_firstorder_results.txt"
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
  num_levels = vLevel  ! Number of vertical levels
  total_points = num_cell * vLevel

  ! Allocate arrays for atmospheric profiles
  ALLOCATE (altitude(num_levels), temperature(num_levels), pressure(num_levels))
  ALLOCATE (wind_speed(num_levels), humidity(num_levels))

  ! Allocate arrays for adjoint verification
  ALLOCATE (A(vLevel, num_cell), parA_sigma(vLevel, num_cell))
  ALLOCATE (A_adj(vLevel, num_cell), parA_sigma_adj(vLevel, num_cell))

  ! Generate atmospheric profiles
  CALL generate_atmospheric_profiles(altitude, temperature, pressure, wind_speed, humidity)

  ! Initialize metrics
  max_diff = 0.0_R_KIND
  sum_error = 0.0_R_KIND
  sum_square_error = 0.0_R_KIND
  sum_relative_error = 0.0_R_KIND
  sum_square_relative_error = 0.0_R_KIND
  count_pass = 0  ! Counter for passing grid points
  n = 0  ! Counter for grid points with significant error

  global_pass = .TRUE.

  A_adj = 0.0_R_KIND
  parA_sigma_adj = 1.0_R_KIND
  ! Generate test data
  CALL generate_test_data(A, temperature, pressure, wind_speed, humidity)

  ! Forward calculation
  CALL calverder%FirstOrder(A, parA_sigma)

  ! Initialize adjoint model
  calverder_adj = CalVerDer_adj_t(sg)

  ! Adjoint calculation (outside the loop)
  CALL calverder_adj%FirstOrder_AD(A_adj, parA_sigma_adj)

  ! Loop over all cells and vertical levels
  DO i = 1, num_cell
    DO k = 1, vLevel

      ! Compute numeric gradient for each grid point
      CALL compute_numeric_gradient(calverder, A, parA_sigma, epsilon, k, i, numeric_grad, "FirstOrder")

      ! Verify adjoint for each grid point
      CALL verify_adjoint(A_adj(k, i), numeric_grad, "FirstOrder", k, i, max_diff, sum_error, sum_square_error, &
                          sum_relative_error, sum_square_relative_error, n, count_pass, global_pass, unit_num)

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
  WRITE (unit_num, *) "Final Metrics for First Order:"
  CALL print_final_metrics("FirstOrder", max_diff, mean_error, std_dev_error, mean_relative_error, std_dev_relative_error, global_pass, tolerance, unit_num)

  CLOSE (unit_num)

  ! Deallocate arrays
  DEALLOCATE (A, parA_sigma, A_adj, parA_sigma_adj)
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
    END IF

    A(k, i) = A(k, i) + epsilon

    ! Compute numeric gradient using central difference
    numeric_grad = (parA_sigma_plus - parA_sigma_minus) / (2.0_R_KIND * epsilon)
  END SUBROUTINE compute_numeric_gradient

  ! Subroutine to verify adjoint gradients
  SUBROUTINE verify_adjoint(adjoint_grad, numeric_grad, order, k, i, max_diff, sum_error, sum_square_error, &
                            sum_relative_error, sum_square_relative_error, n, count_pass, global_pass, unit_num)
    REAL(r_kind), INTENT(IN) :: adjoint_grad, numeric_grad
    CHARACTER(LEN=*), INTENT(IN) :: order
    INTEGER, INTENT(IN) :: k, i, unit_num
    REAL(r_kind), INTENT(INOUT) :: max_diff, sum_error, sum_square_error, sum_relative_error, sum_square_relative_error
    INTEGER, INTENT(INOUT) :: n, count_pass
    LOGICAL, INTENT(INOUT) :: global_pass
    REAL(r_kind) :: diff, rel_error

    diff = ABS(numeric_grad - adjoint_grad)
    IF (diff > max_diff) max_diff = diff

    rel_error = diff / MAX(ABS(numeric_grad), 1.0_R_KIND * 1.0E-8)

    ! Accumulate errors and increment counters for significant errors
    IF (diff > tolerance) THEN
      n = n + 1
      sum_error = sum_error + diff
      sum_square_error = sum_square_error + diff**2
      sum_relative_error = sum_relative_error + rel_error
      sum_square_relative_error = sum_square_relative_error + rel_error**2
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
    WRITE (unit_num, *) "Mean relative error: ", mean_relative_error
    WRITE (unit_num, *) "Standard deviation of relative error: ", std_dev_relative_error

    IF (global_pass) THEN
      WRITE (unit_num, *) order, " Adjoint Test: PASSED"
    ELSE
      WRITE (unit_num, *) order, " Adjoint Test: FAILED"
    END IF
  END SUBROUTINE print_final_metrics

END PROGRAM verify_adj_calverder_firstorder
