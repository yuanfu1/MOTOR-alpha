PROGRAM verify_adj_calverder_without_tl
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

  REAL(r_kind), POINTER :: A(:, :), parA_sigma(:, :), A_adj(:, :), parA_sigma_adj(:, :)
  REAL(r_kind), POINTER :: parA_sigma_half(:, :), parA_sigma_2nd(:, :)
  REAL(r_kind), POINTER :: parA_sigma_half_adj(:, :), parA_sigma_2nd_adj(:, :)
  REAL(r_kind) :: epsilon, diff, max_diff, numeric_grad, adjoint_grad
  INTEGER(i_kind) :: num_cell, vLevel, i, k, blockSize, startCell, endCell
  CHARACTER(LEN=100) :: file_name
  INTEGER :: unit_num
  REAL(r_kind) :: tolerance
  LOGICAL :: PASS
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Set parameters
  epsilon = 1.0_R_KIND * 1.0E-6
  tolerance = 1.0_R_KIND * 1.0E-5
  PASS = .TRUE.
  blockSize = 10  ! 设置每次处理的单元数块大小

  file_name = "adjoint_calverder_verification_results.txt"
  unit_num = 10

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

  ! Allocate arrays
  ALLOCATE (A(vLevel, num_cell), parA_sigma(vLevel, num_cell))
  ALLOCATE (A_adj(vLevel, num_cell), parA_sigma_adj(vLevel, num_cell))
  ALLOCATE (parA_sigma_half(vLevel, num_cell), parA_sigma_2nd(vLevel, num_cell))
  ALLOCATE (parA_sigma_half_adj(vLevel, num_cell), parA_sigma_2nd_adj(vLevel, num_cell))

  ! Generate test data
  CALL generate_test_data(A)
  A_adj = 0.0_R_KIND
  parA_sigma_adj = 1.0_R_KIND
  parA_sigma_half_adj = 1.0_R_KIND
  parA_sigma_2nd_adj = 1.0_R_KIND

  ! Perform forward mode calculations
  CALL calverder%FirstOrder(A, parA_sigma)
  CALL calverder%SecondOrder(A, parA_sigma_2nd)
  CALL calverder%FirstOrder_half(A, parA_sigma_half)

  ! Initialize adjoint model
  calverder_adj = CalVerDer_adj_t(sg)
  CALL calverder_adj%FirstOrder_AD(A_adj, parA_sigma_adj)
  CALL calverder_adj%SecondOrder_AD(A_adj, parA_sigma_2nd_adj)
  CALL calverder_adj%FirstOrder_half_AD(A_adj, parA_sigma_half_adj)

  max_diff = 0.0_R_KIND

  WRITE (unit_num, *) "Adjoint Verification:"

  ! 分块处理循环
  DO startCell = 1, num_cell, blockSize
    endCell = MIN(startCell + blockSize - 1, num_cell)
    DO i = startCell, endCell
      DO k = 1, vLevel
        ! 检验第一阶伴随
        CALL compute_numeric_gradient(calverder, A, parA_sigma, epsilon, k, i, numeric_grad)
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "FirstOrder", k, i, max_diff, PASS, unit_num)

        ! 检验二阶伴随
        CALL compute_numeric_gradient(calverder, A, parA_sigma_2nd, epsilon, k, i, numeric_grad)
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "SecondOrder", k, i, max_diff, PASS, unit_num)

        ! 检验半阶伴随
        CALL compute_numeric_gradient(calverder, A, parA_sigma_half, epsilon, k, i, numeric_grad)
        CALL verify_adjoint(A_adj(k, i), numeric_grad, "FirstOrderHalf", k, i, max_diff, PASS, unit_num)
      END DO
    END DO
  END DO

  WRITE (unit_num, *) "Maximum difference: ", max_diff
  IF (max_diff > tolerance) THEN
    WRITE (unit_num, *) "Adjoint Test: FAILED"
  ELSE
    WRITE (unit_num, *) "Adjoint Test: PASSED"
  END IF
  WRITE (unit_num, *)

  CLOSE (unit_num)
  DEALLOCATE (A, parA_sigma, A_adj, parA_sigma_adj)
  DEALLOCATE (parA_sigma_half, parA_sigma_2nd, parA_sigma_half_adj, parA_sigma_2nd_adj)

CONTAINS

  ! Subroutine to generate test data
  SUBROUTINE generate_test_data(A)
    REAL(r_kind), INTENT(OUT) :: A(:, :)
    INTEGER :: i, j, vLevel, num_cell

    vLevel = SIZE(A, 1)
    num_cell = SIZE(A, 2)

    DO i = 1, vLevel
      DO j = 1, num_cell
        A(i, j) = SIN(2.0_R_KIND * 3.14159_R_KIND * REAL(i) / REAL(vLevel)) * &
                  COS(2.0_R_KIND * 3.14159_R_KIND * REAL(j) / REAL(num_cell))
      END DO
    END DO
  END SUBROUTINE generate_test_data

  ! Subroutine to compute the numeric gradient
  SUBROUTINE compute_numeric_gradient(calverder, A, parA_sigma, epsilon, k, i, numeric_grad)
    TYPE(CalVerDer_t), INTENT(IN) :: calverder
    REAL(r_kind), INTENT(INOUT) :: A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :)
    REAL(r_kind), INTENT(IN) :: epsilon
    INTEGER, INTENT(IN) :: k, i
    REAL(r_kind), INTENT(OUT) :: numeric_grad

    REAL(r_kind) :: parA_sigma_plus, parA_sigma_minus

    ! Perturb A(k, i) by epsilon
    A(k, i) = A(k, i) + epsilon
    CALL calverder%FirstOrder(A, parA_sigma)
    parA_sigma_plus = SUM(parA_sigma)

    ! Perturb A(k, i) by -epsilon
    A(k, i) = A(k, i) - 2.0_R_KIND * epsilon
    CALL calverder%FirstOrder(A, parA_sigma)
    parA_sigma_minus = SUM(parA_sigma)

    ! Restore A(k, i)
    A(k, i) = A(k, i) + epsilon

    ! Compute numeric gradient using central difference
    numeric_grad = (parA_sigma_plus - parA_sigma_minus) / (2.0_R_KIND * epsilon)
  END SUBROUTINE compute_numeric_gradient

  ! Subroutine to verify adjoint gradients
  SUBROUTINE verify_adjoint(adjoint_grad, numeric_grad, order, k, i, max_diff, PASS, unit_num)
    REAL(r_kind), INTENT(IN) :: adjoint_grad, numeric_grad
    CHARACTER(LEN=*), INTENT(IN) :: order
    INTEGER, INTENT(IN) :: k, i, unit_num
    REAL(r_kind), INTENT(INOUT) :: max_diff
    LOGICAL, INTENT(INOUT) :: PASS
    REAL(r_kind) :: diff

    diff = ABS(numeric_grad - adjoint_grad)
    IF (diff > max_diff) max_diff = diff
    IF (diff > 1.0_R_KIND * 1.0E-5) PASS = .FALSE.

    WRITE (unit_num, *) order, " Adjoint Verification -- Cell: ", i, " Level: ", k, &
      " Numeric Grad: ", numeric_grad, " Adjoint Grad: ", adjoint_grad, " Diff: ", diff
  END SUBROUTINE verify_adjoint

END PROGRAM verify_adj_calverder_without_tl
