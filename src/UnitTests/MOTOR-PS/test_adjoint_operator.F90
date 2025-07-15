PROGRAM test_adjoint_operator
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_tlm, Jacobian_tlm, Laplacia_tlm
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES

  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :), pert_oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), adj_oprand(:, :), d_oprand(:, :)

  REAL(r_kind) :: inner_product_tl, inner_product_ad, error_norm
  REAL(r_kind) :: tolerance, epsilon
  CHARACTER(LEN=1024) :: configFile, ncOutputFile
  CHARACTER(LEN=256) :: VERIFY_AD_result_file, ADJ_VALUE_CORRELATION_AD
  INTEGER(i_kind) :: i, j, k_test, l
  LOGICAL :: ieee_support
  REAL(r_kind) :: norm_diff, norm_tlm, norm_x, norm_dx
  REAL(r_kind), ALLOCATABLE :: lambda_values(:)
  INTEGER :: num_lambdas

  REAL(r_kind) :: lambda, ratio, rel_error, abs_error, dot_product, angle
  REAL(r_kind) :: f_base, f_perturb
  REAL(r_kind), ALLOCATABLE :: x1(:, :), x2(:, :), x3(:, :), dx1(:, :), dx2(:, :), dx3(:, :)
  REAL(r_kind), ALLOCATABLE :: M_x(:, :), M_x_perturb(:, :), M_tlm(:, :)

  REAL(r_kind), ALLOCATABLE :: random_perturbation(:, :)

  INTEGER(i_kind) :: vLevel, num_icell

  ! configure file path
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)
  ! associate SingleGrid
  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  gzm%sg => sg
  gzm_tlm%sg => sg
  gzm_adj%sg => sg

  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))
  ALLOCATE (random_perturbation(SIZE(opr_left, 1), SIZE(opr_left, 2)))
  ! Step 1: Generate initial test data
  CALL generate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, sg%vLevel, sg%num_icell)

! ! case 2
!   CALL random_number(opr_left)
!   CALL random_number(opr_right)
!   CALL random_number(oprand)
!   CALL random_number(d_opr_left)
!   CALL random_number(d_opr_right)
!   CALL random_number(d_oprand)
!   CALL random_number(adj_value)
!   CALL random_number(random_perturbation)

!   opr_left = opr_left * 10.0_r_kind
!   opr_right = opr_right * 10.0_r_kind
!   oprand = oprand * 12.6_r_kind
!   d_opr_left = opr_left + random_perturbation * 0.1_r_kind
!   d_opr_right = opr_right + random_perturbation * 0.1_r_kind
!   d_oprand = d_oprand + random_perturbation * 0.1_r_kind
! ! case 3
!    DO i = 1, sg%vLevel
!     d_opr_left(i,:) = (DSIN(sg%cell_cntr(1,:))+DCOS(sg%cell_cntr(2,:)))*10.0D0
!     d_opr_right(i,:) = (DCOS(sg%cell_cntr(1,:))+DSIN(sg%cell_cntr(2,:)))*10.0D0
!     d_oprand(i,:) = DSIN(sg%cell_cntr(1,:)) *1.0D1
!     adj_opr_left(i,:) = 0.0D0
!     adj_opr_right(i,:) = 0.0D0
!     adj_oprand(i,:) = 0.0D0
!   END DO

!   epsilon = 10.0_r_kind

!   opr_left = d_opr_left + 1.56E3
!   opr_right = d_opr_right + 1.66E3 - 10000D0
!   oprand = d_oprand + 1.03E5 + 80000D0

  VERIFY_AD_result_file = "verifying_ajoint.txt"
  ADJ_VALUE_CORRELATION_AD = "verify_adjoint_adj_value.txt"
  ! Step 2: Run the adjoint operator test for Divergen
  PRINT *, 'Testing Divergen_AD...'
  CALL test_adjoint_operator_for_subroutine('Divergen', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, d_value, adj_value, VERIFY_AD_result_file)

  ! Step 3: Run the adjoint operator test for Jacobian
!    PRINT *, 'Testing Jacobian_AD...'
!    CALL test_adjoint_operator_for_subroutine('Jacobian', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, value, d_value, adj_value, VERIFY_AD_result_file)

!    ! Step 4: Run the adjoint operator test for Laplacia
!    PRINT *, 'Testing Laplacia_AD...'
!    CALL test_adjoint_operator_for_subroutine('Laplacia', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, value, d_value, adj_value, VERIFY_AD_result_file)
! ! !!!!! for all kinds of adj_value.

!   PRINT *, 'Testing Divergen_AD...'
!    CALL test_adjoint_operator_alongwith_adj_value('Divergen', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, value, d_value, adj_value, ADJ_VALUE_CORRELATION_AD)

!    ! Step 3: Run the adjoint operator test for Jacobian
!    PRINT *, 'Testing Jacobian_AD...'
!    CALL test_adjoint_operator_alongwith_adj_value('Jacobian', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, value, d_value, adj_value, ADJ_VALUE_CORRELATION_AD)

!    ! Step 4: Run the adjoint operator test for Laplacia
!    PRINT *, 'Testing Laplacia_AD...'
!    CALL test_adjoint_operator_alongwith_adj_value('Laplacia', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, value, d_value, adj_value, ADJ_VALUE_CORRELATION_AD)

  ! Deallocate all arrays
  DEALLOCATE (opr_left, opr_right, oprand, adj_value, d_opr_left, d_opr_right, d_oprand)

CONTAINS
! ---------------------------------------------------------------------------
  SUBROUTINE test_adjoint_operator_for_subroutine(sub_name, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, d_value, adj_value, file_name)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :), d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
    REAL(r_kind), INTENT(OUT)  :: VALUE(:, :), d_value(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_value(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    REAL(r_kind) :: inner_prod_tl, inner_prod_ad
    INTEGER :: vLevel, num_icell, unit, k, i

    REAL(r_kind) :: max_d_value, min_d_value, avg_d_value
    REAL(r_kind) :: max_adj_value, min_adj_value, avg_adj_value

    ! Get the size of the input arrays
    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(opr_left, 2)

    ! Zero the adjoint variables
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    adj_oprand = 0.0_R_KIND

    ! Step 1: Run tangent-linear model
    IF (sub_name == 'Divergen') THEN
      CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (sub_name == 'Jacobian') THEN
      CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (sub_name == 'Laplacia') THEN
      CALL gzm_tlm%Laplacia_TLM(oprand, d_oprand, d_value)
    END IF

    ! Calculate statistics for d_value and adj_value
    max_d_value = MAXVAL(d_value)
    min_d_value = MINVAL(d_value)
    avg_d_value = SUM(d_value) / (sg%vLevel * sg%num_icell)

    max_adj_value = MAXVAL(adj_value)
    min_adj_value = MINVAL(adj_value)
    avg_adj_value = SUM(adj_value) / (vLevel * num_icell)

! Print the statistics
    PRINT *, 'Statistics for d_value: max = ', max_d_value, ', min = ', min_d_value, ', avg = ', avg_d_value
    PRINT *, 'Statistics for adj_value: max = ', max_adj_value, ', min = ', min_adj_value, ', avg = ', avg_adj_value
    adj_value = 100.0_R_KIND
! Calculate the inner product
    inner_prod_tl = 0.0_R_KIND
    DO k = 1, vLevel
      DO i = 1, num_icell
        inner_prod_tl = inner_prod_tl + d_value(k, i) * adj_value(k, i)
      END DO
    END DO
    PRINT *, 'Inner Product TL = ', inner_prod_tl
!stop 1000

    ! Initialize adj_value properly (e.g., random or based on another method)
    !CALL random_number(adj_value)
    !adj_value = 1.0_r_kind
    ! Step 2: Run adjoint model
    IF (sub_name == 'Divergen') THEN
      CALL gzm_adj%Divergen_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    ELSE IF (sub_name == 'Jacobian') THEN
      CALL gzm_adj%Jacobian_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    ELSE IF (sub_name == 'Laplacia') THEN
      CALL gzm_adj%Laplacia_AD(oprand, adj_oprand, adj_value)
    END IF

    inner_prod_tl = 0.0_R_KIND
    inner_prod_ad = 0.0_R_KIND

    IF (sub_name == 'Divergen' .OR. sub_name == 'Jacobian') THEN
      DO k = 1, sg%vLevel
        DO i = 1, sg%num_icell
          !PRINT *, 'd_value(', k, ',', i, ') = ', d_value(k, i)
          !PRINT *, 'adj_value(', k, ',', i, ') = ', adj_value(k, i)
          inner_prod_tl = inner_prod_tl + d_value(k, i) * adj_value(k, i)
          !PRINT *, 'Current inner_prod_tl = ', inner_prod_tl
          inner_prod_ad = inner_prod_ad + d_opr_left(k, i) * adj_opr_left(k, i) + &
                          d_opr_right(k, i) * adj_opr_right(k, i)
          !qPRINT *, 'Current inner_prod_ad = ', inner_prod_ad
        END DO
      END DO
    ELSE IF (sub_name == 'Laplacia') THEN
      DO k = 1, sg%vLevel
        DO i = 1, sg%num_icell
          inner_prod_tl = inner_prod_tl + d_value(k, i) * adj_value(k, i)
          inner_prod_ad = inner_prod_ad + d_oprand(k, i) * adj_oprand(k, i)
        END DO
      END DO
    END IF

    ! Step 5: Compare inner products
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) sub_name//' adjoint'
    WRITE (unit, *) 'Inner Product TL:', inner_prod_tl
    WRITE (unit, *) 'Inner Product AD:', inner_prod_ad
    IF (ABS(inner_prod_tl - inner_prod_ad) < 1.0E-14_R_KIND) THEN
      WRITE (unit, *) 'The adjoint operator test has passed for subroutine ', sub_name
    ELSE
      WRITE (unit, *) 'WARNING: The adjoint operator test has failed for subroutine ', sub_name
    END IF
    CLOSE (unit)
  END SUBROUTINE test_adjoint_operator_for_subroutine

  SUBROUTINE test_adjoint_operator_alongwith_adj_value(sub_name, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, d_value, adj_value, file_name)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :), d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
    REAL(r_kind), INTENT(OUT)  :: VALUE(:, :), d_value(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_value(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: vLevel, num_icell, unit
    REAL(r_kind) :: inner_product_tl, inner_product_adj
    REAL(r_kind) :: max_val

    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(opr_left, 2)

    ! Zero the adjoint variables
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    adj_oprand = 0.0_R_KIND

    ! Run tangent-linear model
    IF (sub_name == 'Divergen') THEN
      CALL gzm%Divergen(opr_left, opr_right, VALUE)
      CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (sub_name == 'Jacobian') THEN
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
      CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (sub_name == 'Laplacia') THEN
      CALL gzm%Laplacia(oprand, VALUE)
      CALL gzm_tlm%Laplacia_TLM(oprand, d_oprand, d_value)
    END IF

    ! Calculate different adj_value scenarios
    ! Scenario 1: Total sum
    adj_value = 1.0_R_KIND
    CALL calculate_and_write_inner_product(sub_name, 'Total Sum', d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)

    ! Scenario 2: Mean
    adj_value = 1.0_R_KIND / REAL(SIZE(VALUE), r_kind)
    CALL calculate_and_write_inner_product(sub_name, 'Mean', d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)

    ! Scenario 3: Sum of Squares
    adj_value = 2.0_R_KIND * VALUE
    CALL calculate_and_write_inner_product(sub_name, 'Sum of Squares', d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)

    ! Scenario 4: Max Value
    ! Scenario 4: Max Value
    adj_value = 0.0_R_KIND
    max_val = MAXVAL(VALUE)

    WHERE (VALUE == max_val)
      adj_value = 1.0_R_KIND
    END WHERE
    CALL calculate_and_write_inner_product(sub_name, 'Max Value', d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)

    ! Scenario 5: Log Sum
    adj_value = 1.0_R_KIND / VALUE
    CALL calculate_and_write_inner_product(sub_name, 'Log Sum', d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)

  END SUBROUTINE test_adjoint_operator_alongwith_adj_value

  SUBROUTINE calculate_and_write_inner_product(sub_name, scenario, d_value, adj_value, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, file_name)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name, scenario
    REAL(r_kind), INTENT(IN) :: d_value(:, :)
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :), adj_value(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit
    REAL(r_kind) :: inner_product_tl, inner_product_adj

    ! Calculate inner product <M' * u, v>
    inner_product_tl = SUM(d_value * adj_value)

    ! Run adjoint operator on adj_value to get the adjoint outputs
    IF (sub_name == 'Divergen') THEN
      CALL gzm_adj%Divergen_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    ELSE IF (sub_name == 'Jacobian') THEN
      CALL gzm_adj%Jacobian_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
    ELSE IF (sub_name == 'Laplacia') THEN
      CALL gzm_adj%Laplacia_AD(oprand, adj_oprand, adj_value)
    END IF

    ! Calculate inner product <u, M* * v>
    IF (sub_name == 'Divergen' .OR. sub_name == 'Jacobian') THEN
      inner_product_adj = SUM(d_opr_left * adj_opr_left + d_opr_right * adj_opr_right)
    ELSE IF (sub_name == 'Laplacia') THEN
      inner_product_adj = SUM(d_oprand * adj_oprand)
    END IF

    ! Compare inner products and write results to file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) sub_name//' '//scenario//' adjoint'
    WRITE (unit, *) 'Inner Product TL:', inner_product_tl
    WRITE (unit, *) 'Inner Product AD:', inner_product_adj
    IF (ABS(inner_product_tl - inner_product_adj) < 1.0E-6_R_KIND) THEN
      WRITE (unit, *) 'The adjoint operator test has passed for scenario ', scenario
    ELSE
      WRITE (unit, *) 'WARNING: The adjoint operator test has failed for scenario ', scenario
    END IF
    CLOSE (unit)
  END SUBROUTINE calculate_and_write_inner_product

! ---------------------------------------------------------------------------
! Subroutine to generate test data

  SUBROUTINE generate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, vLevel, num_icell)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    REAL(r_kind), INTENT(OUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    INTEGER, INTENT(IN) :: vLevel, num_icell
    REAL(r_kind) :: epsilon
    INTEGER :: i, j

    epsilon = 1.0E-6_R_KIND

    ! Initialize operator data with a base pattern (e.g., combining sine, cosine, and a linear gradient)
    DO i = 1, vLevel
      DO j = 1, num_icell
        opr_left(i, j) = SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                         EXP(-REAL(i) / REAL(vLevel)) + 0.1_R_KIND * REAL(j) / REAL(num_icell)

        opr_right(i, j) = COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                          EXP(-REAL(i) / REAL(vLevel)) + 0.1_R_KIND * REAL(i) / REAL(vLevel)

        oprand(i, j) = SIN(3.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel)) * &
                       COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell))
        d_opr_left(i, j) = epsilon * (2.0_R_KIND * RAND() - 1.0_R_KIND)
        d_opr_right(i, j) = epsilon * (2.0_R_KIND * RAND() - 1.0_R_KIND)
        d_oprand(i, j) = epsilon * (2.0_R_KIND * RAND() - 1.0_R_KIND)
      END DO
    END DO
  END SUBROUTINE generate_test_data

END PROGRAM test_adjoint_operator
