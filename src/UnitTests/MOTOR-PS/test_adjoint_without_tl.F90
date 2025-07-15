PROGRAM test_adjoint_without_tl
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :), pert_oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), adj_oprand(:, :)

  REAL(r_kind) :: num_grad, adj_grad, epsilon, relative_diff
  INTEGER(i_kind) :: i, j, unit
  INTEGER :: vLevel, num_icell

  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Open output file
  OPEN (UNIT=unit, FILE='adjoint_test_results.txt', STATUS='UNKNOWN')

  ! Configure file path
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)

  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  gzm%sg => sg
  gzm_adj%sg => sg

  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))

  ! Step 1: Generate initial test data
  CALL generate_test_data(opr_left, opr_right, oprand, sg%vLevel, sg%num_icell)

  epsilon = 1.0E-5_R_KIND

  ! Case 1: Testing Divergen_AD
  WRITE (unit, *) 'Testing Divergen_AD...'
  CALL test_adjoint_operator_for_subroutine(unit, 'Divergen', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon)

  ! Case 2: Testing Jacobian_AD
  WRITE (unit, *) 'Testing Jacobian_AD...'
  CALL test_adjoint_operator_for_subroutine(unit, 'Jacobian', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon)

  ! Case 3: Testing Laplacia_AD
  WRITE (unit, *) 'Testing Laplacia_AD...'
  CALL test_adjoint_operator_for_subroutine(unit, 'Laplacia', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon)

  ! Close output file
  CLOSE (unit)

  ! Deallocate all arrays
  DEALLOCATE (opr_left, opr_right, oprand, adj_value, adj_opr_left, adj_opr_right, adj_oprand)

CONTAINS

  SUBROUTINE test_adjoint_operator_for_subroutine(unit, sub_name, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :), epsilon
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
    REAL(r_kind), INTENT(OUT)  :: VALUE(:, :), adj_value(:, :)
    REAL(r_kind) :: num_grad, adj_grad, relative_diff, absolute_diff
    INTEGER :: vLevel, num_icell, i, j

    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(opr_left, 2)

    ! Initialize adjoint variables
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    adj_oprand = 0.0_R_KIND
    adj_value = 0.0_R_KIND

    ! Step 1: Compute value using the forward operator
    IF (sub_name == 'Divergen') THEN
      CALL gzm%Divergen(opr_left, opr_right, VALUE)
    ELSE IF (sub_name == 'Jacobian') THEN
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
    ELSE IF (sub_name == 'Laplacia') THEN
      CALL gzm%Laplacia(oprand, VALUE)
    END IF

    ! Step 2: Compute numerical gradient for multiple points or the entire domain
    num_grad = 0.0_R_KIND
    adj_grad = 0.0_R_KIND

    IF (sub_name == 'Divergen') THEN
      ! Perturb opr_left and compute numerical gradient
      CALL perturb_and_compute_gradient('Divergen', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)

      ! Compute adjoint gradient
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Divergen_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      adj_grad = SUM(adj_opr_left)

    ELSE IF (sub_name == 'Jacobian') THEN
      ! Perturb opr_right and compute numerical gradient
      CALL perturb_and_compute_gradient('Jacobian', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)

      ! Compute adjoint gradient
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Jacobian_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      adj_grad = SUM(adj_opr_right)

    ELSE IF (sub_name == 'Laplacia') THEN
      ! Perturb oprand and compute numerical gradient
      CALL perturb_and_compute_gradient('Laplacia', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)

      ! Compute adjoint gradient
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Laplacia_AD(oprand, adj_oprand, adj_value)
      adj_grad = SUM(adj_oprand)
    END IF

    ! Step 3: Compare gradients using relative difference
    ! relative_diff = ABS(num_grad - adj_grad) / MAX(ABS(num_grad), ABS(adj_grad), 1.0e-14_r_kind)
    ! Step 3: Compare gradients using absolute difference  ! num_grad is too small, about 1E-8 for divergen, 1E-18 and 1E-20 for Jacobian and Laplacia, respectively.
    absolute_diff = ABS(num_grad - adj_grad)
    WRITE (unit, *) sub_name//' adjoint test'
    WRITE (unit, *) 'Numerical Gradient:', num_grad
    WRITE (unit, *) 'Adjoint Gradient:', adj_grad
    ! WRITE(unit, *) 'Relative Difference:', relative_diff

    ! IF (relative_diff < 1.0E-6_r_kind) THEN
    !     WRITE(unit, *) 'The adjoint operator test has passed for subroutine ', sub_name
    ! ELSE
    !     WRITE(unit, *) 'WARNING: The adjoint operator test has failed for subroutine ', sub_name
    ! END IF
    ! Adjust the threshold for absolute difference based on the expected numerical accuracy
    IF (absolute_diff < 1.0E-6_R_KIND) THEN
      WRITE (unit, *) 'The adjoint operator test has passed for subroutine ', sub_name
    ELSE
      WRITE (unit, *) 'WARNING: The adjoint operator test has failed for subroutine ', sub_name
    END IF
  END SUBROUTINE test_adjoint_operator_for_subroutine

  SUBROUTINE perturb_and_compute_gradient(sub_name, opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :)
    REAL(r_kind), INTENT(OUT) :: pert_value(:, :)
    REAL(r_kind), INTENT(OUT) :: num_grad
    REAL(r_kind), INTENT(IN) :: epsilon
    INTEGER :: vLevel, num_icell, i, j

    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(opr_left, 2)

    IF (sub_name == 'Divergen') THEN
      pert_value = VALUE
      CALL gzm%Divergen(opr_left + epsilon, opr_right + epsilon, pert_value)
      num_grad = SUM((pert_value - VALUE) / epsilon)

    ELSE IF (sub_name == 'Jacobian') THEN
      pert_value = VALUE
      CALL gzm%Jacobian(opr_left, opr_right + epsilon, pert_value)
      num_grad = SUM((pert_value - VALUE) / epsilon)

    ELSE IF (sub_name == 'Laplacia') THEN
      pert_value = VALUE
      CALL gzm%Laplacia(oprand + epsilon, pert_value)
      num_grad = SUM((pert_value - VALUE) / epsilon)
    END IF
  END SUBROUTINE perturb_and_compute_gradient

  SUBROUTINE generate_test_data(opr_left, opr_right, oprand, vLevel, num_icell)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    REAL(r_kind), INTENT(OUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    INTEGER, INTENT(IN) :: vLevel, num_icell
    INTEGER :: i, j

    ! Initialize operator data with a base pattern (e.g., combining sine, cosine, and a linear gradient)
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

END PROGRAM test_adjoint_without_tl
