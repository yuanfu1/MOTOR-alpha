PROGRAM test_grid_point_adjoint_without_tl
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

  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :), pert_oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), adj_oprand(:, :)

  REAL(r_kind) :: epsilon
  INTEGER(i_kind) :: i, unit, test_points
  INTEGER :: vLevel, num_icell

  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Open output file
  WRITE (*, *) 'Attempting to open output file...'
  OPEN (UNIT=unit, FILE='adjoint_test_multiple_grid_point_results_for_gzm.txt', STATUS='UNKNOWN')
  WRITE (*, *) 'Output file opened successfully.'

  ! Configure file paths
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"
  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialization
  WRITE (*, *) 'Initializing mpddGlob and geometry...'
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)

  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  gzm%sg => sg
  gzm_adj%sg => sg

  WRITE (*, *) 'Allocating arrays...'
  ! Allocate arrays
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))

  epsilon = 1.0E-5_R_KIND
  WRITE (*, *) 'Arrays allocated.'

  ! Define the number of test points as you wish
  test_points = MIN(1000, sg%num_icell)

  ! Multiple point test for Divergen_AD
  WRITE (unit, *) 'Testing Divergen_AD for multiple points...'
  DO i = 1, test_points
    CALL test_adjoint_operator_for_point(unit, 'Divergen', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon, i)
  END DO

  ! Multiple point test for Jacobian_AD
  WRITE (unit, *) 'Testing Jacobian_AD for multiple points...'
  DO i = 1, test_points
    CALL test_adjoint_operator_for_point(unit, 'Jacobian', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon, i)
  END DO

  ! Multiple point test for Laplacia_AD
  WRITE (unit, *) 'Testing Laplacia_AD for multiple points...'
  DO i = 1, test_points
    CALL test_adjoint_operator_for_point(unit, 'Laplacia', opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon, i)
  END DO

  ! Close output file
  WRITE (*, *) 'Closing output file...'
  CLOSE (unit)
  WRITE (*, *) 'Output file closed.'

  ! Deallocate memory
  DEALLOCATE (opr_left, opr_right, oprand, adj_value, adj_opr_left, adj_opr_right, adj_oprand)

CONTAINS

  SUBROUTINE test_adjoint_operator_for_point(unit, sub_name, opr_left, opr_right, oprand, adj_opr_left, adj_opr_right, adj_oprand, VALUE, adj_value, epsilon, point_index)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit, point_index
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :), epsilon
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :)
    REAL(r_kind), INTENT(OUT)  :: VALUE(:, :), adj_value(:, :)
    REAL(r_kind) :: num_grad, adj_grad, absolute_diff
    REAL(r_kind), ALLOCATABLE :: pert_value(:, :)

    ALLOCATE (pert_value(SIZE(VALUE, 1), SIZE(VALUE, 2)))

    ! Debugging output
    WRITE (unit, *) 'Entering test_adjoint_operator_for_point for ', sub_name
    WRITE (unit, *) 'Point Index:', point_index

    ! Initialize adjoint variables
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    adj_oprand = 0.0_R_KIND
    adj_value = 0.0_R_KIND

    ! Step 1: Compute value using the forward operator
    IF (sub_name == 'Divergen') THEN
      WRITE (unit, *) 'Calling Divergen...'
      CALL gzm%Divergen(opr_left, opr_right, VALUE)
    ELSE IF (sub_name == 'Jacobian') THEN
      WRITE (unit, *) 'Calling Jacobian...'
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
    ELSE IF (sub_name == 'Laplacia') THEN
      WRITE (unit, *) 'Calling Laplacia...'
      CALL gzm%Laplacia(oprand, VALUE)
    END IF

    WRITE (unit, *) 'Finished computing forward operator'

    ! Step 2: Compute numerical gradient
    num_grad = 0.0_R_KIND
    adj_grad = 0.0_R_KIND

    IF (sub_name == 'Divergen') THEN
      CALL perturb_and_compute_gradient('Divergen', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Divergen_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      adj_grad = SUM(adj_opr_left)

    ELSE IF (sub_name == 'Jacobian') THEN
      CALL perturb_and_compute_gradient('Jacobian', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Jacobian_AD(opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      adj_grad = SUM(adj_opr_right)

    ELSE IF (sub_name == 'Laplacia') THEN
      CALL perturb_and_compute_gradient('Laplacia', opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)
      adj_value = 1.0_R_KIND
      CALL gzm_adj%Laplacia_AD(oprand, adj_oprand, adj_value)
      adj_grad = SUM(adj_oprand)
    END IF

    WRITE (unit, *) 'Finished computing gradients'

    ! Step 3: Compare gradients using absolute difference
    absolute_diff = ABS(num_grad - adj_grad)
    WRITE (unit, *) 'Point Index:', point_index, ' ', sub_name//' adjoint test'
    WRITE (unit, *) 'Numerical Gradient:', num_grad
    WRITE (unit, *) 'Adjoint Gradient:', adj_grad
    WRITE (unit, *) 'Absolute Difference:', absolute_diff

    IF (absolute_diff < 1.0E-6_R_KIND) THEN
      WRITE (unit, *) 'The adjoint operator test has passed for point ', point_index, ' in subroutine ', sub_name
    ELSE
      WRITE (unit, *) 'WARNING: The adjoint operator test has failed for point ', point_index, ' in subroutine ', sub_name
    END IF

    DEALLOCATE (pert_value)
  END SUBROUTINE test_adjoint_operator_for_point

  SUBROUTINE perturb_and_compute_gradient(sub_name, opr_left, opr_right, oprand, VALUE, pert_value, num_grad, epsilon)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(IN) :: VALUE(:, :)
    REAL(r_kind), INTENT(OUT) :: pert_value(:, :)
    REAL(r_kind), INTENT(OUT) :: num_grad
    REAL(r_kind), INTENT(IN) :: epsilon
    INTEGER :: vLevel, num_icell

    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(VALUE, 2)

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

END PROGRAM test_grid_point_adjoint_without_tl
