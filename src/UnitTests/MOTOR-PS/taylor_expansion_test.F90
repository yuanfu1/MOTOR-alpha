PROGRAM test_taylor_expansion
  USE kinds_m, ONLY: r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: VALUE(:, :), value_perturbed(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
  REAL(r_kind), ALLOCATABLE :: d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: random_perturbation(:, :)
  REAL(r_kind) :: epsilon, ratio, lambda, diff_norm, tl_norm
  INTEGER :: vLevel, num_icell, unit
  INTEGER :: i
  REAL(r_kind), PARAMETER :: epsilons(11) = (/1.0E-9_R_KIND, 1.0E-8_R_KIND, 1.0E-7_R_KIND, 1.0E-6_R_KIND, 1.0E-5_R_KIND, &
                                              1.0E-4_R_KIND, 1.0E-3_R_KIND, 1.0E-2_R_KIND, 1.0E-1_R_KIND, 1.0E+0_R_KIND, 1.0E+1_R_KIND/)
  CHARACTER(LEN=1024) :: configFile, ncOutputFile, output_file

  ! Configure file path
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)

  ! Associate SingleGrid
  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  gzm%sg => sg
  gzm_tlm%sg => sg
  gzm_adj%sg => sg

  output_file = 'taylor_expansion_results.txt'

  ! Open file for writing results
  OPEN (UNIT=unit, FILE=output_file, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')

  ! Allocate arrays and initialize
  vLevel = sg%vLevel
  num_icell = sg%num_icell
  ALLOCATE (opr_left(vLevel, num_icell), opr_right(vLevel, num_icell), VALUE(vLevel, num_icell))
  ALLOCATE (d_opr_left(vLevel, num_icell), d_opr_right(vLevel, num_icell), d_value(vLevel, num_icell))
  ALLOCATE (oprand(vLevel, num_icell), d_oprand(vLevel, num_icell))
  ALLOCATE (random_perturbation(SIZE(opr_left, 1), SIZE(opr_left, 2)))

  ! Step 1: Generate initial test data
  CALL RANDOM_NUMBER(opr_left)
  CALL RANDOM_NUMBER(opr_right)
  CALL RANDOM_NUMBER(oprand)
  CALL RANDOM_NUMBER(d_opr_left)
  CALL RANDOM_NUMBER(d_opr_right)
  CALL RANDOM_NUMBER(d_oprand)
  CALL RANDOM_NUMBER(random_perturbation)

  opr_left = opr_left * 10.0_R_KIND
  opr_right = opr_right * 10.0_R_KIND
  oprand = oprand * 12.6_R_KIND
  d_opr_left = opr_left + random_perturbation * 0.1_R_KIND
  d_opr_right = opr_right + random_perturbation * 0.1_R_KIND
  d_oprand = d_oprand + random_perturbation * 0.1_R_KIND

!   CALL simulate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, vLevel, num_icell)

  ! Step 2: Run tests for Divergen
  CALL gzm%Divergen(opr_left, opr_right, VALUE)
  WRITE (unit, *) 'Testing Divergen...'
  CALL test_taylor_for_subroutine('Divergen', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, VALUE, unit)

  ! Step 3: Run tests for Jacobian
  CALL gzm%Jacobian(opr_left, opr_right, VALUE)
  WRITE (unit, *) 'Testing Jacobian...'
  CALL test_taylor_for_subroutine('Jacobian', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, VALUE, unit)

  ! Step 4: Run tests for Laplacia
  CALL gzm%Laplacia(oprand, VALUE)
  WRITE (unit, *) 'Testing Laplacia...'
  CALL test_taylor_for_subroutine('Laplacia', opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, VALUE, unit)

  ! Deallocate all arrays
  DEALLOCATE (opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)

  ! Close the file
  CLOSE (unit)

CONTAINS

  SUBROUTINE test_taylor_for_subroutine(sub_name, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, VALUE, unit)
    USE, INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: sub_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: VALUE(:, :), d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    INTEGER, INTENT(IN) :: unit
    REAL(r_kind) :: lambda, ratio, diff_norm, tl_norm
    INTEGER :: i, vLevel, num_icell
    REAL(r_kind), ALLOCATABLE :: value_perturbed(:, :)

    vLevel = SIZE(opr_left, 1)
    num_icell = SIZE(opr_left, 2)
    ALLOCATE (value_perturbed(vLevel, num_icell))

    WRITE (unit, *) 'Lambda', 'Ratio'
    DO i = 1, SIZE(epsilons)
      lambda = epsilons(i)

      ! Perturb inputs and run forward model with perturbed input
      IF (sub_name == 'Divergen') THEN
        CALL gzm%Divergen(opr_left + lambda * d_opr_left, opr_right + lambda * d_opr_right, value_perturbed)
      ELSE IF (sub_name == 'Jacobian') THEN
        CALL gzm%Jacobian(opr_left + lambda * d_opr_left, opr_right + lambda * d_opr_right, value_perturbed)
      ELSE IF (sub_name == 'Laplacia') THEN
        CALL gzm%Laplacia(oprand + lambda * d_oprand, value_perturbed)
      END IF

      ! Run tangent-linear model
      IF (sub_name == 'Divergen') THEN
        CALL gzm_tlm%Divergen_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      ELSE IF (sub_name == 'Jacobian') THEN
        CALL gzm_tlm%Jacobian_TLM(opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      ELSE IF (sub_name == 'Laplacia') THEN
        CALL gzm_tlm%Laplacia_TLM(oprand, d_oprand, d_value)
      END IF

      ! Compute norms and check for invalid values
      diff_norm = NORM2(value_perturbed - VALUE)
      tl_norm = NORM2(lambda * d_value)

      IF (tl_norm > 0.0_R_KIND) THEN
        ! Compute ratio
        ratio = diff_norm / tl_norm
        WRITE (unit, '(F20.17, 2X, F20.4)') lambda, ratio
      ELSE
        WRITE (unit, *) "Warning: Skipping ratio calculation due to invalid tl_norm."
      END IF
    END DO

    DEALLOCATE (value_perturbed)

  END SUBROUTINE test_taylor_for_subroutine

  SUBROUTINE simulate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, vLevel, num_icell)
    IMPLICIT NONE
    REAL(r_kind), INTENT(OUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    INTEGER, INTENT(IN) :: vLevel, num_icell

    INTEGER :: i, j

    ! Simulate operator data with controlled patterns
    DO i = 1, vLevel
      DO j = 1, num_icell
        opr_left(i, j) = SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                         COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel))

        opr_right(i, j) = COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                          SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel))

        oprand(i, j) = SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel)) * &
                       SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell))
      END DO
    END DO

    ! Generate tangent linear perturbations as small deviations
    d_opr_left = 0.1_R_KIND * opr_left
    d_opr_right = 0.1_R_KIND * opr_right
    d_oprand = 0.1_R_KIND * oprand

  END SUBROUTINE simulate_test_data

  FUNCTION NORM2(array) RESULT(norm_value)
    USE, INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: array(:, :)
    REAL(r_kind) :: norm_value
    norm_value = SQRT(SUM(array**2))

    ! Check for invalid norm values
    IF (norm_value == 0.0_R_KIND) THEN
      PRINT *, "Warning: Norm value is zero. Check your inputs."
    END IF
  END FUNCTION norm2

END PROGRAM test_taylor_expansion
