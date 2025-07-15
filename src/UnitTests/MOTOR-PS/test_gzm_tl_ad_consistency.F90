PROGRAM test_gzm_tl_ad_consistency
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
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
  CHARACTER(LEN=256) :: TL_AD_result_file, TL_result_file, Taylor_result_file
  INTEGER(i_kind) :: i, j, k_test, l
  LOGICAL :: ieee_support
  ! Variables and parameters for Taylor expension verification
  REAL(r_kind) :: norm_diff, norm_tlm, norm_x, norm_dx
  REAL(r_kind), ALLOCATABLE :: lambda_values(:)
  INTEGER :: num_lambdas

  REAL(r_kind) :: lambda, ratio, rel_error, abs_error, dot_product, angle
  REAL(r_kind) :: f_base, f_perturb
  REAL(r_kind), ALLOCATABLE :: x1(:, :), x2(:, :), x3(:, :), dx1(:, :), dx2(:, :), dx3(:, :)
  REAL(r_kind), ALLOCATABLE :: M_x(:, :), M_x_perturb(:, :), M_tlm(:, :)

  REAL(r_kind), ALLOCATABLE :: random_perturbation(:, :)

  INTEGER(i_kind) :: vLevel, num_icell

  ! Taylor expansion coefficients (from the provided diagram)
!   num_lambdas = 11
!   ALLOCATE(lambda_values(num_lambdas))
!   lambda_values = (/ 1.0E-09_r_kind, 1.0E-08_r_kind, 1.0E-07_r_kind, 1.0E-06_r_kind, &
!                     1.0E-05_r_kind, 1.0E-04_r_kind, 1.0E-03_r_kind, 1.0E-02_r_kind, &
!                     1.0E-01_r_kind, 1.0E+00_r_kind, 1.0E+01_r_kind /)

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

  k_test = 10
  tolerance = 1.0E-5_R_KIND
  TL_AD_result_file = "tl_ad_unit_test_result.txt"
  TL_result_file = "tl_system_test_result.txt"
  Taylor_result_file = "taylor_expansion_test_result.txt"
  ! allocate arrays and initialize
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))
  ALLOCATE (random_perturbation(SIZE(opr_left, 1), SIZE(opr_left, 2)))

  ALLOCATE (x1(sg%vLevel, sg%num_icell), x2(sg%vLevel, sg%num_icell), x3(sg%vLevel, sg%num_icell))
  ALLOCATE (dx1(sg%vLevel, sg%num_icell), dx2(sg%vLevel, sg%num_icell), dx3(sg%vLevel, sg%num_icell))
  ALLOCATE (M_x(sg%vLevel, sg%num_icell), M_x_perturb(sg%vLevel, sg%num_icell), M_tlm(sg%vLevel, sg%num_icell))

!    CALL RANDOM_NUMBER(opr_left)
!    CALL RANDOM_NUMBER(opr_right)
  !CALL RANDOM_NUMBER(d_opr_left)
  !CALL RANDOM_NUMBER(d_opr_right)
  !CALL RANDOM_NUMBER(adj_value)

! case 1

  DO i = 1, sg%vLevel
    d_opr_left(i, :) = (DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))) * 10.0D0
    d_opr_right(i, :) = (DCOS(sg%cell_cntr(1, :)) + DSIN(sg%cell_cntr(2, :))) * 10.0D0
    d_oprand(i, :) = DSIN(sg%cell_cntr(1, :)) * 1.0D1
    adj_opr_left(i, :) = 0.0D0
    adj_opr_right(i, :) = 0.0D0
    adj_oprand(i, :) = 0.0D0
  END DO

  epsilon = 10.0_R_KIND

  opr_left = d_opr_left + 1.56E3
  opr_right = d_opr_right + 1.66E3 - 10000D0
  oprand = d_oprand + 1.03E5 + 80000D0
! case 2
  CALL RANDOM_NUMBER(opr_left)
  CALL RANDOM_NUMBER(opr_right)
  CALL RANDOM_NUMBER(oprand)
  CALL RANDOM_NUMBER(d_opr_left)
  CALL RANDOM_NUMBER(d_opr_right)
  CALL RANDOM_NUMBER(d_oprand)
  CALL RANDOM_NUMBER(adj_value)
  CALL RANDOM_NUMBER(random_perturbation)

  opr_left = opr_left * 10.0_R_KIND
  opr_right = opr_right * 10.0_R_KIND
  oprand = oprand * 12.6_R_KIND
  d_opr_left = opr_left + random_perturbation * 0.1_R_KIND
  d_opr_right = opr_right + random_perturbation * 0.1_R_KIND
  d_oprand = d_oprand + random_perturbation * 0.1_R_KIND
!   opr_left = EXP(opr_left)
!   opr_right = EXP(opr_right)
!   opr_left = opr_left + 5.0_r_kind
!   opr_right = opr_right + 5.0_r_kind
! CASE 3
  CALL generate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, sg%vLevel, sg%num_icell)

  ! Initialize x1, x2, x3, dx1, dx2, dx3
  ! These would typically be initialized based on the model's inputs and perturbations
  ! Example initialization (this should be replaced with actual model data):
  x1 = opr_left
  x2 = opr_right
  x3 = oprand
  dx1 = d_opr_left
  dx2 = d_opr_right
  dx3 = d_oprand

  ! set error threshold
  tolerance = 1.0E-10_R_KIND

  ! Perform TL-AD tests
  CALL TL_AD_Test("Divergen", k_test, tolerance, TL_AD_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand)
  CALL TL_AD_Test("Jacobian", k_test, tolerance, TL_AD_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand)
  CALL TL_AD_Test("Laplacia", k_test, tolerance, TL_AD_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand)

  ! Run TL-specific tests
!   CALL TL_RelativeError_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value )
!   CALL TL_AbsoluteError_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!  CALL TL_DotProduct_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_Angle_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_SensitivityCoefficient_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)

!   CALL TL_RelativeError_Test("Jacobian", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value )
!   CALL TL_AbsoluteError_Test("Jacobian", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_DotProduct_Test("Jacobian", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_Angle_Test("Jacobian", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_SensitivityCoefficient_Test("Jacobian", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)

!   CALL TL_RelativeError_Test("Laplacia", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value )
!   CALL TL_AbsoluteError_Test("Laplacia", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_DotProduct_Test("Laplacia", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_Angle_Test("Laplacia", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)
!   CALL TL_SensitivityCoefficient_Test("Laplacia", tolerance, TL_result_file, opr_left, opr_right, oprand, value, d_opr_left, d_opr_right, d_oprand, d_value)

!   CALL TL_RelativeError_Test("Divergen", tolerance, TL_result_file, opr_left, opr_right, value, d_opr_left, d_opr_right, d_value)

  CALL Verify_TLM("Divergen", tolerance, TL_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_opr_left, pert_opr_right, pert_oprand)
  CALL Verify_TLM("Jacobian", tolerance, TL_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_opr_left, pert_opr_right, pert_oprand)
  CALL Verify_TLM("Laplacia", tolerance, TL_result_file, opr_left, opr_right, VALUE, &
                  d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_opr_left, pert_opr_right, pert_oprand)
  ! Perform Taylor expansion test
  ! Perform Taylor expansion test for Divergen_TLM, Jacobian_TLM, and Laplacia_TLM
  ! Set up lambda values for Taylor expansion test
  num_lambdas = 11
  ALLOCATE (lambda_values(num_lambdas))

  lambda_values(1) = 0.1E-09_R_KIND
  lambda_values(2) = 0.1E-08_R_KIND
  lambda_values(3) = 0.1E-07_R_KIND
  lambda_values(4) = 0.1E-06_R_KIND
  lambda_values(5) = 0.1E-05_R_KIND
  lambda_values(6) = 0.1E-04_R_KIND
  lambda_values(7) = 0.1E-03_R_KIND
  lambda_values(8) = 0.1E-02_R_KIND
  lambda_values(9) = 0.1E-01_R_KIND
  lambda_values(10) = 0.1E+00_R_KIND
  lambda_values(11) = 0.1E+01_R_KIND

!   CALL Test_TLM_Taylor("Divergen", gzm, gzm_tlm, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, k_test, lambda_values, num_lambdas, tolerance, Taylor_result_file)
!   CALL Test_TLM_Taylor("Jacobian", gzm, gzm_tlm, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, k_test, lambda_values, num_lambdas, tolerance, Taylor_result_file)
!   CALL Test_TLM_Taylor("Laplacia", gzm, gzm_tlm, opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, k_test, lambda_values, num_lambdas, tolerance, Taylor_result_file)

! !   CALL Test_TLM_Taylor("Divergen", gzm, gzm_tlm, opr_left(:,:), opr_right(:,:), oprand(:,:), &
! !                        d_opr_left(k_test,:), d_opr_right(k_test,:), d_oprand(k_test,:), lambda_values, num_lambdas, tolerance, Taylor_result_file)
! !   CALL Test_TLM_Taylor("Jacobian", gzm, gzm_tlm, opr_left(k_test,:), opr_right(k_test,:), oprand(k_test,:), &
! !                        d_opr_left(k_test,:), d_opr_right(k_test,:), d_oprand(k_test,:), lambda_values, num_lambdas, tolerance, Taylor_result_file)
! !   CALL Test_TLM_Taylor("Laplacia", gzm, gzm_tlm, opr_left(k_test,:), opr_right(k_test,:), oprand(k_test,:), &
! !                        d_opr_left(k_test,:), d_opr_right(k_test,:), d_oprand(k_test,:), lambda_values, num_lambdas, tolerance, Taylor_result_file)

!   ! Deallocate memory after usage
! ! Deallocate memory after usage
! IF (ALLOCATED(x1)) DEALLOCATE(x1)
! IF (ALLOCATED(x2)) DEALLOCATE(x2)
! IF (ALLOCATED(x3)) DEALLOCATE(x3)
! IF (ALLOCATED(dx1)) DEALLOCATE(dx1)
! IF (ALLOCATED(dx2)) DEALLOCATE(dx2)
! IF (ALLOCATED(dx3)) DEALLOCATE(dx3)
! IF (ALLOCATED(M_x)) DEALLOCATE(M_x)
! IF (ALLOCATED(M_x_perturb)) DEALLOCATE(M_x_perturb)
! IF (ALLOCATED(M_tlm)) DEALLOCATE(M_tlm)
  DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, d_oprand, adj_oprand)

CONTAINS

  SUBROUTINE generate_test_data(opr_left, opr_right, oprand, d_opr_left, d_opr_right, d_oprand, vLevel, num_icell)
    USE kinds_m, ONLY: r_kind
    IMPLICIT NONE
    REAL(r_kind), INTENT(OUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
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

        d_opr_left(i, j) = 0.05_R_KIND * COS(2.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel)) * &
                           EXP(-REAL(j) / REAL(num_icell))

        d_opr_right(i, j) = 0.05_R_KIND * SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell)) * &
                            EXP(-REAL(i) / REAL(vLevel))

        d_oprand(i, j) = 0.05_R_KIND * SIN(2.0_R_KIND * 3.14159265359_R_KIND * REAL(i) / REAL(vLevel)) * &
                         COS(3.0_R_KIND * 3.14159265359_R_KIND * REAL(j) / REAL(num_icell))
      END DO
    END DO
  END SUBROUTINE generate_test_data

  SUBROUTINE TL_AD_Test(operator_name, k, tolerance, file_name, opr_left, opr_right, VALUE, &
                        d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand)
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    INTEGER(i_kind), INTENT(IN) :: k
    REAL(r_kind), INTENT(IN) :: tolerance
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :), d_value(:, :)
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :), adj_value(:, :)
    REAL(r_kind) :: error_norm, inner_product_tl, inner_product_ad
    INTEGER :: unit

    ! initializing inner product
    inner_product_tl = 0.0_R_KIND
    inner_product_ad = 0.0_R_KIND

    !adj_value(k, :) = 0.0D0  !0.1E-8_r_kind
    adj_opr_left(k, :) = 0.0_R_KIND
    adj_opr_right(k, :) = 0.0_R_KIND

!   call forward, tangent linear and adjoint model
    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)

      adj_value(k, :) = 1.0_R_KIND  !1/2*SQRT(value(k,:)*value(k,:))
      ! forwardmodel = value
      !CALL Divergen(gzm, opr_left + d_opr_left, opr_right + d_opr_right, value)
      ! pert_forwardmodel = value
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      CALL Divergen_AD(gzm_adj, opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_value(k, j) * adj_value(k, j)
!       print *, "D_VALUE(K,J) after call tlm and adj procedure: ", k, d_value(k,j), adj_value(k,j)
        inner_product_ad = inner_product_ad + opr_left(k, j) * adj_opr_left(k, j) + opr_right(k, j) * adj_opr_right(k, j)
      END DO

      ! calculating error according to L2
      IF (ABS(inner_product_tl) + ABS(inner_product_ad) > 1.0E-12_R_KIND) THEN
        error_norm = ABS(inner_product_tl - inner_product_ad) / (ABS(inner_product_tl) + ABS(inner_product_ad))
      ELSE
        error_norm = 0.0_R_KIND
      END IF

      ! output to a file
      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL-AD Consistency Check at level '//TRIM(ADJUSTL(CHAR(k)))
      WRITE (unit, *) 'Inner Product TL:', inner_product_tl
      WRITE (unit, *) 'Inner Product AD:', inner_product_ad
      WRITE (unit, *) 'Error Norm:', error_norm
      IF (error_norm < tolerance) THEN
        WRITE (unit, *) 'Consistency Check PASSED.'
      ELSE
        WRITE (unit, *) 'Consistency Check FAILED.'
      END IF
      CLOSE (unit)
    END IF

    IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      adj_value(k, :) = 1.0_R_KIND  !1/2*SQRT(value(k,:)*value(k,:))
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      !adj_value = 0.1E-14_r_kind
      CALL Jacobian_AD(gzm_adj, opr_left, opr_right, adj_opr_left, adj_opr_right, adj_value)
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_value(k, j) * adj_value(k, j)
!      print *, "D_VALUE(K,J) after call tlm and adj procedure: ", k, d_value(k,j), adj_value(k,j)
        inner_product_ad = inner_product_ad + opr_left(k, j) * adj_opr_left(k, j) + opr_right(k, j) * adj_opr_right(k, j)
      END DO

      ! calculating error according to L2
      IF (ABS(inner_product_tl) + ABS(inner_product_ad) > 1.0E-12_R_KIND) THEN
        error_norm = ABS(inner_product_tl - inner_product_ad) / (ABS(inner_product_tl) + ABS(inner_product_ad))
      ELSE
        error_norm = 0.0_R_KIND
      END IF

      ! output to a file
      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL-AD Consistency Check at level '//TRIM(ADJUSTL(CHAR(k)))
      WRITE (unit, *) 'Inner Product TL:', inner_product_tl
      WRITE (unit, *) 'Inner Product AD:', inner_product_ad
      WRITE (unit, *) 'Error Norm:', error_norm
      IF (error_norm < tolerance) THEN
        WRITE (unit, *) 'Consistency Check PASSED.'
      ELSE
        WRITE (unit, *) 'Consistency Check FAILED.'
      END IF
      CLOSE (unit)
    END IF

    IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      adj_value(k, :) = 1.0_R_KIND  !1/2*SQRT(value(k,:)*value(k,:))
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
      CALL Laplacia_AD(gzm_adj, oprand, adj_oprand, adj_value)
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_value(k, j) * adj_value(k, j)
        inner_product_ad = inner_product_ad + oprand(k, j) * adj_oprand(k, j)
      END DO
      IF (ABS(inner_product_tl) + ABS(inner_product_ad) > 1.0E-12_R_KIND) THEN
        error_norm = ABS(inner_product_tl - inner_product_ad) / (ABS(inner_product_tl) + ABS(inner_product_ad))
      ELSE
        error_norm = 0.0_R_KIND
      END IF

      !CALL Write_TL_AD_Result(TRIM(operator_name) // " at level " // TRIM(ADJUSTL(CHAR(k))), inner_product_tl, inner_product_ad, error_norm, tolerance, file_name)
      ! output to a file
      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL-AD Consistency Check at level '//TRIM(ADJUSTL(CHAR(k)))
      WRITE (unit, *) 'Inner Product TL:', inner_product_tl
      WRITE (unit, *) 'Inner Product AD:', inner_product_ad
      WRITE (unit, *) 'Error Norm:', error_norm
      IF (error_norm < tolerance) THEN
        WRITE (unit, *) 'Consistency Check PASSED.'
      ELSE
        WRITE (unit, *) 'Consistency Check FAILED.'
      END IF
      CLOSE (unit)
    END IF

  END SUBROUTINE TL_AD_Test

! SUBROUTINE Write_TL_AD_Result(operator_name, inner_tl, inner_ad, error_norm, tolerance, file_name)
!     CHARACTER(LEN=*), INTENT(IN) :: operator_name
!     REAL(r_kind), INTENT(IN) :: inner_tl, inner_ad, error_norm, tolerance
!     CHARACTER(LEN=*), INTENT(IN) :: file_name
!     INTEGER :: unit

!     ! open a file with append mode
!     OPEN(UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
!     WRITE(unit, *) operator_name // ' TL-AD Consistency Check'
!     WRITE(unit, *) 'Inner Product TL:', inner_tl
!     WRITE(unit, *) 'Inner Product AD:', inner_ad
!     WRITE(unit, *) 'Error Norm:', error_norm
!     IF (error_norm < tolerance) THEN
!       WRITE(unit, *) 'Consistency Check PASSED.'
!     ELSE
!       WRITE(unit, *) 'Consistency Check FAILED.'
!     END IF
!     CLOSE(unit)
!   END SUBROUTINE Write_TL_AD_Result

  SUBROUTINE Verify_TLM(operator_name, tolerance, file_name, opr_left, opr_right, VALUE, &
                        d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_opr_left, pert_opr_right, pert_oprand)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    INTEGER(i_kind) :: k, j, n
    REAL(r_kind), INTENT(IN) :: tolerance
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), oprand(:, :)
    REAL(r_kind), INTENT(OUT) :: VALUE(:, :), d_value(:, :)
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :)
    REAL(r_kind), INTENT(INOUT) :: adj_opr_left(:, :), adj_opr_right(:, :), adj_oprand(:, :), adj_value(:, :)
    REAL(r_kind), INTENT(INOUT) :: pert_opr_left(:, :), pert_opr_right(:, :), pert_oprand(:, :)
    INTEGER :: unit

    REAL(r_kind) :: max_error, avg_error

    epsilon = 1.0E-5
    max_error = 0.0
    avg_error = 0.0
    n = 0
    IF (operator_name == "Divergen") THEN
      CALL gzm%Divergen(opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      pert_opr_left = opr_left + epsilon * d_opr_left
      pert_opr_right = opr_right + epsilon * d_opr_right
      CALL Divergen(gzm, pert_opr_left, pert_opr_right, pert_value)

      DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell

        max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
        avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
        n = n + 1

      END DO
      END DO
      avg_error = avg_error / n

      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL Unit Check  '

      WRITE (unit, *) operator_name//' Max Error: ', max_error
      WRITE (unit, *) operator_name//' Avg Error :', avg_error

      IF (max_error < 1E-5) THEN
        WRITE (unit, *) operator_name//' TLM verification passed.'
      ELSE
        WRITE (unit, *) operator_name//' TLM verification failed.'
      END IF
      CLOSE (unit)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL gzm%Jacobian(opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
      pert_opr_left = opr_left + epsilon * d_opr_left
      pert_opr_right = opr_right + epsilon * d_opr_right
      CALL Jacobian(gzm, pert_opr_left, pert_opr_right, pert_value)

      DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell

        max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
        avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
        n = n + 1

      END DO
      END DO
      avg_error = avg_error / n

      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL Unit Check  '

      WRITE (unit, *) operator_name//' Max Error: ', max_error
      WRITE (unit, *) operator_name//' Avg Error :', avg_error

      IF (max_error < 1E-5) THEN
        WRITE (unit, *) operator_name//' TLM verification passed.'
      ELSE
        WRITE (unit, *) operator_name//' TLM verification failed.'
      END IF
      CLOSE (unit)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL gzm%Laplacia(oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
      pert_oprand = oprand + epsilon * d_oprand

      CALL Laplacia(gzm, pert_oprand, pert_value)

      DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell

        max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
        avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
        n = n + 1

      END DO
      END DO
      avg_error = avg_error / n

      OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
      WRITE (unit, *) operator_name//' TL Unit Check  '

      WRITE (unit, *) operator_name//' Max Error: ', max_error
      WRITE (unit, *) operator_name//' Avg Error :', avg_error

      IF (max_error < 1E-5) THEN
        WRITE (unit, *) operator_name//' TLM verification passed.'
      ELSE
        WRITE (unit, *) operator_name//' TLM verification failed.'
      END IF
      CLOSE (unit)
    END IF

  END SUBROUTINE Verify_TLM

  SUBROUTINE TL_RelativeError_Test(operator_name, tolerance, file_name, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)
    ! Purpose: The relative error test is used to assess the relative difference between value and d_value to ensure the accuracy of the tangent linear model
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    REAL(r_kind), INTENT(IN) :: tolerance
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :), VALUE(:, :)
    REAL(r_kind), INTENT(INOUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), d_value(:, :)
    REAL(r_kind) :: relative_error, norm_value, norm_d_value
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit

    ! Call the tangent linear model
    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
    END IF

    ! Calculate the norms of value and d_value
    norm_value = SQRT(SUM(VALUE * VALUE))
    norm_d_value = SQRT(SUM(d_value * d_value))

    ! Calculate the relative error
    IF (norm_value > 0.0_R_KIND) THEN
      relative_error = ABS(norm_d_value - norm_value) / norm_value
    ELSE
      relative_error = 0.0_R_KIND
    END IF

    ! Write the relative error to a file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) operator_name//' Relative Error Test'
    WRITE (unit, *) 'Relative Error:', relative_error
    IF (relative_error < tolerance) THEN
      WRITE (unit, *) 'Test passed with tolerance: ', tolerance
    ELSE
      WRITE (unit, *) 'Test failed with tolerance: ', tolerance
    END IF
    CLOSE (unit)

  END SUBROUTINE TL_RelativeError_Test

  SUBROUTINE TL_AbsoluteError_Test(operator_name, tolerance, file_name, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)
    ! Test Purpose: The absolute error test measures the absolute difference between value and d_value to ensure the accuracy of the tangent linear model.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    REAL(r_kind), INTENT(IN) :: tolerance
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :), VALUE(:, :)
    REAL(r_kind), INTENT(INOUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), d_value(:, :)
    REAL(r_kind) :: absolute_error
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit

    ! Call the tangent linear model
    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
    END IF

    ! Calculate the absolute error
    absolute_error = SUM(ABS(VALUE - d_value))

    ! Write the absolute error to a file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) operator_name//' Absolute Error Test'
    WRITE (unit, *) 'Absolute Error:', absolute_error
    IF (absolute_error < tolerance) THEN
      WRITE (unit, *) 'Test passed with tolerance: ', tolerance
    ELSE
      WRITE (unit, *) 'Test failed with tolerance: ', tolerance
    END IF
    CLOSE (unit)

  END SUBROUTINE TL_AbsoluteError_Test

  SUBROUTINE TL_DotProduct_Test(operator_name, tolerance, file_name, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)
    ! Test Purpose: The dot product test is used to quantify the similarity between value and d_value to verify if the tangent linear model correctly captures directional changes.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    !CLASS(gzm_tlm_t), INTENT(INOUT) :: gzm_tlm
    REAL(r_kind), INTENT(IN) :: tolerance
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :), VALUE(:, :)
    REAL(r_kind), INTENT(INOUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), d_value(:, :)
    REAL(r_kind) :: dot_product
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit

    ! Call the tangent linear model
    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
    END IF

    ! Calculate the dot product between value and d_value
    dot_product = SUM(VALUE * d_value)

    ! Write the dot product to a file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) operator_name//' Dot Product Test'
    WRITE (unit, *) 'Dot Product:', dot_product
    IF (dot_product < tolerance) THEN
      WRITE (unit, *) 'Test passed with tolerance: ', tolerance
    ELSE
      WRITE (unit, *) 'Test failed with tolerance: ', tolerance
    END IF
    CLOSE (unit)

  END SUBROUTINE TL_DotProduct_Test

  SUBROUTINE TL_Angle_Test(operator_name, tolerance, file_name, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)
    ! Test Purpose: The angle test measures the geometric angle between value and d_value to evaluate if the tangent linear model accurately reflects their directional consistency.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    !CLASS(gzm_tlm_t), INTENT(INOUT) :: gzm_tlm
    REAL(r_kind), INTENT(IN) :: tolerance
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :), VALUE(:, :)
    REAL(r_kind), INTENT(INOUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), d_value(:, :)
    REAL(r_kind) :: angle, dot_product, norm_value, norm_d_value
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit

    ! Call the tangent linear model
    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
    END IF
    ! Calculate the dot product, norms of value and d_value, and the angle between them
    dot_product = SUM(VALUE * d_value)
    norm_value = SQRT(SUM(VALUE * VALUE))
    norm_d_value = SQRT(SUM(d_value * d_value))

    IF (norm_value > 0.0_R_KIND .AND. norm_d_value > 0.0_R_KIND) THEN
      angle = ACOS(dot_product / (norm_value * norm_d_value))
    ELSE
      angle = 0.0_R_KIND
    END IF

    ! Write the angle to a file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) operator_name//' Angle Test'
    WRITE (unit, *) 'Angle (radians):', angle
    IF (angle < tolerance) THEN
      WRITE (unit, *) 'Test passed with tolerance: ', tolerance
    ELSE
      WRITE (unit, *) 'Test failed with tolerance: ', tolerance
    END IF
    CLOSE (unit)

  END SUBROUTINE TL_Angle_Test

  SUBROUTINE TL_SensitivityCoefficient_Test(operator_name, tolerance, file_name, opr_left, opr_right, oprand, VALUE, d_opr_left, d_opr_right, d_oprand, d_value)
    ! Test Purpose: The sensitivity coefficient test evaluates the system's response to small perturbations to verify if the tangent linear model correctly reflects the actual trend of changes in the system.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    !CLASS(gzm_tlm_t), INTENT(INOUT) :: gzm_tlm
    REAL(r_kind), INTENT(IN) :: tolerance
    REAL(r_kind), INTENT(INOUT) :: opr_left(:, :), opr_right(:, :), oprand(:, :), VALUE(:, :)
    REAL(r_kind), INTENT(INOUT) :: d_opr_left(:, :), d_opr_right(:, :), d_oprand(:, :), d_value(:, :)
    REAL(r_kind) :: sensitivity_coefficient, perturbation_size
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit

    ! Set the size of the small perturbation
    perturbation_size = 1.0E-10_R_KIND

    IF (operator_name == "Divergen") THEN
      CALL Divergen(gzm, opr_left, opr_right, VALUE)
      CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Jacobian") THEN
      CALL Jacobian(gzm, opr_left, opr_right, VALUE)
      CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, d_opr_left, d_opr_right, d_value)
    ELSE IF (operator_name == "Laplacia") THEN
      CALL Laplacia(gzm, oprand, VALUE)
      CALL Laplacia_TLM(gzm_tlm, oprand, d_oprand, d_value)
    END IF

    ! Calculate the sensitivity coefficient
    sensitivity_coefficient = SUM((d_value - VALUE) / perturbation_size)

    ! Write the sensitivity coefficient to a file
    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
    WRITE (unit, *) operator_name//' Sensitivity Coefficient Test'
    WRITE (unit, *) 'Sensitivity Coefficient:', sensitivity_coefficient
    IF (ABS(sensitivity_coefficient) < tolerance) THEN
      WRITE (unit, *) 'Test passed with tolerance: ', tolerance
    ELSE
      WRITE (unit, *) 'Test failed with tolerance: ', tolerance
    END IF
    CLOSE (unit)

  END SUBROUTINE TL_SensitivityCoefficient_Test

! SUBROUTINE Test_TLM_Taylor(operator_name, model, model_tlm, x1, x2, x3, dx1, dx2, dx3, k_layer, lambda_values, num_lambdas, tolerance, file_name)
!     USE, INTRINSIC :: IEEE_ARITHMETIC
!     CHARACTER(LEN=*), INTENT(IN) :: operator_name
!     TYPE(gzm_t), INTENT(INOUT) :: model
!     TYPE(gzm_tlm_t), INTENT(INOUT) :: model_tlm
!     REAL(r_kind), INTENT(IN) :: x1(:,:), x2(:,:), x3(:,:), dx1(:,:), dx2(:,:), dx3(:,:)
!     INTEGER, INTENT(IN) :: k_layer  ! 特定层数
!     REAL(r_kind), INTENT(IN) :: lambda_values(:)
!     INTEGER, INTENT(IN) :: num_lambdas
!     REAL(r_kind), INTENT(IN) :: tolerance
!     CHARACTER(LEN=*), INTENT(IN) :: file_name

!     REAL(r_kind) :: norm_diff, norm_tlm, taylor_ratio
!     REAL(r_kind) :: lambda
!     INTEGER :: l, unit
!     REAL(r_kind), DIMENSION(1, SIZE(x1,2)) :: M_x, M_x_perturb, M_tlm

!     !
!     OPEN(UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
!     WRITE(unit, *) operator_name // " Taylor Expansion Test"
!     WRITE(unit, *) "Lambda Values and Corresponding Taylor Ratios:"
!     WRITE(unit, *) "Lambda", "Taylor Ratio"

!     !  M(x)
!     IF (operator_name == "Divergen") THEN
!         CALL Divergen(model, RESHAPE(x1(k_layer,:), (/1, SIZE(x1, 2)/)), RESHAPE(x2(k_layer,:), (/1, SIZE(x2, 2)/)), M_x)
!     ELSE IF (operator_name == "Jacobian") THEN
!         CALL Jacobian(model, RESHAPE(x1(k_layer,:), (/1, SIZE(x1, 2)/)), RESHAPE(x2(k_layer,:), (/1, SIZE(x2, 2)/)), M_x)
!     ELSE IF (operator_name == "Laplacia") THEN
!         CALL Laplacia(model, RESHAPE(x3(k_layer,:), (/1, SIZE(x3, 2)/)), M_x)
!     END IF

!     ! Taylor Expansion Test
!     DO l = 1, num_lambdas
!         lambda = lambda_values(l)

!         !  M(x + lambda * dx)
!         IF (operator_name == "Divergen") THEN
!             CALL Divergen(model, RESHAPE(x1(k_layer,:) + lambda * dx1(k_layer,:), (/1, SIZE(x1, 2)/)), RESHAPE(x2(k_layer,:) + lambda * dx2(k_layer,:), (/1, SIZE(x2, 2)/)), M_x_perturb)
!             CALL Divergen_TLM(model_tlm, RESHAPE(dx1(k_layer,:), (/1, SIZE(dx1, 2)/)), RESHAPE(dx2(k_layer,:), (/1, SIZE(dx2, 2)/)), M_tlm)
!         ELSE IF (operator_name == "Jacobian") THEN
!             CALL Jacobian(model, RESHAPE(x1(k_layer,:) + lambda * dx1(k_layer,:), (/1, SIZE(x1, 2)/)), RESHAPE(x2(k_layer,:) + lambda * dx2(k_layer,:), (/1, SIZE(x2, 2)/)), M_x_perturb)
!             ! CALL Jacobian_TLM(model_tlm, RESHAPE(dx1(k_layer,:), (/1, SIZE(dx1, 2)/)), RESHAPE(dx2(k_layer,:), (/1, SIZE(dx2, 2)/)), M_tlm)
!         ELSE IF (operator_name == "Laplacia") THEN
!             CALL Laplacia(model, RESHAPE(x3(k_layer,:) + lambda * dx3(k_layer,:), (/1, SIZE(x3, 2)/)), M_x_perturb)
!             CALL Laplacia_TLM(model_tlm, RESHAPE(dx3(k_layer,:), (/1, SIZE(dx3, 2)/)), M_tlm)
!         END IF

!         ! Calculate norms
!         norm_diff = NORM2_1D(RESHAPE(M_x_perturb - M_x, [SIZE(M_x_perturb)]))
!         norm_tlm = NORM2_1D(RESHAPE(lambda * M_tlm, [SIZE(M_tlm)]))
!         print *, norm_tlm, norm_diff
!         !  Taylor Expansion Ratio
!         IF (norm_tlm > 0.0_r_kind) THEN
!             taylor_ratio = norm_diff / norm_tlm
!         ELSE
!             taylor_ratio = IEEE_VALUE(taylor_ratio, IEEE_QUIET_NAN)  ! 如果 norm_tlm 为零，则分配 NaN
!         END IF

!
!         print *,  lambda, taylor_ratio
!         WRITE(unit, '(E15.6, 2X, E15.6)') lambda, taylor_ratio
!     END DO

!
!     CLOSE(unit)
! END SUBROUTINE Test_TLM_Taylor

  FUNCTION NORM2_1D(array) RESULT(norm)
    REAL(r_kind), INTENT(IN) :: array(:)
    REAL(r_kind) :: norm

    ! Compute the L2 norm of the 1D array
    norm = SQRT(SUM(array**2))
  END FUNCTION NORM2_1D

  FUNCTION NORM2_2D(array) RESULT(norm)
    REAL(r_kind), INTENT(IN) :: array(:, :)
    REAL(r_kind) :: norm

    ! Compute the L2 norm of the 2D array
    norm = SQRT(SUM(array**2))
  END FUNCTION NORM2_2D

END PROGRAM test_gzm_tl_ad_consistency
