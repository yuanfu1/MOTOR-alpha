PROGRAM test_gzm_tlm_numericalgradient
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES

  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: num_grad(:, :)

  REAL(r_kind) :: epsilon, error_norm
  REAL(r_kind) :: lower_bound, upper_bound
  INTEGER(i_kind) :: i, j
  CHARACTER(LEN=1024) :: configFile, ncOutputFile
  LOGICAL :: ieee_support, invalid_flag, overflow_flag, divide_by_zero_flag, underflow_flag

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ! Associate SingleGrid with gzm and gzm_tlm
  sg => geometry%mg%sg(6)
  sg%bdy_type = sg%cell_type

  gzm%sg => sg
  gzm_tlm%sg => sg

  ! Ensure arrays are allocated and initialized properly
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (num_grad(SIZE(VALUE, 1), SIZE(VALUE, 2)))

  ! small perturbation
  epsilon = 1.0E-6

  lower_bound = 0.1_R_KIND
  upper_bound = 1.0_R_KIND

  CALL RANDOM_NUMBER(opr_left)
  CALL RANDOM_NUMBER(opr_right)
  opr_left = lower_bound + (upper_bound - lower_bound) * opr_left
  opr_right = lower_bound + (upper_bound - lower_bound) * opr_right

  PRINT *, 'Min value of opr_left:', MINVAL(opr_left)
  PRINT *, 'Max value of opr_left:', MAXVAL(opr_left)
  PRINT *, 'Mean value of opr_left:', SUM(opr_left) / SIZE(opr_left)

  PRINT *, 'Min value of opr_right:', MINVAL(opr_right)
  PRINT *, 'Max value of opr_right:', MAXVAL(opr_right)
  PRINT *, 'Mean value of opr_right:', SUM(opr_right) / SIZE(opr_right)

  ! Initialize perturbations
  lower_bound = 0.001_R_KIND
  upper_bound = 0.01_R_KIND
  CALL RANDOM_NUMBER(d_opr_left)
  CALL RANDOM_NUMBER(d_opr_right)
  d_opr_left = lower_bound + (upper_bound - lower_bound) * d_opr_left
  d_opr_right = lower_bound + (upper_bound - lower_bound) * d_opr_right

  PRINT *, 'Min value of d_opr_left:', MINVAL(d_opr_left)
  PRINT *, 'Max value of d_opr_left:', MAXVAL(d_opr_left)
  PRINT *, 'Mean value of d_opr_left:', SUM(d_opr_left) / SIZE(d_opr_left)

  PRINT *, 'Min value of d_opr_right:', MINVAL(d_opr_right)
  PRINT *, 'Max value of d_opr_right:', MAXVAL(d_opr_right)
  PRINT *, 'Mean value of d_opr_right:', SUM(d_opr_right) / SIZE(d_opr_right)

  ! Check for IEEE floating-point exception support
  ieee_support = IEEE_SUPPORT_DATATYPE(opr_left)

  IF (ieee_support) THEN
    PRINT *, "IEEE floating-point exception support is available."

    ! Disable specific exceptions
    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_UNDERFLOW, .FALSE.)
  ELSE
    PRINT *, "No IEEE floating-point exception support available."
  END IF

  ! Call original model (Forward)
  CALL Divergen(gzm, opr_left, opr_right, VALUE)
  PRINT *, 'Forward Model Output (value):'
  PRINT *, VALUE

  ! Check for any IEEE floating-point exceptions after the computation
  IF (ieee_support) THEN
    CALL IEEE_GET_FLAG(IEEE_INVALID, invalid_flag)
    CALL IEEE_GET_FLAG(IEEE_OVERFLOW, overflow_flag)
    CALL IEEE_GET_FLAG(IEEE_DIVIDE_BY_ZERO, divide_by_zero_flag)
    CALL IEEE_GET_FLAG(IEEE_UNDERFLOW, underflow_flag)

    IF (invalid_flag) THEN
      PRINT *, 'IEEE_INVALID_FLAG detected during forward model calculation.'
    END IF
    IF (overflow_flag) THEN
      PRINT *, 'IEEE_OVERFLOW detected during forward model calculation.'
    END IF
    IF (divide_by_zero_flag) THEN
      PRINT *, 'IEEE_DIVIDE_BY_ZERO detected during forward model calculation.'
    END IF
    IF (underflow_flag) THEN
      PRINT *, 'IEEE_UNDERFLOW detected during forward model calculation.'
    END IF
  END IF

  ! Call tangent linear model (TLM)
  CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)

  ! Test the numerical gradient for Divergen operator
  CALL Numerical_Gradient(gzm, opr_left, opr_right, VALUE, epsilon, num_grad, "Divergen")
  error_norm = NORM2(d_value - num_grad)
  CALL Write_Result("Divergen", error_norm, "tlm_unit_test_result.txt")

  ! Test Jacobian operator
  CALL Jacobian(gzm, opr_left, opr_right, VALUE)
  CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
  CALL Numerical_Gradient(gzm, opr_left, opr_right, VALUE, epsilon, num_grad, "Jacobian")
  error_norm = NORM2(d_value - num_grad)
  CALL Write_Result("Jacobian", error_norm, "tlm_unit_test_result.txt")

  ! Test Laplacia operator
  CALL Laplacia(gzm, opr_left, VALUE)
  CALL Laplacia_TLM(gzm_tlm, opr_left, VALUE, d_opr_left, d_value)
  CALL Numerical_Gradient(gzm, opr_left, opr_right, VALUE, epsilon, num_grad, "Laplacia")
  error_norm = NORM2(d_value - num_grad)
  CALL Write_Result("Laplacia", error_norm, "tlm_unit_test_result.txt")

  DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, num_grad)

CONTAINS

  SUBROUTINE Numerical_Gradient(gzm, opr_left, opr_right, VALUE, epsilon, num_grad, operator_name)
    CLASS(gzm_t) :: gzm
    REAL(r_kind), INTENT(IN) :: opr_left(:, :), opr_right(:, :), VALUE(:, :), epsilon
    REAL(r_kind), INTENT(OUT) :: num_grad(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_value(:, :)
    INTEGER :: i, j

    ! 分配内存
    ALLOCATE (pert_opr_left(SIZE(opr_left, 1), SIZE(opr_left, 2)), pert_value(SIZE(VALUE, 1), SIZE(VALUE, 2)))

    ! 初始化 num_grad
    num_grad = 0.0_R_KIND

    ! 计算数值梯度
    DO i = 1, SIZE(opr_left, 1)
      DO j = 1, SIZE(opr_left, 2)
        pert_opr_left = opr_left
        pert_opr_left(i, j) = opr_left(i, j) + epsilon

        IF (operator_name == "Divergen") THEN
          CALL Divergen(gzm, pert_opr_left, opr_right, pert_value)
        ELSE IF (operator_name == "Jacobian") THEN
          CALL Jacobian(gzm, pert_opr_left, opr_right, pert_value)
        ELSE IF (operator_name == "Laplacia") THEN
          CALL Laplacia(gzm, pert_opr_left, pert_value)
        END IF

        num_grad(i, j) = (pert_value(i, j) - VALUE(i, j)) / epsilon

        ! 打印当前计算状态
        IF (MOD(i, 10) == 0 .AND. MOD(j, 1000) == 0) THEN
          PRINT *, 'Progress: i =', i, 'j =', j, 'Current num_grad(i,j) =', num_grad(i, j)
        END IF
      END DO
    END DO

    ! 清理内存
    DEALLOCATE (pert_opr_left, pert_value)
  END SUBROUTINE Numerical_Gradient

  FUNCTION NORM2(arr) RESULT(norm)
    REAL(r_kind), INTENT(IN) :: arr(:, :)
    REAL(r_kind) :: norm
    norm = SQRT(SUM(arr**2))
  END FUNCTION NORM2

  SUBROUTINE Write_Result(operator_name, error_norm, file_name)
    CHARACTER(LEN=*), INTENT(IN) :: operator_name
    REAL(r_kind), INTENT(IN) :: error_norm
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER :: unit
    unit = 10

    OPEN (UNIT=unit, FILE=file_name, STATUS='UNKNOWN', ACTION='WRITE')
    WRITE (unit, *) operator_name//' TLM Unit Test'
    WRITE (unit, *) 'Error Norm:', error_norm
    IF (error_norm < 1.0E-5) THEN
      WRITE (unit, *) 'Test PASSED.'
    ELSE
      WRITE (unit, *) 'Test FAILED.'
    END IF
    CLOSE (unit)
  END SUBROUTINE Write_Result

END PROGRAM test_gzm_tlm_numericalgradient
