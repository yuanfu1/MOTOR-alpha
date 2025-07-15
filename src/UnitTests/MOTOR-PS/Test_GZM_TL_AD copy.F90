PROGRAM Test_GZM_TL_AD
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_tlm, ONLY: gzm_tl_t
  USE gzm_adj, ONLY: gzm_ad_t

  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(gzm_tl_t) :: gzm_tl
  TYPE(gzm_ad_t) :: gzm_ad

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  INTEGER(i_kind) :: i, j, k, size_3
  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind) :: inner_product_1, inner_product_2, epsilon
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGZM_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testGZM.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ! Associate SingleGrid with gzm_tl and gzm_ad
  ASSOCIATE (sg => geometry%mg%sg(5))

    sg%bdy_type = sg%cell_type

    gzm_tl = gzm_tl_t(sg)
    gzm_ad = gzm_ad_t(sg)

    ! Ensure arrays are allocated and initialized properly
    ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), &
              VALUE(sg%vLevel, sg%num_icell), d_opr_left(sg%vLevel, sg%num_icell), &
              d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell), &
              adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))

    ! Initialize input and perturbations with random numbers
    opr_left = 0.05D0 !CALL RANDOM_NUMBER(opr_left)
    opr_right = 0.05D0 !CALL RANDOM_NUMBER(opr_right)
    !  d_opr_left = 0.05D0 !CALL RANDOM_NUMBER(opr_left)
    !  d_opr_right = 0.05D0 !CALL RANDOM_NUMBER(opr_right)
    d_value = 0.0D0
    DO i = 1, sg%vLevel
      d_opr_left(i, :) = DSIN(sg%cell_cntr(1, :)) * 1000
      d_opr_right(i, :) = DSIN(sg%cell_cntr(1, :)) * 1000
    END DO

    ! CALL RANDOM_NUMBER(d_opr_left)
    ! CALL RANDOM_NUMBER(d_opr_right)

    ! Initialize value with zeros
    VALUE = 0.0_R_KIND

    ! Call TL model
    CALL gzm_tl%Divergen_tl(opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)

    ! Initialize adjoint variables
    adj_value = d_value
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND

    ! Call AD model
    CALL gzm_ad%Divergen_ad(opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)

    ! Compute inner products
    inner_product_1 = SUM(d_value * adj_value)
    inner_product_2 = SUM(d_opr_right * adj_opr_right) + SUM(d_opr_left * adj_opr_left)

    ! inner_product_1 = SUM(d_opr_left) + SUM(d_opr_right)
    ! inner_product_2 = SUM(adj_opr_left)+SUM(adj_opr_right)

    PRINT *, "Inner product d_opr . adj_opr:", inner_product_1
    PRINT *, "Inner product d_value . adj_value:", inner_product_2

    epsilon = ABS(inner_product_1 - inner_product_2) / MAX(ABS(inner_product_1), ABS(inner_product_2))
    PRINT *, "Epsilon (should be close to zero):", epsilon

    ! Debug output
    ! PRINT *, "Debug: Test_GZM_TL_AD"
    ! PRINT *, "opr_left:", opr_left
    ! PRINT *, "opr_right:", opr_right
    ! PRINT *, "value:", value
    ! PRINT *, "d_opr_left:", d_opr_left
    ! PRINT *, "d_opr_right:", d_opr_right
    ! PRINT *, "d_value:", d_value
    ! PRINT *, "adj_opr_left:", adj_opr_left
    ! PRINT *, "adj_opr_right:", adj_opr_right
    ! PRINT *, "adj_value:", adj_value

    IF (epsilon < 1.0E-10) THEN
      PRINT *, "TL-AD consistency check passed."
    ELSE
      PRINT *, "TL-AD consistency check failed."
    END IF

    DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, adj_opr_left, adj_opr_right, adj_value)
  END ASSOCIATE

  CALL mpddGlob%finalize

END PROGRAM Test_GZM_TL_AD
