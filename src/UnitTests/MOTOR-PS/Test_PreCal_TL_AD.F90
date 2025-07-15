PROGRAM Test_PreCal_TL_AD
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE PreCal_TLM_m, ONLY: PreCal_TLM_t
  USE PreCal_AD_m, ONLY: PreCal_AD_t
  USE State_m, ONLY: State_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t

  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(PreCal_TLM_t) :: precal_tlm
  TYPE(PreCal_AD_t) :: precal_ad
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: X, X_tl, X_adj

  INTEGER(i_kind) :: i, j, k, size_3
  REAL(r_kind), ALLOCATABLE :: x(:, :), y(:, :), tl_y(:, :), adj_x(:, :)
  REAL(r_kind) :: inner_product_1, inner_product_2, epsilon
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testPreCal_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testPreCal.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(5))
    CALL precal_tlm%initialize(sg)
    CALL precal_ad%initialize(sg)
    CALL X%initialize(configFile, sg)
    CALL X_tl%initialize(configFile, sg)
    CALL X_adj%initialize(configFile, sg)

    ALLOCATE (x(sg%vLevel, sg%num_cell), y(sg%vLevel, sg%num_cell), &
              tl_y(sg%vLevel, sg%num_cell), adj_x(sg%vLevel, sg%num_cell))

    x = 0.005D0
    y = 0.005D0

    ! Call TL model
    CALL precal_tlm%PreCals_TLM(X_tl, X)

    ! Call AD model
    CALL precal_ad%PreCals_AD(X_adj, X)

    ! Compute inner products
    inner_product_1 = SUM(x * adj_x)
    inner_product_2 = SUM(y * tl_y)

    PRINT *, "Inner product x . adj_x:", inner_product_1
    PRINT *, "Inner product y . tl_y:", inner_product_2

    epsilon = ABS(inner_product_1 / inner_product_2)
    PRINT *, "Epsilon (should be close to one):", epsilon

    IF ((epsilon - 1.0D0) < 1.0E-10) THEN
      PRINT *, "TL-AD consistency check passed."
    ELSE
      PRINT *, "TL-AD consistency check failed."
    END IF

    DEALLOCATE (x, y, tl_y, adj_x)
  END ASSOCIATE

  CALL mpddGlob%finalize

END PROGRAM Test_PreCal_TL_AD
