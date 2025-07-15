PROGRAM Test_CalVerDer_TL_AD
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE CalVerDer_TL_m, ONLY: CalVerDer_TL_t
  USE CalVerDer_AD_m, ONLY: CalVerDer_AD_t
  USE Verification_m

  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(CalVerDer_TL_t) :: calverder_tl
  TYPE(CalVerDer_AD_t) :: calverder_ad

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  INTEGER(i_kind) :: i, j, k, size_3
  REAL(r_kind), ALLOCATABLE :: x(:, :), y(:, :), tl_y(:, :), adj_x(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_x(:, :), pert_tl_y(:, :), pert_adj_x(:, :)
  REAL(r_kind) :: inner_product_1, inner_product_2, epsilon
  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testCalVerDer_TLAD.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testCalVerDer.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)  ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(5))
    calverder_tl%sg => sg
    calverder_ad%sg => sg

    ALLOCATE (x(sg%vLevel, sg%num_cell), y(sg%vLevel, sg%num_cell), &
              tl_y(sg%vLevel, sg%num_cell), adj_x(sg%vLevel, sg%num_cell), &
              pert_x(sg%vLevel, sg%num_cell), pert_tl_y(sg%vLevel, sg%num_cell), pert_adj_x(sg%vLevel, sg%num_cell))

    x = 0.005D0
    y = 0.005D0

    ! Call TL model
    CALL calverder_tl%FirstOrder_TL(x, x, tl_y)

    ! Call AD model
    CALL calverder_ad%FirstOrder_AD(y, y, adj_x)

    ! Compute inner products
    inner_product_1 = SUM(x * adj_x)
    inner_product_2 = SUM(y * tl_y)

    PRINT *, "Inner product x . adj_x:", inner_product_1
    PRINT *, "Inner product y . tl_y:", inner_product_2

    epsilon = ABS(inner_product_1 / inner_product_2)
    PRINT *, "Epsilon (should be close to one):", epsilon

    IF (ABS(epsilon - 1.0D0) < 1.0E-10) THEN
      PRINT *, "TL-AD consistency check passed."
    ELSE
      PRINT *, "TL-AD consistency check failed."
    END IF

    ! Perturbation for double differencing method
    epsilon = 1.0E-60_R_KIND
    pert_x = x + epsilon * tl_y

    ! Call TL model with perturbation
    CALL calverder_tl%FirstOrder_TL(pert_x, pert_x, pert_tl_y)

    ! Verify TL using double differencing method
    CALL Verify_TL_DoubleDifferencing(tl_y, pert_tl_y, tl_y, epsilon, sg)

    ! Perturbation for adjoint double differencing method
    pert_adj_x = y + epsilon * adj_x

    ! Call AD model with perturbation
    CALL calverder_ad%FirstOrder_AD(pert_adj_x, pert_adj_x, pert_adj_x)

    ! Verify AD using double differencing method
    CALL Verify_AD_DoubleDifferencing(adj_x, pert_adj_x, adj_x, epsilon, sg)

    ! Verify TL-AD consistency
    CALL Verify_TL_AD_Consistency(x, adj_x, y, tl_y, sg)

    DEALLOCATE (x, y, tl_y, adj_x, pert_x, pert_tl_y, pert_adj_x)
  END ASSOCIATE

  CALL mpddGlob%finalize

END PROGRAM Test_CalVerDer_TL_AD
