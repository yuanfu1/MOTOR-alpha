PROGRAM unit_test_gzm_tlm_ad_points
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  USE Verification_m
  IMPLICIT NONE

  TYPE(SingleGrid_t), POINTER :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  INTEGER(i_kind) :: i, j
  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), d_oprand(:, :), adj_oprand(:, :)
  REAL(r_kind) :: epsilon

  CHARACTER(LEN=1024) :: configFile, ncOutputFile

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

  ! Initialize gzm, gzm_tlm and gzm_adj
  gzm%sg => sg
  CALL gzm_tlm%initialize(sg)
  gzm_adj%sg => sg

  ! Ensure arrays are allocated and initialized properly
  ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
  ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
  ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
  ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
  ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell))

  ! Initialize input data
  opr_left = 0.05D0
  opr_right = 0.05D0
  oprand = 0.05D0
  d_value = 0.0D0
  VALUE = 0.0D0
  adj_value = 0.0D0
  DO i = 1, sg%vLevel
    d_opr_left(i, :) = (DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))) / 1.0D5 / 2.0D0
    d_opr_right(i, :) = (DCOS(sg%cell_cntr(1, :)) + DSIN(sg%cell_cntr(2, :))) / 1.0D5 / 2.0D0
    d_oprand(i, :) = DSIN(sg%cell_cntr(1, :)) / 1000D0
    adj_opr_left(i, :) = 0.0D0
    adj_opr_right(i, :) = 0.0D0
    adj_oprand(i, :) = 0.0D0
  END DO

  opr_left = d_opr_left * 1.0D5
  opr_right = d_opr_right * 1.0D5
  oprand = d_oprand

  epsilon = 1E-6_R_KIND  ! Small perturbation

  ! Test Divergence for each grid point
  DO j = 1, sg%num_icell
    CALL gzm%Divergen(opr_left(:, j), opr_right(:, j), VALUE(:, j))
    CALL gzm_tlm%Divergen_TLM(opr_left(:, j), opr_right(:, j), VALUE(:, j), d_opr_left(:, j), d_opr_right(:, j), d_value(:, j))
    pert_opr_left(:, j) = opr_left(:, j) + epsilon * d_opr_left(:, j)
    pert_opr_right(:, j) = opr_right(:, j) + epsilon * d_opr_right(:, j)
    CALL gzm%Divergen(pert_opr_left(:, j), pert_opr_right(:, j), pert_value(:, j))
    CALL Verify_TL_FiniteDifference(VALUE(:, j), d_value(:, j), pert_value(:, j), epsilon, sg, "Divergence")

    adj_value(:, j) = d_value(:, j)
    CALL gzm_adj%Divergen_AD(opr_left(:, j), opr_right(:, j), VALUE(:, j), adj_opr_left(:, j), adj_opr_right(:, j), adj_value(:, j))
    CALL Verify_AD_FiniteDifference(VALUE(:, j), adj_value(:, j), pert_value(:, j), epsilon, sg, "Divergence")

    CALL Verify_TLM_AD(d_opr_left(:, j), adj_opr_left(:, j), d_opr_right(:, j), adj_opr_right(:, j), sg, "Divergence")
    CALL Verify_Model_Equations(d_value(:, j), adj_value(:, j), "Divergence")
    CALL Verify_Spectral_Radius(d_value(:, j), "Divergence")
  END DO

  ! Test Jacobian for each grid point
  DO j = 1, sg%num_icell
    CALL gzm%Jacobian(opr_left(:, j), opr_right(:, j), VALUE(:, j))
    CALL gzm_tlm%Jacobian_TLM(opr_left(:, j), opr_right(:, j), VALUE(:, j), d_opr_left(:, j), d_opr_right(:, j), d_value(:, j))
    pert_opr_left(:, j) = opr_left(:, j) + epsilon * d_opr_left(:, j)
    pert_opr_right(:, j) = opr_right(:, j) + epsilon * d_opr_right(:, j)
    CALL gzm%Jacobian(pert_opr_left(:, j), pert_opr_right(:, j), pert_value(:, j))
    CALL Verify_TL_FiniteDifference(VALUE(:, j), d_value(:, j), pert_value(:, j), epsilon, sg, "Jacobian")

    adj_value(:, j) = d_value(:, j)
    CALL gzm_adj%Jacobian_AD(opr_left(:, j), opr_right(:, j), VALUE(:, j), adj_opr_left(:, j), adj_opr_right(:, j), adj_value(:, j))
    CALL Verify_AD_FiniteDifference(VALUE(:, j), adj_value(:, j), pert_value(:, j), epsilon, sg, "Jacobian")

    CALL Verify_TLM_AD(d_opr_left(:, j), adj_opr_left(:, j), d_opr_right(:, j), adj_opr_right(:, j), sg, "Jacobian")
    CALL Verify_Model_Equations(d_value(:, j), adj_value(:, j), "Jacobian")
    CALL Verify_Spectral_Radius(d_value(:, j), "Jacobian")
  END DO

  ! Test Laplacia for each grid point
  DO j = 1, sg%num_icell
    CALL gzm%Laplacia(oprand(:, j), VALUE(:, j))
    CALL gzm_tlm%Laplacia_TLM(oprand(:, j), VALUE(:, j), d_oprand(:, j), d_value(:, j))
    pert_oprand(:, j) = oprand(:, j) + epsilon * d_oprand(:, j)
    CALL gzm%Laplacia(pert_oprand(:, j), pert_value(:, j))
    CALL Verify_TL_FiniteDifference(VALUE(:, j), d_value(:, j), pert_value(:, j), epsilon, sg, "Laplacia")

    adj_value(:, j) = d_value(:, j)
    CALL gzm_adj%Laplacia_AD(oprand(:, j), VALUE(:, j), adj_oprand(:, j), adj_value(:, j))
    CALL Verify_AD_FiniteDifference(VALUE(:, j), adj_value(:, j), pert_value(:, j), epsilon, sg, "Laplacia")

    CALL Verify_TLM_AD(d_oprand(:, j), adj_oprand(:, j), d_oprand(:, j), adj_oprand(:, j), sg, "Laplacia")
    CALL Verify_Model_Equations(d_value(:, j), adj_value(:, j), "Laplacia")
    CALL Verify_Spectral_Radius(d_value(:, j), "Laplacia")
  END DO

  ! Clean up
  DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, pert_opr_left, pert_opr_right, pert_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand)

  CALL mpddGlob%finalize
END PROGRAM unit_test_gzm_tlm_ad_points
