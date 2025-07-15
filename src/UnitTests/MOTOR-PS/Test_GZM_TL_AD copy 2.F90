PROGRAM Test_GZM_TLM_ADJ
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t, Divergen, Jacobian, Laplacia
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD

  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(gzm_t) :: gzm
  TYPE(gzm_tlm_t) :: gzm_tlm
  TYPE(gzm_adj_t) :: gzm_adj

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  INTEGER(i_kind) :: i, j, k
  REAL(r_kind), ALLOCATABLE :: opr_left(:, :), opr_right(:, :), VALUE(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_opr_left(:, :), pert_opr_right(:, :), pert_value(:, :)
  REAL(r_kind), ALLOCATABLE :: d_opr_left(:, :), d_opr_right(:, :), d_value(:, :)
  REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_value(:, :)
  REAL(r_kind), ALLOCATABLE :: oprand(:, :), d_oprand(:, :), adj_oprand(:, :), pert_oprand(:, :)
  REAL(r_kind) :: inner_product_tl, inner_product_ad, epsilon

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
  ASSOCIATE (sg => geometry%mg%sg(6))

    sg%bdy_type = sg%cell_type

    ! Initialize gzm, gzm_tlm and gzm_adj
    gzm%sg => sg
    gzm_tlm%sg => sg
    gzm_adj%sg => sg

    ! Ensure arrays are allocated and initialized properly
    ALLOCATE (opr_left(sg%vLevel, sg%num_icell), opr_right(sg%vLevel, sg%num_icell), VALUE(sg%vLevel, sg%num_icell))
    ALLOCATE (d_opr_left(sg%vLevel, sg%num_icell), d_opr_right(sg%vLevel, sg%num_icell), d_value(sg%vLevel, sg%num_icell))
    ALLOCATE (pert_opr_left(sg%vLevel, sg%num_icell), pert_opr_right(sg%vLevel, sg%num_icell), pert_value(sg%vLevel, sg%num_icell))
    ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_value(sg%vLevel, sg%num_icell))
    ALLOCATE (oprand(sg%vLevel, sg%num_icell), d_oprand(sg%vLevel, sg%num_icell), adj_oprand(sg%vLevel, sg%num_icell), pert_oprand(sg%vLevel, sg%num_icell))

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

    epsilon = 0.5E0_R_KIND  ! Small perturbation

    ! Test Divergen
    CALL gzm%Divergen(opr_left, opr_right, VALUE)
    CALL Divergen_TLM(gzm_tlm, opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
    pert_opr_left = opr_left + epsilon * d_opr_left
    pert_opr_right = opr_right + epsilon * d_opr_right
    CALL Divergen(gzm, pert_opr_left, pert_opr_right, pert_value)
    CALL Verify_Divergen_TLM(VALUE, d_value, pert_value, epsilon, sg)

    ! Initialize adjoint variables for Divergen
    adj_value = d_value
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    CALL Divergen_AD(gzm_adj, opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)

    ! Verify Divergen (TL-AD consistency)
    CALL Verify_Divergen_TLM_AD(d_opr_left, d_opr_right, adj_opr_left, adj_opr_right, sg)

    ! Test Jacobian
    CALL Jacobian(gzm, opr_left, opr_right, VALUE)
    PRINT *, 'test gzm data: ', MAXVAL(VALUE)
    CALL Jacobian_TLM(gzm_tlm, opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value)
    PRINT *, 'test gzm_tlm data: ', MAXVAL(d_value)
    pert_opr_left = opr_left + epsilon * d_opr_left
    pert_opr_right = opr_right + epsilon * d_opr_right
    CALL Jacobian(gzm, pert_opr_left, pert_opr_right, pert_value)
    PRINT *, 'pert value is: ', MAXVAL(pert_value)
    CALL Verify_Jacobian_TLM(VALUE, d_value, pert_value, epsilon, sg)

    ! Initialize adjoint variables for Jacobian
    adj_value = d_value
    adj_opr_left = 0.0_R_KIND
    adj_opr_right = 0.0_R_KIND
    CALL Jacobian_AD(gzm_adj, opr_left, opr_right, VALUE, adj_opr_left, adj_opr_right, adj_value)

    ! Verify Jacobian (TL-AD consistency)
    CALL Verify_Jacobian_TLM_AD(d_opr_left, d_opr_right, adj_opr_left, adj_opr_right, sg)

    ! Test Laplacia
    CALL Laplacia(gzm, oprand, VALUE)
    CALL Laplacia_TLM(gzm_tlm, oprand, VALUE, d_oprand, d_value)
    pert_oprand = oprand + epsilon * d_oprand
    CALL Laplacia(gzm, pert_oprand, pert_value)
    CALL Verify_Laplacia_TLM(VALUE, d_value, pert_value, epsilon, sg)

    ! Initialize adjoint variables for Laplacia
    adj_value = d_value
    adj_oprand = 0.0_R_KIND
    CALL Laplacia_AD(gzm_adj, oprand, VALUE, adj_oprand, adj_value)

    ! Verify Laplacia (TL-AD consistency)
    CALL Verify_Laplacia_TLM_AD(d_oprand, adj_oprand, sg)

    ! Clean up
    DEALLOCATE (opr_left, opr_right, VALUE, d_opr_left, d_opr_right, d_value, pert_opr_left, pert_opr_right, pert_value, adj_opr_left, adj_opr_right, adj_value, oprand, d_oprand, adj_oprand, pert_oprand)
  END ASSOCIATE

  CALL mpddGlob%finalize

CONTAINS

  SUBROUTINE Verify_Divergen_TLM(VALUE, d_value, pert_value, epsilon, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), d_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
          avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (TLM Divergen):', max_error
    PRINT *, 'Avg Error (TLM Divergen):', avg_error

    IF (max_error < 1E-5) THEN
      PRINT *, 'TLM Divergen verification passed.'
    ELSE
      PRINT *, 'TLM Divergen verification failed.'
    END IF
  END SUBROUTINE Verify_Divergen_TLM

  SUBROUTINE Verify_Jacobian_TLM(VALUE, d_value, pert_value, epsilon, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), d_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
          avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (TLM Jacobian):', max_error
    PRINT *, 'Avg Error (TLM Jacobian):', avg_error

    IF (max_error < 1E-5) THEN
      PRINT *, 'TLM Jacobian verification passed.'
    ELSE
      PRINT *, 'TLM Jacobian verification failed.'
    END IF
  END SUBROUTINE Verify_Jacobian_TLM

  SUBROUTINE Verify_Laplacia_TLM(VALUE, d_value, pert_value, epsilon, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: VALUE(:, :), d_value(:, :), pert_value(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon))
          avg_error = avg_error + ABS(d_value(k, j) - (pert_value(k, j) - VALUE(k, j)) / epsilon)
          n = n + 1
        END IF
      END DO
    END DO
    avg_error = avg_error / n

    PRINT *, 'Max Error (TLM Laplacia):', max_error
    PRINT *, 'Avg Error (TLM Laplacia):', avg_error

    IF (max_error < 1E-5) THEN
      PRINT *, 'TLM Laplacia verification passed.'
    ELSE
      PRINT *, 'TLM Laplacia verification failed.'
    END IF
  END SUBROUTINE Verify_Laplacia_TLM

  SUBROUTINE Verify_Divergen_TLM_AD(d_opr_left, d_opr_right, adj_opr_left, adj_opr_right, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), d_opr_right(:, :), adj_opr_left(:, :), adj_opr_right(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: inner_product_tl, inner_product_ad
    INTEGER :: k, j

    inner_product_tl = 0.0
    inner_product_ad = 0.0

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_opr_left(k, j) * adj_opr_left(k, j) + d_opr_right(k, j) * adj_opr_right(k, j)
        inner_product_ad = inner_product_ad + adj_opr_left(k, j) * d_opr_left(k, j) + adj_opr_right(k, j) * d_opr_right(k, j)
      END DO
    END DO

    PRINT *, 'Inner Product (TL Divergen):', inner_product_tl
    PRINT *, 'Inner Product (AD Divergen):', inner_product_ad

    IF (ABS(inner_product_tl - inner_product_ad) < 1E-5) THEN
      PRINT *, 'TL-AD Divergen verification passed.'
    ELSE
      PRINT *, 'TL-AD Divergen verification failed.'
    END IF
  END SUBROUTINE Verify_Divergen_TLM_AD

  SUBROUTINE Verify_Jacobian_TLM_AD(d_opr_left, d_opr_right, adj_opr_left, adj_opr_right, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: d_opr_left(:, :), d_opr_right(:, :), adj_opr_left(:, :), adj_opr_right(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: inner_product_tl, inner_product_ad
    INTEGER :: k, j

    inner_product_tl = 0.0
    inner_product_ad = 0.0

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_opr_left(k, j) * adj_opr_left(k, j) + d_opr_right(k, j) * adj_opr_right(k, j)
        inner_product_ad = inner_product_ad + adj_opr_left(k, j) * d_opr_left(k, j) + adj_opr_right(k, j) * d_opr_right(k, j)
      END DO
    END DO

    PRINT *, 'Inner Product (TL Jacobian):', inner_product_tl
    PRINT *, 'Inner Product (AD Jacobian):', inner_product_ad

    IF (ABS(inner_product_tl - inner_product_ad) < 1E-5) THEN
      PRINT *, 'TL-AD Jacobian verification passed.'
    ELSE
      PRINT *, 'TL-AD Jacobian verification failed.'
    END IF
  END SUBROUTINE Verify_Jacobian_TLM_AD

  SUBROUTINE Verify_Laplacia_TLM_AD(d_oprand, adj_oprand, sg)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: d_oprand(:, :), adj_oprand(:, :)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: inner_product_tl, inner_product_ad
    INTEGER :: k, j

    inner_product_tl = 0.0
    inner_product_ad = 0.0

    DO k = 1, sg%vLevel
      DO j = 1, sg%num_icell
        inner_product_tl = inner_product_tl + d_oprand(k, j) * adj_oprand(k, j)
        inner_product_ad = inner_product_ad + adj_oprand(k, j) * d_oprand(k, j)
      END DO
    END DO

    PRINT *, 'Inner Product (TL Laplacia):', inner_product_tl
    PRINT *, 'Inner Product (AD Laplacia):', inner_product_ad

    IF (ABS(inner_product_tl - inner_product_ad) < 1E-5) THEN
      PRINT *, 'TL-AD Laplacia verification passed.'
    ELSE
      PRINT *, 'TL-AD Laplacia verification failed.'
    END IF
  END SUBROUTINE Verify_Laplacia_TLM_AD

END PROGRAM Test_GZM_TLM_ADJ
