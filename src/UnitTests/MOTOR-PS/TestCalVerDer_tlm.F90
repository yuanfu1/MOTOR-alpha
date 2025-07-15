PROGRAM TestCalVerDer_tlm
  USE YAMLRead_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE CalVerDer_tlm_m, ONLY: CalVerDer_tlm_t
  USE mpddGlob_m, ONLY: mpddGlob_t  ! 确保引用了 mpddGlob_t 类型

  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(CalVerDer_t) :: calVerDer
  TYPE(CalVerDer_tlm_t) :: calVerDer_tlm
  TYPE(mpddGlob_t), TARGET :: mpddGlob

  INTEGER(i_kind) :: i, j, k
  REAL(r_kind), ALLOCATABLE :: A(:, :), parA_sigma(:, :)
  REAL(r_kind), ALLOCATABLE :: dA(:, :), dparA_sigma(:, :)
  REAL(r_kind), ALLOCATABLE :: pert_A(:, :), pert_parA_sigma(:, :)

  REAL(r_kind) :: epsilon
  CHARACTER(LEN=1024) :: configFile
  INTEGER(i_kind) :: gLevel, vLevel, tSlots, rank, group
  INTEGER(i_kind), DIMENSION(2) :: proc_layout
  REAL(r_kind) :: start_time, end_time

  ! Initialize
  configFile = 'path_to_config_file.yaml'
  CALL mpddGlob%initialize()

  gLevel = 1
  vLevel = 10
  tSlots = 1
  rank = 0
  group = 0
  proc_layout = (/1, 1/)
  start_time = 0.0_R_KIND
  end_time = 1.0_R_KIND

  CALL sg%initialize(mpddGlob, gLevel, vLevel, group, rank, proc_layout, tSlots, start_time, end_time, configFile)

  calVerDer%sg => sg
  calVerDer_tlm%sg => sg

  ! Ensure arrays are allocated and initialized properly
  ALLOCATE (A(sg%vLevel, sg%num_cell))
  ALLOCATE (parA_sigma(sg%vLevel, sg%num_cell))
  ALLOCATE (dA(sg%vLevel, sg%num_cell))
  ALLOCATE (dparA_sigma(sg%vLevel, sg%num_cell))
  ALLOCATE (pert_A(sg%vLevel, sg%num_cell))
  ALLOCATE (pert_parA_sigma(sg%vLevel, sg%num_cell))

  ! Initialize input data
  A = 0.05D0
  DO i = 1, sg%vLevel
    dA(i, :) = DSIN(sg%cell_cntr(1, :)) * 1000.0_R_KIND
  END DO
  epsilon = 1.0E-6_R_KIND  ! Small perturbation

  ! Test FirstOrder
  CALL calVerDer%FirstOrder(A, parA_sigma)
  CALL calVerDer_tlm%FirstOrder_tlm(A, parA_sigma, dA, dparA_sigma)
  pert_A = A + epsilon * dA
  CALL calVerDer%FirstOrder(pert_A, pert_parA_sigma)
  CALL Verify_TLM(parA_sigma, dparA_sigma, pert_parA_sigma, epsilon, sg)

  ! Clean up
  DEALLOCATE (A, parA_sigma, dA, dparA_sigma, pert_A, pert_parA_sigma)

CONTAINS

  SUBROUTINE Verify_TLM(parA_sigma, dparA_sigma, pert_parA_sigma, epsilon, sg)
    REAL(r_kind), INTENT(IN) :: parA_sigma(:, :), dparA_sigma(:, :), pert_parA_sigma(:, :), epsilon
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind) :: max_error, avg_error
    INTEGER :: k, j, n

    max_error = 0.0
    avg_error = 0.0
    n = 0
    DO k = 1, sg%vLevel
      DO j = 1, sg%num_cell
        IF (sg%bdy_type(j) .LT. 1) THEN
          max_error = MAX(max_error, ABS(dparA_sigma(k, j) - (pert_parA_sigma(k, j) - parA_sigma(k, j)) / epsilon))
          avg_error = avg_error + ABS(dparA_sigma(k, j) - (pert_parA_sigma(k, j) - parA_sigma(k, j)) / epsilon)
          n =
