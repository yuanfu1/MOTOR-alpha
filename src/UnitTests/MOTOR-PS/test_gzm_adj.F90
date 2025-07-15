MODULE test_gzm_adj_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_adj_m, ONLY: gzm_adj_t, adj_constructor

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: run_tests

CONTAINS

  SUBROUTINE run_tests()
    IMPLICIT NONE

    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(gzm_adj_t) :: gzm_adj

    REAL(r_kind), ALLOCATABLE :: adj_opr_left(:, :), adj_opr_right(:, :), adj_VALUE(:, :)

    CALL init_single_grid(sg)

    gzm_adj = adj_constructor(sg)

    ALLOCATE (adj_opr_left(sg%vLevel, sg%num_icell), adj_opr_right(sg%vLevel, sg%num_icell), adj_VALUE(sg%vLevel, sg%num_icell))

    CALL init_test_data(adj_opr_left, adj_opr_right)

    CALL gzm_adj%adj_Divergen(adj_opr_left, adj_opr_right, adj_VALUE)
    PRINT *, "Adjoint Divergen results:"
    PRINT *, adj_VALUE

    CALL gzm_adj%adj_Jacobian(adj_opr_left, adj_opr_right, adj_VALUE)
    PRINT *, "Adjoint Jacobian results:"
    PRINT *, adj_VALUE

    CALL gzm_adj%adj_Laplacia(adj_opr_left, adj_VALUE)
    PRINT *, "Adjoint Laplacia results:"
    PRINT *, adj_VALUE

    DEALLOCATE (adj_opr_left, adj_opr_right, adj_VALUE)
  END SUBROUTINE run_tests

  SUBROUTINE init_single_grid(sg)
    TYPE(SingleGrid_t), INTENT(OUT) :: sg

    sg%vLevel = 10
    sg%num_icell = 100
    sg%num_cell = 100
    sg%num_edge = 4
    sg%numQuadPerEdge = 3
    ALLOCATE (sg%edge_stcl(10, 4, 100), sg%coef_func(10, 3, 4, 100), sg%coef_norm(10, 3, 4, 100), sg%coef_tang(10, 3, 4, 100))
    ALLOCATE (sg%coef_gl(3), sg%edge_leng(4, 100), sg%cell_area(100), sg%bdy_type(100))
    ALLOCATE (sg%zHght(10, 100), sg%sigma(10), sg%Hz(10, 100), sg%parz_parsigma(10, 100), sg%parHz_parsigma(10, 100))
    ALLOCATE (sg%F_z_z(10, 100), sg%Lz(10, 100), sg%F_Hz_z(10, 100), sg%F_invHz_z(10, 100))
    sg%edge_stcl = 1.0
    sg%coef_func = 1.0
    sg%coef_norm = 1.0
    sg%coef_tang = 1.0
    sg%coef_gl = 1.0
    sg%edge_leng = 1.0
    sg%cell_area = 1.0
    sg%bdy_type = 0
    sg%zHght = 1.0
    sg%sigma = 1.0
    sg%Hz = 1.0
    sg%parz_parsigma = 1.0
    sg%parHz_parsigma = 1.0
    sg%F_z_z = 1.0
    sg%Lz = 1.0
    sg%F_Hz_z = 1.0
    sg%F_invHz_z = 1.0
  END SUBROUTINE init_single_grid

  SUBROUTINE init_test_data(adj_opr_left, adj_opr_right)
    REAL(r_kind), INTENT(OUT) :: adj_opr_left(:, :), adj_opr_right(:, :)
    adj_opr_left = 1.0
    adj_opr_right = 1.0
  END SUBROUTINE init_test_data

END MODULE test_gzm_adj_m

PROGRAM test_runner
  USE test_gzm_adj_m
  CALL run_tests()
END PROGRAM test_runner
