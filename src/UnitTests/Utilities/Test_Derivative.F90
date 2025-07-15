PROGRAM Derivative_test
  USE Derivative_m, ONLY: FirstOrderHalfLayerShortestDist_t, FirstOrderHalfLayerMinNorm_t, &
                          FirstOrder_t, SecOrderShortestDist_t, SecOrderMinNorm_t
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE
  TYPE(FirstOrderHalfLayerShortestDist_t) :: HalfLayerShortestDist
  TYPE(FirstOrderHalfLayerMinNorm_t) :: HalfLayerMinNorm
  TYPE(FirstOrder_t) :: FirstOrder
  TYPE(SecOrderShortestDist_t) :: SecOrderShortestDist
  TYPE(SecOrderMinNorm_t) :: SecOrderMinNorm
  REAL(r_kind) :: z(4), f(4), f_n(4), z2(5), f2(5), z3(3), f3(3), f2_n(5), z_out
  INTEGER(i_kind) :: ks
  REAL(r_kind) :: At_frst_ana, At_scd_ana, At_first, At_second
  REAL(r_kind) :: A3(3), A2(5)

  ks = 3
  z = [10.0D0, 11.0D0, 13.0D0, 13.5D0] - 10.0D0
  z_out = 11.5D0 - 10.0D0
  ! z2 = ([10.0D0, 11.0D0, 13.0D0, 15.0D0, 16.0D0] - 10.0D0)/10.0D0
  ! z3 = ([10.0D0, 11.0D0, 13.0D0] - 10.0D0)/10.0D0

  z2 = ([11.5D0, 12.0D0, 13.0D0, 14.0D0, 14.5D0] - 10.0D0) / 10.0D0
  z3 = ([10.5D0, 11.0D0, 12.0D0] - 10.0D0) / 10.0D0

  A3 = SIN(z3 / 2.0D0); 
  A2 = SIN(z2 / 2.0D0)
  At_frst_ana = 1.0D0 / 2.0D0 * COS(z3(2) / 2.0D0)
  At_scd_ana = -1.0D0 / 4.0D0 * SIN(z2(3) / 2.0D0)
  CALL HalfLayerShortestDist%FirstOrder(z, z_out, f)
  CALL HalfLayerMinNorm%FirstOrder(z, z_out, f_n)
  CALL SecOrderShortestDist%SecondOrder(f2, z2, ks)
  CALL SecOrderMinNorm%SecondOrder(f2_n, z2, ks)
  ks = 2
  CALL FirstOrder%FirstOrder(f3, z3, ks)
  At_first = SUM(f3 * A3)
  At_second = SUM(f2_n * A2)

  PRINT *, 'coef is ', f
  PRINT *, 'coef_n is ', f_n
  PRINT *, 'coef 2 is', f2
  PRINT *, 'coef 2_n is', f2_n
  PRINT *, 'coef 3 is', f3
  PRINT *, 'velues are ', At_first, At_frst_ana, At_second, At_scd_ana
  PRINT *, 'error first is', DABS(At_first - At_frst_ana), DABS(At_second - At_scd_ana)

  ! PRINT *, 'index is', id2
END PROGRAM

