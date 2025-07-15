!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.utility.test code for Tendency
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2022/02/23, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

PROGRAM Tendency_test
  USE Tendency_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  REAL(r_kind) :: uwnd1(3, 9), vwnd1(3, 9), wwnd1(3, 9), T1(3, 9), &
                  uwnd2(5, 9), vwnd2(5, 9), wwnd2(5, 9), T2(5, 9)
  REAL(r_kind) :: x1(9), y1(9), x2(9), y2(9), z1(3), sigma1(3), z2(5), sigma2(5), &
                  wwnd_s1(3), rightHand1(3), Ten1(3), Ten_test1(3), &
                  wwnd_s2(5), rightHand2(5), Ten2(5), Ten_test2(5), &
                  dist_we1, dist_sn1, dist_we2, dist_sn2, topo(9), eps1, eps2
  INTEGER(i_kind) :: i, j, k

  dist_we1 = 0.8D0
  dist_sn1 = 0.6D0
  dist_we2 = 0.8D0
  dist_sn2 = 0.6D0

  x1 = (/1.0D0, 2.0D0, 3.0D0, 2.8D0, 2.0D0, 1.2D0, 1.0D0, 2.0D0, 3.0D0/)
  y1 = (/1.0D0, 1.4D0, 1.0D0, 2.0D0, 2.0D0, 2.0D0, 3.0D0, 2.6D0, 3.0D0/)
  z1 = (/5.0D0, 7.0D0, 1.0D01/)
  x2 = x1!(/1.5D0, 2.0D0, 2.5D0, 2.4D0, 2.0D0, 1.6D0, 1.5D0, 2.0D0, 2.5D0/)
  y2 = y1!(/1.5D0, 1.7D0, 1.5D0, 2.0D0, 2.0D0, 2.0D0, 2.5D0, 2.3D0, 2.5D0/)
  z2 = (/5.0D0, 6.0D0, 7.0D0, 8.5D0, 1.0D01/)
  topo = 1.0D0
  sigma1 = z1(3) * (z1 - topo(5)) / (z1(3) - topo(5))
  sigma2 = z2(5) * (z2 - topo(5)) / (z2(5) - topo(5))

  DO k = 1, 3
    uwnd1(k, :) = SIN(x1 + y1 + z1(k))
    vwnd1(k, :) = COS(x1 + y1 + z1(k))
    wwnd1(k, :) = z1(k) / 10.0D0; 
    T1(k, :) = EXP(-z1(k) / 10.0D0) * (x1 + y1)
    ! wwnd_s1(k) = wwnd1(k, 5)*z1(3)/(z1(3) - topo)
    rightHand1(k) = 0.0D0!z1(k)/10.0D0

    Ten_test1(k) = -uwnd1(k, 5) * EXP(-z1(k) / 1.0D01) &
                   - vwnd1(k, 5) * EXP(-z1(k) / 1.0D01) &
                   - wwnd1(k, 5) * (-0.1D0 * EXP(-z1(k) / 1.0D01)) * (x1(5) + y1(5)) &
                   + rightHand1(k)

  END DO

  DO k = 1, 5
    uwnd2(k, :) = SIN(x2 + y2 + z2(k))
    vwnd2(k, :) = COS(x2 + y2 + z2(k))
    wwnd2(k, :) = z2(k) / 10.0D0; 
    T2(k, :) = EXP(-z2(k) / 10.0D0) * (x2 + y2)
    ! wwnd_s1(k) = wwnd1(k, 5)*z1(3)/(z1(3) - topo)
    rightHand2(k) = 0.0D0!z2(k)/10.0D0
    Ten_test2(k) = -uwnd2(k, 5) * EXP(-z2(k) / 1.0D01) &
                   - vwnd2(k, 5) * EXP(-z2(k) / 1.0D01) &
                   - wwnd2(k, 5) * (-0.1D0 * EXP(-z2(k) / 1.0D01)) * (x2(5) + y2(5)) &
                   + rightHand2(k)
  END DO

  CALL Tendency(uwnd1(:, 5), vwnd1(:, 5), wwnd1(:, 5), T1(:, 5), T1(:, 6), T1(:, 4), T1(:, 2), T1(:, 8), &
                sigma1, z1(3), topo(5), topo(5), topo(5), topo(5), topo(5), dist_we1, dist_sn1, rightHand1, Ten1)
  CALL Tendency(uwnd2(:, 5), vwnd2(:, 5), wwnd2(:, 5), T2(:, 5), T2(:, 6), T2(:, 4), T2(:, 2), T2(:, 8), &
                sigma2, z2(5), topo(5), topo(5), topo(5), topo(5), topo(5), dist_we2, dist_sn2, rightHand2, Ten2)
  eps1 = dabs(Ten1(2) - Ten_test1(2))
  eps2 = dabs(Ten2(3) - Ten_test2(3))

  PRINT *, "eps1:", eps1
  PRINT *, "eps1:", eps2
  PRINT *, "eps ratio:", eps1 / eps2

  IF (eps1 / eps2 > 3.9D0) THEN
    PRINT *, 'tested passed'
  ELSE
    PRINT *, 'tested failed'
  END IF

END PROGRAM Tendency_test
