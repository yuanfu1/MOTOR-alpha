!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.utility.test code for vertical interpolation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2021/10/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

PROGRAM interp1d_test
  USE Interp1D_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double

  IMPLICIT NONE
  INTEGER(i_kind) :: i
  REAL(r_kind) :: dT
  REAL(r_kind), DIMENSION(6) :: zin = [(i, i=2, 12, 2)], Tin = [(i, i=12, 2, -2)]
  REAL(r_kind), DIMENSION(9) :: zout = [(i, i=-1, 15, 2)], Ttest = [(i, i=15, -1, -2)], Tout(1:9)

  PRINT *, "zin:", zin
  PRINT *, "Tin:", Tin
  PRINT *, "zout:", zout
  ! PRINT *, "Tout:", Tout

  CALL interp1d(zin, Tin, zout, Tout)
  PRINT *, "Tout:", Tout
  dT = SUM(ABS(Ttest - Tout))
  PRINT *, "dT:", dT
  IF (DT < 1D-15) THEN
    PRINT *, 'Test passed'
  ELSE
    PRINT *, 'Test failed'
  END IF

END PROGRAM interp1d_test
