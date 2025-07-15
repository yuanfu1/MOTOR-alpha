!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

MODULE Tendency_m

  USE kinds_m, ONLY: i_kind, r_kind

CONTAINS
  SUBROUTINE Tendency(uwnd, vwnd, wwnd, Y1, Y_w, Y_e, Y_s, Y_n, &
                      sigma, ztop, topo_c, topo_w, topo_e, &
                      topo_s, topo_n, dist_we, dist_sn, &
                      rightHand, Ten)

    IMPLICIT NONE

    REAL(r_kind), DIMENSION(:), INTENT(IN) :: uwnd, vwnd, wwnd, &
                                              Y1, Y_e, Y_w, Y_n, Y_s, sigma
    REAL(r_kind), INTENT(IN) :: ztop, topo_c, topo_w, topo_n, topo_s, &
                                topo_e, dist_we, dist_sn
    REAL(r_kind), INTENT(IN), OPTIONAL :: rightHand(:)
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: wwnd_s, H_0, H_n1, H_p1, par_Ysigma
    INTEGER(i_kind) :: vLevel
    REAL(r_kind), INTENT(INOUT) :: Ten(:)
    REAL(r_kind) :: par_temp, par_sigmax, par_sigmay

    vLevel = SIZE(sigma)

    ALLOCATE (par_Ysigma(vLevel), H_0(vLevel - 2), H_n1(vLevel - 2), &
              H_p1(vLevel - 2), wwnd_s(vLevel))

    wwnd_s = wwnd * ztop / (ztop - topo_c)
    par_temp = 1 / (ztop - topo_c)**2.0D0
    par_sigmax = par_temp * (topo_e - topo_w) / 2.0D0 / dist_we
    par_sigmay = par_temp * (topo_n - topo_s) / 2.0D0 / dist_sn

    H_0 = (2.0D0 * sigma(2:vLevel - 1) - sigma(1:vLevel - 2) - sigma(3:vLevel)) / &
          (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / (sigma(3:vLevel) - sigma(2:vLevel - 1))
    H_n1 = (sigma(3:vLevel) - sigma(2:vLevel - 1)) / (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / &
           (sigma(3:vLevel) - sigma(1:vLevel - 2)) ! n represents negative
    H_p1 = (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / (sigma(3:vLevel) - sigma(2:vLevel - 1)) / &
           (sigma(1:vLevel - 2) - sigma(3:vLevel)) ! p represents positive

    par_Ysigma(2:vLevel - 1) = (Y1(2:vLevel - 1) * H_0 + &
                                Y1(3:vLevel) * H_p1 + &
                                Y1(1:vLevel - 2) * H_n1)
    par_Ysigma(vLevel) = 0.0D0

    IF (PRESENT(rightHand)) THEN
      Ten(2:vLevel - 1) = -1.0D0 * uwnd(2:vLevel - 1) * ((Y_e(2:vLevel - 1) &
                                                          - Y_w(2:vLevel - 1)) / dist_we / 2.0D0 &
                                                         - par_Ysigma(2:vLevel - 1) * par_sigmax) &
                          - vwnd(2:vLevel - 1) * ((Y_n(2:vLevel - 1) &
                                                   - Y_s(2:vLevel - 1)) / dist_sn / 2.0D0 &
                                                  - par_Ysigma(2:vLevel - 1) * par_sigmay) &
                          - wwnd_s(2:vLevel - 1) * par_Ysigma(2:vLevel - 1) &
                          + rightHand(2:vLevel - 1)
      Ten(vLevel) = -uwnd(vLevel) * (Y_e(vLevel) &
                                     - Y_w(vLevel)) / dist_we / 2.0D0 &
                    - vwnd(vLevel) * (Y_n(vLevel) &
                                      - Y_s(vLevel)) / dist_sn / 2.0D0 &
                    + rightHand(vLevel)
    ELSE
      Ten(2:vLevel - 1) = -1.0D0 * uwnd(2:vLevel - 1) * ((Y_e(2:vLevel - 1) &
                                                          - Y_w(2:vLevel - 1)) / dist_we / 2.0D0 &
                                                         - par_Ysigma(2:vLevel - 1) * par_sigmax) &
                          - vwnd(2:vLevel - 1) * ((Y_n(2:vLevel - 1) &
                                                   - Y_s(2:vLevel - 1)) / dist_sn / 2.0D0 &
                                                  - par_Ysigma(2:vLevel - 1) * par_sigmay) &
                          - wwnd_s(2:vLevel - 1) * par_Ysigma(2:vLevel - 1)
      Ten(vLevel) = -uwnd(vLevel) * (Y_e(vLevel) &
                                     - Y_w(vLevel)) / dist_we / 2.0D0 &
                    - vwnd(vLevel) * (Y_n(vLevel) &
                                      - Y_s(vLevel)) / dist_sn / 2.0D0
    END IF

    DEALLOCATE (par_Ysigma, H_0, H_n1, H_p1, wwnd_s)

  END SUBROUTINE Tendency
END MODULE Tendency_m
