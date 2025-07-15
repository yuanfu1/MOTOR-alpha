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

MODULE divergence_m

  USE kinds_m, ONLY: i_kind, r_kind

CONTAINS
  SUBROUTINE divergence(fu_c, fu_w, fu_e, &
                        fv_c, fv_s, fv_n, &
                        fw, sigma, ztop, &
                        topo_c, topo_w, topo_e, &
                        topo_s, topo_n, &
                        dist_we, dist_sn, &
                        divgc)

    IMPLICIT NONE

    REAL(r_kind), DIMENSION(:), INTENT(IN) :: fu_c, fu_e, fu_w, &
                                              fv_c, fv_n, fv_s, fw, sigma
    REAL(r_kind), INTENT(IN) :: ztop, topo_c, topo_e, topo_w, topo_n, topo_s, &
                                dist_we, dist_sn
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: H_0, H_n1, H_p1, par_fusigma, &
                                               par_fvsigma, par_fwsigma
    INTEGER(i_kind) :: vLevel
    REAL(r_kind), INTENT(INOUT) :: divgc(:)
    REAL(r_kind) :: rat, par_temp, par_sigmax, par_sigmay

    rat = ztop / (ztop - topo_c)
    vLevel = SIZE(sigma)

    ALLOCATE (H_0(vLevel - 2), H_p1(vLevel - 2), H_n1(vLevel - 2), &
              par_fusigma(vLevel), par_fvsigma(vLevel), &
              par_fwsigma(vLevel))

    par_temp = 1 / (ztop - topo_c)**2.0D0
    par_sigmax = par_temp * (topo_e - topo_w) / 2.0D0 / dist_we
    par_sigmay = par_temp * (topo_n - topo_s) / 2.0D0 / dist_sn

    H_0 = (2.0D0 * sigma(2:vLevel - 1) - sigma(1:vLevel - 2) - sigma(3:vLevel)) / &
          (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / (sigma(3:vLevel) - sigma(2:vLevel - 1))
    H_n1 = (sigma(3:vLevel) - sigma(2:vLevel - 1)) / (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / &
           (sigma(3:vLevel) - sigma(1:vLevel - 2)) ! n represents negative
    H_p1 = (sigma(1:vLevel - 2) - sigma(2:vLevel - 1)) / (sigma(3:vLevel) - sigma(2:vLevel - 1)) / &
           (sigma(1:vLevel - 2) - sigma(3:vLevel)) ! p represents positive

    par_fusigma(1) = 0.0D0
    par_fusigma(2:vLevel - 1) = (fu_c(2:vLevel - 1) * H_0 + &
                                 fu_c(3:vLevel) * H_p1 + &
                                 fu_c(1:vLevel - 2) * H_n1)
    par_fusigma(vLevel) = 0.0D0

    par_fvsigma(1) = 0.0D0
    par_fvsigma(2:vLevel - 1) = (fv_c(2:vLevel - 1) * H_0 + &
                                 fv_c(3:vLevel) * H_p1 + &
                                 fv_c(1:vLevel - 2) * H_n1)
    par_fvsigma(vLevel) = 0.0D0

    par_fwsigma(1) = 0.0D0
    par_fwsigma(2:vLevel - 1) = (fw(2:vLevel - 1) * H_0 + &
                                 fw(3:vLevel) * H_p1 + &
                                 fw(1:vLevel - 2) * H_n1)
    par_fwsigma(vLevel) = 0.0D0

    divgc(1) = (fu_e(1) - fu_w(1)) / dist_we / 2.0D0 &
               + (fv_n(1) - fv_s(1)) / dist_sn / 2.0D0
    divgc(2:vLevel - 1) = (fu_e(2:vLevel - 1) &
                           - fu_w(2:vLevel - 1)) / dist_we / 2.0D0 &
                          + (fv_n(2:vLevel - 1) &
                             - fv_s(2:vLevel - 1)) / dist_sn / 2.0D0 &
                          - par_fusigma(2:vLevel - 1) * par_sigmax &
                          - par_fvsigma(2:vLevel - 1) * par_sigmay &
                          + rat * par_fwsigma(2:vLevel - 1)
    divgc(vLevel) = (fu_e(vLevel) &
                     - fu_w(vLevel)) / dist_we / 2.0D0 &
                    + (fv_n(vLevel) &
                       - fv_s(vLevel)) / dist_sn / 2.0D0

  END SUBROUTINE divergence
END MODULE divergence_m
