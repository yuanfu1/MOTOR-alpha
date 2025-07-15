! Description:
!> @file
!!   Nullify/zero a PC-RTTOV coefficients structure.
!
!> @brief
!!   Nullify/zero a PC-RTTOV coefficients structure.
!!
!!
!! @param[in,out]  coef_pccomp     the PC coefficients structure to nullify/zero
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_nullify_coef_pccomp(coef_pccomp)
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_pccomp
  IMPLICIT NONE
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
!INTF_END

  coef_pccomp%fmv_pc_comp_pc    = 0_jpim
  coef_pccomp%fmv_pc_cld        = 0_jpim
  coef_pccomp%fmv_pc_aer        = 0_jpim
  coef_pccomp%fmv_pc_nlte       = 0_jpim
  coef_pccomp%fmv_pc_msets      = 0_jpim
  coef_pccomp%fmv_pc_bands      = 0_jpim
  coef_pccomp%fmv_pc_mnum       = 0_jpim
  coef_pccomp%fmv_pc_mchn       = 0_jpim
  coef_pccomp%fmv_pc_nchn       = 0_jpim
  coef_pccomp%fmv_pc_nchn_noise = 0_jpim
  coef_pccomp%fmv_pc_nche       = 0_jpim
  coef_pccomp%fmv_pc_gas        = 0_jpim
  coef_pccomp%fmv_pc_gas_lim    = 0_jpim
  coef_pccomp%fmv_pc_nlev       = 0_jpim

  coef_pccomp%lim_pc_prfl_pmin  = 0._jprb
  coef_pccomp%lim_pc_prfl_pmax  = 0._jprb
  coef_pccomp%lim_pc_prfl_tsmin = 0._jprb
  coef_pccomp%lim_pc_prfl_tsmax = 0._jprb
  coef_pccomp%lim_pc_prfl_skmin = 0._jprb
  coef_pccomp%lim_pc_prfl_skmax = 0._jprb
  coef_pccomp%lim_pc_prfl_wsmin = 0._jprb
  coef_pccomp%lim_pc_prfl_wsmax = 0._jprb

  NULLIFY (coef_pccomp%fmv_pc_sets)
  NULLIFY (coef_pccomp%eigen)
  NULLIFY (coef_pccomp%emiss_chn)
  NULLIFY (coef_pccomp%emiss_c1)
  NULLIFY (coef_pccomp%emiss_c2)
  NULLIFY (coef_pccomp%emiss_c3)
  NULLIFY (coef_pccomp%emiss_c4)
  NULLIFY (coef_pccomp%emiss_c5)
  NULLIFY (coef_pccomp%emiss_c6)
  NULLIFY (coef_pccomp%emiss_c7)
  NULLIFY (coef_pccomp%emiss_c8)
  NULLIFY (coef_pccomp%emiss_c9)
  NULLIFY (coef_pccomp%ref_pc_prfl_p)
  NULLIFY (coef_pccomp%ref_pc_prfl_mr)
  NULLIFY (coef_pccomp%lim_pc_prfl_tmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_tmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_qmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_qmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_ozmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_ozmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_gasmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_gasmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_aermin)
  NULLIFY (coef_pccomp%lim_pc_prfl_aermax)
  NULLIFY (coef_pccomp%co2_pc_ref)
  NULLIFY (coef_pccomp%n2o_pc_ref)
  NULLIFY (coef_pccomp%co_pc_ref)
  NULLIFY (coef_pccomp%ch4_pc_ref)
  NULLIFY (coef_pccomp%co2_pc_min)
  NULLIFY (coef_pccomp%n2o_pc_min)
  NULLIFY (coef_pccomp%co_pc_min)
  NULLIFY (coef_pccomp%ch4_pc_min)
  NULLIFY (coef_pccomp%co2_pc_max)
  NULLIFY (coef_pccomp%n2o_pc_max)
  NULLIFY (coef_pccomp%co_pc_max)
  NULLIFY (coef_pccomp%ch4_pc_max)
  NULLIFY (coef_pccomp%noise_in)
  NULLIFY (coef_pccomp%noise)
  NULLIFY (coef_pccomp%noise_r)
  NULLIFY (coef_pccomp%ff_ori_chn_in)
  NULLIFY (coef_pccomp%ff_cwn_in)
  NULLIFY (coef_pccomp%ff_bco_in)
  NULLIFY (coef_pccomp%ff_bcs_in)
  NULLIFY (coef_pccomp%planck1_in)
  NULLIFY (coef_pccomp%planck2_in)
  NULLIFY (coef_pccomp%pcreg)
END SUBROUTINE
