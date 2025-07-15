! Description:
!> @file
!!   Nullify/zero the clear-sky optical depth coefficients structure.
!
!> @brief
!!   Nullify/zero the clear-sky optical depth coefficients structure.
!!
!! @param[in,out]  coef     the coefficient structure to nullify/zero
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_nullify_coef(coef)

!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  IMPLICIT NONE
  TYPE(rttov_coef), INTENT(INOUT) :: coef
!INTF_END
  coef%id_platform      = 0_jpim
  coef%id_sat           = 0_jpim
  coef%id_inst          = 0_jpim
  coef%id_sensor        = 0_jpim
  coef%id_comp_lvl      = 0_jpim
  coef%id_creation_date = (/0_jpim, 0_jpim, 0_jpim/)
  coef%line_by_line     = 'xxxx'
  coef%readme_srf       = 'xxxx'
  coef%id_creation      = 'xxxx'
  coef%id_common_name   = 'xxxx'
  coef%fmv_model_def    = 'xxxx'
  coef%fmv_model_ver    = 0_jpim
  coef%id_comp_pc       = 0_jpim
  coef%inczeeman        = .FALSE.
  coef%fc_planck_c1     = 0._jprb
  coef%fc_planck_c2     = 0._jprb
  coef%fc_sat_height    = 0._jprb
  coef%fmv_ori_nchn     = 0_jpim
  coef%fmv_chn          = 0_jpim
  coef%fmv_gas          = 0_jpim
  coef%nmixed           = 0_jpim
  coef%nwater           = 0_jpim
  coef%nozone           = 0_jpim
  coef%nwvcont          = 0_jpim
  coef%nco2             = 0_jpim
  coef%nn2o             = 0_jpim
  coef%nco              = 0_jpim
  coef%nch4             = 0_jpim
  coef%nso2             = 0_jpim
  coef%nlevels          = 0_jpim
  coef%nlayers          = 0_jpim
  coef%ncmixed          = 0_jpim
  coef%ncwater          = 0_jpim
  coef%ncozone          = 0_jpim
  coef%ncwvcont         = 0_jpim
  coef%ncco2            = 0_jpim
  coef%ncn2o            = 0_jpim
  coef%ncco             = 0_jpim
  coef%ncch4            = 0_jpim
  coef%ncso2            = 0_jpim
  coef%nccmixed         = 0_jpim
  coef%nccwater         = 0_jpim
  coef%nccozone         = 0_jpim
  coef%nccwvcont        = 0_jpim
  coef%nccco2           = 0_jpim
  coef%nccn2o           = 0_jpim
  coef%nccco            = 0_jpim
  coef%nccch4           = 0_jpim
  coef%nccso2           = 0_jpim
  coef%ws_nomega        = 1_jpim  ! used in automatic array declaration
  coef%ssirem_ver       = 0_jpim
  coef%iremis_version   = 0_jpim
  coef%ratoe            = 0._jprb
  coef%solarcoef        = .FALSE.
  coef%nltecoef         = .FALSE.
  coef%pmc_shift        = .FALSE.
  coef%pmc_lengthcell   = 0._jprb
  coef%pmc_tempcell     = 0._jprb
  coef%pmc_betaplus1    = 0._jprb
  coef%pmc_nlay         = 0_jpim
  coef%pmc_nvar         = 0_jpim
  NULLIFY (coef%fmv_gas_id)
  NULLIFY (coef%fmv_gas_pos)
  NULLIFY (coef%fmv_var)
  NULLIFY (coef%fmv_coe)
  NULLIFY (coef%fmv_ncorr)
  NULLIFY (coef%fmv_lvl)
  NULLIFY (coef%ff_ori_chn)
  NULLIFY (coef%ff_val_chn)
  NULLIFY (coef%ff_cwn)
  NULLIFY (coef%ff_bco)
  NULLIFY (coef%ff_bcs)
  NULLIFY (coef%ff_gam)
  NULLIFY (coef%tt_val_chn)
  NULLIFY (coef%tt_a0)
  NULLIFY (coef%tt_a1)
  NULLIFY (coef%pw_val_chn)
  NULLIFY (coef%ss_val_chn)
  NULLIFY (coef%ss_solar_spectrum)
  NULLIFY (coef%ss_rayleigh_ext)
  NULLIFY (coef%refl_visnir_ow)
  NULLIFY (coef%refl_visnir_fw)
  NULLIFY (coef%woc_waopc_ow)
  NULLIFY (coef%woc_waopc_fw)
  NULLIFY (coef%ws_npoint)
  NULLIFY (coef%ws_k_omega)
  NULLIFY (coef%fastem_polar)
  NULLIFY (coef%pol_phi)
  NULLIFY (coef%pol_fac_v)
  NULLIFY (coef%pol_fac_h)
  NULLIFY (coef%ssirem_a0)
  NULLIFY (coef%ssirem_a1)
  NULLIFY (coef%ssirem_a2)
  NULLIFY (coef%ssirem_xzn1)
  NULLIFY (coef%ssirem_xzn2)
  NULLIFY (coef%iremis_coef)
  NULLIFY (coef%ref_prfl_p)
  NULLIFY (coef%ref_prfl_t)
  NULLIFY (coef%ref_prfl_mr)
  NULLIFY (coef%bkg_prfl_mr)
  NULLIFY (coef%lim_prfl_p)
  NULLIFY (coef%env_prfl_tmax)
  NULLIFY (coef%env_prfl_tmin)
  NULLIFY (coef%env_prfl_gmax)
  NULLIFY (coef%env_prfl_gmin)
  NULLIFY (coef%lim_prfl_tmax)
  NULLIFY (coef%lim_prfl_tmin)
  NULLIFY (coef%lim_prfl_gmax)
  NULLIFY (coef%lim_prfl_gmin)
  NULLIFY (coef%thermal)
  NULLIFY (coef%thermal_corr)
  NULLIFY (coef%solar)
  NULLIFY (coef%solar_corr)
  NULLIFY (coef%planck1)
  NULLIFY (coef%planck2)
  NULLIFY (coef%frequency_ghz)
  NULLIFY (coef%dp)
  NULLIFY (coef%dpp)
  NULLIFY (coef%tstar, coef%tstar_r, coef%tstar_wsum_r, coef%tstarmod_wsum_r, coef%tstar_uwsum_r)
  NULLIFY (coef%to3star, coef%to3star_r)
  NULLIFY (coef%wstar,coef%wstar_r,coef%wstar_wsum_r,coef%wtstar_wsum_r)
  NULLIFY (coef%ostar,coef%ostar_r,coef%ostar_wsum_r)
  NULLIFY (coef%co2star,coef%co2star_r,coef%co2star_wsum_r)
  NULLIFY (coef%n2ostar,coef%n2ostar_r,coef%n2ostar_wsum_r,coef%n2otstar_wsum_r)
  NULLIFY (coef%costar,coef%costar_r,coef%costar_wsum_r,coef%cotstar_wsum_r)
  NULLIFY (coef%ch4star,coef%ch4star_r,coef%ch4star_wsum_r,coef%ch4tstar_wsum_r)
  NULLIFY (coef%so2star,coef%so2star_r,coef%so2star_wsum_r,coef%so2tstar_wsum_r)
  NULLIFY (coef%nlte_coef)
  NULLIFY (coef%pmc_coef)
  NULLIFY (coef%pmc_pnominal)
  NULLIFY (coef%pmc_ppmc)
  NULLIFY (coef%bounds)
END SUBROUTINE 
