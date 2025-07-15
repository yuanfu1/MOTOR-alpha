! Description:
!> @file
!!   Deallocate a clear-sky optical depth coefficients structure.
!!   This should usually be called via rttov_dealloc_coefs.
!
!> @brief
!!   Deallocate a clear-sky optical depth coefficients structure.
!!   This should usually be called via rttov_dealloc_coefs.
!!
!! @param[out]     err      status on exit
!! @param[in,out]  coef     the coefficient structure to deallocate
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
SUBROUTINE rttov_dealloc_coef(err, coef)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)   :: err
  TYPE(rttov_coef),   INTENT(INOUT) :: coef
!INTF_END
#include "rttov_nullify_coef.interface"
#include "rttov_errorreport.interface"
!- End of header --------------------------------------------------------

  TRY

  IF (ASSOCIATED(coef%fmv_gas_id)) DEALLOCATE (coef%fmv_gas_id, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%fmv_gas_pos)) DEALLOCATE (coef%fmv_gas_pos, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%fmv_var)) DEALLOCATE (coef%fmv_var, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%fmv_coe)) DEALLOCATE (coef%fmv_coe, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%fmv_ncorr)) DEALLOCATE (coef%fmv_ncorr, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%fmv_lvl)) DEALLOCATE (coef%fmv_lvl, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ff_ori_chn)) DEALLOCATE (coef%ff_ori_chn, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ff_val_chn)) DEALLOCATE (coef%ff_val_chn, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ff_cwn)) DEALLOCATE (coef%ff_cwn, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ff_bco)) DEALLOCATE (coef%ff_bco, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ff_bcs)) DEALLOCATE (coef%ff_bcs, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ff_gam)) DEALLOCATE (coef%ff_gam, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%fastem_polar)) DEALLOCATE (coef%fastem_polar, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%pol_phi)) DEALLOCATE (coef%pol_phi, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%pol_fac_v)) DEALLOCATE (coef%pol_fac_v, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%pol_fac_h)) DEALLOCATE (coef%pol_fac_h, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ssirem_a0)) DEALLOCATE (coef%ssirem_a0, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ssirem_a1)) DEALLOCATE (coef%ssirem_a1, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ssirem_a2)) DEALLOCATE (coef%ssirem_a2, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ssirem_xzn1)) DEALLOCATE (coef%ssirem_xzn1, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ssirem_xzn2)) DEALLOCATE (coef%ssirem_xzn2, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%iremis_coef)) DEALLOCATE (coef%iremis_coef, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ref_prfl_p)) DEALLOCATE (coef%ref_prfl_p, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ref_prfl_t)) DEALLOCATE (coef%ref_prfl_t, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ref_prfl_mr)) DEALLOCATE (coef%ref_prfl_mr, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%bkg_prfl_mr)) DEALLOCATE (coef%bkg_prfl_mr, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%lim_prfl_p)) DEALLOCATE (coef%lim_prfl_p, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%env_prfl_tmax)) DEALLOCATE (coef%env_prfl_tmax, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%env_prfl_tmin)) DEALLOCATE (coef%env_prfl_tmin, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%env_prfl_gmin)) DEALLOCATE (coef%env_prfl_gmin, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%env_prfl_gmax)) DEALLOCATE (coef%env_prfl_gmax, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%lim_prfl_tmax)) DEALLOCATE (coef%lim_prfl_tmax, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%lim_prfl_tmin)) DEALLOCATE (coef%lim_prfl_tmin, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%lim_prfl_gmin)) DEALLOCATE (coef%lim_prfl_gmin, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%lim_prfl_gmax)) DEALLOCATE (coef%lim_prfl_gmax, STAT = err)
  THROW(err.NE.0)


  CALL dealloc_fast_coefs(err, coef%thermal)
  THROW(err.NE.0)

  CALL dealloc_fast_coefs(err, coef%thermal_corr)
  THROW(err.NE.0)

  IF (coef%solarcoef) THEN
    CALL dealloc_fast_coefs(err, coef%solar)
    THROW(err.NE.0)

    CALL dealloc_fast_coefs(err, coef%solar_corr)
    THROW(err.NE.0)
  ENDIF

  IF (ASSOCIATED(coef%bounds)) DEALLOCATE (coef%bounds, STAT = err)
  THROW(err.NE.0)


  IF (coef%nltecoef) THEN
    IF (ASSOCIATED(coef%nlte_coef%coef)) DEALLOCATE(coef%nlte_coef%coef, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%nlte_coef%sol_zen_angle)) DEALLOCATE(coef%nlte_coef%sol_zen_angle, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%nlte_coef%sat_zen_angle)) DEALLOCATE(coef%nlte_coef%sat_zen_angle, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%nlte_coef%cos_sol)) DEALLOCATE(coef%nlte_coef%cos_sol, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%nlte_coef%sec_sat)) DEALLOCATE(coef%nlte_coef%sec_sat, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%nlte_coef)) DEALLOCATE(coef%nlte_coef, STAT = err)
    THROW(err.NE.0)
  ENDIF


  IF (coef%pmc_shift) THEN
    IF (ASSOCIATED(coef%pmc_ppmc)) DEALLOCATE(coef%pmc_ppmc, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%pmc_coef)) DEALLOCATE(coef%pmc_coef, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef%pmc_pnominal)) DEALLOCATE(coef%pmc_pnominal, STAT = err)
    THROW(err.NE.0)
  ENDIF


  IF (ASSOCIATED(coef%planck1)) DEALLOCATE (coef%planck1, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%planck2)) DEALLOCATE (coef%planck2, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%frequency_ghz)) DEALLOCATE (coef%frequency_ghz, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%dp)) DEALLOCATE (coef%dp, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%dpp)) DEALLOCATE (coef%dpp, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%tstar)) &
    DEALLOCATE (coef%tstar, coef%tstar_r, coef%tstar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%tstarmod_wsum_r)) &
    DEALLOCATE (coef%tstarmod_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%tstar_uwsum_r)) &
    DEALLOCATE (coef%tstar_uwsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%to3star)) &
    DEALLOCATE (coef%to3star, coef%to3star_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%wstar)) &
    DEALLOCATE (coef%wstar, coef%wstar_r, coef%wstar_wsum_r, coef%wtstar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ostar)) &
    DEALLOCATE (coef%ostar, coef%ostar_r, coef%ostar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%co2star)) &
    DEALLOCATE (coef%co2star, coef%co2star_r, coef%co2star_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%n2ostar)) &
    DEALLOCATE (coef%n2ostar, coef%n2ostar_r, coef%n2ostar_wsum_r, coef%n2otstar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%costar)) &
    DEALLOCATE (coef%costar, coef%costar_r, coef%costar_wsum_r, coef%cotstar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ch4star)) &
    DEALLOCATE (coef%ch4star, coef%ch4star_r, coef%ch4star_wsum_r, coef%ch4tstar_wsum_r, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%so2star)) &
    DEALLOCATE (coef%so2star, coef%so2star_r, coef%so2star_wsum_r, coef%so2tstar_wsum_r, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%tt_val_chn)) DEALLOCATE (coef%tt_val_chn, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%tt_a0)) DEALLOCATE (coef%tt_a0, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%tt_a1)) DEALLOCATE (coef%tt_a1, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%pw_val_chn)) DEALLOCATE (coef%pw_val_chn, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ss_val_chn)) DEALLOCATE (coef%ss_val_chn, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ss_solar_spectrum)) DEALLOCATE (coef%ss_solar_spectrum, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ss_rayleigh_ext)) DEALLOCATE (coef%ss_rayleigh_ext, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%refl_visnir_ow)) DEALLOCATE (coef%refl_visnir_ow, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%refl_visnir_fw)) DEALLOCATE (coef%refl_visnir_fw, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%woc_waopc_ow)) DEALLOCATE (coef%woc_waopc_ow, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%woc_waopc_fw)) DEALLOCATE (coef%woc_waopc_fw, STAT = err)
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ws_k_omega)) DEALLOCATE (coef%ws_k_omega, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef%ws_npoint)) DEALLOCATE (coef%ws_npoint, STAT = err)
  THROW(err.NE.0)


  CALL rttov_nullify_coef(coef)

  CATCH

CONTAINS

  SUBROUTINE dealloc_fast_coefs(err, fast_coefs)
    USE rttov_types, ONLY : rttov_fast_coef
    INTEGER(jpim),                  INTENT(OUT)   :: err
    TYPE(rttov_fast_coef), POINTER, INTENT(INOUT) :: fast_coefs(:)
    INTEGER(jpim) :: ichan, igas
    TRY
    IF (ASSOCIATED(fast_coefs)) THEN
      DO ichan = 1, SIZE(fast_coefs)
        IF (ASSOCIATED(fast_coefs(ichan)%gasarray)) THEN
          DO igas = 1, coef%fmv_gas
            IF (ASSOCIATED(fast_coefs(ichan)%gasarray(igas)%coef)) &
              DEALLOCATE(fast_coefs(ichan)%gasarray(igas)%coef, STAT = err)
              THROW(err.NE.0)
          ENDDO
          DEALLOCATE(fast_coefs(ichan)%gasarray, STAT = err)
          THROW(err.NE.0)
        ENDIF
      ENDDO
      DEALLOCATE(fast_coefs, STAT = err)
      THROW(err.NE.0)
    ENDIF
    CATCH
  END SUBROUTINE dealloc_fast_coefs

END SUBROUTINE rttov_dealloc_coef
