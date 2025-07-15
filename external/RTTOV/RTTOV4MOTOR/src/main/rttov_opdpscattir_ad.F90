! Description:
!> @file
!!   AD of combined optical properties of aerosols and/or clouds.
!
!> @brief
!!   AD of combined optical properties of aerosols and/or clouds.
!!
!! @details
!!   AD of combined optical properties of aerosols and/or clouds.
!!
!! @param[in]     nlayers                number of layers in input profile
!! @param[in]     chanprof               specifies channels and profiles to simulate
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     dom_nstr               number of DOM streams
!! @param[in]     aux                    additional internal profile variables
!! @param[in,out] aux_ad                 additional internal profile variable increments
!! @param[in]     ircld                  computed cloud column data
!! @param[in]     profiles               input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad            atmospheric profiles and surface variable increments
!! @param[in]     profiles_int           profiles in internal units
!! @param[in,out] profiles_int_ad        profile increments in internal units
!! @param[in]     aer_opt_param          explicit aerosol optical properties per channel (optional)
!! @param[in,out] aer_opt_param_ad       explicit aerosol optical property increments (optional)
!! @param[in]     cld_opt_param          explicit cloud optical properties per channel (optional)
!! @param[in,out] cld_opt_param_ad       explicit cloud optical property increments (optional)
!! @param[in]     do_thermal             flag to indicate if any thermal (emissive) simulations are being performed
!! @param[in]     thermal                per-channel flag to indicate if thermal (emissive) simulations are being performed
!! @param[in]     do_solar               flag to indicate if any solar simulations are being performed
!! @param[in]     solar                  per-channel flag to indicate if any solar simulations are being performed
!! @param[in]     do_rayleigh_dom        flag to indicate DOM Rayleigh simulation
!! @param[in]     coef                   optical depth coefficients structure
!! @param[in]     coef_scatt             visible/IR scattering coefficients structure
!! @param[in]     coef_mfasis_cld        MFASIS cloud coefficients structure
!! @param[in]     angles                 geometry structure
!! @param[in]     raytracing             raytracing structure
!! @param[in,out] raytracing_ad          AD of raytracing structure
!! @param[in]     trans_scatt_ir         computed optical depths
!! @param[in,out] trans_scatt_ir_ad      computed optical depth increments
!! @param[in]     trans_scatt_ir_dyn     computed optical depths per cloud column and phase functions
!! @param[in,out] trans_scatt_ir_dyn_ad  computed optical depth increments per cloud column and phase function increments
!!
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
SUBROUTINE rttov_opdpscattir_ad( &
              nlayers,             &
              chanprof,            &
              opts,                &
              dom_nstr,            &
              aux,                 &
              aux_ad,              &
              ircld,               &
              profiles,            &
              profiles_ad,         &
              profiles_int,        &
              profiles_int_ad,     &
              aer_opt_param,       &
              aer_opt_param_ad,    &
              cld_opt_param,       &
              cld_opt_param_ad,    &
              do_thermal,          &
              thermal,             &
              do_solar,            &
              solar,               &
              do_rayleigh_dom,     &
              coef,                &
              coef_scatt,          &
              coef_mfasis_cld,     &
              angles,              &
              raytracing,          &
              raytracing_ad,       &
              trans_scatt_ir,      &
              trans_scatt_ir_ad,   &
              trans_scatt_ir_dyn,  &
              trans_scatt_ir_dyn_ad)

  USE rttov_types, ONLY :  &
       rttov_chanprof,              &
       rttov_options,               &
       rttov_coef,                  &
       rttov_profile,               &
       rttov_opt_param,             &
       rttov_geometry,              &
       rttov_raytracing,            &
       rttov_transmission_scatt_ir, &
       rttov_coef_scatt,            &
       rttov_coef_mfasis,           &
       rttov_profile_aux,           &
       rttov_ircld
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       realtol,            &
       deg2rad,            &
       rad2deg,            &
       ncldtyp,            &
       nwcl_max,           &
       vis_scatt_dom,      &
       vis_scatt_single,   &
       vis_scatt_mfasis,   &
       ir_scatt_dom,       &
       ir_scatt_chou,      &
       nphangle_lores,     &
       phangle_lores,      &
       nphangle_hires,     &
       phangle_hires,      &
       baran_ngauss,       &
       clw_scheme_deff,    &
       ice_scheme_baum,    &
       ice_scheme_baran2018
  USE rttov_types, ONLY : rttov_scatt_ir_aercld, rttov_optp_data
  USE rttov_scattering_mod, ONLY :  &
       spline_interp,               &
       normalise,                   &
       calc_legendre_coef_gauss,    &
       spline_interp_ad,            &
       normalise_ad,                &
       calc_legendre_coef_gauss_ad
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options),               INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                INTENT(IN)    :: dom_nstr
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux_ad
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),               INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),               INTENT(INOUT) :: profiles_int_ad(SIZE(profiles))
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(INOUT) :: aer_opt_param_ad
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: cld_opt_param
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(INOUT) :: cld_opt_param_ad
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_thermal
  LOGICAL(KIND=jplm),                INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_solar
  LOGICAL(KIND=jplm),                INTENT(IN)    :: solar(SIZE(chanprof))
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_rayleigh_dom
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_coef_scatt),            INTENT(IN)    :: coef_scatt
  TYPE(rttov_coef_mfasis),           INTENT(IN)    :: coef_mfasis_cld
  TYPE(rttov_geometry),              INTENT(IN)    :: angles(:)
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_ad
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_ad
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_dyn_ad
!INTF_END

#include "rttov_baran2014_calc_optpar.interface"
#include "rttov_baran2014_calc_optpar_ad.interface"
#include "rttov_baran2018_calc_optpar.interface"
#include "rttov_baran2018_calc_optpar_ad.interface"
#include "rttov_baran_calc_phase.interface"
#include "rttov_baran_calc_phase_ad.interface"

!   LOGICAL(KIND=jplm) :: do_dom, do_chou_scaling, do_single_scatt
  LOGICAL(KIND=jplm) :: do_mfasis, do_aer_or_ray_dom
  LOGICAL(KIND=jplm) :: do_dom_chan, do_chou_scaling_chan, do_single_scatt_chan, do_mfasis_chan, do_chan
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: j, k, k1, k2
  INTEGER(KIND=jpim) :: prof, chan, phchan
  INTEGER(KIND=jpim) :: clw_scheme, ice_scheme, col, coli, iaer, icld
  INTEGER(KIND=jpim) :: lev, lay
  REAL(KIND=jprb)    :: opd_ad, opdsun_ad
  REAL(KIND=jprb)    :: dgfrac, dgfrac_ad
  REAL(KIND=jprb)    :: abso, abso_ad
  REAL(KIND=jprb)    :: sca, sca_ad
  REAL(KIND=jprb)    :: sumsca, sumsca_ad
  REAL(KIND=jprb)    :: bpr, bpr_ad
  REAL(KIND=jprb)    :: asym, asym_ad
  REAL(KIND=jprb)    :: clw, clw_ad, tmpval, tmpval2, ray_phup, cosscata
  REAL(KIND=jprb)    :: afac, sfac, gfac, frach, frach_ad
  REAL(KIND=jprb)    :: pfac1(0:dom_nstr)
  REAL(KIND=jprb)    :: pfac2(coef_scatt%optp_aer%nphangle)

  INTEGER(KIND=jpim) :: nmom
  REAL(KIND=jprb)    :: legcoef(0:dom_nstr), thislegcoef(0:dom_nstr)
  REAL(KIND=jprb)    :: legcoefcld(0:dom_nstr)
  REAL(KIND=jprb)    :: legcoef_ad(0:dom_nstr), thislegcoef_ad(0:dom_nstr)
  REAL(KIND=jprb)    :: legcoefcld_ad(0:dom_nstr)
  REAL(KIND=jprb)    :: aer_phfn(coef_scatt%optp_aer%nphangle)
  REAL(KIND=jprb)    :: aer_phfn_ad(coef_scatt%optp_aer%nphangle)
  REAL(KIND=jprb)    :: wcldeff_phfn(coef_scatt%optp_wcl_deff%nphangle)
  REAL(KIND=jprb)    :: wcldeff_phfn_ad(coef_scatt%optp_wcl_deff%nphangle)
  REAL(KIND=jprb)    :: icl_phfn(coef_scatt%optp_icl_baum%nphangle)
  REAL(KIND=jprb)    :: icl_phfn_ad(coef_scatt%optp_icl_baum%nphangle)
  REAL(KIND=jprb)    :: zminphadiff
  REAL(KIND=jprb)    :: relazi, musat, musun
  REAL(KIND=jprb)    :: musat_ad, musun_ad, phasint_ad
  INTEGER(KIND=jpim) :: thisnphangle
  REAL(KIND=jprb)    :: thisphangle(nphangle_hires), thiscosphangle(nphangle_hires)
  REAL(KIND=jprb)    :: baran_phfn(1:nphangle_hires), baran_phfn_ad(1:nphangle_hires)
  REAL(KIND=jprb)    :: baran_phfn_interp(baran_ngauss), baran_phfn_interp_unnorm(baran_ngauss)
  REAL(KIND=jprb)    :: baran_phfn_interp_ad(baran_ngauss)

  TYPE(rttov_scatt_ir_aercld), POINTER :: aer, aer_ad, cld, cld_ad
  TYPE(rttov_optp_data),       POINTER :: optp_data

  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!-----End of header-----------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_AD', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)

!   do_dom = (do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
!            (do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
!   do_chou_scaling = do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
!   do_single_scatt = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
  do_mfasis = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
  do_aer_or_ray_dom = opts%rt_ir%addaerosl .OR. do_rayleigh_dom

  ! trans_scatt_ir_ad and trans_scatt_ir_dyn_ad are initialised in rttov_ad

  IF (do_aer_or_ray_dom) THEN
    aer    => trans_scatt_ir%aer
    aer_ad => trans_scatt_ir_ad%aer
  ENDIF
  IF (opts%rt_ir%addclouds) THEN
    cld    => trans_scatt_ir%cld
    cld_ad => trans_scatt_ir_ad%cld
  ENDIF

  IF (do_thermal .OR. .NOT. do_mfasis) THEN
    !----------------------------------------------------------------------------
    ! CALCULATE TOTAL OPTICAL DEPTHS AND PARAMETERS FOR EACH CLOUD COLUMN
    !----------------------------------------------------------------------------
    DO j = 1, nchanprof
      IF (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis) CYCLE

      chan = chanprof(j)%chan
      prof = chanprof(j)%prof

      do_dom_chan = (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                    (thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
      do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single

      DO col = ircld%ncolumn(prof), 0, -1
        opd_ad = 0._jprb
        opdsun_ad = 0._jprb
        IF (col == 0) THEN
          IF (do_aer_or_ray_dom) THEN
            DO lay = nlayers, 1, -1
              lev = lay + 1
              IF (solar(j)) THEN
                opdsun_ad = opdsun_ad + trans_scatt_ir_dyn_ad%opdpacsun(lev,col,j)
                trans_scatt_ir_ad%opdpaclsun(0,lay,j) = &
                    trans_scatt_ir_ad%opdpaclsun(0,lay,j) + opdsun_ad
              ENDIF
              IF (thermal(j)) THEN
                opd_ad = opd_ad + trans_scatt_ir_dyn_ad%opdpac(lev,col,j)
                trans_scatt_ir_ad%opdpacl(0,lay,j) = &
                    trans_scatt_ir_ad%opdpacl(0,lay,j) + opd_ad
              ENDIF
            ENDDO ! layers
          ELSE
            IF (solar(j)) trans_scatt_ir_dyn_ad%opdpacsun(:,col,j) = 0._jprb
            IF (thermal(j)) trans_scatt_ir_dyn_ad%opdpac(:,col,j) = 0._jprb
          ENDIF
        ELSE
          DO lay = nlayers, 1, -1
            lev = lay + 1
            coli = ircld%icldarr(col,lay,prof) ! This is 0 or 1 for clear(+aer)/cloud(+aer)
            IF (solar(j)) THEN
              opdsun_ad = opdsun_ad + trans_scatt_ir_dyn_ad%opdpacsun(lev,col,j)
              trans_scatt_ir_ad%opdpaclsun(coli,lay,j) = &
                  trans_scatt_ir_ad%opdpaclsun(coli,lay,j) + opdsun_ad
            ENDIF
            IF (thermal(j)) THEN
              opd_ad = opd_ad + trans_scatt_ir_dyn_ad%opdpac(lev,col,j)
              trans_scatt_ir_ad%opdpacl(coli,lay,j) = &
                  trans_scatt_ir_ad%opdpacl(coli,lay,j) + opd_ad
            ENDIF
          ENDDO ! layers
        ENDIF ! col == 0

!         IF (solar(j)) THEN
!           opdsun_ad = 0._jprb
!           trans_scatt_ir_dyn_ad%opdpacsun(1,col,j) = 0._jprb
!         ENDIF
!         IF (thermal(j)) THEN
!           opd_ad = 0._jprb
!           trans_scatt_ir_dyn_ad%opdpac(1,col,j) = 0._jprb
!         ENDIF
      ENDDO ! col

      IF (do_dom_chan .OR. do_single_scatt_chan) THEN

        DO lay = nlayers, 1, -1

          IF (solar(j)) THEN
            IF (opts%rt_ir%addclouds .AND. do_aer_or_ray_dom) THEN
              IF (ABS(cld%opdpsca(lay,j) + aer%opdpsca(lay,j)) > realtol) THEN
                tmpval = 1._jprb / (cld%opdpsca(lay,j) + aer%opdpsca(lay,j))

                IF (do_single_scatt_chan) THEN
                  aer_ad%phtotdo(lay,j) = aer_ad%phtotdo(lay,j) + trans_scatt_ir_ad%phdo(1,lay,j) * tmpval
                  cld_ad%phtotdo(lay,j) = cld_ad%phtotdo(lay,j) + trans_scatt_ir_ad%phdo(1,lay,j) * tmpval

                  tmpval2 = (aer%phtotdo(lay,j) + cld%phtotdo(lay,j)) * tmpval**2
                  cld_ad%opdpsca(lay,j) = cld_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phdo(1,lay,j) * tmpval2
                  aer_ad%opdpsca(lay,j) = aer_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phdo(1,lay,j) * tmpval2
                ENDIF

                aer_ad%phtotup(lay,j) = aer_ad%phtotup(lay,j) + trans_scatt_ir_ad%phup(1,lay,j) * tmpval
                cld_ad%phtotup(lay,j) = cld_ad%phtotup(lay,j) + trans_scatt_ir_ad%phup(1,lay,j) * tmpval

                tmpval2 = (aer%phtotup(lay,j) + cld%phtotup(lay,j)) * tmpval**2
                cld_ad%opdpsca(lay,j) = cld_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phup(1,lay,j) * tmpval2
                aer_ad%opdpsca(lay,j) = aer_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phup(1,lay,j) * tmpval2

              ENDIF
            ELSEIF (opts%rt_ir%addclouds) THEN
              IF (ABS(cld%opdpsca(lay,j)) > realtol) THEN
                tmpval = 1._jprb / cld%opdpsca(lay,j)

                IF (do_single_scatt_chan) THEN
                  cld_ad%phtotdo(lay,j) = cld_ad%phtotdo(lay,j) + trans_scatt_ir_ad%phdo(1,lay,j) * tmpval

                  tmpval2 = cld%phtotdo(lay,j) * tmpval**2
                  cld_ad%opdpsca(lay,j) = cld_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phdo(1,lay,j) * tmpval2
                ENDIF

                cld_ad%phtotup(lay,j) = cld_ad%phtotup(lay,j) + trans_scatt_ir_ad%phup(1,lay,j) * tmpval

                tmpval2 = cld%phtotup(lay,j) * tmpval**2
                cld_ad%opdpsca(lay,j) = cld_ad%opdpsca(lay,j) - trans_scatt_ir_ad%phup(1,lay,j) * tmpval2
              ENDIF
            ENDIF
          ENDIF

          ! Determine combined phase functions for aer+cld case
          IF (do_dom_chan) THEN
            IF (do_aer_or_ray_dom .AND. opts%rt_ir%addclouds) THEN
              IF (aer%sca(lay,j) + cld%sca(lay,j) > 0._jprb) THEN

                ! NB In the direct model trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) initially contains
                !    just the cloud phase fn. If aerosols are also present then it is overwritten
                !    with the combined phase function. We need to be careful about this in the AD/K.
                !    We don't want to store the cloud phase fn in this case because the phase fns take
                !    a lot of memory.
                IF (cld%sca(lay,j) > 0._jprb) THEN
                  legcoefcld = (trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) * &
                                (aer%sca(lay,j) + cld%sca(lay,j)) - &
                                trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay) * &
                                aer%sca(lay,j)) / cld%sca(lay,j)
                ELSE
                  legcoefcld = 0._jprb
                ENDIF
                legcoefcld_ad = 0._jprb

                trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,0,lay) = &
                    trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,0,lay) + &
                    trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * aer%sca(lay,j) / &
                    (aer%sca(lay,j) + cld%sca(lay,j))

                aer_ad%sca(lay,j) = aer_ad%sca(lay,j) + &
                    SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * &
                        trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay)) / &
                    (aer%sca(lay,j) + cld%sca(lay,j))

                legcoefcld_ad(:) = legcoefcld_ad(:) + &
                    trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * cld%sca(lay,j) / &
                    (aer%sca(lay,j) + cld%sca(lay,j))

                cld_ad%sca(lay,j) = cld_ad%sca(lay,j) + &
                    SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * legcoefcld(:)) / &
                    (aer%sca(lay,j) + cld%sca(lay,j))

                tmpval = SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * &
                             trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay)) / &
                         (aer%sca(lay,j) + cld%sca(lay,j))
                aer_ad%sca(lay,j) = aer_ad%sca(lay,j) - tmpval
                cld_ad%sca(lay,j) = cld_ad%sca(lay,j) - tmpval

                ! Overwrite the combined phase fn AD with just the cloud phase fn AD
                trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) = legcoefcld_ad(:)
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! layers

        ! Determine final layer *nadir* absorption and scattering optical depths
        IF (opts%rt_ir%addclouds) THEN
          trans_scatt_ir_ad%opdpabs(0,:,j) = trans_scatt_ir_ad%opdpabs(0,:,j) + &
              trans_scatt_ir_ad%opdpabs(1,:,j)
          cld_ad%opdpabs(:,j) = cld_ad%opdpabs(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpabs(1,:,j)
          trans_scatt_ir_ad%opdpsca(0,:,j) = trans_scatt_ir_ad%opdpsca(0,:,j) + &
              trans_scatt_ir_ad%opdpsca(1,:,j)
          cld_ad%opdpsca(:,j) = cld_ad%opdpsca(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpsca(1,:,j)
        ENDIF
        IF (do_aer_or_ray_dom) THEN
          aer_ad%opdpsca(:,j) = aer_ad%opdpsca(:,j) + coef%ff_gam(chan) * trans_scatt_ir_ad%opdpsca(0,:,j)
          aer_ad%opdpabs(:,j) = aer_ad%opdpabs(:,j) + coef%ff_gam(chan) * trans_scatt_ir_ad%opdpabs(0,:,j)
        ELSE
          trans_scatt_ir_ad%opdpabs(0,:,j) = 0._jprb
          trans_scatt_ir_ad%opdpsca(0,:,j) = 0._jprb
        ENDIF
      ENDIF ! do dom or do single_scatt


      ! For layer-specific quantities store just the non-cloudy and cloudy values
      ! When used later the code looks up the appropriate value for each cloud column

      ! Determine layer total aerosol/cloud optical depths
      IF (solar(j)) THEN
        IF (opts%rt_ir%addclouds) THEN
          cld_ad%opdpsun(:,j) = cld_ad%opdpsun(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpaclsun(1,:,j) * raytracing%patheff(:,prof)
          raytracing_ad%patheff(:,prof) = raytracing_ad%patheff(:,prof) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpaclsun(1,:,j) * cld%opdpsun(:,j)
          trans_scatt_ir_ad%opdpaclsun(0,:,j) = trans_scatt_ir_ad%opdpaclsun(0,:,j) + &
              trans_scatt_ir_ad%opdpaclsun(1,:,j)
        ENDIF
        IF (do_aer_or_ray_dom) THEN
          aer_ad%opdpsun(:,j) = aer_ad%opdpsun(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpaclsun(0,:,j) * raytracing%patheff(:,prof)
          raytracing_ad%patheff(:,prof) = raytracing_ad%patheff(:,prof) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpaclsun(0,:,j) * aer%opdpsun(:,j)
        ENDIF
      ENDIF
      IF (thermal(j)) THEN
        IF (opts%rt_ir%addclouds) THEN
          cld_ad%opdp(:,j) = cld_ad%opdp(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpacl(1,:,j) * raytracing%pathsat(:,prof)
          raytracing_ad%pathsat(:,prof) = raytracing_ad%pathsat(:,prof) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpacl(1,:,j) * cld%opdp(:,j)
          trans_scatt_ir_ad%opdpacl(0,:,j) = trans_scatt_ir_ad%opdpacl(0,:,j) + &
              trans_scatt_ir_ad%opdpacl(1,:,j)
        ENDIF
        IF (opts%rt_ir%addaerosl) THEN
          aer_ad%opdp(:,j) = aer_ad%opdp(:,j) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpacl(0,:,j) * raytracing%pathsat(:,prof)
          raytracing_ad%pathsat(:,prof) = raytracing_ad%pathsat(:,prof) + &
              coef%ff_gam(chan) * trans_scatt_ir_ad%opdpacl(0,:,j) * aer%opdp(:,j)
        ENDIF
      ENDIF
    ENDDO ! channels
  ENDIF ! .NOT. MFASIS


  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF CLOUDS
  !----------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    DO j = 1, nchanprof
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou

      !------------------------------------------------------------------------
      ! Calculate total cloud optical depths
      !------------------------------------------------------------------------
      IF (solar(j) .AND. .NOT. do_mfasis) THEN
        ! Full optical depth for solar channels
        cld_ad%opdpsca(:,j) = cld_ad%opdpsca(:,j) + cld_ad%opdpsun(:,j)
        cld_ad%opdpabs(:,j) = cld_ad%opdpabs(:,j) + cld_ad%opdpsun(:,j)
      ENDIF

      IF (thermal(j)) THEN
        IF (do_chou_scaling_chan) THEN
          ! Chou-scaled optical depth for thermal channels
          cld_ad%opdpabs(:,j) = cld_ad%opdpabs(:,j) + cld_ad%opdp(:,j)
          cld_ad%opdpscabpr(:,j) = cld_ad%opdpscabpr(:,j) + cld_ad%opdp(:,j)
        ELSE
          ! Full optical depth for thermal channels
          cld_ad%opdpsca(:,j) = cld_ad%opdpsca(:,j) + cld_ad%opdp(:,j)
          cld_ad%opdpabs(:,j) = cld_ad%opdpabs(:,j) + cld_ad%opdp(:,j)
        ENDIF
      ENDIF
    ENDDO

    IF (opts%rt_ir%user_cld_opt_param) THEN

      CALL calc_user_opt_param_ad(cld_opt_param, cld_opt_param_ad, cld, cld_ad, 1_jpim, .FALSE._jplm)

    ELSE

      DO j = 1, nchanprof
        chan = chanprof(j)%chan
        prof = chanprof(j)%prof
        relazi = profiles(prof)%azangle - profiles(prof)%sunazangle

        do_dom_chan = (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                      (thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
        do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
        do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
        do_mfasis_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
        do_chan = thermal(j) .OR. solar(j)

        IF (do_mfasis_chan) THEN
          clw_scheme = coef_mfasis_cld%clw_scheme
          ice_scheme = coef_mfasis_cld%ice_scheme
        ELSE
          clw_scheme = profiles(prof)%clw_scheme
          ice_scheme = profiles(prof)%ice_scheme
        ENDIF

        DO lay = nlayers, 1, -1
          !------------------------------------------------------------------
          ! Calculate final phase function Leg. coefs for all cloud types
          !------------------------------------------------------------------
          IF (do_dom_chan) THEN

            ! Retrieve direct model results (avoid recomputing)
            sumsca = cld%sca(lay,j)
            IF (do_aer_or_ray_dom) THEN
              ! We want just the cloud phase fn not cloud+aerosol combined
              legcoef(:) = trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) * &
                          (aer%sca(lay,j) + cld%sca(lay,j)) - &
                          trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay) * aer%sca(lay,j)
            ELSE
              legcoef(:) = sumsca * trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay)
            ENDIF

            sumsca_ad = 0._jprb
            legcoef_ad = 0._jprb

            IF (ABS(sumsca) > realtol) THEN
              sumsca_ad = sumsca_ad - &
                  SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) * legcoef(:)) / sumsca ** 2_jpim
              legcoef_ad(:) = & !legcoef_ad(:) + &
                  trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,1,lay) / sumsca
              sumsca_ad = sumsca_ad + cld_ad%sca(lay,j)
            ENDIF
          ENDIF

          DO icld = 1, ncldtyp

            IF (do_dom_chan) THEN
              thislegcoef(:) = 0._jprb
              thislegcoef_ad(:) = 0._jprb
            ENDIF

            sca_ad = 0._jprb

            IF (icld <= nwcl_max) THEN
              !--------------------------------------------------------------
              ! Water clouds
              !--------------------------------------------------------------

              IF (clw_scheme == clw_scheme_deff) THEN
                !--------------------------------------------------------------
                ! Water clouds - Deff scheme
                !--------------------------------------------------------------

                ! Combine all input CLW values
                IF (icld > 1) CYCLE
                clw = SUM(profiles_int(prof)%cloud(1:nwcl_max,lay))

                IF (.NOT. clw > 0._jprb .OR. .NOT. aux%clw_dg(lay,prof) > 0._jprb) CYCLE

                ! Direct model results

                ! Interpolate optical parameters based on aux%clw_dg

                optp_data => coef_scatt%optp_wcl_deff%data(1)

                IF (aux%clw_dg(lay,prof) >= optp_data%deff(1) .AND. &
                    aux%clw_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                  ! Find deff index below this channel
                  DO k1 = 1, optp_data%ndeff - 1
                    IF (optp_data%deff(k1+1) > aux%clw_dg(lay,prof)) EXIT
                  ENDDO
                  k2 = k1 + 1
                  dgfrac = (aux%clw_dg(lay,prof) - optp_data%deff(k1)) / &
                           (optp_data%deff(k2) - optp_data%deff(k1))
                ELSE
                  ! Take first or last value if deff lies beyond data range
                  IF (aux%clw_dg(lay,prof) < optp_data%deff(1)) THEN
                    k1 = 1
                    k2 = 1
                  ELSE
                    k1 = optp_data%ndeff
                    k2 = optp_data%ndeff
                  ENDIF
                  dgfrac = 0._jprb
                ENDIF

                clw_ad = 0._jprb
                dgfrac_ad = 0._jprb

                IF (do_mfasis_chan) THEN
                  clw_ad = clw_ad + &
                    coef%ff_gam(chan) * raytracing%ltick(lay,prof) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    coef%ff_gam(chan) * clw * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                  dgfrac_ad = & !dgfrac_ad + &
                    coef%ff_gam(chan) * clw * raytracing%ltick(lay,prof) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan) + &
                     optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))

                ELSEIF (do_chan) THEN

                  IF (do_dom_chan) THEN
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                      optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                      (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))
                  ENDIF

                  IF (solar(j)) THEN
                    phchan = coef_scatt%optp_wcl_deff%chan_pha_index(chan)
                    wcldeff_phfn(:) = &
                      optp_data%pha(:,1,k1,phchan) + dgfrac * &
                      (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
                  ENDIF

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                    (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * clw


                  ! Adjoint calculations

                  dgfrac_ad = 0._jprb

                  IF (solar(j)) THEN
                    musat = -1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_wcl_deff%phfn_int%zminphadiff * rad2deg

                    musat_ad = 0._jprb
                    musun_ad = 0._jprb
                    wcldeff_phfn_ad = 0._jprb

                    phchan = coef_scatt%optp_wcl_deff%chan_pha_index(chan)

                    IF (do_single_scatt_chan) THEN
                      phasint_ad = cld_ad%phtotdo(lay,j) * cld%partsca(icld,lay,j)
                      cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld_ad%phtotdo(lay,j) * cld%phintdo(icld,lay,j)

                      CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                           wcldeff_phfn, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                           coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint_ad, wcldeff_phfn_ad)

                      musat_ad = -musat_ad
                    ENDIF

                    musat = -musat

                    phasint_ad = cld_ad%phtotup(lay,j) * cld%partsca(icld,lay,j)
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                      cld_ad%phtotup(lay,j) * cld%phintup(icld,lay,j)

                    CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                         wcldeff_phfn, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                         coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint_ad, wcldeff_phfn_ad)

                    dgfrac_ad = dgfrac_ad + SUM(wcldeff_phfn_ad * &
                        (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan)))

                    raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
                    raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
                  ENDIF

                  IF (do_dom_chan) THEN
                    sca_ad = sca_ad + SUM(thislegcoef(:) * legcoef_ad(:))
                    thislegcoef_ad(:) = thislegcoef_ad(:) + legcoef_ad(:) * sca
                    dgfrac_ad = dgfrac_ad + SUM(thislegcoef_ad(0:nmom) * &
                        (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan)))
                    sca_ad = sca_ad + sumsca_ad
                  ENDIF

                  IF (do_chou_scaling_chan) THEN
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                      cld%partbpr(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                    cld_ad%partbpr(icld,lay,j) = cld_ad%partbpr(icld,lay,j) + &
                      cld%partsca(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                    dgfrac_ad = dgfrac_ad + cld_ad%partbpr(icld,lay,j) * &
                        (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                  ENDIF

                  cld_ad%partsca(icld,lay,j) = &
                      cld_ad%partsca(icld,lay,j) + cld_ad%opdpsca(lay,j)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      sca * cld_ad%partsca(icld,lay,j)
                  sca_ad = sca_ad + raytracing%ltick(lay,prof) * cld_ad%partsca(icld,lay,j)
                  clw_ad = clw_ad + &
                      sca_ad * (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
                  dgfrac_ad = dgfrac_ad + sca_ad * (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)) * &
                      clw

                  dgfrac_ad = dgfrac_ad + cld_ad%opdpabs(lay,j) * &
                      clw * raytracing%ltick(lay,prof) * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      cld_ad%opdpabs(lay,j) * clw * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
                  clw_ad = clw_ad + &
                      cld_ad%opdpabs(lay,j) * raytracing%ltick(lay,prof) * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
                ENDIF

                IF (aux%clw_dg(lay,prof) >= optp_data%deff(1) .AND. &
                    aux%clw_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                  aux_ad%clw_dg(lay,prof) = aux_ad%clw_dg(lay,prof) + &
                      dgfrac_ad / (optp_data%deff(k2) - optp_data%deff(k1))
!                   ELSE
!                     dgfrac_ad = 0._jprb
                ENDIF

                profiles_int_ad(prof)%cloud(1:nwcl_max,lay) = &
                  profiles_int_ad(prof)%cloud(1:nwcl_max,lay) + clw_ad

                clw_ad = 0._jprb

              ELSE
                !--------------------------------------------------------------
                ! Water clouds - OPAC optical parameters
                !--------------------------------------------------------------

                IF (.NOT. profiles_int(prof)%cloud(icld,lay) > 0._jprb) CYCLE

                optp_data => coef_scatt%optp_wcl_opac%data(icld)

                IF (do_mfasis_chan) THEN
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                    coef%ff_gam(chan) * optp_data%confac * &
                    (optp_data%abs(1,1,chan) + optp_data%sca(1,1,chan)) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * raytracing%ltick(lay,prof)

                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    coef%ff_gam(chan) * optp_data%confac * &
                    (optp_data%abs(1,1,chan) + optp_data%sca(1,1,chan)) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * profiles_int(prof)%cloud(icld,lay)

                ELSEIF (do_chan) THEN

                  ! Direct model results for this cloud type

                  IF (do_dom_chan) THEN
                    sca = optp_data%sca(1,1,chan) * &
                          profiles_int(prof)%cloud(icld,lay) * optp_data%confac
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
!                     legcoef(0:nmom) = legcoef(0:nmom) + optp_data%legcoef(1:nmom+1,1,1,chan) * sca
                  ENDIF


                  ! Adjoint calculation

                  IF (solar(j)) THEN
                    musat = -1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_wcl_opac%phfn_int%zminphadiff * rad2deg

                    musat_ad = 0._jprb
                    musun_ad = 0._jprb

                    phchan = coef_scatt%optp_wcl_opac%chan_pha_index(chan)

                    IF (do_single_scatt_chan) THEN
                      phasint_ad = cld_ad%phtotdo(lay,j) * cld%partsca(icld,lay,j)
                      cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                          cld_ad%phtotdo(lay,j) * cld%phintdo(icld,lay,j)

                      CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                           optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                           coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint_ad)

                      musat_ad = -musat_ad
                    ENDIF

                    musat = -musat

                    phasint_ad = cld_ad%phtotup(lay,j) * cld%partsca(icld,lay,j)
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld_ad%phtotup(lay,j) * cld%phintup(icld,lay,j)

                    CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                         optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                         coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint_ad)

                    raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
                    raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
                  ENDIF

                  IF (do_dom_chan) THEN
                    sca_ad = sca_ad + SUM(legcoef_ad(0:nmom) * optp_data%legcoef(1:nmom+1,1,1,chan))
                    sca_ad = sca_ad + sumsca_ad
                    profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                        optp_data%sca(1,1,chan) * sca_ad * optp_data%confac
                  ENDIF

                  IF (do_chou_scaling_chan) THEN
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld%partbpr(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                  ENDIF

                  cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + cld_ad%opdpsca(lay,j)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      optp_data%confac * optp_data%sca(1,1,chan) * &
                      profiles_int(prof)%cloud(icld,lay) * cld_ad%partsca(icld,lay,j)
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                      optp_data%confac * optp_data%sca(1,1,chan) * &
                      raytracing%ltick(lay,prof) * cld_ad%partsca(icld,lay,j)

                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      optp_data%confac * optp_data%abs(1,1,chan) * &
                      profiles_int(prof)%cloud(icld,lay) * cld_ad%opdpabs(lay,j)
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                      optp_data%confac * optp_data%abs(1,1,chan) * &
                      raytracing%ltick(lay,prof) * cld_ad%opdpabs(lay,j)
                ENDIF
              ENDIF

            ELSE
              !--------------------------------------------------------------
              ! Ice clouds
              !--------------------------------------------------------------
              IF (.NOT. profiles_int(prof)%cloud(icld,lay) > 0._jprb) CYCLE

              IF (ice_scheme == ice_scheme_baum) THEN
                !--------------------------------------------------------------
                ! Optical parameters from Baum database
                !--------------------------------------------------------------

                ! Direct model results

                ! Interpolate optical parameters based on aux%ice_dg

                optp_data => coef_scatt%optp_icl_baum%data(1)

                IF (aux%ice_dg(lay,prof) >= optp_data%deff(1) .AND. &
                    aux%ice_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                  ! Find deff index below this channel
                  DO k1 = 1, optp_data%ndeff - 1
                    IF (optp_data%deff(k1+1) > aux%ice_dg(lay,prof)) EXIT
                  ENDDO
                  k2 = k1 + 1
                  dgfrac = (aux%ice_dg(lay,prof) - optp_data%deff(k1)) / &
                           (optp_data%deff(k2) - optp_data%deff(k1))
                ELSE
                  ! Take first or last value if deff lies beyond data range
                  IF (aux%ice_dg(lay,prof) < optp_data%deff(1)) THEN
                    k1 = 1
                    k2 = 1
                  ELSE
                    k1 = optp_data%ndeff
                    k2 = optp_data%ndeff
                  ENDIF
                  dgfrac = 0._jprb
                ENDIF

                dgfrac_ad = 0._jprb

                IF (do_mfasis_chan) THEN
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                    coef%ff_gam(chan) * raytracing%ltick(lay,prof) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    coef%ff_gam(chan) * profiles_int(prof)%cloud(icld,lay) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                  dgfrac_ad = & !dgfrac_ad + &
                    coef%ff_gam(chan) * profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                    trans_scatt_ir_ad%opdpext(icld,lay,j) * &
                    (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan) + &
                     optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))

                ELSEIF (do_chan) THEN
                  IF (do_dom_chan) THEN
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                        optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                        (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))
                  ENDIF

                  IF (solar(j)) THEN
                    phchan = coef_scatt%optp_icl_baum%chan_pha_index(chan)
                    icl_phfn(:) = &
                        optp_data%pha(:,1,k1,phchan) + dgfrac * &
                        (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
                  ENDIF

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * profiles_int(prof)%cloud(icld,lay)


                  ! Adjoint calculations

                  dgfrac_ad = 0._jprb

                  IF (solar(j)) THEN
                    musat = -1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_icl_baum%phfn_int%zminphadiff * rad2deg

                    musat_ad = 0._jprb
                    musun_ad = 0._jprb
                    icl_phfn_ad = 0._jprb

                    phchan = coef_scatt%optp_icl_baum%chan_pha_index(chan)

                    IF (do_single_scatt_chan) THEN
                      phasint_ad = cld_ad%phtotdo(lay,j) * cld%partsca(icld,lay,j)
                      cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                          cld_ad%phtotdo(lay,j) * cld%phintdo(icld,lay,j)

                      CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                           icl_phfn, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                           coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint_ad, icl_phfn_ad)

                      musat_ad = -musat_ad
                    ENDIF

                    musat = -musat

                    phasint_ad = cld_ad%phtotup(lay,j) * cld%partsca(icld,lay,j)
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld_ad%phtotup(lay,j) * cld%phintup(icld,lay,j)

                    CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                         icl_phfn, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                         coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint_ad, icl_phfn_ad)

                    dgfrac_ad = dgfrac_ad + SUM(icl_phfn_ad * &
                        (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan)))

                    raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
                    raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
                  ENDIF

                  IF (do_dom_chan) THEN
                    sca_ad = sca_ad + SUM(thislegcoef(:) * legcoef_ad(:))
                    thislegcoef_ad(:) = thislegcoef_ad(:) + legcoef_ad(:) * sca
                    dgfrac_ad = dgfrac_ad + SUM(thislegcoef_ad(0:nmom) * &
                        (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan)))
                    sca_ad = sca_ad + sumsca_ad
                  ENDIF

                  IF (do_chou_scaling_chan) THEN
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld%partbpr(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                    bpr_ad = & !bpr_ad + &
                        cld%partsca(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                    dgfrac_ad = dgfrac_ad + bpr_ad * (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                  ENDIF

                  cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + cld_ad%opdpsca(lay,j)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + sca * cld_ad%partsca(icld,lay,j)
                  sca_ad = sca_ad + raytracing%ltick(lay,prof) * cld_ad%partsca(icld,lay,j)
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                      sca_ad * (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
                  dgfrac_ad = dgfrac_ad + sca_ad * (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)) * &
                      profiles_int(prof)%cloud(icld,lay)

                  dgfrac_ad = dgfrac_ad + cld_ad%opdpabs(lay,j) * &
                      profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      cld_ad%opdpabs(lay,j) * profiles_int(prof)%cloud(icld,lay) * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
                  profiles_int_ad(prof)%cloud(icld,lay) = profiles_int_ad(prof)%cloud(icld,lay) + &
                      cld_ad%opdpabs(lay,j) * raytracing%ltick(lay,prof) * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
                ENDIF

                IF (do_chan) THEN
                  IF (aux%ice_dg(lay,prof) >= optp_data%deff(1) .AND. &
                      aux%ice_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                    aux_ad%ice_dg(lay,prof) = aux_ad%ice_dg(lay,prof) + &
                        dgfrac_ad / (optp_data%deff(k2) - optp_data%deff(k1))
!                     ELSE
!                       dgfrac_ad = 0._jprb
                  ENDIF
                ENDIF

              ELSE
                !--------------------------------------------------------------
                ! Optical parameters computed using Baran scheme
                !--------------------------------------------------------------

                ! Direct model results

                IF (ice_scheme == ice_scheme_baran2018) THEN
                  CALL rttov_baran2018_calc_optpar(coef_scatt%optp_icl_baran2018, chan, profiles(prof)%t(lay), &
                      profiles_int(prof)%cloud(icld,lay), abso, sca, bpr, asym)
                ELSE
                  CALL rttov_baran2014_calc_optpar(coef_scatt%optp_icl_baran2014, chan, profiles(prof)%t(lay), &
                      profiles_int(prof)%cloud(icld,lay), abso, sca, bpr, asym)
                ENDIF

                ! Baran phase function
                IF (do_single_scatt_chan .OR. do_dom_chan) THEN

                  ! Set up angular grid for Baran phase fn
                  IF (solar(j)) THEN
                    thisnphangle = nphangle_hires
                    thisphangle = phangle_hires
                    thiscosphangle = coef_scatt%optp_icl_baran2018%phfn_int%cosphangle
                  ELSE
                    thisnphangle = nphangle_lores
                    thisphangle(1:nphangle_lores) = phangle_lores
                    thiscosphangle(1:nphangle_lores) = COS(phangle_lores * deg2rad)
                  ENDIF

                  ! Compute Baran phase fn
                  CALL rttov_baran_calc_phase(asym, thisphangle(1:thisnphangle), baran_phfn(1:thisnphangle))

                  ! Compute Legendre coefficients
                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    nmom = dom_nstr

                    CALL spline_interp(thisnphangle, thiscosphangle(thisnphangle:1:-1), &
                                       baran_phfn(thisnphangle:1:-1), baran_ngauss, &
                                       coef_scatt%optp_icl_baran2018%q, baran_phfn_interp)
                    baran_phfn_interp_unnorm = baran_phfn_interp
                    CALL normalise(baran_ngauss, coef_scatt%optp_icl_baran2018%w, baran_phfn_interp)
                    CALL calc_legendre_coef_gauss(coef_scatt%optp_icl_baran2018%q, coef_scatt%optp_icl_baran2018%w, &
                                                  baran_phfn_interp, dom_nstr, dom_nstr, nmom, thislegcoef(:))

!                       legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                  ENDIF
                ENDIF

                ! Adjoint calculations

                abso_ad = 0._jprb
                bpr_ad = 0._jprb
                asym_ad = 0._jprb

                ! Baran phase function
                IF (do_single_scatt_chan .OR. do_dom_chan) THEN

                  baran_phfn_ad = 0._jprb

                  IF (solar(j)) THEN
                    musat = -1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_icl_baran2018%phfn_int%zminphadiff * rad2deg

                    musat_ad = 0._jprb
                    musun_ad = 0._jprb

                    IF (do_single_scatt_chan) THEN
                      phasint_ad = cld_ad%phtotdo(lay,j) * cld%partsca(icld,lay,j)
                      cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                          cld_ad%phtotdo(lay,j) * cld%phintdo(icld,lay,j)

                      CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                           baran_phfn, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                           coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint_ad, baran_phfn_ad)

                      musat_ad = -musat_ad
                    ENDIF

                    musat = -musat

                    phasint_ad = cld_ad%phtotup(lay,j) * cld%partsca(icld,lay,j)
                    cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                        cld_ad%phtotup(lay,j) * cld%phintup(icld,lay,j)

                    CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180.0_jprb - relazi, zminphadiff, &
                                         baran_phfn, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                         coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint_ad, baran_phfn_ad)

                    raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
                    raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
                  ENDIF

                  ! Compute Legendre coefficients
                  IF (do_dom_chan) THEN
                    sca_ad = sca_ad + SUM(thislegcoef(:) * legcoef_ad(:))
                    thislegcoef_ad(:) = thislegcoef_ad(:) + legcoef_ad(:) * sca

                    baran_phfn_interp_ad = 0._jprb
                    CALL calc_legendre_coef_gauss_ad(coef_scatt%optp_icl_baran2018%q, coef_scatt%optp_icl_baran2018%w, &
                                                     baran_phfn_interp, baran_phfn_interp_ad, dom_nstr, &
                                                     thislegcoef_ad(:))

                    CALL normalise_ad(baran_ngauss, coef_scatt%optp_icl_baran2018%w, &
                                      baran_phfn_interp_unnorm, baran_phfn_interp_ad)

                    CALL spline_interp_ad(thisnphangle, thiscosphangle(thisnphangle:1:-1), &
                                          baran_phfn(thisnphangle:1:-1), baran_phfn_ad(thisnphangle:1:-1), &
                                          baran_ngauss, coef_scatt%optp_icl_baran2018%q, baran_phfn_interp_ad)

                    sca_ad = sca_ad + sumsca_ad
                  ENDIF

                  ! Compute Baran phase fn
                  CALL rttov_baran_calc_phase_ad(asym, asym_ad, thisphangle(1:thisnphangle), &
                        baran_phfn_ad(1:thisnphangle))
                ENDIF

                IF (do_chou_scaling_chan) THEN
                  cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + &
                      cld%partbpr(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                  bpr_ad = bpr_ad + &
                      cld%partsca(icld,lay,j) * cld_ad%opdpscabpr(lay,j)
                ENDIF

                cld_ad%partsca(icld,lay,j) = cld_ad%partsca(icld,lay,j) + cld_ad%opdpsca(lay,j)
                raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    sca * cld_ad%partsca(icld,lay,j)
                sca_ad = sca_ad + raytracing%ltick(lay,prof) * cld_ad%partsca(icld,lay,j)

                raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + abso * cld_ad%opdpabs(lay,j)
                abso_ad = abso_ad + raytracing%ltick(lay,prof) * cld_ad%opdpabs(lay,j)

                IF (ice_scheme == ice_scheme_baran2018) THEN
                  CALL rttov_baran2018_calc_optpar_ad(coef_scatt%optp_icl_baran2018, chan, &
                          profiles(prof)%t(lay), profiles_int(prof)%cloud(icld,lay), &
                          profiles_ad(prof)%t(lay), profiles_int_ad(prof)%cloud(icld,lay), &
                          abso_ad, sca_ad, bpr_ad, asym_ad)
                ELSE
                  CALL rttov_baran2014_calc_optpar_ad(coef_scatt%optp_icl_baran2014, chan, &
                          profiles(prof)%t(lay), profiles_int(prof)%cloud(icld,lay), &
                          profiles_ad(prof)%t(lay), profiles_int_ad(prof)%cloud(icld,lay), &
                          abso_ad, sca_ad, bpr_ad, asym_ad)
                ENDIF
              ENDIF ! ice_scheme
            ENDIF ! water or ice
          ENDDO ! cloud types

!           IF (do_dom_chan) THEN
!             sumsca_ad = 0._jprb
!             legcoef_ad = 0._jprb
!           ENDIF

        ENDDO ! layers
      ENDDO ! channels
    ENDIF ! sccld coef file

!     IF (do_dom) cld_ad%sca = 0._jprb
!     IF (do_solar .AND. .NOT. do_mfasis) cld_ad%phtotup = 0._jprb
!     IF (do_single_scatt) cld_ad%phtotdo = 0._jprb
!     IF (do_chou_scaling) cld_ad%opdpscabpr = 0._jprb
!     IF (do_thermal .OR. .NOT. do_mfasis) THEN
!       cld_ad%opdpabs = 0._jprb
!       cld_ad%opdpsca = 0._jprb
!     ENDIF
  ENDIF ! addclouds


  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF AEROSOLS
  !----------------------------------------------------------------------------
  IF (do_aer_or_ray_dom) THEN
    DO j = 1, nchanprof
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
      do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
      do_mfasis_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
      do_chan = thermal(j) .OR. solar(j)

      !------------------------------------------------------------------------
      ! Calculate total aerosol optical depths
      !------------------------------------------------------------------------
      IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
        ! Full optical depth for solar channels
        aer_ad%opdpsca(:,j) = aer_ad%opdpsca(:,j) + aer_ad%opdpsun(:,j)
        aer_ad%opdpabs(:,j) = aer_ad%opdpabs(:,j) + aer_ad%opdpsun(:,j)

        IF (do_single_scatt_chan) THEN
          WHERE (ABS(aer%opdpsca(:,j)) > realtol)
            aer_ad%phtotdo(:,j) = aer_ad%phtotdo(:,j) + &
                trans_scatt_ir_ad%phdo(0,:,j) / aer%opdpsca(:,j)
            aer_ad%opdpsca(:,j) = aer_ad%opdpsca(:,j) - &
                trans_scatt_ir_ad%phdo(0,:,j) * aer%phtotdo(:,j) / aer%opdpsca(:,j)**2
          ENDWHERE
        ENDIF
        WHERE (ABS(aer%opdpsca(:,j)) > realtol)
          aer_ad%phtotup(:,j) = aer_ad%phtotup(:,j) + &
              trans_scatt_ir_ad%phup(0,:,j) / aer%opdpsca(:,j)
          aer_ad%opdpsca(:,j) = aer_ad%opdpsca(:,j) - &
              trans_scatt_ir_ad%phup(0,:,j) * aer%phtotup(:,j) / aer%opdpsca(:,j)**2
        ENDWHERE
      ENDIF

      IF (thermal(j)) THEN
        IF (do_chou_scaling_chan) THEN
          ! Chou-scaled optical depth for thermal channels
          aer_ad%opdpabs(:,j) = aer_ad%opdpabs(:,j) + aer_ad%opdp(:,j)
          aer_ad%opdpscabpr(:,j) = aer_ad%opdpscabpr(:,j) + aer_ad%opdp(:,j)
        ELSE
          ! Full optical depth for thermal channels
          aer_ad%opdpsca(:,j) = aer_ad%opdpsca(:,j) + aer_ad%opdp(:,j)
          aer_ad%opdpabs(:,j) = aer_ad%opdpabs(:,j) + aer_ad%opdp(:,j)
        ENDIF
      ENDIF
    ENDDO

    IF (opts%rt_ir%user_aer_opt_param) THEN

      CALL calc_user_opt_param_ad(aer_opt_param, aer_opt_param_ad, aer, aer_ad, 0_jpim, do_rayleigh_dom)

    ELSE

      DO j = 1, nchanprof
        chan = chanprof(j)%chan
        prof = chanprof(j)%prof
        relazi = profiles(prof)%azangle - profiles(prof)%sunazangle

        do_dom_chan = (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                      (thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
        do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
        do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
        do_mfasis_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
        do_chan = thermal(j) .OR. solar(j)

        IF (do_rayleigh_dom .AND. solar(j) .AND. &
            10000._jprb / coef%ff_cwn(chan) <= opts%rt_ir%rayleigh_max_wavelength) THEN
          ! Due to strict plane-parallel geometry, ray_phup (used below) is fixed for all
          ! layers so can be calculated outside layer loop, and has zero TL/AD/K
          cosscata = - angles(prof)%coszen * angles(prof)%coszen_sun - &
                       angles(prof)%sinzen * angles(prof)%sinzen_sun * COS(relazi * deg2rad)
          ray_phup = 0.75_jprb * (1._jprb + cosscata * cosscata)
        ENDIF

        DO lay = nlayers, 1, -1
          !------------------------------------------------------------------
          ! Calculate final phase function Leg. coefs for all aerosol types
          !------------------------------------------------------------------
          IF (do_dom_chan) THEN

            ! Retrieve direct model results (avoid recomputing)
            sumsca = aer%sca(lay,j)
            legcoef(:) = sumsca * trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay)

            sumsca_ad = 0._jprb
            legcoef_ad = 0._jprb

            IF (ABS(sumsca) > realtol) THEN
              sumsca_ad = sumsca_ad - &
                  SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,0,lay) * legcoef(:)) / sumsca ** 2_jpim
              legcoef_ad(:) = &!legcoef_ad(:) + &
                  trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,0,lay) / sumsca
              sumsca_ad = sumsca_ad + aer_ad%sca(lay,j)
            ENDIF
          ENDIF


          DO iaer = 1, coef_scatt%optp_aer%ntypes
            IF (profiles_int(prof)%aerosols(iaer,lay) <= 0._jprb) CYCLE

            ! Direct model results for this aerosol type

            optp_data => coef_scatt%optp_aer%data(iaer)

            IF (do_dom_chan .OR. do_single_scatt_chan) thislegcoef(:) = 0._jprb

            k = optp_data%nrelhum
            IF (k /= 1 .AND. aux%relhum(lay,prof) <= optp_data%relhum(k)) THEN
              ! Interpolate scattering parameters to actual value of relative humidity
              DO k = 1, optp_data%nrelhum - 1
                IF (aux%relhum(lay,prof) >= optp_data%relhum(k) .AND. &
                    aux%relhum(lay,prof) <= optp_data%relhum(k+1)) THEN

                  frach = (aux%relhum(lay,prof) - optp_data%relhum(k)) / &
                          (optp_data%relhum(k+1) - optp_data%relhum(k))
                  afac  = (optp_data%abs(k+1,1,chan) - optp_data%abs(k,1,chan))
                  abso  = optp_data%abs(k,1,chan) + afac * frach
                  sfac  = (optp_data%sca(k+1,1,chan) - optp_data%sca(k,1,chan))
                  sca   = optp_data%sca(k,1,chan) + sfac * frach

                  IF (do_chou_scaling_chan) THEN
                    gfac = (optp_data%bpr(k+1,1,chan) - optp_data%bpr(k,1,chan))
                    bpr  = optp_data%bpr(k,1,chan) + gfac * frach
                  ENDIF

                  IF (do_dom_chan) THEN
                    nmom = MIN(MAX(optp_data%nmom(k,chan), optp_data%nmom(k+1,chan)), dom_nstr)
                    pfac1(0:nmom) = (optp_data%legcoef(1:nmom+1,k+1,1,chan) - &
                                     optp_data%legcoef(1:nmom+1,k,1,chan))
                    thislegcoef(0:nmom) = optp_data%legcoef(1:nmom+1,k,1,chan) + pfac1(0:nmom) * frach
                  ENDIF

                  IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
                    phchan = coef_scatt%optp_aer%chan_pha_index(chan)
                    pfac2(:) = (optp_data%pha(:,k+1,1,phchan) - &
                                optp_data%pha(:,k,1,phchan))
                    aer_phfn(:) = optp_data%pha(:,k,1,phchan) + pfac2(:) * frach
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ELSE
              ! Particle doesn't change with rel. hum. (k=1) or rel. hum. exceeds max (k=max_rh_index)
              abso = optp_data%abs(k,1,chan)
              sca  = optp_data%sca(k,1,chan)
              IF (do_chou_scaling_chan) THEN
                bpr = optp_data%bpr(k,1,chan)
              ENDIF
              IF (do_dom_chan) THEN
                nmom = MIN(optp_data%nmom(k,chan), dom_nstr)
                thislegcoef(0:nmom) = optp_data%legcoef(1:nmom+1,k,1,chan)
              ENDIF
              IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
                phchan = coef_scatt%optp_aer%chan_pha_index(chan)
                aer_phfn(:) = optp_data%pha(:,k,1,phchan)
              ENDIF
            ENDIF


            ! Adjoint calculation

            IF (do_mfasis_chan) THEN
              profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                coef%ff_gam(chan) * (abso + sca) * raytracing%ltick(lay,prof) * &
                trans_scatt_ir_ad%opdpext(iaer,lay,j)

              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                coef%ff_gam(chan) * (abso + sca) * profiles_int(prof)%aerosols(iaer,lay) * &
                trans_scatt_ir_ad%opdpext(iaer,lay,j)

              abso_ad = & !abso_ad + &
                coef%ff_gam(chan) * profiles_int(prof)%aerosols(iaer,lay) * raytracing%ltick(lay,prof) * &
                trans_scatt_ir_ad%opdpext(iaer,lay,j)

              sca_ad = & !sca_ad + &
                coef%ff_gam(chan) * profiles_int(prof)%aerosols(iaer,lay) * raytracing%ltick(lay,prof) * &
                trans_scatt_ir_ad%opdpext(iaer,lay,j)
            ELSEIF (do_chan) THEN
              !------------------------------------------------------------------
              ! Accumulate total aerosol optical parameters (all particles)
              !------------------------------------------------------------------
              ! For pre-defined particle types OD calculated as (aerosol amount in layer) *
              !   (aerosol ext coeff for given layer RH) * (layer thickness)
              ! OD is accumulated for each aerosol type - arrays are zeroed above.

              sca_ad = 0._jprb

              IF (solar(j)) THEN
                musat = -1._jprb / raytracing%pathsat(lay,prof)
                musun = -1._jprb / raytracing%pathsun(lay,prof)
                zminphadiff = coef_scatt%optp_aer%phfn_int%zminphadiff * rad2deg

                musat_ad = 0._jprb
                musun_ad = 0._jprb
                aer_phfn_ad = 0._jprb

                IF (do_single_scatt_chan) THEN
                  profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                      aer_ad%phtotdo(lay,j) * aer%phintdo(iaer,lay,j) * &
                      sca * raytracing%ltick(lay,prof)
                  phasint_ad = &
                      aer_ad%phtotdo(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                      sca * raytracing%ltick(lay,prof)
                  sca_ad = sca_ad + &
                      aer_ad%phtotdo(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                      aer%phintdo(iaer,lay,j) * raytracing%ltick(lay,prof)
                  raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                      aer_ad%phtotdo(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                      aer%phintdo(iaer,lay,j) * sca

                  CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                       aer_phfn, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                       coef_scatt%optp_aer%phfn_int%iphangle, phasint_ad, aer_phfn_ad)

                  musat_ad = -musat_ad
                ENDIF

                musat = -musat

                profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                    aer_ad%phtotup(lay,j) * aer%phintup(iaer,lay,j) * &
                    sca * raytracing%ltick(lay,prof)
                phasint_ad = &
                    aer_ad%phtotup(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                    sca * raytracing%ltick(lay,prof)
                sca_ad = sca_ad + &
                    aer_ad%phtotup(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                    aer%phintup(iaer,lay,j) * raytracing%ltick(lay,prof)
                raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    aer_ad%phtotup(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                    aer%phintup(iaer,lay,j) * sca

                CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180._jprb - relazi, zminphadiff, &
                                     aer_phfn, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                     coef_scatt%optp_aer%phfn_int%iphangle, phasint_ad, aer_phfn_ad)

                raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
                raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
              ENDIF

              IF (do_dom_chan) THEN
                thislegcoef_ad = 0._jprb
                sca_ad = sca_ad + &
                    SUM(legcoef_ad(:) * thislegcoef(:)) * profiles_int(prof)%aerosols(iaer,lay)
                profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                    SUM(legcoef_ad(:) * thislegcoef(:)) * sca
                thislegcoef_ad(:) = thislegcoef_ad(:) + &
                    legcoef_ad(:) * profiles_int(prof)%aerosols(iaer,lay) * sca

                sca_ad = sca_ad + sumsca_ad * profiles_int(prof)%aerosols(iaer,lay)
                profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + sumsca_ad * sca
              ENDIF

              IF (do_chou_scaling_chan) THEN
  !               bpr_ad = 0._jprb
                raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                    aer_ad%opdpscabpr(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * sca * bpr
                bpr_ad = & !bpr_ad + &
                    aer_ad%opdpscabpr(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                    sca * raytracing%ltick(lay,prof)
                sca_ad = sca_ad + &
                    aer_ad%opdpscabpr(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                    bpr * raytracing%ltick(lay,prof)
                profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                    aer_ad%opdpscabpr(lay,j) * sca * bpr * raytracing%ltick(lay,prof)
              ENDIF

              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                  aer_ad%opdpsca(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * sca
              sca_ad = sca_ad + &
                  aer_ad%opdpsca(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                  raytracing%ltick(lay,prof)
              profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                  aer_ad%opdpsca(lay,j) * sca * raytracing%ltick(lay,prof)

              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                  aer_ad%opdpabs(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * abso
              abso_ad = & !abso_ad + &
                  aer_ad%opdpabs(lay,j) * profiles_int(prof)%aerosols(iaer,lay) * &
                  raytracing%ltick(lay,prof)
              profiles_int_ad(prof)%aerosols(iaer,lay) = profiles_int_ad(prof)%aerosols(iaer,lay) + &
                  aer_ad%opdpabs(lay,j) * abso * raytracing%ltick(lay,prof)
            ENDIF

            IF (do_chan) THEN
              k = optp_data%nrelhum
              IF (k /= 1 .AND. aux%relhum(lay,prof) <= optp_data%relhum(k)) THEN
                ! Interpolate scattering parameters to actual value of relative humidity
                DO k = 1, optp_data%nrelhum - 1
                  IF (aux%relhum(lay,prof) >= optp_data%relhum(k) .AND. &
                      aux%relhum(lay,prof) <= optp_data%relhum(k+1)) THEN
                    afac  = (optp_data%abs(k+1,1,chan) - optp_data%abs(k,1,chan))
                    sfac  = (optp_data%sca(k+1,1,chan) - optp_data%sca(k,1,chan))

                    frach_ad = 0._jprb

                    IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
                      phchan = coef_scatt%optp_aer%chan_pha_index(chan)
                      pfac2(:) = (optp_data%pha(:,k+1,1,phchan) - &
                                  optp_data%pha(:,k,1,phchan))
                      frach_ad = frach_ad + SUM(pfac2(:) * aer_phfn_ad(:))
                    ENDIF

                    IF (do_dom_chan) THEN
                      nmom = MIN(MAX(optp_data%nmom(k,chan), optp_data%nmom(k+1,chan)), dom_nstr)
                      pfac1(0:nmom) = (optp_data%legcoef(1:nmom+1,k+1,1,chan) - &
                                       optp_data%legcoef(1:nmom+1,k,1,chan))
                      frach_ad = frach_ad + SUM(pfac1(0:nmom) * thislegcoef_ad(0:nmom))
                    ENDIF

                    IF (do_chou_scaling_chan) THEN
                      gfac  = (optp_data%bpr(k+1,1,chan) - optp_data%bpr(k,1,chan))
                      frach_ad = frach_ad + gfac * bpr_ad
                    ENDIF

                    frach_ad = frach_ad + sfac * sca_ad
                    frach_ad = frach_ad + afac * abso_ad
                    aux_ad%relhum(lay,prof) = aux_ad%relhum(lay,prof) + frach_ad / &
                                              (optp_data%relhum(k+1) - optp_data%relhum(k))
                    EXIT
                  ENDIF
                ENDDO
!               ELSE
                ! Particle doesn't change with rel. hum. (k=1) or rel. hum. exceeds max (k=max_rh_index)
                !   - not strictly necessary to zero the variables
              ENDIF
            ENDIF

          ENDDO ! aer types

          IF (do_dom_chan) THEN
            IF (do_rayleigh_dom .AND. solar(j) .AND. &
                10000._jprb / coef%ff_cwn(chan) <= opts%rt_ir%rayleigh_max_wavelength .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)

              sca_ad = 1._jprb * legcoef_ad(0) + 0.5_jprb * legcoef_ad(2) + sumsca_ad
              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) - &
                                              sca_ad * sca / raytracing%ltick(lay,prof)
              trans_scatt_ir_ad%ray_sca(lay,j) = trans_scatt_ir_ad%ray_sca(lay,j) + &
                                                 sca_ad / raytracing%ltick(lay,prof) + &
                                                 ray_phup * aer_ad%phtotup(lay,j) + &
                                                 aer_ad%opdpsca(lay,j)
            ENDIF
!             sumsca_ad = 0._jprb
!             legcoef_ad(:) = 0._jprb
          ENDIF

        ENDDO ! layers
      ENDDO ! channels
    ENDIF ! scaer coef file

!     aer_ad%opdp(:,:) = 0._jprb

!     IF (do_dom) aer_ad%sca = 0._jprb
!     IF (do_solar .AND. .NOT. do_mfasis) aer_ad%phtotup = 0._jprb
!     IF (do_single_scatt) aer_ad%phtotdo = 0._jprb
!     IF (do_chou_scaling) aer_ad%opdpscabpr = 0._jprb
!     IF (do_thermal .OR. .NOT. do_mfasis) THEN
!       aer_ad%opdpabs = 0._jprb
!       aer_ad%opdpsca = 0._jprb
!     ENDIF
  ENDIF ! do_aer_or_ray_dom

!   IF (do_mfasis) trans_scatt_ir_ad%opdpext = 0._jprb
!   IF (do_dom) THEN
!     DO j = 1, nchanprof
!       IF (ASSOCIATED(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef)) &
!           trans_scatt_ir_dyn_ad%phasefn(j)%legcoef = 0._jprb
!     ENDDO
!   ENDIF
!   IF (do_single_scatt) trans_scatt_ir_ad%phdo = 0._jprb
  IF (do_solar .AND. .NOT. do_mfasis) THEN
!     trans_scatt_ir_ad%phup = 0._jprb
!     trans_scatt_ir_ad%opdpaclsun = 0._jprb
  ENDIF
  IF (do_thermal) THEN
!     trans_scatt_ir_ad%opdpacl = 0._jprb
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_AD', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_user_opt_param_ad(opt_param, opt_param_ad, aercld, aercld_ad, iaercld, do_rayleigh_dom)
    ! Process aer/cld user optical property inputs
    TYPE(rttov_opt_param),       INTENT(IN)              :: opt_param
    TYPE(rttov_opt_param),       INTENT(INOUT), OPTIONAL :: opt_param_ad
    TYPE(rttov_scatt_ir_aercld), INTENT(IN)              :: aercld
    TYPE(rttov_scatt_ir_aercld), INTENT(INOUT)           :: aercld_ad
    INTEGER(jpim),               INTENT(IN)              :: iaercld ! aer=>0, cld=>1
    LOGICAL(jplm),               INTENT(IN)              :: do_rayleigh_dom

    INTEGER(jpim) :: j, lay, prof, chan
    LOGICAL(jplm) :: do_chou_scaling_chan, do_dom_chan, do_single_scatt_chan
    LOGICAL(jplm) :: do_opt_param_lay, do_rayleigh_dom_chan, do_opt_param_ad
    REAL(jprb)    :: musat, musun, zminphadiff, phasint_ad
    REAL(jprb)    :: ray_lcoef(0:dom_nstr), opdpsca, opdpsca_ad, sca, sca_ad

    ! This is useful for the parallel interface which always passes opt_param_ad
    ! for convenience even when not allocated
    do_opt_param_ad = .FALSE.
    IF (PRESENT(opt_param_ad) .AND. .NOT. opts%dev%no_opt_param_tladk) &
      do_opt_param_ad = (ASSOCIATED(opt_param_ad%sca))

    IF (do_rayleigh_dom) THEN
      ray_lcoef = 0.
      ray_lcoef(0) = 1._jprb
      ray_lcoef(2) = 0.5_jprb
    ENDIF

    DO j = 1, nchanprof
      prof = chanprof(j)%prof
      chan = chanprof(j)%chan

      do_dom_chan = (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                    (thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
      do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
      do_rayleigh_dom_chan = do_rayleigh_dom .AND. solar(j) .AND. &
          10000._jprb / coef%ff_cwn(chan) <= opts%rt_ir%rayleigh_max_wavelength

      relazi = profiles(prof)%azangle - profiles(prof)%sunazangle

      IF (do_rayleigh_dom) THEN
        ! Due to strict plane-parallel geometry, ray_phup (used below) is fixed for all
        ! layers so can be calculated outside layer loop, and has zero TL/AD/K
        cosscata = - angles(prof)%coszen * angles(prof)%coszen_sun - &
                     angles(prof)%sinzen * angles(prof)%sinzen_sun * COS(relazi * deg2rad)
        ray_phup = 0.75_jprb * (1._jprb + cosscata * cosscata)
      ENDIF
      
      DO lay = 1, nlayers

        do_opt_param_lay = opt_param%abs(lay,j) > 0._jprb .OR. opt_param%sca(lay,j) > 0._jprb

        opdpsca = opt_param%sca(lay,j) * raytracing%ltick(lay,prof)
        opdpsca_ad = 0._jprb

        IF (do_opt_param_ad) THEN
          IF (do_opt_param_lay .AND. solar(j)) THEN
            musat = -1._jprb / raytracing%pathsat(lay,prof)
            musun = -1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = opt_param%phasefn_int%zminphadiff * rad2deg

            musat_ad = 0._jprb
            musun_ad = 0._jprb

            IF (do_single_scatt_chan) THEN
              phasint_ad = aercld_ad%phtotdo(lay,j) * opdpsca
              opdpsca_ad = opdpsca_ad + aercld_ad%phtotdo(lay,j) * aercld%phintdo(1,lay,j)

              CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                   opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                   opt_param%phasefn_int%iphangle, phasint_ad, opt_param_ad%pha(:,lay,j))

              musat_ad = -musat_ad
            ENDIF

            musat = -musat

            phasint_ad = aercld_ad%phtotup(lay,j) * opdpsca
            opdpsca_ad = opdpsca_ad + aercld_ad%phtotup(lay,j) * aercld%phintup(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180._jprb - relazi, zminphadiff, &
                                 opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                 opt_param%phasefn_int%iphangle, phasint_ad, opt_param_ad%pha(:,lay,j))

            raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
            raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
          ENDIF

          IF (do_dom_chan) THEN
            IF (do_rayleigh_dom_chan .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)

              sca_ad = aercld_ad%sca(lay,j)
              IF (do_opt_param_lay) THEN
                opt_param_ad%sca(lay,j) = opt_param_ad%sca(lay,j) + aercld_ad%sca(lay,j)

                sca_ad = sca_ad + &
                  SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,iaercld,lay) * &
                      (ray_lcoef(:) - trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay))) / &
                  (sca + opt_param%sca(lay,j))
                opt_param_ad%legcoef(1:dom_nstr+1,lay,j) = opt_param_ad%legcoef(1:dom_nstr+1,lay,j) + &
                  trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,iaercld,lay) * opt_param%sca(lay,j) / &
                  (sca + opt_param%sca(lay,j))
                opt_param_ad%sca(lay,j) = opt_param_ad%sca(lay,j) + &
                  SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,iaercld,lay) * &
                      (opt_param%legcoef(1:dom_nstr+1,lay,j) - trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay))) / &
                  (sca + opt_param%sca(lay,j))
              ENDIF

              trans_scatt_ir_ad%ray_sca(lay,j) = trans_scatt_ir_ad%ray_sca(lay,j) + &
                                                 sca_ad / raytracing%ltick(lay,prof)
              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) - &
                                              sca_ad * sca / raytracing%ltick(lay,prof)
            ELSEIF (do_opt_param_lay) THEN
              opt_param_ad%sca(lay,j) = opt_param_ad%sca(lay,j) + aercld_ad%sca(lay,j)
              opt_param_ad%legcoef(1:dom_nstr+1,lay,j) = opt_param_ad%legcoef(1:dom_nstr+1,lay,j) + &
                trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,iaercld,lay)
            ENDIF
          ENDIF

          IF (do_opt_param_lay .AND. do_chou_scaling_chan) THEN
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                opt_param%sca(lay,j) * opt_param%bpr(lay,j) * aercld_ad%opdpscabpr(lay,j)
            opt_param_ad%sca(lay,j)       = opt_param_ad%sca(lay,j) + &
                raytracing%ltick(lay,prof) * opt_param%bpr(lay,j) * aercld_ad%opdpscabpr(lay,j)
            opt_param_ad%bpr(lay,j)       = opt_param_ad%bpr(lay,j) + &
                raytracing%ltick(lay,prof) * opt_param%sca(lay,j) * aercld_ad%opdpscabpr(lay,j)
          ENDIF

          IF (do_rayleigh_dom_chan .AND. &
              profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

            trans_scatt_ir_ad%ray_sca(lay,j) = trans_scatt_ir_ad%ray_sca(lay,j) + &
                                               ray_phup * aercld_ad%phtotup(lay,j) + &
                                               aercld_ad%opdpsca(lay,j)
            opdpsca_ad = opdpsca_ad + aercld_ad%opdpsca(lay,j)
          ELSEIF (do_opt_param_lay) THEN
            opdpsca_ad = opdpsca_ad + aercld_ad%opdpsca(lay,j)
          ENDIF

          IF (do_opt_param_lay) THEN
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                            opt_param%sca(lay,j) * opdpsca_ad
            opt_param_ad%sca(lay,j)       = opt_param_ad%sca(lay,j) + &
                                            raytracing%ltick(lay,prof) * opdpsca_ad

            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                            opt_param%abs(lay,j) * aercld_ad%opdpabs(lay,j)
            opt_param_ad%abs(lay,j)       = opt_param_ad%abs(lay,j) + &
                                            raytracing%ltick(lay,prof) * aercld_ad%opdpabs(lay,j)
          ENDIF
        ELSE
          ! No opt_param_ad
          IF (do_opt_param_lay .AND. solar(j)) THEN
            musat = -1._jprb / raytracing%pathsat(lay,prof)
            musun = -1._jprb / raytracing%pathsun(lay,prof)
            zminphadiff = opt_param%phasefn_int%zminphadiff * rad2deg

            musat_ad = 0._jprb
            musun_ad = 0._jprb

            IF (do_single_scatt_chan) THEN
              phasint_ad = aercld_ad%phtotdo(lay,j) * opdpsca
              opdpsca_ad = opdpsca_ad + aercld_ad%phtotdo(lay,j) * aercld%phintdo(1,lay,j)

              CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                                   opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                   opt_param%phasefn_int%iphangle, phasint_ad)

              musat_ad = -musat_ad
            ENDIF

            musat = -musat

            phasint_ad = aercld_ad%phtotup(lay,j) * opdpsca
            opdpsca_ad = opdpsca_ad  + aercld_ad%phtotup(lay,j) * aercld%phintup(1,lay,j)

            CALL int_phase_fn_ad(musat, musat_ad, musun, musun_ad, 180._jprb - relazi, zminphadiff, &
                                 opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                 opt_param%phasefn_int%iphangle, phasint_ad)

            raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) + musun_ad * musun**2
            raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) - musat_ad * musat**2
          ENDIF

          IF (do_dom_chan) THEN
            IF (do_rayleigh_dom_chan .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)

              sca_ad = aercld_ad%sca(lay,j)
              IF (do_opt_param_lay) THEN
                sca_ad = sca_ad + &
                  SUM(trans_scatt_ir_dyn_ad%phasefn(j)%legcoef(:,iaercld,lay) * &
                      (ray_lcoef(:) - trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay))) / &
                  (sca + opt_param%sca(lay,j))
              ENDIF

              trans_scatt_ir_ad%ray_sca(lay,j) = trans_scatt_ir_ad%ray_sca(lay,j) + &
                                                 sca_ad / raytracing%ltick(lay,prof)
              raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) - &
                                              sca_ad * sca / raytracing%ltick(lay,prof)
            ELSEIF (do_opt_param_lay) THEN
              aercld_ad%sca(lay,j) = 0._jprb
            ENDIF
          ENDIF

          IF (do_opt_param_lay .AND. do_chou_scaling_chan) THEN
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                opt_param%sca(lay,j) * opt_param%bpr(lay,j) * aercld_ad%opdpscabpr(lay,j)
          ENDIF

          IF (do_rayleigh_dom_chan .AND. &
              profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

            trans_scatt_ir_ad%ray_sca(lay,j) = trans_scatt_ir_ad%ray_sca(lay,j) + &
                                               ray_phup * aercld_ad%phtotup(lay,j) + &
                                               aercld_ad%opdpsca(lay,j)
            opdpsca_ad = opdpsca_ad + aercld_ad%opdpsca(lay,j)
          ELSEIF (do_opt_param_lay) THEN
            opdpsca_ad = opdpsca_ad + aercld_ad%opdpsca(lay,j)
          ENDIF

          IF (do_opt_param_lay) THEN
            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                            opt_param%sca(lay,j) * opdpsca_ad

            raytracing_ad%ltick(lay,prof) = raytracing_ad%ltick(lay,prof) + &
                                            opt_param%abs(lay,j) * aercld_ad%opdpabs(lay,j)
          ENDIF
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE calc_user_opt_param_ad

  SUBROUTINE int_phase_fn_ad(musat, musat_ad, musun, musun_ad, relazi, zminphadiff, &
                             pha, cospha, ipha, phasint_ad, pha_ad)
    ! Interpolate phase function to scattering angle
    ! pha_ad may be omitted if the phase function is static
    REAL(KIND=jprb),           INTENT(IN)    :: musat, musun, relazi
    REAL(KIND=jprb),           INTENT(INOUT) :: musat_ad, musun_ad
    REAL(KIND=jprb),           INTENT(IN)    :: zminphadiff
    REAL(KIND=jprb),           INTENT(IN)    :: pha(:)
    REAL(KIND=jprb),           INTENT(IN)    :: cospha(:)
    INTEGER(KIND=jpim),        INTENT(IN)    :: ipha(:)
    REAL(KIND=jprb),           INTENT(IN)    :: phasint_ad
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: pha_ad(:)

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, ztmpx_ad
    REAL(KIND=jprb)    :: scattangle, scattangle_ad, deltap, deltap_ad, delta

    ztmpx = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))

    ztmpx_ad      = 0._jprb
    scattangle_ad = 0._jprb
    deltap_ad     = 0._jprb

    scattangle    = musat * musun + ztmpx * COS(relazi * deg2rad)
    ikk           = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk            = ipha(ikk) - 1_jpim
    deltap        = pha(kk + 1) - pha(kk)
    delta         = cospha(kk) - cospha(kk + 1)

    IF (PRESENT(pha_ad)) THEN
      pha_ad(kk) = pha_ad(kk) + phasint_ad
      deltap_ad = deltap_ad + phasint_ad * (cospha(kk) - scattangle) / delta
      scattangle_ad = scattangle_ad - deltap * phasint_ad / delta
      pha_ad(kk + 1) = pha_ad(kk + 1) + deltap_ad
      pha_ad(kk) = pha_ad(kk) - deltap_ad
    ELSE
      scattangle_ad = scattangle_ad - deltap * phasint_ad / delta
    ENDIF

    musat_ad = musat_ad + scattangle_ad * musun
    musun_ad = musun_ad + scattangle_ad * musat
    ztmpx_ad = ztmpx_ad + scattangle_ad * COS(relazi * deg2rad)

    IF (ABS(ztmpx) < realtol) THEN
      ztmpx_ad = 0._jprb
    ELSE
      musat_ad = musat_ad - (1._jprb - musun ** 2) * musat * ztmpx_ad / ztmpx
      musun_ad = musun_ad - (1._jprb - musat ** 2) * musun * ztmpx_ad / ztmpx
    ENDIF

  END SUBROUTINE int_phase_fn_ad

END SUBROUTINE rttov_opdpscattir_ad
