! Description:
!> @file
!!   TL of combined optical properties of aerosols and/or clouds.
!
!> @brief
!!   TL of combined optical properties of aerosols and/or clouds.
!!
!! @details
!!   TL of combined optical properties of aerosols and/or clouds.
!!
!! @param[in]     nlayers                number of layers in input profile
!! @param[in]     chanprof               specifies channels and profiles to simulate
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     dom_nstr               number of DOM streams
!! @param[in]     aux                    additional internal profile variables
!! @param[in,out] aux_tl                 additional internal profile variable perturbations
!! @param[in]     ircld                  computed cloud column data
!! @param[in]     profiles               input atmospheric profiles and surface variables
!! @param[in]     profiles_tl            atmospheric profiles and surface variable perturbations
!! @param[in]     profiles_int           profiles in internal units
!! @param[in]     profiles_int_tl        profile perturbations in internal units
!! @param[in]     aer_opt_param          explicit aerosol optical properties per channel (optional)
!! @param[in]     aer_opt_param_tl       explicit aerosol optical property perturbations (optional)
!! @param[in]     cld_opt_param          explicit cloud optical properties per channel (optional)
!! @param[in]     cld_opt_param_tl       explicit cloud optical property perturbations (optional)
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
!! @param[in]     raytracing_tl          TL of raytracing structure
!! @param[in]     trans_scatt_ir         computed optical depths
!! @param[in,out] trans_scatt_ir_tl      computed optical depth perturbations
!! @param[in]     trans_scatt_ir_dyn     computed optical depths per cloud column and phase functions
!! @param[in,out] trans_scatt_ir_dyn_tl  computed optical depth perturbations per cloud column and phase function perturbations
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
SUBROUTINE rttov_opdpscattir_tl( &
              nlayers,             &
              chanprof,            &
              opts,                &
              dom_nstr,            &
              aux,                 &
              aux_tl,              &
              ircld,               &
              profiles,            &
              profiles_tl,         &
              profiles_int,        &
              profiles_int_tl,     &
              aer_opt_param,       &
              aer_opt_param_tl,    &
              cld_opt_param,       &
              cld_opt_param_tl,    &
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
              raytracing_tl,       &
              trans_scatt_ir,      &
              trans_scatt_ir_tl,   &
              trans_scatt_ir_dyn,  &
              trans_scatt_ir_dyn_tl)

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
       spline_interp_tl,            &
       normalise_tl,                &
       calc_legendre_coef_gauss_tl
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options),               INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                INTENT(IN)    :: dom_nstr
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux_tl
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_int_tl(SIZE(profiles))
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: aer_opt_param_tl
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: cld_opt_param
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: cld_opt_param_tl
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
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_dyn_tl
!INTF_END

#include "rttov_baran2014_calc_optpar_tl.interface"
#include "rttov_baran2018_calc_optpar_tl.interface"
#include "rttov_baran_calc_phase_tl.interface"

  LOGICAL(KIND=jplm) :: do_dom, do_chou_scaling, do_single_scatt, do_mfasis, do_aer_or_ray_dom
  LOGICAL(KIND=jplm) :: do_dom_chan, do_chou_scaling_chan, do_single_scatt_chan, do_mfasis_chan, do_chan
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: j, k, k1, k2
  INTEGER(KIND=jpim) :: prof, chan, phchan
  INTEGER(KIND=jpim) :: clw_scheme, ice_scheme, col, coli, iaer, icld
  INTEGER(KIND=jpim) :: lev, lay
  REAL(KIND=jprb)    :: opd_tl, opdsun_tl
  REAL(KIND=jprb)    :: dgfrac, dgfrac_tl
  REAL(KIND=jprb)    :: abso, abso_tl
  REAL(KIND=jprb)    :: sca, sca_tl
  REAL(KIND=jprb)    :: sumsca, sumsca_tl
  REAL(KIND=jprb)    :: bpr, bpr_tl
  REAL(KIND=jprb)    :: asym, asym_tl
  REAL(KIND=jprb)    :: clw, clw_tl, tmpval, ray_phup, cosscata
  REAL(KIND=jprb)    :: afac, sfac, gfac, frach, frach_tl
  REAL(KIND=jprb)    :: pfac1(0:dom_nstr)
  REAL(KIND=jprb)    :: pfac2(coef_scatt%optp_aer%nphangle)

  INTEGER(KIND=jpim) :: nmom
  REAL(KIND=jprb)    :: legcoef(0:dom_nstr), thislegcoef(0:dom_nstr)
  REAL(KIND=jprb)    :: legcoefcld(0:dom_nstr)
  REAL(KIND=jprb)    :: legcoef_tl(0:dom_nstr), thislegcoef_tl(0:dom_nstr)
  REAL(KIND=jprb)    :: aer_pha(coef_scatt%optp_aer%nphangle)
  REAL(KIND=jprb)    :: aer_pha_tl(coef_scatt%optp_aer%nphangle)
  REAL(KIND=jprb)    :: wcldeff_pha(coef_scatt%optp_wcl_deff%nphangle)
  REAL(KIND=jprb)    :: wcldeff_pha_tl(coef_scatt%optp_wcl_deff%nphangle)
  REAL(KIND=jprb)    :: icl_pha(coef_scatt%optp_icl_baum%nphangle)
  REAL(KIND=jprb)    :: icl_pha_tl(coef_scatt%optp_icl_baum%nphangle)
  REAL(KIND=jprb)    :: zminphadiff
  REAL(KIND=jprb)    :: relazi, musat, musun
  REAL(KIND=jprb)    :: musat_tl, musun_tl, phasint_tl
  INTEGER(KIND=jpim) :: thisnphangle
  REAL(KIND=jprb)    :: thisphangle(nphangle_hires), thiscosphangle(nphangle_hires)
  REAL(KIND=jprb)    :: baran_pha(1:nphangle_hires), baran_pha_tl(1:nphangle_hires)
  REAL(KIND=jprb)    :: baran_pha_interp(baran_ngauss), baran_pha_interp_tl(baran_ngauss)

  TYPE(rttov_scatt_ir_aercld), POINTER :: aer, aer_tl, cld, cld_tl
  TYPE(rttov_optp_data),       POINTER :: optp_data

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)

  do_dom = (do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
           (do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
  do_chou_scaling = do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
  do_single_scatt = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
  do_mfasis = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
  do_aer_or_ray_dom = opts%rt_ir%addaerosl .OR. do_rayleigh_dom

  IF (do_thermal) trans_scatt_ir_tl%opdpacl = 0._jprb
  IF (do_solar .AND. .NOT. do_mfasis) THEN
    trans_scatt_ir_tl%opdpaclsun = 0._jprb
    trans_scatt_ir_tl%phup = 0._jprb
  ENDIF
  IF (do_single_scatt) trans_scatt_ir_tl%phdo = 0._jprb
  IF (do_dom) THEN
    DO j = 1, nchanprof
      IF (ASSOCIATED(trans_scatt_ir_dyn_tl%phasefn(j)%legcoef)) &
          trans_scatt_ir_dyn_tl%phasefn(j)%legcoef = 0._jprb
    ENDDO
  ENDIF
  IF (do_mfasis) trans_scatt_ir_tl%opdpext = 0._jprb

  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF AEROSOLS
  !----------------------------------------------------------------------------

  IF (do_aer_or_ray_dom) THEN

    aer    => trans_scatt_ir%aer
    aer_tl => trans_scatt_ir_tl%aer

    IF (do_thermal .OR. .NOT. do_mfasis) THEN
      aer_tl%opdpabs = 0._jprb
      aer_tl%opdpsca = 0._jprb
    ENDIF
    IF (do_chou_scaling) aer_tl%opdpscabpr = 0._jprb
    IF (do_solar .AND. .NOT. do_mfasis) aer_tl%phtotup = 0._jprb
    IF (do_single_scatt) aer_tl%phtotdo = 0._jprb
    IF (do_dom) aer_tl%sca = 0._jprb

    IF (opts%rt_ir%user_aer_opt_param) THEN

      CALL calc_user_opt_param_tl(aer_opt_param, aer_opt_param_tl, aer, aer_tl, 0_jpim, do_rayleigh_dom)

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

        DO lay = 1, nlayers

          !----------------------------------------------------------------------
          ! Calculate combined aerosol parameters for layer
          !----------------------------------------------------------------------
          IF (do_dom_chan) THEN
            sumsca = 0._jprb
            sumsca_tl = 0._jprb
            legcoef(:) = 0._jprb
            legcoef_tl(:) = 0._jprb

            IF (do_rayleigh_dom .AND. solar(j) .AND. &
                10000._jprb / coef%ff_cwn(chan) <= opts%rt_ir%rayleigh_max_wavelength .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN
              aer_tl%opdpsca(lay,j) = trans_scatt_ir_tl%ray_sca(lay,j)
              aer_tl%phtotup(lay,j) = ray_phup * trans_scatt_ir_tl%ray_sca(lay,j)
              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)
              sca_tl = (trans_scatt_ir_tl%ray_sca(lay,j) - raytracing_tl%ltick(lay,prof) * sca) / &
                       raytracing%ltick(lay,prof)
              sumsca = sca
              sumsca_tl = sca_tl
              legcoef(0) = 1._jprb * sca
              legcoef(2) = 0.5_jprb * sca
              legcoef_tl(0) = 1._jprb * sca_tl
              legcoef_tl(2) = 0.5_jprb * sca_tl
            ENDIF
          ENDIF

          DO iaer = 1, coef_scatt%optp_aer%ntypes
            IF (profiles_int(prof)%aerosols(iaer,lay) <= 0._jprb) CYCLE

            optp_data => coef_scatt%optp_aer%data(iaer)

            IF (do_dom_chan) THEN
              thislegcoef(:) = 0._jprb
              thislegcoef_tl(:) = 0._jprb
            ENDIF

            k = optp_data%nrelhum
            IF (k /= 1 .AND. aux%relhum(lay,prof) <= optp_data%relhum(k)) THEN
              ! Interpolate scattering parameters to actual value of relative humidity
              DO k = 1, coef_scatt%optp_aer%data(iaer)%nrelhum - 1
                IF (aux%relhum(lay,prof) >= optp_data%relhum(k) .AND. &
                    aux%relhum(lay,prof) <= optp_data%relhum(k+1)) THEN
                  frach = (aux%relhum(lay,prof) - optp_data%relhum(k)) / &
                          (optp_data%relhum(k+1) - optp_data%relhum(k))
                  afac  = (optp_data%abs(k+1,1,chan) - optp_data%abs(k,1,chan))
                  abso  = optp_data%abs(k,1,chan) + afac * frach
                  sfac  = (optp_data%sca(k+1,1,chan) - optp_data%sca(k,1,chan))
                  sca   = optp_data%sca(k,1,chan) + sfac * frach

                  frach_tl = aux_tl%relhum(lay,prof) / &
                             (optp_data%relhum(k+1) - optp_data%relhum(k))
                  abso_tl  = afac * frach_tl
                  sca_tl   = sfac * frach_tl

                  IF (do_chou_scaling_chan) THEN
                    gfac  = (optp_data%bpr(k+1,1,chan) - optp_data%bpr(k,1,chan))
                    bpr   = optp_data%bpr(k,1,chan) + gfac * frach
                    bpr_tl = gfac * frach_tl
                  ENDIF

                  IF (do_dom_chan) THEN
                    nmom = MIN(MAX(optp_data%nmom(k,chan), optp_data%nmom(k+1,chan)), dom_nstr)
                    pfac1(0:nmom) = (optp_data%legcoef(1:nmom+1,k+1,1,chan) - &
                                     optp_data%legcoef(1:nmom+1,k,1,chan))
                    thislegcoef(0:nmom) = optp_data%legcoef(1:nmom+1,k,1,chan) + pfac1(0:nmom) * frach
                    thislegcoef_tl(0:nmom) = pfac1(0:nmom) * frach_tl
                  ENDIF

                  IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
                    phchan = coef_scatt%optp_aer%chan_pha_index(chan)
                    pfac2(:) = (optp_data%pha(:,k+1,1,phchan) - optp_data%pha(:,k,1,phchan))
                    aer_pha(:) = optp_data%pha(:,k,1,phchan) + pfac2(:) * frach
                    aer_pha_tl(:) = pfac2(:) * frach_tl
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ELSE
              ! Particle doesn't change with rel. hum. (k=1) or rel. hum. exceeds max (k=max_rh_index)
              abso = optp_data%abs(k,1,chan)
              sca  = optp_data%sca(k,1,chan)

              abso_tl = 0._jprb
              sca_tl  = 0._jprb

              IF (do_chou_scaling_chan) THEN
                bpr = optp_data%bpr(k,1,chan)
                bpr_tl = 0._jprb
              ENDIF

              IF (do_dom_chan) THEN
                nmom = MIN(optp_data%nmom(k,chan), dom_nstr)
                thislegcoef(0:nmom) = optp_data%legcoef(1:nmom+1,k,1,chan)
                thislegcoef_tl(0:nmom) = 0._jprb
              ENDIF

              IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
                phchan = coef_scatt%optp_aer%chan_pha_index(chan)
                aer_pha(:) = optp_data%pha(:,k,1,phchan)
                aer_pha_tl(:) = 0._jprb
              ENDIF
            ENDIF

            IF (do_mfasis_chan) THEN
              trans_scatt_ir_tl%opdpext(iaer,lay,j) = coef%ff_gam(chan) * &
                (profiles_int_tl(prof)%aerosols(iaer,lay) * raytracing%ltick(lay,prof) * (abso + sca) + &
                 profiles_int(prof)%aerosols(iaer,lay) * raytracing_tl%ltick(lay,prof) * (abso + sca) + &
                 profiles_int(prof)%aerosols(iaer,lay) * raytracing%ltick(lay,prof) * (abso_tl + sca_tl))
            ELSEIF (do_chan) THEN
              !------------------------------------------------------------------
              ! Accumulate total aerosol optical parameters (all particles)
              !------------------------------------------------------------------
              ! For pre-defined particle types OD calculated as (aerosol amount in layer) *
              !   (aerosol ext coeff for given layer RH) * (layer thickness)
              ! OD is accumulated for each aerosol type - arrays are zeroed above.
              aer_tl%opdpabs(lay,j) = aer_tl%opdpabs(lay,j) + &
                profiles_int_tl(prof)%aerosols(iaer,lay) * abso * raytracing%ltick(lay,prof) + &
                profiles_int(prof)%aerosols(iaer,lay) * abso_tl * raytracing%ltick(lay,prof) + &
                profiles_int(prof)%aerosols(iaer,lay) * abso * raytracing_tl%ltick(lay,prof)

              aer_tl%opdpsca(lay,j) = aer_tl%opdpsca(lay,j) + &
                profiles_int_tl(prof)%aerosols(iaer,lay) * sca * raytracing%ltick(lay,prof) + &
                profiles_int(prof)%aerosols(iaer,lay) * sca_tl * raytracing%ltick(lay,prof) + &
                profiles_int(prof)%aerosols(iaer,lay) * sca * raytracing_tl%ltick(lay,prof)

              IF (do_chou_scaling_chan) THEN
                aer_tl%opdpscabpr(lay,j) = aer_tl%opdpscabpr(lay,j) + &
                  profiles_int_tl(prof)%aerosols(iaer,lay) * sca * bpr * raytracing%ltick(lay,prof) + &
                  profiles_int(prof)%aerosols(iaer,lay) * sca_tl * bpr * raytracing%ltick(lay,prof) + &
                  profiles_int(prof)%aerosols(iaer,lay) * sca * bpr_tl * raytracing%ltick(lay,prof) + &
                  profiles_int(prof)%aerosols(iaer,lay) * sca * bpr * raytracing_tl%ltick(lay,prof)
              ENDIF

              IF (do_dom_chan) THEN
                sumsca = sumsca + profiles_int(prof)%aerosols(iaer,lay) * sca
                sumsca_tl = sumsca_tl + profiles_int_tl(prof)%aerosols(iaer,lay) * sca + &
                                        profiles_int(prof)%aerosols(iaer,lay) * sca_tl
                legcoef(:) = legcoef(:) + thislegcoef(:) * profiles_int(prof)%aerosols(iaer,lay) * sca
                legcoef_tl(:) = legcoef_tl(:) + &
                  thislegcoef_tl(:) * profiles_int(prof)%aerosols(iaer,lay) * sca + &
                  thislegcoef(:) * profiles_int_tl(prof)%aerosols(iaer,lay) * sca + &
                  thislegcoef(:) * profiles_int(prof)%aerosols(iaer,lay) * sca_tl
              ENDIF

              IF (solar(j)) THEN
                musat =  1._jprb / raytracing%pathsat(lay,prof)
                musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
                musun = -1._jprb / raytracing%pathsun(lay,prof)
                musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
                zminphadiff = coef_scatt%optp_aer%phfn_int%zminphadiff * rad2deg

                CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                     aer_pha, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                     coef_scatt%optp_aer%phfn_int%iphangle, phasint_tl, aer_pha_tl)

                aer_tl%phtotup(lay,j) = aer_tl%phtotup(lay,j) + &
                    profiles_int_tl(prof)%aerosols(iaer,lay) * &
                    aer%phintup(iaer,lay,j) * sca * raytracing%ltick(lay,prof) + &
                    profiles_int(prof)%aerosols(iaer,lay) * &
                    phasint_tl * sca * raytracing%ltick(lay,prof) + &
                    profiles_int(prof)%aerosols(iaer,lay) * &
                    aer%phintup(iaer,lay,j) * sca_tl * raytracing%ltick(lay,prof) + &
                    profiles_int(prof)%aerosols(iaer,lay) * &
                    aer%phintup(iaer,lay,j) * sca * raytracing_tl%ltick(lay,prof)

                IF (do_single_scatt_chan) THEN
                  musat = -musat
                  musat_tl = -musat_tl

                  CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                     aer_pha, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                     coef_scatt%optp_aer%phfn_int%iphangle, phasint_tl, aer_pha_tl)

                  aer_tl%phtotdo(lay,j) = aer_tl%phtotdo(lay,j) + &
                      profiles_int_tl(prof)%aerosols(iaer,lay) * &
                      aer%phintdo(iaer,lay,j) * sca * raytracing%ltick(lay,prof) + &
                      profiles_int(prof)%aerosols(iaer,lay) * &
                      phasint_tl * sca * raytracing%ltick(lay,prof) + &
                      profiles_int(prof)%aerosols(iaer,lay) * &
                      aer%phintdo(iaer,lay,j) * sca_tl * raytracing%ltick(lay,prof) + &
                      profiles_int(prof)%aerosols(iaer,lay) * &
                      aer%phintdo(iaer,lay,j) * sca * raytracing_tl%ltick(lay,prof)
                ENDIF
              ENDIF
            ENDIF
          ENDDO ! aer types

          !------------------------------------------------------------------
          ! Calculate final phase function Leg. coefs for all aerosol types
          !------------------------------------------------------------------
          IF (do_dom_chan) THEN
            IF (ABS(sumsca) > realtol) THEN
              aer_tl%sca(lay,j) = sumsca_tl
              trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,0,lay) = &
                  legcoef_tl(:) / sumsca - sumsca_tl * legcoef(:) / sumsca ** 2_jpim
            ENDIF
          ENDIF
        ENDDO ! layers
      ENDDO ! chanprof
    ENDIF ! scaer coef file

    DO j = 1, nchanprof
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
      do_single_scatt_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
      do_mfasis_chan = solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis

      !------------------------------------------------------------------------
      ! Calculate total aerosol optical depths
      !------------------------------------------------------------------------
      IF (thermal(j)) THEN
        IF (do_chou_scaling_chan) THEN
          ! Chou-scaled optical depth for thermal channels
          aer_tl%opdp(:,j) = aer_tl%opdpabs(:,j) + aer_tl%opdpscabpr(:,j)
        ELSE
          ! Full optical depth for thermal channels
          aer_tl%opdp(:,j) = aer_tl%opdpabs(:,j) + aer_tl%opdpsca(:,j)
        ENDIF
      ENDIF

      IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
        WHERE (ABS(aer%opdpsca(:,j)) > realtol)
          trans_scatt_ir_tl%phup(0,:,j) = aer_tl%phtotup(:,j) / aer%opdpsca(:,j) - &
              aer_tl%opdpsca(:,j) * aer%phtotup(:,j) / aer%opdpsca(:,j)**2
        ENDWHERE
        IF (do_single_scatt_chan) THEN
          WHERE (ABS(aer%opdpsca(:,j)) > realtol)
            trans_scatt_ir_tl%phdo(0,:,j) = aer_tl%phtotdo(:,j) / aer%opdpsca(:,j) - &
                aer_tl%opdpsca(:,j) * aer%phtotdo(:,j) / aer%opdpsca(:,j)**2
          ENDWHERE
        ENDIF

        ! Full optical depth for solar channels
        aer_tl%opdpsun(:,j) = aer_tl%opdpabs(:,j) + aer_tl%opdpsca(:,j)
      ENDIF
    ENDDO ! chanprof
  ENDIF ! do_aer_or_ray_dom


  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF CLOUDS
  !----------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN

    cld    => trans_scatt_ir%cld
    cld_tl => trans_scatt_ir_tl%cld

    IF (do_thermal .OR. .NOT. do_mfasis) THEN
      cld_tl%opdpabs = 0._jprb
      cld_tl%opdpsca = 0._jprb
    ENDIF
    IF (do_chou_scaling) cld_tl%opdpscabpr = 0._jprb
    IF (do_solar .AND. .NOT. do_mfasis) cld_tl%phtotup = 0._jprb
    IF (do_single_scatt) cld_tl%phtotdo = 0._jprb
    IF (do_dom) cld_tl%sca = 0._jprb

    IF (opts%rt_ir%user_cld_opt_param) THEN

      CALL calc_user_opt_param_tl(cld_opt_param, cld_opt_param_tl, cld, cld_tl, 1_jpim, .FALSE._jplm)

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

        DO lay = 1, nlayers
          IF (do_dom_chan) THEN
            sumsca = 0._jprb
            sumsca_tl = 0._jprb
            legcoef = 0._jprb
            legcoef_tl = 0._jprb
          ENDIF

          DO icld = 1, ncldtyp

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

                clw_tl = SUM(profiles_int_tl(prof)%cloud(1:nwcl_max,lay))

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
                  dgfrac_tl = aux_tl%clw_dg(lay,prof) / &
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
                  dgfrac_tl = 0._jprb
                ENDIF

                IF (do_mfasis_chan) THEN
                  trans_scatt_ir_tl%opdpext(icld,lay,j) = &
                    coef%ff_gam(chan) * &
                    ((clw_tl * raytracing%ltick(lay,prof) + &
                      clw * raytracing_tl%ltick(lay,prof)) * &
                     (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                      optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) + &
                    clw * raytracing%ltick(lay,prof) * &
                    dgfrac_tl * (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan) + &
                                 optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                ELSEIF (do_chan) THEN
                  cld_tl%opdpabs(lay,j) = cld_tl%opdpabs(lay,j) + &
                    (clw_tl * raytracing%ltick(lay,prof) + &
                     clw * raytracing_tl%ltick(lay,prof)) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                    (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))) + &
                    clw * raytracing%ltick(lay,prof) * &
                    dgfrac_tl * (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                    (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * clw
                  sca_tl = &
                    dgfrac_tl * (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)) * clw + &
                    (optp_data%sca(1,k1,chan) + dgfrac * &
                    (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * clw_tl
                  cld_tl%partsca(icld,lay,j) = sca_tl * raytracing%ltick(lay,prof) + &
                    sca * raytracing_tl%ltick(lay,prof)
                  cld_tl%opdpsca(lay,j) = cld_tl%opdpsca(lay,j) + cld_tl%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
                    cld_tl%partbpr(icld,lay,j) = &
                      dgfrac_tl * (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                    cld_tl%opdpscabpr(lay,j) = cld_tl%opdpscabpr(lay,j) + &
                      cld_tl%partbpr(icld,lay,j) * cld%partsca(icld,lay,j) + &
                      cld%partbpr(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    sumsca_tl = sumsca_tl + sca_tl
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                      optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                      (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))
                    thislegcoef_tl(0:nmom) = &
                      dgfrac_tl * (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                    legcoef_tl(:) = legcoef_tl(:) + thislegcoef_tl(:) * sca + thislegcoef(:) * sca_tl
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
                    zminphadiff = coef_scatt%optp_wcl_deff%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_wcl_deff%chan_pha_index(chan)
                    wcldeff_pha(:) = &
                      optp_data%pha(:,1,k1,phchan) + dgfrac * &
                      (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
                    wcldeff_pha_tl(:) = &
                      dgfrac_tl * (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))

                    CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                         wcldeff_pha, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                         coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint_tl, wcldeff_pha_tl)

                    cld_tl%phtotup(lay,j) = cld_tl%phtotup(lay,j) + &
                        phasint_tl * cld%partsca(icld,lay,j) + &
                        cld%phintup(icld,lay,j) * cld_tl%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      musat_tl = -musat_tl

                      CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                           wcldeff_pha, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                           coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint_tl, wcldeff_pha_tl)

                    cld_tl%phtotdo(lay,j) = cld_tl%phtotdo(lay,j) + &
                        phasint_tl * cld%partsca(icld,lay,j) + &
                        cld%phintdo(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                    ENDIF
                  ENDIF
                ENDIF

              ELSE
                !--------------------------------------------------------------
                ! Water clouds - OPAC optical parameters
                !--------------------------------------------------------------
                IF (.NOT. profiles_int(prof)%cloud(icld,lay) > 0._jprb) CYCLE

                optp_data => coef_scatt%optp_wcl_opac%data(icld)

                IF (do_mfasis_chan) THEN
                  trans_scatt_ir_tl%opdpext(icld,lay,j) = &
                    coef%ff_gam(chan) * optp_data%confac * &
                    (optp_data%abs(1,1,chan) + optp_data%sca(1,1,chan)) * &
                    (profiles_int_tl(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) + &
                     profiles_int(prof)%cloud(icld,lay) * raytracing_tl%ltick(lay,prof))

                ELSEIF (do_chan) THEN
                  cld_tl%opdpabs(lay,j) = cld_tl%opdpabs(lay,j) + &
                      optp_data%confac * optp_data%abs(1,1,chan) * &
                      (profiles_int_tl(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) + &
                       profiles_int(prof)%cloud(icld,lay) * raytracing_tl%ltick(lay,prof))
                  cld_tl%partsca(icld,lay,j) = &
                      optp_data%confac * optp_data%sca(1,1,chan) * &
                      (profiles_int_tl(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) + &
                       profiles_int(prof)%cloud(icld,lay) * raytracing_tl%ltick(lay,prof))
                  cld_tl%opdpsca(lay,j) = cld_tl%opdpsca(lay,j) + cld_tl%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
!                     bpr_tl = 0._jprb !optp_data%bpr(1,1,chan)
                    cld_tl%opdpscabpr(lay,j) = cld_tl%opdpscabpr(lay,j) + &
                        cld%partbpr(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    sca = optp_data%sca(1,1,chan) * &
                          profiles_int(prof)%cloud(icld,lay) * optp_data%confac
                    sca_tl = optp_data%sca(1,1,chan) * &
                             profiles_int_tl(prof)%cloud(icld,lay) * optp_data%confac
                    sumsca = sumsca + sca
                    sumsca_tl = sumsca_tl + sca_tl
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    legcoef(0:nmom) = legcoef(0:nmom) + optp_data%legcoef(1:nmom+1,1,1,chan) * sca
                    legcoef_tl(0:nmom) = legcoef_tl(0:nmom) + optp_data%legcoef(1:nmom+1,1,1,chan) * sca_tl
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
                    zminphadiff = coef_scatt%optp_wcl_opac%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_wcl_opac%chan_pha_index(chan)
                    CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                         optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                         coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint_tl)

                    cld_tl%phtotup(lay,j) = cld_tl%phtotup(lay,j) + &
                        phasint_tl * cld%partsca(icld,lay,j) + &
                        cld%phintup(icld,lay,j) * cld_tl%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      musat_tl = -musat_tl

                      CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                         optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                         coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint_tl)

                      cld_tl%phtotdo(lay,j) = cld_tl%phtotdo(lay,j) + &
                          phasint_tl * cld%partsca(icld,lay,j) + &
                          cld%phintdo(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                    ENDIF
                  ENDIF
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
                  dgfrac_tl = aux_tl%ice_dg(lay,prof) / &
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
                  dgfrac_tl = 0._jprb
                ENDIF

                IF (do_mfasis_chan) THEN
                  trans_scatt_ir_tl%opdpext(icld,lay,j) = &
                    coef%ff_gam(chan) * &
                    ((profiles_int_tl(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) + &
                      profiles_int(prof)%cloud(icld,lay) * raytracing_tl%ltick(lay,prof)) * &
                     (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                      optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) + &
                    profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                    dgfrac_tl * (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan) + &
                                 optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))

                ELSEIF (do_chan) THEN
                  cld_tl%opdpabs(lay,j) = cld_tl%opdpabs(lay,j) + &
                      (profiles_int_tl(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) + &
                       profiles_int(prof)%cloud(icld,lay) * raytracing_tl%ltick(lay,prof)) * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))) + &
                      profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                      dgfrac_tl * (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan))

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * profiles_int(prof)%cloud(icld,lay)
                  sca_tl = &
                      dgfrac_tl * (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)) * &
                      profiles_int(prof)%cloud(icld,lay) + &
                      (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * &
                      profiles_int_tl(prof)%cloud(icld,lay)
                  cld_tl%partsca(icld,lay,j) = &
                      sca_tl * raytracing%ltick(lay,prof) + sca * raytracing_tl%ltick(lay,prof)
                  cld_tl%opdpsca(lay,j) = cld_tl%opdpsca(lay,j) + cld_tl%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
                    bpr_tl = dgfrac_tl * (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                    cld_tl%opdpscabpr(lay,j) = cld_tl%opdpscabpr(lay,j) + &
                        bpr_tl * cld%partsca(icld,lay,j) + &
                        cld%partbpr(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    sumsca_tl = sumsca_tl + sca_tl
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                        optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                        (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))
                    thislegcoef_tl(0:nmom) = &
                        dgfrac_tl * (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                    legcoef_tl(:) = legcoef_tl(:) + thislegcoef_tl(:) * sca + thislegcoef(:) * sca_tl
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
                    zminphadiff = coef_scatt%optp_icl_baum%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_icl_baum%chan_pha_index(chan)
                    icl_pha(:) = &
                        optp_data%pha(:,1,k1,phchan) + dgfrac * &
                        (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
                    icl_pha_tl(:) = &
                        dgfrac_tl * (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))

                    CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                         icl_pha, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                         coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint_tl, icl_pha_tl)

                    cld_tl%phtotup(lay,j) = cld_tl%phtotup(lay,j) + &
                        phasint_tl * cld%partsca(icld,lay,j) + &
                        cld%phintup(icld,lay,j) * cld_tl%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      musat_tl = -musat_tl

                      CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                           icl_pha, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                           coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint_tl, icl_pha_tl)

                      cld_tl%phtotdo(lay,j) = cld_tl%phtotdo(lay,j) + &
                          phasint_tl * cld%partsca(icld,lay,j) + &
                          cld%phintdo(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                    ENDIF
                  ENDIF
                ENDIF
              ELSE
                !--------------------------------------------------------------
                ! Optical parameters computed using Baran scheme
                !--------------------------------------------------------------
                IF (ice_scheme == ice_scheme_baran2018) THEN
                  CALL rttov_baran2018_calc_optpar_tl(coef_scatt%optp_icl_baran2018, chan, &
                        profiles(prof)%t(lay), profiles_int(prof)%cloud(icld,lay), &
                        profiles_tl(prof)%t(lay), profiles_int_tl(prof)%cloud(icld,lay), &
                        abso, sca, bpr, asym, abso_tl, sca_tl, bpr_tl, asym_tl)
                ELSE
                  CALL rttov_baran2014_calc_optpar_tl(coef_scatt%optp_icl_baran2014, chan, &
                        profiles(prof)%t(lay), profiles_int(prof)%cloud(icld,lay), &
                        profiles_tl(prof)%t(lay), profiles_int_tl(prof)%cloud(icld,lay), &
                        abso, sca, bpr, asym, abso_tl, sca_tl, bpr_tl, asym_tl)
                ENDIF

                cld_tl%opdpabs(lay,j) = cld_tl%opdpabs(lay,j) + &
                    abso_tl * raytracing%ltick(lay,prof) + abso * raytracing_tl%ltick(lay,prof)
                cld_tl%partsca(icld,lay,j) = &
                    sca_tl * raytracing%ltick(lay,prof) + sca * raytracing_tl%ltick(lay,prof)
                cld_tl%opdpsca(lay,j) = cld_tl%opdpsca(lay,j) + cld_tl%partsca(icld,lay,j)

                IF (do_chou_scaling_chan) THEN
                  cld_tl%opdpscabpr(lay,j) = cld_tl%opdpscabpr(lay,j) + &
                      bpr_tl * cld%partsca(icld,lay,j) + &
                      cld%partbpr(icld,lay,j) * cld_tl%partsca(icld,lay,j)
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
                  CALL rttov_baran_calc_phase_tl(asym, asym_tl, thisphangle(1:thisnphangle), &
                        baran_pha(1:thisnphangle), baran_pha_tl(1:thisnphangle))

                  ! Compute Legendre coefficients
                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    sumsca_tl = sumsca_tl + sca_tl
                    nmom = dom_nstr

                    CALL spline_interp_tl(thisnphangle, thiscosphangle(thisnphangle:1:-1), &
                                          baran_pha(thisnphangle:1:-1), baran_pha_tl(thisnphangle:1:-1), &
                                          baran_ngauss, coef_scatt%optp_icl_baran2018%q, baran_pha_interp, baran_pha_interp_tl)

                    CALL normalise_tl(baran_ngauss, coef_scatt%optp_icl_baran2018%w, baran_pha_interp, baran_pha_interp_tl)

                    CALL calc_legendre_coef_gauss_tl(coef_scatt%optp_icl_baran2018%q, coef_scatt%optp_icl_baran2018%w, &
                                                     baran_pha_interp, baran_pha_interp_tl, dom_nstr, dom_nstr, &
                                                     nmom, thislegcoef(:), thislegcoef_tl(:))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                    legcoef_tl(:) = legcoef_tl(:) + thislegcoef_tl(:) * sca + thislegcoef(:) * sca_tl
                  ENDIF

                  ! Evaluate phase function for solar scattering
                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
                    zminphadiff = coef_scatt%optp_icl_baran2018%phfn_int%zminphadiff * rad2deg

                    CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                         baran_pha, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                         coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint_tl, baran_pha_tl)

                    cld_tl%phtotup(lay,j) = cld_tl%phtotup(lay,j) + &
                        phasint_tl * cld%partsca(icld,lay,j) + &
                        cld%phintup(icld,lay,j) * cld_tl%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                           baran_pha, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                           coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint_tl, baran_pha_tl)

                      cld_tl%phtotdo(lay,j) = cld_tl%phtotdo(lay,j) + &
                          phasint_tl * cld%partsca(icld,lay,j) + &
                          cld%phintdo(icld,lay,j) * cld_tl%partsca(icld,lay,j)
                    ENDIF
                  ENDIF
                ENDIF

              ENDIF ! ice_scheme
            ENDIF ! water or ice
          ENDDO ! cloud types

          !------------------------------------------------------------------
          ! Calculate final phase function Leg. coefs for all cloud types
          !------------------------------------------------------------------
          IF (do_dom_chan) THEN
            IF (ABS(sumsca) > realtol) THEN
              cld_tl%sca(lay,j) = sumsca_tl
              trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,1,lay) = &
                  legcoef_tl(:) / sumsca - sumsca_tl * legcoef(:) / sumsca ** 2_jpim
            ENDIF
          ENDIF
        ENDDO ! layers
      ENDDO ! chanprof
    ENDIF ! sccld coef file

    DO j = 1, nchanprof
      do_chou_scaling_chan = thermal(j) .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou

      !------------------------------------------------------------------------
      ! Calculate total cloud optical depths
      !------------------------------------------------------------------------
      IF (thermal(j)) THEN
        IF (do_chou_scaling_chan) THEN
          ! Chou-scaled optical depth for thermal channels
          cld_tl%opdp(:,j) = cld_tl%opdpabs(:,j) + cld_tl%opdpscabpr(:,j)
        ELSE
          ! Full optical depth for thermal channels
          cld_tl%opdp(:,j) = cld_tl%opdpabs(:,j) + cld_tl%opdpsca(:,j)
        ENDIF
      ENDIF

      ! Full optical depth for solar channels
      IF (solar(j) .AND. .NOT. do_mfasis) cld_tl%opdpsun(:,j) = cld_tl%opdpabs(:,j) + cld_tl%opdpsca(:,j)
    ENDDO ! chanprof
  ENDIF ! addclouds

  IF (do_mfasis .AND. .NOT. do_thermal) RETURN

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

    ! For layer-specific quantities store just the non-cloudy and cloudy values
    ! When used later the code looks up the appropriate value for each cloud column

    ! Determine layer total aerosol/cloud optical depths
    IF (thermal(j)) THEN
      IF (opts%rt_ir%addaerosl) THEN
        trans_scatt_ir_tl%opdpacl(0,:,j) = coef%ff_gam(chan) * &
            (aer_tl%opdp(:,j) * raytracing%pathsat(:,prof) + &
             aer%opdp(:,j) * raytracing_tl%pathsat(:,prof))
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir_tl%opdpacl(1,:,j) = &
            trans_scatt_ir_tl%opdpacl(0,:,j) + coef%ff_gam(chan) * &
            (cld_tl%opdp(:,j) * raytracing%pathsat(:,prof) + &
             cld%opdp(:,j) * raytracing_tl%pathsat(:,prof))
      ENDIF
    ENDIF
    IF (solar(j)) THEN
      IF (do_aer_or_ray_dom) THEN
        trans_scatt_ir_tl%opdpaclsun(0,:,j) = coef%ff_gam(chan) * &
            (aer_tl%opdpsun(:,j) * raytracing%patheff(:,prof) + &
             aer%opdpsun(:,j) * raytracing_tl%patheff(:,prof))
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir_tl%opdpaclsun(1,:,j) = &
            trans_scatt_ir_tl%opdpaclsun(0,:,j) + coef%ff_gam(chan) * &
            (cld_tl%opdpsun(:,j) * raytracing%patheff(:,prof) + &
             cld%opdpsun(:,j) * raytracing_tl%patheff(:,prof))
      ENDIF
    ENDIF

    IF (do_dom_chan .OR. do_single_scatt_chan) THEN
      ! Determine final layer *nadir* absorption and scattering optical depths
      IF (do_aer_or_ray_dom) THEN
        trans_scatt_ir_tl%opdpabs(0,:,j) = coef%ff_gam(chan) * aer_tl%opdpabs(:,j)
        trans_scatt_ir_tl%opdpsca(0,:,j) = coef%ff_gam(chan) * aer_tl%opdpsca(:,j)
      ELSE
        trans_scatt_ir_tl%opdpabs(0,:,j) = 0._jprb
        trans_scatt_ir_tl%opdpsca(0,:,j) = 0._jprb
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir_tl%opdpabs(1,:,j) = trans_scatt_ir_tl%opdpabs(0,:,j) + &
            coef%ff_gam(chan) * cld_tl%opdpabs(:,j)
        trans_scatt_ir_tl%opdpsca(1,:,j) = trans_scatt_ir_tl%opdpsca(0,:,j) + &
            coef%ff_gam(chan) * cld_tl%opdpsca(:,j)
      ENDIF

      DO lay = 1, nlayers

        IF (do_dom_chan) THEN
          ! Determine combined phase functions for aer+cld case
          IF (do_aer_or_ray_dom .AND. opts%rt_ir%addclouds) THEN
            IF (aer%sca(lay,j) + cld%sca(lay,j) > 0._jprb) THEN

              ! NB In the direct model trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) initially contains
              !    just the cloud phase fn. If aerosols are also present then it is overwritten
              !    with the combined phase function. We need to be careful about this in the TL/AD/K.
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

              trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,1,lay) = &
                  (trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,0,lay) * aer%sca(lay,j) + &
                   trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay) * aer_tl%sca(lay,j) + &
                   trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,1,lay) * cld%sca(lay,j) + &
                   legcoefcld(:) * cld_tl%sca(lay,j) - &
                   trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) * & ! Re-use direct calculation
                   (aer_tl%sca(lay,j) + cld_tl%sca(lay,j))) / &
                  (aer%sca(lay,j) + cld%sca(lay,j))
            ENDIF
          ENDIF
        ENDIF

        IF (solar(j)) THEN
          IF (opts%rt_ir%addclouds .AND. do_aer_or_ray_dom) THEN
            IF (ABS(cld%opdpsca(lay,j) + aer%opdpsca(lay,j)) > realtol) THEN
              tmpval = 1._jprb / (cld%opdpsca(lay,j) + aer%opdpsca(lay,j))
              trans_scatt_ir_tl%phup(1,lay,j) = &
                  (aer_tl%phtotup(lay,j) + cld_tl%phtotup(lay,j)) * tmpval - &
                  (cld_tl%opdpsca(lay,j) + aer_tl%opdpsca(lay,j)) * &
                  (aer%phtotup(lay,j) + cld%phtotup(lay,j)) * tmpval**2
              IF (do_single_scatt_chan) &
                trans_scatt_ir_tl%phdo(1,lay,j) = &
                    (aer_tl%phtotdo(lay,j) + cld_tl%phtotdo(lay,j)) * tmpval - &
                    (cld_tl%opdpsca(lay,j) + aer_tl%opdpsca(lay,j)) * &
                    (aer%phtotdo(lay,j) + cld%phtotdo(lay,j)) * tmpval**2
            ENDIF
          ELSEIF (opts%rt_ir%addclouds) THEN
            IF (ABS(cld%opdpsca(lay,j)) > realtol) THEN
              tmpval = 1._jprb / cld%opdpsca(lay,j)
              trans_scatt_ir_tl%phup(1,lay,j) = &
                 cld_tl%phtotup(lay,j) * tmpval - &
                 cld_tl%opdpsca(lay,j) * cld%phtotup(lay,j) * tmpval**2
              IF (do_single_scatt_chan) &
                trans_scatt_ir_tl%phdo(1,lay,j) = &
                   cld_tl%phtotdo(lay,j) * tmpval - &
                   cld_tl%opdpsca(lay,j) * cld%phtotdo(lay,j) * tmpval**2
            ENDIF
          ENDIF
        ENDIF

      ENDDO ! layers
    ENDIF ! do dom or do single_scatt

    DO col = 0, ircld%ncolumn(prof)
      IF (thermal(j)) THEN
        opd_tl = 0._jprb
        trans_scatt_ir_dyn_tl%opdpac(1,col,j) = 0._jprb
      ENDIF
      IF (solar(j)) THEN
        opdsun_tl = 0._jprb
        trans_scatt_ir_dyn_tl%opdpacsun(1,col,j) = 0._jprb
      ENDIF
      IF (col == 0) THEN
        IF (do_aer_or_ray_dom) THEN
          DO lay = 1, nlayers
            lev = lay + 1
            IF (thermal(j)) THEN
              opd_tl = opd_tl + trans_scatt_ir_tl%opdpacl(0,lay,j)
              trans_scatt_ir_dyn_tl%opdpac(lev,col,j) = opd_tl
            ENDIF
            IF (solar(j)) THEN
              opdsun_tl = opdsun_tl + trans_scatt_ir_tl%opdpaclsun(0,lay,j)
              trans_scatt_ir_dyn_tl%opdpacsun(lev,col,j) = opdsun_tl
            ENDIF
          ENDDO ! layers
        ELSE
          IF (thermal(j)) trans_scatt_ir_dyn_tl%opdpac(:,col,j) = 0._jprb
          IF (solar(j)) trans_scatt_ir_dyn_tl%opdpacsun(:,col,j) = 0._jprb
        ENDIF
      ELSE
        DO lay = 1, nlayers
          lev = lay + 1
          coli = ircld%icldarr(col,lay,prof) ! This is 0 or 1 for clear(+aer)/cloud(+aer)
          IF (thermal(j)) THEN
            opd_tl = opd_tl + trans_scatt_ir_tl%opdpacl(coli,lay,j)
            trans_scatt_ir_dyn_tl%opdpac(lev,col,j) = opd_tl
          ENDIF
          IF (solar(j)) THEN
            opdsun_tl = opdsun_tl + trans_scatt_ir_tl%opdpaclsun(coli,lay,j)
            trans_scatt_ir_dyn_tl%opdpacsun(lev,col,j) = opdsun_tl
          ENDIF
        ENDDO ! layers
      ENDIF ! col == 0
    ENDDO ! col
  ENDDO ! channels

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_user_opt_param_tl(opt_param, opt_param_tl, aercld, aercld_tl, iaercld, do_rayleigh_dom)
    ! Process aer/cld user optical property inputs
    TYPE(rttov_opt_param),       INTENT(IN)           :: opt_param
    TYPE(rttov_opt_param),       INTENT(IN), OPTIONAL :: opt_param_tl
    TYPE(rttov_scatt_ir_aercld), INTENT(IN)           :: aercld
    TYPE(rttov_scatt_ir_aercld), INTENT(INOUT)        :: aercld_tl
    INTEGER(jpim),               INTENT(IN)           :: iaercld ! aer=>0, cld=>1
    LOGICAL(jplm),               INTENT(IN)           :: do_rayleigh_dom

    INTEGER(jpim) :: j, lay, prof, chan
    LOGICAL(jplm) :: do_chou_scaling_chan, do_dom_chan, do_single_scatt_chan
    LOGICAL(jplm) :: do_opt_param_lay, do_rayleigh_dom_chan, do_opt_param_tl
    REAL(jprb)    :: musat, musun, zminphadiff, phasint_tl
    REAL(jprb)    :: ray_lcoef(0:dom_nstr), opdpsca, opdpsca_tl, sca, sca_tl

    ! This is useful for the parallel interface which always passes opt_param_tl
    ! for convenience even when not allocated
    do_opt_param_tl = .FALSE.
    IF (PRESENT(opt_param_tl) .AND. .NOT. opts%dev%no_opt_param_tladk) &
      do_opt_param_tl = (ASSOCIATED(opt_param_tl%sca))

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

        IF (do_opt_param_tl) THEN

          IF (do_opt_param_lay) THEN
            aercld_tl%opdpabs(lay,j) = opt_param%abs(lay,j) * raytracing_tl%ltick(lay,prof) + &
                                       opt_param_tl%abs(lay,j) * raytracing%ltick(lay,prof)

            opdpsca = opt_param%sca(lay,j) * raytracing%ltick(lay,prof)
            opdpsca_tl = opt_param%sca(lay,j) * raytracing_tl%ltick(lay,prof) + &
                         opt_param_tl%sca(lay,j) * raytracing%ltick(lay,prof)
          ELSE
            opdpsca = 0._jprb
            opdpsca_tl = 0._jprb
          ENDIF

          IF (do_rayleigh_dom_chan .AND. &
              profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

            aercld_tl%opdpsca(lay,j) = trans_scatt_ir_tl%ray_sca(lay,j) + opdpsca_tl
            aercld_tl%phtotup(lay,j) = ray_phup * trans_scatt_ir_tl%ray_sca(lay,j)
          ELSEIF (do_opt_param_lay) THEN
            aercld_tl%opdpsca(lay,j) = opdpsca_tl
          ENDIF

          IF (do_opt_param_lay .AND. do_chou_scaling_chan) THEN
            aercld_tl%opdpscabpr(lay,j) = &
                opt_param%sca(lay,j) * opt_param%bpr(lay,j) * raytracing_tl%ltick(lay,prof) + &
                opt_param%sca(lay,j) * opt_param_tl%bpr(lay,j) * raytracing%ltick(lay,prof) + &
                opt_param_tl%sca(lay,j) * opt_param%bpr(lay,j) * raytracing%ltick(lay,prof)
          ENDIF

          IF (do_dom_chan) THEN
            IF (do_rayleigh_dom_chan .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)
              sca_tl = (trans_scatt_ir_tl%ray_sca(lay,j) - &
                        raytracing_tl%ltick(lay,prof) * sca) / raytracing%ltick(lay,prof)

              IF (do_opt_param_lay) THEN
                trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,iaercld,lay) = &
                  (ray_lcoef(:) * sca_tl + &
                   opt_param_tl%legcoef(1:dom_nstr+1,lay,j) * opt_param%sca(lay,j) + &
                   opt_param%legcoef(1:dom_nstr+1,lay,j) * opt_param_tl%sca(lay,j) - &
                   trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay) * &
                   (sca_tl + opt_param_tl%sca(lay,j))) / &
                  (sca + opt_param%sca(lay,j))
                aercld_tl%sca(lay,j) = sca_tl + opt_param_tl%sca(lay,j)
              ELSE
                aercld_tl%sca(lay,j) = sca_tl
              ENDIF

            ELSEIF (do_opt_param_lay) THEN
              aercld_tl%sca(lay,j) = opt_param_tl%sca(lay,j)
              trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,iaercld,lay) = opt_param_tl%legcoef(1:dom_nstr+1,lay,j)
            ENDIF
          ENDIF

          IF (do_opt_param_lay .AND. solar(j)) THEN
            musat =  1._jprb / raytracing%pathsat(lay,prof)
            musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
            musun = -1._jprb / raytracing%pathsun(lay,prof)
            musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
            zminphadiff = opt_param%phasefn_int%zminphadiff * rad2deg

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                 opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                 opt_param%phasefn_int%iphangle, phasint_tl, opt_param_tl%pha(:,lay,j))

            aercld_tl%phtotup(lay,j) = aercld_tl%phtotup(lay,j) + phasint_tl * opdpsca + &
                                       aercld%phintup(1,lay,j) * opdpsca_tl
            IF (do_single_scatt_chan) THEN
              musat = -musat
              musat_tl = -musat_tl

              CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                   opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                   opt_param%phasefn_int%iphangle, phasint_tl, opt_param_tl%pha(:,lay,j))

              aercld_tl%phtotdo(lay,j) = aercld_tl%phtotdo(lay,j) + phasint_tl * opdpsca + &
                                         aercld%phintdo(1,lay,j) * opdpsca_tl
            ENDIF
          ENDIF
        ELSE
          ! No opt_param_tl

          IF (do_opt_param_lay) THEN
            aercld_tl%opdpabs(lay,j) = opt_param%abs(lay,j) * raytracing_tl%ltick(lay,prof)
            opdpsca = opt_param%sca(lay,j) * raytracing%ltick(lay,prof)
            opdpsca_tl = opt_param%sca(lay,j) * raytracing_tl%ltick(lay,prof)
          ELSE
            opdpsca = 0._jprb
            opdpsca_tl = 0._jprb
          ENDIF

          IF (do_rayleigh_dom_chan .AND. &
              profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

            aercld_tl%opdpsca(lay,j) = trans_scatt_ir_tl%ray_sca(lay,j) + opdpsca_tl
            aercld_tl%phtotup(lay,j) = ray_phup * trans_scatt_ir_tl%ray_sca(lay,j)
          ELSEIF (do_opt_param_lay) THEN
            aercld_tl%opdpsca(lay,j) = opdpsca_tl
          ENDIF

          IF (do_opt_param_lay .AND. do_chou_scaling_chan) THEN
            aercld_tl%opdpscabpr(lay,j) = opt_param%sca(lay,j) * opt_param%bpr(lay,j) * raytracing_tl%ltick(lay,prof)
          ENDIF

          IF (do_dom_chan) THEN
            IF (do_rayleigh_dom_chan .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)
              sca_tl = (trans_scatt_ir_tl%ray_sca(lay,j) - &
                        raytracing_tl%ltick(lay,prof) * sca) / raytracing%ltick(lay,prof)

              IF (do_opt_param_lay) THEN
                trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,iaercld,lay) = &
                  (ray_lcoef(:) - trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay)) * &
                  sca_tl / (sca + opt_param%sca(lay,j))
              ENDIF
              aercld_tl%sca(lay,j) = sca_tl

            ELSEIF (do_opt_param_lay) THEN
              aercld_tl%sca(lay,j) = 0._jprb
              trans_scatt_ir_dyn_tl%phasefn(j)%legcoef(:,iaercld,lay) = 0._jprb
            ENDIF
          ENDIF

          IF (do_opt_param_lay .AND. solar(j)) THEN
            musat =  1._jprb / raytracing%pathsat(lay,prof)
            musat_tl = -raytracing_tl%pathsat(lay,prof) * musat**2
            musun = -1._jprb / raytracing%pathsun(lay,prof)
            musun_tl = raytracing_tl%pathsun(lay,prof) * musun**2
            zminphadiff = opt_param%phasefn_int%zminphadiff * rad2deg

            CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, 180.0_jprb - relazi, zminphadiff, &
                                 opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                 opt_param%phasefn_int%iphangle, phasint_tl)

            aercld_tl%phtotup(lay,j) = aercld_tl%phtotup(lay,j) + phasint_tl * opdpsca + &
                                       aercld%phintup(1,lay,j) * opdpsca_tl

            IF (do_single_scatt_chan) THEN
              musat = -musat
              musat_tl = -musat_tl

              CALL int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                                   opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                                   opt_param%phasefn_int%iphangle, phasint_tl)

              aercld_tl%phtotdo(lay,j) = aercld_tl%phtotdo(lay,j) + phasint_tl * opdpsca + &
                                         aercld%phintdo(1,lay,j) * opdpsca_tl
            ENDIF
          ENDIF
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE calc_user_opt_param_tl

  SUBROUTINE int_phase_fn_tl(musat, musat_tl, musun, musun_tl, relazi, zminphadiff, &
                             pha, cospha, ipha, phasint_tl, pha_tl)
    ! Interpolate phase function to scattering angle
    ! pha_tl may be omitted if the phase function is static
    REAL(KIND=jprb),           INTENT(IN)  :: musat, musat_tl, musun, musun_tl, relazi
    REAL(KIND=jprb),           INTENT(IN)  :: zminphadiff
    REAL(KIND=jprb),           INTENT(IN)  :: pha(:)
    REAL(KIND=jprb),           INTENT(IN)  :: cospha(:)
    INTEGER(KIND=jpim),        INTENT(IN)  :: ipha(:)
    REAL(KIND=jprb),           INTENT(OUT) :: phasint_tl
    REAL(KIND=jprb), OPTIONAL, INTENT(IN)  :: pha_tl(:)

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, ztmpx_tl
    REAL(KIND=jprb)    :: scattangle, scattangle_tl, deltap, deltap_tl, delta

    phasint_tl = 0._jprb

    ztmpx = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))
    IF (ABS(ztmpx) < realtol) THEN
      ztmpx_tl = 0._jprb
    ELSE
      ztmpx_tl = - ((1._jprb - musun ** 2) * musat * musat_tl + &
                    (1._jprb - musat ** 2) * musun * musun_tl) / ztmpx
    ENDIF

    scattangle = musat * musun + ztmpx * COS(relazi * deg2rad)
    scattangle_tl = musat_tl * musun + musat * musun_tl + ztmpx_tl * COS(relazi * deg2rad)
    ikk        = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk         = ipha(ikk) - 1_jpim
    deltap     = pha(kk + 1) - pha(kk)
    delta      = cospha(kk) - cospha(kk + 1)
    IF (PRESENT(pha_tl)) THEN
      deltap_tl   = pha_tl(kk + 1) - pha_tl(kk)
      phasint_tl  = pha_tl(kk) + deltap_tl * (cospha(kk) - scattangle) / delta - &
                                 deltap * scattangle_tl / delta
    ELSE
      phasint_tl  = - deltap * scattangle_tl / delta
    ENDIF

  END SUBROUTINE int_phase_fn_tl

END SUBROUTINE rttov_opdpscattir_tl
