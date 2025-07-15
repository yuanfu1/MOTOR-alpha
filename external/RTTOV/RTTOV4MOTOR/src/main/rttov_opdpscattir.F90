! Description:
!> @file
!!   Calculate combined optical properties of aerosols and/or clouds.
!
!> @brief
!!   Calculate combined optical properties of aerosols and/or clouds.
!!
!! @details
!!   This subroutine is complicated by the fact that each scattering
!!   option requires different quantities to be calculated and where
!!   multiple particle types are present in a layer the combined optical
!!   properties of all types must be calculated.
!!
!!   Aerosols: some predefined aerosol properties vary with relative
!!     humidity and are interpolated accordingly.
!!
!!   Clouds: the OPAC water cloud types are fairly straight-forward.
!!     The "Deff" properties are interpolated according to CLW Deff.
!!     Ice cloud properties are either taken from the sccld coefficient
!!     file and interpolated according to Deff or computed from the
!!     Baran parameterisation.
!!
!!   Phase fn Legendre coefficients are stored in the coefficient files
!!   for all channels and the explicit phase fns are stored for solar
!!   channels. For the Baran ice scheme the phase fn is computed from
!!   the asymmetry parameter via a parameterisation and the Legendre
!!   coefficients are computed from the phase fn via Gaussian quadrature.
!!
!!   For aerosols and/or clouds the user may supply optical properties
!!   explicitly: for these cases the code is much simpler as there is no
!!   need to combine optical properties for multiple particle types.
!!
!!   The output layer optical depths in opdpacl and opdpaclsun and the
!!   accumulated optical depths in opdpac and opdpacsun take the local
!!   path angle into account. These are used for Chou-scaling and solar
!!   single-scattering.
!!
!!   The output layer absorption and scattering optical depths in opdpabs
!!   and opdpsca are the *nadir* optical depths. These are used for DOM
!!   calculations and for the SSA calculation for solar single-scattering.
!!   In the latter case the local path angle is taken into account in
!!   transmit_solar.
!!
!!   Chou-scaling (IR only):
!!     outputs required are the total optical depth (with Chou-scaled
!!     scattering op dep) for the local path. The back-scattering parameter
!!     (bpr) is required and this is used to scale the scattering op deps.
!!
!!   Solar single-scattering (VIS/NIR only):
!!     outputs required are total optical depths for the local path, the
!!     absorption and scattering optical depths (these are computed at nadir
!!     for compatibility with DOM and the local path angle is introduced in
!!     transmit_solar), and the phase fn evaluated at scattering angles.
!!     Scattering optical depths are not Chou-scaled (bpr is not required).
!!     The phase fn is required, but not the phase fn Legendre expansion. The
!!     phase fn is interpolated for the upward- and downward-scattered beams.
!!
!!   DOM - IR:
!!     outputs are the nadir absorption and scattering optical depths and the
!!     first dom_nstreams Legendre coefficients for the combined phase function
!!     for all scattering particles.
!!
!!   DOM - solar:
!!     outputs are same as DOM IR, but in addition the explicit phase fn is
!!     required (as for solar single-scattering). In this case the phase fn
!!     is interpolated only for the upward-scattered beam. The total optical
!!     depths are also computed (as for visible single-scattering) as these
!!     are used for Rayleigh single-scattering calculations where relevant.
!!
!!   MFASIS - solar-only:
!!     the only output is the nadir extinction optical depth.
!!
!!   All layer outputs are computed for the clear (aerosol-only) and (when
!!   relevant) the cloudy(+aerosol) cases. The appropriate values are used
!!   later in each cloud column according to whether the layer is clear or
!!   cloudy in that particular column.
!!
!!   If Rayleigh multiple-scattering is enabled via the solar DOM solver, then
!!   it is treated as if it were an additional aerosol. The optical properties
!!   are calculated at the start of the aerosol section, and are thus
!!   included with aerosol (if present), and then combined later
!!   with cloud properties (if present). For user aerosol optical properties,
!!   the Rayleigh parameters are combined with the input aerosol properties.
!!
!!
!! @param[in]     nlayers             number of layers in input profile
!! @param[in]     chanprof            specifies channels and profiles to simulate
!! @param[in]     opts                options to configure the simulations
!! @param[in]     dom_nstr            number of DOM streams
!! @param[in,out] aux                 additional internal profile variables
!! @param[in,out] ircld               computed cloud column data
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in]     profiles_int        profiles in internal units
!! @param[in]     aer_opt_param       explicit aerosol optical properties per channel (optional)
!! @param[in]     cld_opt_param       explicit cloud optical properties per channel (optional)
!! @param[in]     do_thermal          flag to indicate if any thermal (emissive) simulations are being performed
!! @param[in]     thermal             per-channel flag to indicate if thermal (emissive) simulations are being performed
!! @param[in]     do_solar            flag to indicate if any solar simulations are being performed
!! @param[in]     solar               per-channel flag to indicate if any solar simulations are being performed
!! @param[in]     do_rayleigh_dom     flag to indicate DOM Rayleigh simulation
!! @param[in]     coef                optical depth coefficients structure
!! @param[in]     coef_scatt          visible/IR scattering coefficients structure
!! @param[in]     coef_mfasis_cld     MFASIS cloud coefficients structure
!! @param[in]     angles              geometry structure
!! @param[in]     raytracing          raytracing structure
!! @param[in,out] trans_scatt_ir      computed optical depths
!! @param[in,out] trans_scatt_ir_dyn  computed optical depths per cloud column and phase functions
!! @param[in,out] transmission        output transmission structure
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
SUBROUTINE rttov_opdpscattir( &
              nlayers,            &
              chanprof,           &
              opts,               &
              dom_nstr,           &
              aux,                &
              ircld,              &
              profiles,           &
              profiles_int,       &
              aer_opt_param,      &
              cld_opt_param,      &
              do_thermal,         &
              thermal,            &
              do_solar,           &
              solar,              &
              do_rayleigh_dom,    &
              coef,               &
              coef_scatt,         &
              coef_mfasis_cld,    &
              angles,             &
              raytracing,         &
              trans_scatt_ir,     &
              trans_scatt_ir_dyn, &
              transmission)

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
       rttov_ircld,                 &
       rttov_transmission
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       realtol,            &
       deg2rad,            &
       rad2deg,            &
       max_optical_depth,  &
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
       calc_legendre_coef_gauss
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options),               INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                INTENT(IN)    :: dom_nstr
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux
  TYPE(rttov_ircld),                 INTENT(INOUT) :: ircld
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param),   OPTIONAL, INTENT(IN)    :: cld_opt_param
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
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_dyn
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission
!INTF_END

#include "rttov_baran2014_calc_optpar.interface"
#include "rttov_baran2018_calc_optpar.interface"
#include "rttov_baran_calc_phase.interface"

  LOGICAL(KIND=jplm) :: do_dom, do_chou_scaling, do_single_scatt, do_mfasis, do_aer_or_ray_dom
  LOGICAL(KIND=jplm) :: do_dom_chan, do_chou_scaling_chan, do_single_scatt_chan, do_mfasis_chan, do_chan
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: j, k, k1, k2
  INTEGER(KIND=jpim) :: prof, chan, phchan
  INTEGER(KIND=jpim) :: clw_scheme, ice_scheme, col, coli, iaer, icld
  INTEGER(KIND=jpim) :: lev, lay, levsurf
  REAL(KIND=jprb)    :: opd, opdsun
  REAL(KIND=jprb)    :: dgfrac, abso, sca, sumsca, bpr, asym
  REAL(KIND=jprb)    :: clw, tmpval, ray_phup, cosscata
  REAL(KIND=jprb)    :: afac, sfac, gfac, frach
  REAL(KIND=jprb)    :: pfac1(0:dom_nstr)
  REAL(KIND=jprb)    :: pfac2(coef_scatt%optp_aer%nphangle)

  INTEGER(KIND=jpim) :: nmom
  REAL(KIND=jprb)    :: legcoef(0:dom_nstr), thislegcoef(0:dom_nstr)
  REAL(KIND=jprb)    :: aer_pha(coef_scatt%optp_aer%nphangle)
  REAL(KIND=jprb)    :: wcldeff_pha(coef_scatt%optp_wcl_deff%nphangle)
  REAL(KIND=jprb)    :: icl_pha(coef_scatt%optp_icl_baum%nphangle)
  REAL(KIND=jprb)    :: zminphadiff
  REAL(KIND=jprb)    :: relazi, musat, musun, phasint
  INTEGER(KIND=jpim) :: thisnphangle
  REAL(KIND=jprb)    :: thisphangle(nphangle_hires), thiscosphangle(nphangle_hires)
  REAL(KIND=jprb)    :: baran_pha(1:nphangle_hires)
  REAL(KIND=jprb)    :: baran_pha_interp(baran_ngauss)

  TYPE(rttov_scatt_ir_aercld), POINTER :: aer, cld
  TYPE(rttov_optp_data),       POINTER :: optp_data

  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!-----End of header------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)

  do_dom = (do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
           (do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom)
  do_chou_scaling = do_thermal .AND. opts%rt_ir%ir_scatt_model == ir_scatt_chou
  do_single_scatt = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
  do_mfasis = do_solar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
  do_aer_or_ray_dom = opts%rt_ir%addaerosl .OR. do_rayleigh_dom

  IF (do_thermal) trans_scatt_ir%opdpacl = 0._jprb
  IF (do_solar .AND. .NOT. do_mfasis) THEN
    trans_scatt_ir%opdpaclsun = 0._jprb
    trans_scatt_ir%phup = 0._jprb
  ENDIF
  IF (do_single_scatt) trans_scatt_ir%phdo = 0._jprb
  IF (do_dom) THEN
    DO j = 1, nchanprof
      IF (ASSOCIATED(trans_scatt_ir_dyn%phasefn(j)%legcoef)) trans_scatt_ir_dyn%phasefn(j)%legcoef = 0._jprb
    ENDDO
  ENDIF
  IF (do_mfasis) trans_scatt_ir%opdpext = 0._jprb

  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF AEROSOLS
  !----------------------------------------------------------------------------

  IF (do_aer_or_ray_dom) THEN

    aer => trans_scatt_ir%aer

    IF (do_thermal .OR. .NOT. do_mfasis) THEN
      aer%opdpabs = 0._jprb
      aer%opdpsca = 0._jprb
    ENDIF
    IF (do_chou_scaling) aer%opdpscabpr = 0._jprb
    IF (do_solar .AND. .NOT. do_mfasis) aer%phtotup = 0._jprb
    IF (do_single_scatt) aer%phtotdo = 0._jprb
    IF (do_dom) aer%sca = 0._jprb

    IF (opts%rt_ir%user_aer_opt_param) THEN

      !----------------------------------------------------------------------
      ! User-supplied aerosol optical properties
      !----------------------------------------------------------------------
      CALL calc_user_opt_param(aer_opt_param, aer, 0_jpim, do_rayleigh_dom)

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

        ! Initialisation for Rayleigh phase function:
        !   0.75 * (1 + cos(theta)^2) for scattering angle theta
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
            legcoef(:) = 0._jprb

            ! If DOM Rayleigh enabled calculate the Rayleigh properties
            IF (do_rayleigh_dom .AND. solar(j) .AND. &
                10000._jprb / coef%ff_cwn(chan) <= opts%rt_ir%rayleigh_max_wavelength .AND. &
                profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

              aer%opdpsca(lay,j) = trans_scatt_ir%ray_sca(lay,j)
              aer%phtotup(lay,j) = ray_phup * trans_scatt_ir%ray_sca(lay,j)

              ! Ensure Legendre coefficients are weighted consistently with the aerosol ones below
              ! Legendre coefs are: (1, 0, 0.5, 0, 0, 0, ...)
              sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)
              sumsca = sca
              legcoef(0) = 1._jprb * sca
              legcoef(2) = 0.5_jprb * sca
            ENDIF
          ENDIF

          DO iaer = 1, coef_scatt%optp_aer%ntypes
            IF (profiles_int(prof)%aerosols(iaer,lay) <= 0._jprb) CYCLE

            optp_data => coef_scatt%optp_aer%data(iaer)

            IF (do_dom_chan) thislegcoef(:) = 0._jprb

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
                    pfac2(:) = (optp_data%pha(:,k+1,1,phchan) - optp_data%pha(:,k,1,phchan))
                    aer_pha(:) = optp_data%pha(:,k,1,phchan) + pfac2(:) * frach
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
                aer_pha(:) = optp_data%pha(:,k,1,phchan)
              ENDIF
            ENDIF

            IF (do_mfasis_chan) THEN
              trans_scatt_ir%opdpext(iaer,lay,j) = coef%ff_gam(chan) * &
                profiles_int(prof)%aerosols(iaer,lay) * raytracing%ltick(lay,prof) * (abso + sca)
            ELSEIF (do_chan) THEN
              !------------------------------------------------------------------
              ! Accumulate total aerosol optical parameters (all particles)
              !------------------------------------------------------------------
              ! For pre-defined particle types OD calculated as (aerosol amount in layer) *
              !   (aerosol ext coeff for given layer RH) * (layer thickness)
              ! OD is accumulated for each aerosol type - arrays are zeroed above.
              aer%opdpabs(lay,j) = aer%opdpabs(lay,j) + &
                profiles_int(prof)%aerosols(iaer,lay) * abso * raytracing%ltick(lay,prof)
              aer%opdpsca(lay,j) = aer%opdpsca(lay,j) + &
                profiles_int(prof)%aerosols(iaer,lay) * sca * raytracing%ltick(lay,prof)

              IF (do_chou_scaling_chan) THEN
                aer%opdpscabpr(lay,j) = aer%opdpscabpr(lay,j) + &
                  profiles_int(prof)%aerosols(iaer,lay) * sca * bpr * raytracing%ltick(lay,prof)
              ENDIF

              IF (do_dom_chan) THEN
                sumsca = sumsca + profiles_int(prof)%aerosols(iaer,lay) * sca
                legcoef(:) = legcoef(:) + thislegcoef(:) * profiles_int(prof)%aerosols(iaer,lay) * sca
              ENDIF

              IF (solar(j)) THEN
                musat =  1._jprb / raytracing%pathsat(lay,prof)
                musun = -1._jprb / raytracing%pathsun(lay,prof)
                zminphadiff = coef_scatt%optp_aer%phfn_int%zminphadiff * rad2deg

                CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                  aer_pha, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                  coef_scatt%optp_aer%phfn_int%iphangle, phasint)

                aer%phintup(iaer,lay,j) = phasint
                aer%phtotup(lay,j) = aer%phtotup(lay,j) + &
                                     profiles_int(prof)%aerosols(iaer,lay) * &
                                     phasint * sca * raytracing%ltick(lay,prof)

                IF (do_single_scatt_chan) THEN
                  musat = -musat
                  CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                    aer_pha, coef_scatt%optp_aer%phfn_int%cosphangle, &
                                    coef_scatt%optp_aer%phfn_int%iphangle, phasint)

                  aer%phintdo(iaer,lay,j) = phasint
                  aer%phtotdo(lay,j) = aer%phtotdo(lay,j) + &
                                       profiles_int(prof)%aerosols(iaer,lay) * &
                                       phasint * sca * raytracing%ltick(lay,prof)
                ENDIF
              ENDIF
            ENDIF
          ENDDO ! aer types

          !------------------------------------------------------------------
          ! Calculate final phase function Leg. coefs for all aerosol types
          !------------------------------------------------------------------
          IF (do_dom_chan) THEN
            IF (ABS(sumsca) > realtol) THEN
              aer%sca(lay,j) = sumsca
              trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay) = legcoef(:) / sumsca
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
          aer%opdp(:,j) = aer%opdpabs(:,j) + aer%opdpscabpr(:,j)
        ELSE
          ! Full optical depth for thermal channels
          aer%opdp(:,j) = aer%opdpabs(:,j) + aer%opdpsca(:,j)
        ENDIF
      ENDIF

      IF (solar(j) .AND. .NOT. do_mfasis_chan) THEN
        WHERE (ABS(aer%opdpsca(:,j)) > realtol)
          trans_scatt_ir%phup(0,:,j) = aer%phtotup(:,j) / aer%opdpsca(:,j)
        ENDWHERE
        IF (do_single_scatt_chan) THEN
          WHERE (ABS(aer%opdpsca(:,j)) > realtol)
            trans_scatt_ir%phdo(0,:,j) = aer%phtotdo(:,j) / aer%opdpsca(:,j)
          ENDWHERE
        ENDIF

        ! Full optical depth for solar channels
        aer%opdpsun(:,j) = aer%opdpabs(:,j) + aer%opdpsca(:,j)
      ENDIF
    ENDDO ! chanprof
  ENDIF ! do_aer_or_ray_dom


  !----------------------------------------------------------------------------
  ! CALCULATE OPTICAL DEPTHS OF CLOUDS
  !----------------------------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN

    cld => trans_scatt_ir%cld

    IF (do_thermal .OR. .NOT. do_mfasis) THEN
      cld%opdpabs = 0._jprb
      cld%opdpsca = 0._jprb
    ENDIF
    IF (do_chou_scaling) cld%opdpscabpr = 0._jprb
    IF (do_solar .AND. .NOT. do_mfasis) cld%phtotup = 0._jprb
    IF (do_single_scatt) cld%phtotdo = 0._jprb
    IF (do_dom) cld%sca = 0._jprb

    IF (opts%rt_ir%user_cld_opt_param) THEN

      !----------------------------------------------------------------------
      ! User-supplied cloud optical properties
      !----------------------------------------------------------------------
      CALL calc_user_opt_param(cld_opt_param, cld, 1_jpim, .FALSE._jplm)

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
            legcoef = 0._jprb
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

                IF (do_mfasis_chan) THEN
                  trans_scatt_ir%opdpext(icld,lay,j) = coef%ff_gam(chan) * &
                    clw * raytracing%ltick(lay,prof) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
                ELSEIF (do_chan) THEN
                  cld%opdpabs(lay,j) = cld%opdpabs(lay,j) + &
                    clw * raytracing%ltick(lay,prof) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                    (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                    (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * clw
                  cld%partsca(icld,lay,j) = sca * raytracing%ltick(lay,prof)
                  cld%opdpsca(lay,j) = cld%opdpsca(lay,j) + cld%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
                    cld%partbpr(icld,lay,j) = optp_data%bpr(1,k1,chan) + dgfrac * &
                      (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                    cld%opdpscabpr(lay,j) = cld%opdpscabpr(lay,j) + &
                      cld%partbpr(icld,lay,j) * cld%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                      optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                      (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_wcl_deff%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_wcl_deff%chan_pha_index(chan)
                    wcldeff_pha(:) = &
                      optp_data%pha(:,1,k1,phchan) + dgfrac * &
                      (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))

                    CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                      wcldeff_pha, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                      coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint)

                    cld%phintup(icld,lay,j) = phasint
                    cld%phtotup(lay,j) = cld%phtotup(lay,j) + phasint * cld%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                        wcldeff_pha, coef_scatt%optp_wcl_deff%phfn_int%cosphangle, &
                                        coef_scatt%optp_wcl_deff%phfn_int%iphangle, phasint)

                      cld%phintdo(icld,lay,j) = phasint
                      cld%phtotdo(lay,j) = cld%phtotdo(lay,j) + phasint * cld%partsca(icld,lay,j)
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
                  trans_scatt_ir%opdpext(icld,lay,j) = coef%ff_gam(chan) * &
                    profiles_int(prof)%cloud(icld,lay) * optp_data%confac * &
                    raytracing%ltick(lay,prof) * &
                    (optp_data%abs(1,1,chan) + optp_data%sca(1,1,chan))
                ELSEIF (do_chan) THEN
                  cld%opdpabs(lay,j) = cld%opdpabs(lay,j) + &
                      profiles_int(prof)%cloud(icld,lay) * optp_data%confac * &
                      optp_data%abs(1,1,chan) * raytracing%ltick(lay,prof)
                  cld%partsca(icld,lay,j) = &
                      profiles_int(prof)%cloud(icld,lay) * optp_data%confac * &
                      optp_data%sca(1,1,chan) * raytracing%ltick(lay,prof)
                  cld%opdpsca(lay,j) = cld%opdpsca(lay,j) + cld%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
                    cld%partbpr(icld,lay,j)  = optp_data%bpr(1,1,chan)
                    cld%opdpscabpr(lay,j) = cld%opdpscabpr(lay,j) + &
                        cld%partbpr(icld,lay,j) * cld%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    sca = optp_data%sca(1,1,chan) * profiles_int(prof)%cloud(icld,lay) * &
                          optp_data%confac
                    sumsca = sumsca + sca
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    legcoef(0:nmom) = legcoef(0:nmom) + optp_data%legcoef(1:nmom+1,1,1,chan) * sca
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_wcl_opac%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_wcl_opac%chan_pha_index(chan)
                    CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                      optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                      coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint)

                    cld%phintup(icld,lay,j) = phasint
                    cld%phtotup(lay,j) = cld%phtotup(lay,j) + phasint * cld%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                        optp_data%pha(:,1,1,phchan), coef_scatt%optp_wcl_opac%phfn_int%cosphangle, &
                                        coef_scatt%optp_wcl_opac%phfn_int%iphangle, phasint)

                      cld%phintdo(icld,lay,j) = phasint
                      cld%phtotdo(lay,j) = cld%phtotdo(lay,j) + phasint * cld%partsca(icld,lay,j)
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

                IF (do_mfasis_chan) THEN
                  trans_scatt_ir%opdpext(icld,lay,j) = coef%ff_gam(chan) * &
                    profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                    (optp_data%abs(1,k1,chan) + dgfrac * &
                     (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)) + &
                     optp_data%sca(1,k1,chan) + dgfrac * &
                     (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
                ELSEIF (do_chan) THEN
                  cld%opdpabs(lay,j) = cld%opdpabs(lay,j) + &
                      profiles_int(prof)%cloud(icld,lay) * raytracing%ltick(lay,prof) * &
                      (optp_data%abs(1,k1,chan) + dgfrac * &
                      (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))

                  sca = (optp_data%sca(1,k1,chan) + dgfrac * &
                      (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan))) * profiles_int(prof)%cloud(icld,lay)
                  cld%partsca(icld,lay,j) = sca * raytracing%ltick(lay,prof)
                  cld%opdpsca(lay,j) = cld%opdpsca(lay,j) + cld%partsca(icld,lay,j)

                  IF (do_chou_scaling_chan) THEN
                    cld%partbpr(icld,lay,j) = optp_data%bpr(1,k1,chan) + &
                        dgfrac * (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
                    cld%opdpscabpr(lay,j) = cld%opdpscabpr(lay,j) + &
                        cld%partbpr(icld,lay,j) * cld%partsca(icld,lay,j)
                  ENDIF

                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    nmom = MIN(optp_data%nmom(1,chan), dom_nstr)
                    thislegcoef(0:nmom) = &
                      optp_data%legcoef(1:nmom+1,1,k1,chan) + dgfrac * &
                      (optp_data%legcoef(1:nmom+1,1,k2,chan) - optp_data%legcoef(1:nmom+1,1,k1,chan))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                  ENDIF

                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_icl_baum%phfn_int%zminphadiff * rad2deg

                    phchan = coef_scatt%optp_icl_baum%chan_pha_index(chan)
                    icl_pha(:) = &
                      optp_data%pha(:,1,k1,phchan) + dgfrac * &
                      (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))

                    CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                      icl_pha, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                      coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint)

                    cld%phintup(icld,lay,j) = phasint
                    cld%phtotup(lay,j) = cld%phtotup(lay,j) + phasint * cld%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                        icl_pha, coef_scatt%optp_icl_baum%phfn_int%cosphangle, &
                                        coef_scatt%optp_icl_baum%phfn_int%iphangle, phasint)

                      cld%phintdo(icld,lay,j) = phasint
                      cld%phtotdo(lay,j) = cld%phtotdo(lay,j) + phasint * cld%partsca(icld,lay,j)
                    ENDIF
                  ENDIF
                ENDIF

              ELSE
                !--------------------------------------------------------------
                ! Optical parameters computed using Baran scheme
                !--------------------------------------------------------------
                IF (ice_scheme == ice_scheme_baran2018) THEN
                  CALL rttov_baran2018_calc_optpar(coef_scatt%optp_icl_baran2018, chan, profiles(prof)%t(lay), &
                      profiles_int(prof)%cloud(icld,lay), abso, sca, bpr, asym)
                ELSE
                  CALL rttov_baran2014_calc_optpar(coef_scatt%optp_icl_baran2014, chan, profiles(prof)%t(lay), &
                      profiles_int(prof)%cloud(icld,lay), abso, sca, bpr, asym)
                ENDIF

                cld%opdpabs(lay,j) = cld%opdpabs(lay,j) + abso * raytracing%ltick(lay,prof)
                cld%partsca(icld,lay,j) = sca * raytracing%ltick(lay,prof)
                cld%opdpsca(lay,j) = cld%opdpsca(lay,j) + cld%partsca(icld,lay,j)

                IF (do_chou_scaling_chan) THEN
                  cld%partbpr(icld,lay,j) = bpr
                  cld%opdpscabpr(lay,j) = cld%opdpscabpr(lay,j) + &
                      cld%partbpr(icld,lay,j) * cld%partsca(icld,lay,j)
                ENDIF

                ! Baran phase function
                IF (do_single_scatt_chan .OR. do_dom_chan) THEN

                  ! Set up angular grid for Baran phase fn
                  IF (solar(j)) THEN
                    ! For *all* solar channels use higher resolution angle grid
                    ! (for mixed thermal+solar channels with solar scattering we need iphangle and we have this on
                    !  the hi-res grid; saves calculating it for both hi-res and lo-res grids - this could be changed)
                    thisnphangle = nphangle_hires
                    thisphangle = phangle_hires
                    thiscosphangle = coef_scatt%optp_icl_baran2018%phfn_int%cosphangle
                  ELSE
                    thisnphangle = nphangle_lores
                    thisphangle(1:nphangle_lores) = phangle_lores
                    thiscosphangle(1:nphangle_lores) = COS(phangle_lores * deg2rad)
                  ENDIF

                  ! Compute Baran phase fn
                  CALL rttov_baran_calc_phase(asym, thisphangle(1:thisnphangle), baran_pha(1:thisnphangle))

                  ! Compute Legendre coefficients
                  IF (do_dom_chan) THEN
                    thislegcoef(:) = 0._jprb
                    sumsca = sumsca + sca
                    nmom = dom_nstr

                    CALL spline_interp(thisnphangle, thiscosphangle(thisnphangle:1:-1), &
                                       baran_pha(thisnphangle:1:-1), baran_ngauss, &
                                       coef_scatt%optp_icl_baran2018%q, baran_pha_interp)
                    CALL normalise(baran_ngauss, coef_scatt%optp_icl_baran2018%w, baran_pha_interp)
                    CALL calc_legendre_coef_gauss(coef_scatt%optp_icl_baran2018%q, &
                                                  coef_scatt%optp_icl_baran2018%w, &
                                                  baran_pha_interp, dom_nstr, dom_nstr, nmom, thislegcoef(:))

                    legcoef(:) = legcoef(:) + thislegcoef(:) * sca
                  ENDIF

                  ! Evaluate phase function for solar scattering
                  IF (solar(j)) THEN
                    musat =  1._jprb / raytracing%pathsat(lay,prof)
                    musun = -1._jprb / raytracing%pathsun(lay,prof)
                    zminphadiff = coef_scatt%optp_icl_baran2018%phfn_int%zminphadiff * rad2deg

                    CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                                      baran_pha, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                      coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint)

                    cld%phintup(icld,lay,j) = phasint
                    cld%phtotup(lay,j) = cld%phtotup(lay,j) + phasint * cld%partsca(icld,lay,j)

                    IF (do_single_scatt_chan) THEN
                      musat = -musat
                      CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                                        baran_pha, coef_scatt%optp_icl_baran2018%phfn_int%cosphangle, &
                                        coef_scatt%optp_icl_baran2018%phfn_int%iphangle, phasint)

                      cld%phintdo(icld,lay,j) = phasint
                      cld%phtotdo(lay,j) = cld%phtotdo(lay,j) + phasint * cld%partsca(icld,lay,j)
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
              cld%sca(lay,j) = sumsca
              trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) = legcoef(:) / sumsca
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
          cld%opdp(:,j) = cld%opdpabs(:,j) + cld%opdpscabpr(:,j)
        ELSE
          ! Full optical depth for thermal channels
          cld%opdp(:,j) = cld%opdpabs(:,j) + cld%opdpsca(:,j)
        ENDIF
      ENDIF

      ! Full optical depth for solar channels
      IF (solar(j) .AND. .NOT. do_mfasis) cld%opdpsun(:,j) = cld%opdpabs(:,j) + cld%opdpsca(:,j)
    ENDDO ! chanprof

    !------------------------------------------------------------------------
    ! Output cloudy transmittances for all channels - direct model only
    !------------------------------------------------------------------------
    ! This is done here (rather than rttov_transmit* for example) because we can
    ! do the calculations for thermal and solar calculations in one place here
    ! and if we want to output this for other solvers in the future (e.g. MFASIS)
    ! it would need to be done here anyway.
    IF (do_thermal .OR. .NOT. do_mfasis) THEN
      DO j = 1, nchanprof
        IF (solar(j) .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis) CYCLE
        IF (.NOT. (thermal(j) .OR. solar(j))) CYCLE

        prof = chanprof(j)%prof
        chan = chanprof(j)%chan

        DO lay = 1, nlayers
          transmission%tau_levels_cld(lay+1,j) = transmission%tau_levels_cld(lay,j) + &
            coef%ff_gam(chan) * raytracing%pathsat(lay,prof) * &
            (cld%opdpabs(lay,j) + cld%opdpsca(lay,j))
        ENDDO

        levsurf = aux%s(prof)%nearestlev_surf
        transmission%tau_total_cld(j) = transmission%tau_levels_cld(levsurf,j) + &
          aux%s(prof)%pfraction_surf * &
          (transmission%tau_levels_cld(levsurf-1,j) - transmission%tau_levels_cld(levsurf,j))

        transmission%tau_levels_cld(:,j) = MIN(max_optical_depth, transmission%tau_levels_cld(:,j))
        transmission%tau_levels_cld(:,j) = EXP(-transmission%tau_levels_cld(:,j))
        transmission%tau_total_cld(j)    = MIN(max_optical_depth, transmission%tau_total_cld(j))
        transmission%tau_total_cld(j)    = EXP(-transmission%tau_total_cld(j))
      ENDDO
    ENDIF
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
        trans_scatt_ir%opdpacl(0,:,j) = &
          aer%opdp(:,j) * raytracing%pathsat(:,prof) * coef%ff_gam(chan)
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir%opdpacl(1,:,j) = trans_scatt_ir%opdpacl(0,:,j) + &
          cld%opdp(:,j) * raytracing%pathsat(:,prof) * coef%ff_gam(chan)
      ENDIF
    ENDIF
    IF (solar(j)) THEN
      IF (do_aer_or_ray_dom) THEN
        trans_scatt_ir%opdpaclsun(0,:,j) = &
          aer%opdpsun(:,j) * raytracing%patheff(:,prof) * coef%ff_gam(chan)
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir%opdpaclsun(1,:,j) = trans_scatt_ir%opdpaclsun(0,:,j) + &
          cld%opdpsun(:,j) * raytracing%patheff(:,prof) * coef%ff_gam(chan)
      ENDIF
    ENDIF

    IF (do_dom_chan .OR. do_single_scatt_chan) THEN

      ! Determine final layer *nadir* absorption and scattering optical depths
      IF (do_aer_or_ray_dom) THEN
        trans_scatt_ir%opdpabs(0,:,j) = coef%ff_gam(chan) * aer%opdpabs(:,j)
        trans_scatt_ir%opdpsca(0,:,j) = coef%ff_gam(chan) * aer%opdpsca(:,j)
      ELSE
        trans_scatt_ir%opdpabs(0,:,j) = 0._jprb
        trans_scatt_ir%opdpsca(0,:,j) = 0._jprb
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        trans_scatt_ir%opdpabs(1,:,j) = trans_scatt_ir%opdpabs(0,:,j) + &
          coef%ff_gam(chan) * cld%opdpabs(:,j)
        trans_scatt_ir%opdpsca(1,:,j) = trans_scatt_ir%opdpsca(0,:,j) + &
          coef%ff_gam(chan) * cld%opdpsca(:,j)
      ENDIF

      DO lay = 1, nlayers
        IF (do_dom_chan) THEN
          ! Determine combined phase functions for aer+cld case
          IF (do_aer_or_ray_dom .AND. opts%rt_ir%addclouds) THEN
            IF (aer%sca(lay,j) + cld%sca(lay,j) > 0._jprb) THEN
              trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) = &
                  (trans_scatt_ir_dyn%phasefn(j)%legcoef(:,0,lay) * aer%sca(lay,j) + &
                   trans_scatt_ir_dyn%phasefn(j)%legcoef(:,1,lay) * cld%sca(lay,j)) / &
                  (aer%sca(lay,j) + cld%sca(lay,j))
            ENDIF
          ENDIF
        ENDIF

        IF (solar(j)) THEN
          IF (opts%rt_ir%addclouds .AND. do_aer_or_ray_dom) THEN
            IF (ABS(cld%opdpsca(lay,j) + aer%opdpsca(lay,j)) > realtol) THEN
              tmpval = 1._jprb / (cld%opdpsca(lay,j) + aer%opdpsca(lay,j))
              trans_scatt_ir%phup(1,lay,j) = &
                (aer%phtotup(lay,j) + cld%phtotup(lay,j)) * tmpval
              IF (do_single_scatt_chan) &
                trans_scatt_ir%phdo(1,lay,j) = &
                  (aer%phtotdo(lay,j) + cld%phtotdo(lay,j)) * tmpval
            ENDIF
          ELSEIF (opts%rt_ir%addclouds) THEN
            IF (ABS(cld%opdpsca(lay,j)) > realtol) THEN
              tmpval = 1._jprb / cld%opdpsca(lay,j)
              trans_scatt_ir%phup(1,lay,j) = cld%phtotup(lay,j) * tmpval
              IF (do_single_scatt_chan) &
                trans_scatt_ir%phdo(1,lay,j) = cld%phtotdo(lay,j) * tmpval
            ENDIF
          ENDIF
        ENDIF

      ENDDO ! layers
    ENDIF ! do dom or do single_scatt

    DO col = 0, ircld%ncolumn(prof)
      IF (thermal(j)) THEN
        opd = 0._jprb
        trans_scatt_ir_dyn%opdpac(1,col,j) = 0._jprb
      ENDIF
      IF (solar(j)) THEN
        opdsun = 0._jprb
        trans_scatt_ir_dyn%opdpacsun(1,col,j) = 0._jprb
      ENDIF
      IF (col == 0) THEN
        IF (do_aer_or_ray_dom) THEN
          DO lay = 1, nlayers
            lev = lay + 1
            IF (thermal(j)) THEN
              opd = opd + trans_scatt_ir%opdpacl(0,lay,j)
              trans_scatt_ir_dyn%opdpac(lev,col,j) = opd
            ENDIF
            IF (solar(j)) THEN
              opdsun = opdsun + trans_scatt_ir%opdpaclsun(0,lay,j)
              trans_scatt_ir_dyn%opdpacsun(lev,col,j) = opdsun
            ENDIF
          ENDDO ! layers
        ELSE
          IF (thermal(j)) trans_scatt_ir_dyn%opdpac(:,col,j) = 0._jprb
          IF (solar(j)) trans_scatt_ir_dyn%opdpacsun(:,col,j) = 0._jprb
        ENDIF
      ELSE
        DO lay = 1, nlayers
          lev = lay + 1
          coli = ircld%icldarr(col,lay,prof) ! This is 0 or 1 for clear(+aer)/cloud(+aer)
          IF (thermal(j)) THEN
            opd = opd + trans_scatt_ir%opdpacl(coli,lay,j)
            trans_scatt_ir_dyn%opdpac(lev,col,j) = opd
          ENDIF
          IF (solar(j)) THEN
            opdsun = opdsun + trans_scatt_ir%opdpaclsun(coli,lay,j)
            trans_scatt_ir_dyn%opdpacsun(lev,col,j) = opdsun
          ENDIF
        ENDDO ! layers
      ENDIF ! col == 0
    ENDDO ! col
  ENDDO ! channels

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_user_opt_param(opt_param, aercld, iaercld, do_rayleigh_dom)
    ! Process aer/cld user optical property inputs
    TYPE(rttov_opt_param),       INTENT(IN)    :: opt_param
    TYPE(rttov_scatt_ir_aercld), INTENT(INOUT) :: aercld
    INTEGER(jpim),               INTENT(IN)    :: iaercld ! aer=>0, cld=>1
    LOGICAL(jplm),               INTENT(IN)    :: do_rayleigh_dom

    INTEGER(jpim) :: j, lay, prof, chan
    LOGICAL(jplm) :: do_chou_scaling_chan, do_dom_chan, do_single_scatt_chan
    LOGICAL(jplm) :: do_opt_param_lay, do_rayleigh_dom_chan
    REAL(jprb)    :: musat, musun, zminphadiff, phasint
    REAL(jprb)    :: ray_lcoef(0:dom_nstr), opdpsca, sca

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
        cosscata = - angles(prof)%coszen * angles(prof)%coszen_sun - &
                     angles(prof)%sinzen * angles(prof)%sinzen_sun * COS(relazi * deg2rad)
        ray_phup = 0.75_jprb * (1._jprb + cosscata * cosscata)
      ENDIF

      DO lay = 1, nlayers

        do_opt_param_lay = opt_param%abs(lay,j) > 0._jprb .OR. opt_param%sca(lay,j) > 0._jprb

        IF (do_opt_param_lay) THEN
          aercld%opdpabs(lay,j) = opt_param%abs(lay,j) * raytracing%ltick(lay,prof)
          opdpsca = opt_param%sca(lay,j) * raytracing%ltick(lay,prof)
        ELSE
          opdpsca = 0._jprb
        ENDIF

        IF (do_rayleigh_dom_chan .AND. &
            profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

          aercld%opdpsca(lay,j) = trans_scatt_ir%ray_sca(lay,j) + opdpsca
          aercld%phtotup(lay,j) = ray_phup * trans_scatt_ir%ray_sca(lay,j)

        ELSEIF (do_opt_param_lay) THEN
          aercld%opdpsca(lay,j) = opdpsca
        ENDIF

        IF (do_opt_param_lay .AND. do_chou_scaling_chan) THEN
          aercld%opdpscabpr(lay,j) = opt_param%sca(lay,j) * opt_param%bpr(lay,j) * raytracing%ltick(lay,prof)
        ENDIF

        IF (do_dom_chan) THEN
          IF (do_rayleigh_dom_chan .AND. &
              profiles(prof)%p(lay+1) >= opts%rt_ir%rayleigh_min_pressure) THEN

            sca = trans_scatt_ir%ray_sca(lay,j) / raytracing%ltick(lay,prof)

            IF (do_opt_param_lay) THEN
              trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay) = &
                (ray_lcoef(:) * sca + &
                 opt_param%legcoef(1:dom_nstr+1,lay,j) * opt_param%sca(lay,j)) / &
                (sca + opt_param%sca(lay,j))
              aercld%sca(lay,j) = sca + opt_param%sca(lay,j)
            ELSE
              trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay) = ray_lcoef(:)
              aercld%sca(lay,j) = sca
            ENDIF

          ELSEIF (do_opt_param_lay) THEN
            aercld%sca(lay,j) = opt_param%sca(lay,j)
            trans_scatt_ir_dyn%phasefn(j)%legcoef(:,iaercld,lay) = opt_param%legcoef(1:dom_nstr+1,lay,j)
          ENDIF
        ENDIF

        IF (do_opt_param_lay .AND. solar(j)) THEN
          musat =  1._jprb / raytracing%pathsat(lay,prof)
          musun = -1._jprb / raytracing%pathsun(lay,prof)
          zminphadiff = opt_param%phasefn_int%zminphadiff * rad2deg

          CALL int_phase_fn(musat, musun, 180.0_jprb - relazi, zminphadiff, &
                            opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                            opt_param%phasefn_int%iphangle, phasint)

          aercld%phintup(1,lay,j) = phasint
          aercld%phtotup(lay,j) = aercld%phtotup(lay,j) + phasint * opdpsca

          IF (do_single_scatt_chan) THEN
            musat = -musat
            CALL int_phase_fn(musat, musun, relazi, zminphadiff, &
                              opt_param%pha(:,lay,j), opt_param%phasefn_int%cosphangle, &
                              opt_param%phasefn_int%iphangle, phasint)

            aercld%phintdo(1,lay,j) = phasint
            aercld%phtotdo(lay,j) = aercld%phtotdo(lay,j) + phasint * opdpsca
          ENDIF
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE calc_user_opt_param

  SUBROUTINE int_phase_fn(musat, musun, relazi, zminphadiff, pha, cospha, ipha, phasint)
    ! Interpolate phase function to scattering angle
    REAL(KIND=jprb),    INTENT(IN)  :: musat, musun, relazi
    REAL(KIND=jprb),    INTENT(IN)  :: zminphadiff
    REAL(KIND=jprb),    INTENT(IN)  :: pha(:)
    REAL(KIND=jprb),    INTENT(IN)  :: cospha(:)
    INTEGER(KIND=jpim), INTENT(IN)  :: ipha(:)
    REAL(KIND=jprb),    INTENT(OUT) :: phasint

    INTEGER(KIND=jpim) :: ikk, kk
    REAL(KIND=jprb)    :: ztmpx, scattangle, deltap, delta

    phasint = 0._jprb
    ztmpx   = SQRT((1._jprb - musat ** 2) * (1._jprb - musun ** 2))
    scattangle = musat * musun + ztmpx * COS(relazi * deg2rad)
    ikk        = MAX(1_jpim, INT(ACOS(scattangle) * zminphadiff, jpim))
    kk         = ipha(ikk) - 1_jpim
    deltap     = pha(kk + 1) - pha(kk)
    delta      = cospha(kk) - cospha(kk + 1)
    phasint    = pha(kk) + deltap * (cospha(kk) - scattangle) / delta

  END SUBROUTINE int_phase_fn

END SUBROUTINE rttov_opdpscattir
