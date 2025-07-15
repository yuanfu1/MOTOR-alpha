! Description:
!> @file
!!   K of conversion of units of input gas, cloud and aerosol profiles.
!
!> @brief
!!   K of conversion of units of input gas, cloud and aerosol profiles.
!!
!! @details
!!   This does not correspond precisely to the AD for
!!   reasons of efficiency: multiplicative factors which
!!   only vary per-profile are recalculated only when the
!!   profile number changes in the loop over profiles_k.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     chanprof         RTTOV chanprof structure
!! @param[in]     coefs            coefficients structure
!! @param[in]     profiles         profiles structure containing input profiles
!! @param[in,out] profiles_k       profiles structure containing profile increments
!! @param[in]     profiles_int     profiles structure using internal units
!! @param[in,out] profiles_int_k   profiles structure containing converted profile increments
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
SUBROUTINE rttov_convert_profile_units_k(opts, chanprof, coefs, profiles, profiles_k, &
                                         profiles_int, profiles_int_k)

  USE rttov_types, ONLY : &
      rttov_options,      &
      rttov_chanprof,     &
      rttov_coefs,        &
      rttov_profile
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_const, ONLY : &
      gas_unit_specconc,  &
      gas_unit_ppmv,      &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_n2o,         &
      gas_id_co,          &
      gas_id_ch4,         &
      gas_id_so2,         &
      mair,               &
      mh2o,               &
      rgc,                &
      gas_mass

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),  INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coefs),    INTENT(IN)    :: coefs
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),  INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile),  INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),  INTENT(INOUT) :: profiles_int_k(SIZE(chanprof))
!INTF_END

  INTEGER(jpim) :: ichan, nchanprof, prof, lastprof, nlevels, lay, nlayers
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: cldaerfactor(profiles(1)%nlayers)
  REAL(jprb)    :: factor_k(profiles(1)%nlevels), factor2m_k
  REAL(jprb)    :: cldaerfactor_k(profiles(1)%nlayers)
  REAL(jprb)    :: p(profiles(1)%nlayers), t(profiles(1)%nlayers), q(profiles(1)%nlayers)
  REAL(jprb)    :: p_k(profiles(1)%nlayers), t_k(profiles(1)%nlayers), q_k(profiles(1)%nlayers)
  REAL(jprb)    :: qlev(profiles(1)%nlevels)
  REAL(jprb)    :: qlev_k(profiles(1)%nlevels)
  LOGICAL(jplm) :: docld, doaer
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_K', 0_jpim, ZHOOK_HANDLE)

! As noted above multiplicative factors which only vary per-profile are
! recalculated only when the profile number changes in the loop over
! profiles_k. This is much more efficient than re-calculating them in
! every loop, but it means this code differs superficially to the AD.

  nchanprof = SIZE(profiles_k)


  ! Cloud/aerosol unit conversion

  docld = opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param
  doaer = opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param

  IF (docld .OR. doaer) THEN
    nlevels = profiles(1)%nlevels
    nlayers = profiles(1)%nlayers

    lastprof = -1
    DO ichan = 1, nchanprof
      prof = chanprof(ichan)%prof

      IF (profiles(prof)%mmr_cldaer) THEN

        ! Direct

        IF (lastprof /= prof) THEN
          IF (profiles(1)%gas_units == gas_unit_specconc) THEN
            ! Input units were kg/kg moist
            q = 0.5_jprb * (profiles(prof)%q(1:nlayers) + profiles(prof)%q(2:nlevels))
          ELSE
            ! Convert ppmv dry to kg/kg moist
            qlev = profiles_int(prof)%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(prof)%q)
            q = 0.5_jprb * (qlev(1:nlayers) + qlev(2:nlevels))
          ENDIF

          factor(1:nlayers) = (1._jprb + q * (Mair - Mh2o) / Mh2o) * rgc / Mair
          p = 0.5_jprb * 100._jprb * (profiles(prof)%p(1:nlayers) + profiles(prof)%p(2:nlevels))
          t = 0.5_jprb * (profiles(prof)%t(1:nlayers) + profiles(prof)%t(2:nlevels))
          cldaerfactor(:) = p(:) / (t(:) * factor(1:nlayers))

          lastprof = prof
        ENDIF

        ! Adjoint

        cldaerfactor_k = 0._jprb

        IF (doaer) THEN
          DO lay = nlayers, 1, -1
            cldaerfactor_k(lay) = cldaerfactor_k(lay) + SUM(profiles_int_k(ichan)%aerosols(:,lay) * &
                profiles(prof)%aerosols(:,lay) * coefs%coef_scatt%optp_aer%data(:)%confac)
            profiles_k(ichan)%aerosols(:,lay) = profiles_k(ichan)%aerosols(:,lay) + &
                profiles_int_k(ichan)%aerosols(:,lay) *  cldaerfactor(lay) * coefs%coef_scatt%optp_aer%data(:)%confac
          ENDDO
        ENDIF

        IF (docld) THEN
          DO lay = nlayers, 1, -1
            cldaerfactor_k(lay) = cldaerfactor_k(lay) + SUM(profiles_int_k(ichan)%cloud(:,lay) * &
                profiles(prof)%cloud(:,lay))
            profiles_k(ichan)%cloud(:,lay) = profiles_k(ichan)%cloud(:,lay) + &
                profiles_int_k(ichan)%cloud(:,lay) *  cldaerfactor(lay)
          ENDDO
        ENDIF


        IF (opts%interpolation%lgradp) THEN
          p_k(:) = cldaerfactor_k(:) / (t(:) * factor(1:nlayers))
          profiles_k(ichan)%p(1:nlayers) = profiles_k(ichan)%p(1:nlayers) + 0.5_jprb * 100._jprb * p_k
          profiles_k(ichan)%p(2:nlevels) = profiles_k(ichan)%p(2:nlevels) + 0.5_jprb * 100._jprb * p_k
        ENDIF

        t_k(:) = -cldaerfactor_k(:) * p(:) / (t(:)**2 * factor(1:nlayers))
        factor_k(1:nlayers) = -cldaerfactor_k(:) * p(:) / (t(:) * factor(1:nlayers)**2)

        profiles_k(ichan)%t(1:nlayers) = profiles_k(ichan)%t(1:nlayers) + 0.5_jprb * t_k
        profiles_k(ichan)%t(2:nlevels) = profiles_k(ichan)%t(2:nlevels) + 0.5_jprb * t_k

        q_k = (factor_k(1:nlayers) * (Mair - Mh2o) / Mh2o) * rgc / Mair

        IF (profiles(1)%gas_units == gas_unit_specconc) THEN
          ! Input units were kg/kg moist
          profiles_k(ichan)%q(1:nlayers) = profiles_k(ichan)%q(1:nlayers) + 0.5_jprb * q_k
          profiles_k(ichan)%q(2:nlevels) = profiles_k(ichan)%q(2:nlevels) + 0.5_jprb * q_k
        ELSE
          ! Convert ppmv dry to kg/kg moist
          qlev_k = 0._jprb
          qlev_k(1:nlayers) = 0.5_jprb * q_k
          qlev_k(2:nlevels) = qlev_k(2:nlevels) + 0.5_jprb * q_k

          profiles_int_k(ichan)%q = profiles_int_k(ichan)%q + qlev_k * (1._jprb - qlev) / &
                    (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(prof)%q)
        ENDIF

      ELSE
        IF (docld) profiles_k(ichan)%cloud    = profiles_k(ichan)%cloud + profiles_int_k(ichan)%cloud
        IF (doaer) profiles_k(ichan)%aerosols = profiles_k(ichan)%aerosols + profiles_int_k(ichan)%aerosols
      ENDIF
    ENDDO
  ENDIF


  ! Gas units conversion

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!     ppmvdry_gas = kgkgwet_gas / (1._jprb - kgkgwet_h2o) * 1.E+06_jprb * Mair / gas_mass(gas_id)
!     factor = 1.E+06_jprb * Mair / (1._jprb - kgkgwet_h2o)
!     factor_tl = 1.E-06_jprb * kgkgwet_h2o_tl * factor**2 / Mair
!     ppmvdry_gas_tl = (kgkgwet_gas_tl * factor + kgkgwet_gas * factor_tl) / gas_mass(gas_id)

!     kgkgwet_gas_ad = kgkgwet_gas_ad + factor * ppmvdry_gas_ad / gas_mass(gas_id)
!     factor_ad = factor_ad + kgkgwet_gas * ppmvdry_gas_ad / gas_mass(gas_id)
!     kgkgwet_h2o_ad = kgkgwet_h2o_ad + 1.E-06_jprb * factor_ad * factor**2 / Mair

    lastprof = -1
    DO ichan = 1, nchanprof
      prof = chanprof(ichan)%prof

      IF (lastprof /= prof) THEN
        factor = 1.E+06_jprb * Mair / (1._jprb - profiles(prof)%q)
        IF (opts%rt_all%use_q2m) THEN
          factor2m = 1.E+06_jprb * Mair / (1._jprb - profiles(prof)%s2m%q)
        ENDIF
        lastprof = prof
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + &
                            factor * profiles_int_k(ichan)%q / gas_mass(gas_id_watervapour)
      ! NB First use of factor_k
      factor_k = profiles(prof)%q * profiles_int_k(ichan)%q / gas_mass(gas_id_watervapour)

      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + &
                                  factor2m * profiles_int_k(ichan)%s2m%q / gas_mass(gas_id_watervapour)
        ! NB First use of factor2m_k
        factor2m_k = profiles(prof)%s2m%q * profiles_int_k(ichan)%s2m%q / gas_mass(gas_id_watervapour)
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%o3) .AND. ASSOCIATED(profiles_k(ichan)%o3)) THEN
        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + &
                               factor * profiles_int_k(ichan)%o3 / gas_mass(gas_id_ozone)
        factor_k = factor_k + profiles(prof)%o3 * profiles_int_k(ichan)%o3 / gas_mass(gas_id_ozone)

        IF (opts%rt_all%use_q2m) THEN
          profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + &
                                    factor2m * profiles_int_k(ichan)%s2m%o / gas_mass(gas_id_ozone)
          factor2m_k = factor2m_k + profiles(prof)%s2m%o * profiles_int_k(ichan)%s2m%o / gas_mass(gas_id_ozone)
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%co2) .AND. ASSOCIATED(profiles_k(ichan)%co2)) THEN
        profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + &
                                factor * profiles_int_k(ichan)%co2 / gas_mass(gas_id_co2)
        factor_k = factor_k + profiles(prof)%co2 * profiles_int_k(ichan)%co2 / gas_mass(gas_id_co2)
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%co)  .AND. ASSOCIATED(profiles_k(ichan)%co))  THEN
        profiles_k(ichan)%co = profiles_k(ichan)%co + &
                               factor * profiles_int_k(ichan)%co / gas_mass(gas_id_co)
        factor_k = factor_k + profiles(prof)%co * profiles_int_k(ichan)%co / gas_mass(gas_id_co)
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%ch4) .AND. ASSOCIATED(profiles_k(ichan)%ch4)) THEN
        profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + &
                                factor * profiles_int_k(ichan)%ch4 / gas_mass(gas_id_ch4)
        factor_k = factor_k + profiles(prof)%ch4 * profiles_int_k(ichan)%ch4 / gas_mass(gas_id_ch4)
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%so2) .AND. ASSOCIATED(profiles_k(ichan)%so2)) THEN
        profiles_k(ichan)%so2 = profiles_k(ichan)%so2 + &
                                factor * profiles_int_k(ichan)%so2 / gas_mass(gas_id_so2)
        factor_k = factor_k + profiles(prof)%so2 * profiles_int_k(ichan)%so2 / gas_mass(gas_id_so2)
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%n2o) .AND. ASSOCIATED(profiles_k(ichan)%n2o)) THEN
        profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + &
                                factor * profiles_int_k(ichan)%n2o / gas_mass(gas_id_n2o)
        factor_k = factor_k + profiles(prof)%n2o * profiles_int_k(ichan)%n2o / gas_mass(gas_id_n2o)
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + 1.E-06_jprb * factor_k * factor**2 / Mair
      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + 1.E-06_jprb * factor2m_k * factor2m**2 / Mair
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!     ppmvdry_gas = ppmvwet_gas / (1._jprb - ppmvwet_h2o * 1.E-06_jprb)
!     factor = 1. / (1._jprb - ppmvwet_h2o * 1.E-06_jprb)
!     factor_tl = 1.E-06_jprb * ppmvwet_h2o_tl * factor**2
!     ppmvdry_gas_tl = ppmvwet_gas_tl * factor + ppmvwet_gas * factor_tl

!     ppmvwet_gas_ad = ppmvwet_gas_ad + factor * ppmvdry_gas_ad
!     factor_ad = factor_ad + ppmvwet_gas * ppmvdry_gas_ad
!     ppmvwet_h2o_ad = ppmvwet_h2o_ad + 1.E-06_jprb * factor_ad * factor**2

    lastprof = -1
    DO ichan = 1, nchanprof
      prof = chanprof(ichan)%prof

      IF (lastprof /= prof) THEN
        factor = 1._jprb / (1._jprb - profiles(prof)%q * 1.E-06_jprb)
        IF (opts%rt_all%use_q2m) THEN
          factor2m = 1._jprb / (1._jprb - profiles(prof)%s2m%q * 1.E-06_jprb)
        ENDIF
        lastprof = prof
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + factor * profiles_int_k(ichan)%q
      ! NB First use of factor_k
      factor_k = profiles(prof)%q * profiles_int_k(ichan)%q

      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + factor2m * profiles_int_k(ichan)%s2m%q
        ! NB First use of factor2m_k
        factor2m_k = profiles(prof)%s2m%q * profiles_int_k(ichan)%s2m%q
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%o3) .AND. ASSOCIATED(profiles_k(ichan)%o3)) THEN
        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + factor * profiles_int_k(ichan)%o3
        factor_k = factor_k + profiles(prof)%o3 * profiles_int_k(ichan)%o3

        IF (opts%rt_all%use_q2m) THEN
          profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + factor2m * profiles_int_k(ichan)%s2m%o
          factor2m_k = factor2m_k + profiles(prof)%s2m%o * profiles_int_k(ichan)%s2m%o
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%co2) .AND. ASSOCIATED(profiles_k(ichan)%co2)) THEN
        profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + factor * profiles_int_k(ichan)%co2
        factor_k = factor_k + profiles(prof)%co2 * profiles_int_k(ichan)%co2
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%co)  .AND. ASSOCIATED(profiles_k(ichan)%co))  THEN
        profiles_k(ichan)%co = profiles_k(ichan)%co + factor * profiles_int_k(ichan)%co
        factor_k = factor_k + profiles(prof)%co * profiles_int_k(ichan)%co
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%ch4) .AND. ASSOCIATED(profiles_k(ichan)%ch4)) THEN
        profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + factor * profiles_int_k(ichan)%ch4
        factor_k = factor_k + profiles(prof)%ch4 * profiles_int_k(ichan)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%so2) .AND. ASSOCIATED(profiles_k(ichan)%so2)) THEN
        profiles_k(ichan)%so2 = profiles_k(ichan)%so2 + factor * profiles_int_k(ichan)%so2
        factor_k = factor_k + profiles(prof)%so2 * profiles_int_k(ichan)%so2
      ENDIF

      IF (ASSOCIATED(profiles_int_k(ichan)%n2o) .AND. ASSOCIATED(profiles_k(ichan)%n2o)) THEN
        profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + factor * profiles_int_k(ichan)%n2o
        factor_k = factor_k + profiles(prof)%n2o * profiles_int_k(ichan)%n2o
      ENDIF

      profiles_k(ichan)%q = profiles_k(ichan)%q + 1.E-06_jprb * factor_k * factor**2
      IF (opts%rt_all%use_q2m) THEN
        profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + 1.E-06_jprb * factor2m_k * factor2m**2
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. ppmv dry)

    DO ichan = 1, nchanprof
      profiles_k(ichan)%q = profiles_k(ichan)%q + profiles_int_k(ichan)%q
      profiles_k(ichan)%s2m%q = profiles_k(ichan)%s2m%q + profiles_int_k(ichan)%s2m%q
      IF (ASSOCIATED(profiles_int_k(ichan)%o3) .AND. ASSOCIATED(profiles_k(ichan)%o3)) THEN
        profiles_k(ichan)%o3 = profiles_k(ichan)%o3 + profiles_int_k(ichan)%o3
        profiles_k(ichan)%s2m%o = profiles_k(ichan)%s2m%o + profiles_int_k(ichan)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_int_k(ichan)%co2) .AND. ASSOCIATED(profiles_k(ichan)%co2)) &
          profiles_k(ichan)%co2 = profiles_k(ichan)%co2 + profiles_int_k(ichan)%co2
      IF (ASSOCIATED(profiles_int_k(ichan)%co)  .AND. ASSOCIATED(profiles_k(ichan)%co))  &
          profiles_k(ichan)%co  = profiles_k(ichan)%co  + profiles_int_k(ichan)%co
      IF (ASSOCIATED(profiles_int_k(ichan)%ch4) .AND. ASSOCIATED(profiles_k(ichan)%ch4)) &
          profiles_k(ichan)%ch4 = profiles_k(ichan)%ch4 + profiles_int_k(ichan)%ch4
      IF (ASSOCIATED(profiles_int_k(ichan)%so2) .AND. ASSOCIATED(profiles_k(ichan)%so2)) &
          profiles_k(ichan)%so2 = profiles_k(ichan)%so2 + profiles_int_k(ichan)%so2
      IF (ASSOCIATED(profiles_int_k(ichan)%n2o) .AND. ASSOCIATED(profiles_k(ichan)%n2o)) &
          profiles_k(ichan)%n2o = profiles_k(ichan)%n2o + profiles_int_k(ichan)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_profile_units_k
