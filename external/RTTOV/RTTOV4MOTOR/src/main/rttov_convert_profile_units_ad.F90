! Description:
!> @file
!!   AD of conversion of units of input gas, cloud and aerosol profiles.
!
!> @brief
!!   AD of conversion of units of input gas, cloud and aerosol profiles.
!!
!! @details
!!   AD of conversion of units of input gas, cloud and aerosol profiles.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     coefs            coefficients structure
!! @param[in]     profiles         profiles structure containing input profiles
!! @param[in,out] profiles_ad      profiles structure containing profile increments
!! @param[in]     profiles_int     profiles structure using internal units
!! @param[in,out] profiles_int_ad  profiles structure containing converted profile increments
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
SUBROUTINE rttov_convert_profile_units_ad(opts, coefs, profiles, profiles_ad, &
                                          profiles_int, profiles_int_ad)

  USE rttov_types, ONLY : &
      rttov_options,      &
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

  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(rttov_coefs),   INTENT(IN)    :: coefs
  TYPE(rttov_profile), INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile), INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_profile), INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile), INTENT(INOUT) :: profiles_int_ad(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles, nlevels, lay, nlayers
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: cldaerfactor(profiles(1)%nlayers)
  REAL(jprb)    :: factor_ad(profiles(1)%nlevels), factor2m_ad
  REAL(jprb)    :: cldaerfactor_ad(profiles(1)%nlayers)
  REAL(jprb)    :: p(profiles(1)%nlayers), t(profiles(1)%nlayers), q(profiles(1)%nlayers)
  REAL(jprb)    :: p_ad(profiles(1)%nlayers), t_ad(profiles(1)%nlayers), q_ad(profiles(1)%nlayers)
  REAL(jprb)    :: qlev(profiles(1)%nlevels)
  REAL(jprb)    :: qlev_ad(profiles(1)%nlevels)
  LOGICAL(jplm) :: docld, doaer
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_AD', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)


  ! Cloud/aerosol unit conversion

  docld = opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param
  doaer = opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param

  IF (docld .OR. doaer) THEN
    nlevels = profiles(1)%nlevels
    nlayers = profiles(1)%nlayers

    DO iprof = 1, nprofiles

      IF (profiles(iprof)%mmr_cldaer) THEN

        ! Direct

        IF (profiles(1)%gas_units == gas_unit_specconc) THEN
          ! Input units were kg/kg moist
          q = 0.5_jprb * (profiles(iprof)%q(1:nlayers) + profiles(iprof)%q(2:nlevels))
        ELSE
          ! Convert ppmv dry to kg/kg moist
          qlev = profiles_int(iprof)%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(iprof)%q)
          q = 0.5_jprb * (qlev(1:nlayers) + qlev(2:nlevels))
        ENDIF

        factor(1:nlayers) = (1._jprb + q * (Mair - Mh2o) / Mh2o) * rgc / Mair
        p = 0.5_jprb * 100._jprb * (profiles(iprof)%p(1:nlayers) + profiles(iprof)%p(2:nlevels))
        t = 0.5_jprb * (profiles(iprof)%t(1:nlayers) + profiles(iprof)%t(2:nlevels))
        cldaerfactor(:) = p(:) / (t(:) * factor(1:nlayers))


        ! Adjoint

        cldaerfactor_ad = 0._jprb
        qlev_ad = 0._jprb

        IF (doaer) THEN
          DO lay = nlayers, 1, -1
            cldaerfactor_ad(lay) = cldaerfactor_ad(lay) + SUM(profiles_int_ad(iprof)%aerosols(:,lay) * &
                profiles(iprof)%aerosols(:,lay) * coefs%coef_scatt%optp_aer%data(:)%confac)
            profiles_ad(iprof)%aerosols(:,lay) = profiles_ad(iprof)%aerosols(:,lay) + &
                profiles_int_ad(iprof)%aerosols(:,lay) *  cldaerfactor(lay) * coefs%coef_scatt%optp_aer%data(:)%confac
          ENDDO
        ENDIF

        IF (docld) THEN
          DO lay = nlayers, 1, -1
            cldaerfactor_ad(lay) = cldaerfactor_ad(lay) + SUM(profiles_int_ad(iprof)%cloud(:,lay) * &
                profiles(iprof)%cloud(:,lay))
            profiles_ad(iprof)%cloud(:,lay) = profiles_ad(iprof)%cloud(:,lay) + &
                profiles_int_ad(iprof)%cloud(:,lay) *  cldaerfactor(lay)
          ENDDO
        ENDIF


        IF (opts%interpolation%lgradp) THEN
          p_ad(:) = cldaerfactor_ad(:) / (t(:) * factor(1:nlayers))
          profiles_ad(iprof)%p(1:nlayers) = profiles_ad(iprof)%p(1:nlayers) + 0.5_jprb * 100._jprb * p_ad
          profiles_ad(iprof)%p(2:nlevels) = profiles_ad(iprof)%p(2:nlevels) + 0.5_jprb * 100._jprb * p_ad
        ENDIF

        t_ad(:) = -cldaerfactor_ad(:) * p(:) / (t(:)**2 * factor(1:nlayers))
        factor_ad(1:nlayers) = -cldaerfactor_ad(:) * p(:) / (t(:) * factor(1:nlayers)**2)

        profiles_ad(iprof)%t(1:nlayers) = profiles_ad(iprof)%t(1:nlayers) + 0.5_jprb * t_ad
        profiles_ad(iprof)%t(2:nlevels) = profiles_ad(iprof)%t(2:nlevels) + 0.5_jprb * t_ad

        q_ad = (factor_ad(1:nlayers) * (Mair - Mh2o) / Mh2o) * rgc / Mair

        IF (profiles(1)%gas_units == gas_unit_specconc) THEN
          ! Input units were kg/kg moist
          profiles_ad(iprof)%q(1:nlayers) = profiles_ad(iprof)%q(1:nlayers) + 0.5_jprb * q_ad
          profiles_ad(iprof)%q(2:nlevels) = profiles_ad(iprof)%q(2:nlevels) + 0.5_jprb * q_ad
        ELSE
          ! Convert ppmv dry to kg/kg moist
          qlev_ad(1:nlayers) = qlev_ad(1:nlayers) + 0.5_jprb * q_ad
          qlev_ad(2:nlevels) = qlev_ad(2:nlevels) + 0.5_jprb * q_ad

          profiles_int_ad(iprof)%q = profiles_int_ad(iprof)%q + qlev_ad * (1._jprb - qlev) / &
                    (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(iprof)%q)
        ENDIF

      ELSE
        IF (docld) profiles_ad(iprof)%cloud    = profiles_ad(iprof)%cloud + profiles_int_ad(iprof)%cloud
        IF (doaer) profiles_ad(iprof)%aerosols = profiles_ad(iprof)%aerosols + profiles_int_ad(iprof)%aerosols
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

    DO iprof = 1, nprofiles

      factor = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%q)
      IF (opts%rt_all%use_q2m) factor2m = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%s2m%q)

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + &
                             factor * profiles_int_ad(iprof)%q / gas_mass(gas_id_watervapour)
      ! NB First use of factor_ad
      factor_ad = profiles(iprof)%q * profiles_int_ad(iprof)%q / gas_mass(gas_id_watervapour)

      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + &
                                   factor2m * profiles_int_ad(iprof)%s2m%q / gas_mass(gas_id_watervapour)
        ! NB First use of factor2m_ad
        factor2m_ad = profiles(iprof)%s2m%q * profiles_int_ad(iprof)%s2m%q / gas_mass(gas_id_watervapour)
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + &
                                factor * profiles_int_ad(iprof)%o3 / gas_mass(gas_id_ozone)
        factor_ad = factor_ad + profiles(iprof)%o3 * profiles_int_ad(iprof)%o3 / gas_mass(gas_id_ozone)

        IF (opts%rt_all%use_q2m) THEN
          profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + &
                                     factor2m * profiles_int_ad(iprof)%s2m%o / gas_mass(gas_id_ozone)
          factor2m_ad = factor2m_ad + profiles(iprof)%s2m%o * profiles_int_ad(iprof)%s2m%o / gas_mass(gas_id_ozone)
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) THEN
        profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + &
                                 factor * profiles_int_ad(iprof)%co2 / gas_mass(gas_id_co2)
        factor_ad = factor_ad + profiles(iprof)%co2 * profiles_int_ad(iprof)%co2 / gas_mass(gas_id_co2)
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  THEN
        profiles_ad(iprof)%co = profiles_ad(iprof)%co + &
                                factor * profiles_int_ad(iprof)%co / gas_mass(gas_id_co)
        factor_ad = factor_ad + profiles(iprof)%co * profiles_int_ad(iprof)%co / gas_mass(gas_id_co)
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) THEN
        profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + &
                                 factor * profiles_int_ad(iprof)%ch4 / gas_mass(gas_id_ch4)
        factor_ad = factor_ad + profiles(iprof)%ch4 * profiles_int_ad(iprof)%ch4 / gas_mass(gas_id_ch4)
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%so2) .AND. ASSOCIATED(profiles_ad(iprof)%so2)) THEN
        profiles_ad(iprof)%so2 = profiles_ad(iprof)%so2 + &
                                 factor * profiles_int_ad(iprof)%so2 / gas_mass(gas_id_so2)
        factor_ad = factor_ad + profiles(iprof)%so2 * profiles_int_ad(iprof)%so2 / gas_mass(gas_id_so2)
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) THEN
        profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + &
                                 factor * profiles_int_ad(iprof)%n2o / gas_mass(gas_id_n2o)
        factor_ad = factor_ad + profiles(iprof)%n2o * profiles_int_ad(iprof)%n2o / gas_mass(gas_id_n2o)
      ENDIF

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + 1.E-06_jprb * factor_ad * factor**2 / Mair
      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + 1.E-06_jprb * factor2m_ad * factor2m**2 / Mair
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

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      IF (opts%rt_all%use_q2m) THEN
        factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)
      ENDIF

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + factor * profiles_int_ad(iprof)%q
      ! NB First use of factor_ad
      factor_ad = profiles(iprof)%q * profiles_int_ad(iprof)%q

      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + factor2m * profiles_int_ad(iprof)%s2m%q
        ! NB First use of factor2m_ad
        factor2m_ad = profiles(iprof)%s2m%q * profiles_int_ad(iprof)%s2m%q
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + factor * profiles_int_ad(iprof)%o3
        factor_ad = factor_ad + profiles(iprof)%o3 * profiles_int_ad(iprof)%o3

        IF (opts%rt_all%use_q2m) THEN
          profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + factor2m * profiles_int_ad(iprof)%s2m%o
          factor2m_ad = factor2m_ad + profiles(iprof)%s2m%o * profiles_int_ad(iprof)%s2m%o
        ENDIF
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) THEN
        profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + factor * profiles_int_ad(iprof)%co2
        factor_ad = factor_ad + profiles(iprof)%co2 * profiles_int_ad(iprof)%co2
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  THEN
        profiles_ad(iprof)%co = profiles_ad(iprof)%co + factor * profiles_int_ad(iprof)%co
        factor_ad = factor_ad + profiles(iprof)%co * profiles_int_ad(iprof)%co
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) THEN
        profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + factor * profiles_int_ad(iprof)%ch4
        factor_ad = factor_ad + profiles(iprof)%ch4 * profiles_int_ad(iprof)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%so2) .AND. ASSOCIATED(profiles_ad(iprof)%so2)) THEN
        profiles_ad(iprof)%so2 = profiles_ad(iprof)%so2 + factor * profiles_int_ad(iprof)%so2
        factor_ad = factor_ad + profiles(iprof)%so2 * profiles_int_ad(iprof)%so2
      ENDIF

      IF (ASSOCIATED(profiles_int_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) THEN
        profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + factor * profiles_int_ad(iprof)%n2o
        factor_ad = factor_ad + profiles(iprof)%n2o * profiles_int_ad(iprof)%n2o
      ENDIF

      profiles_ad(iprof)%q = profiles_ad(iprof)%q + 1.E-06_jprb * factor_ad * factor**2
      IF (opts%rt_all%use_q2m) THEN
        profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + 1.E-06_jprb * factor2m_ad * factor2m**2
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. ppmv dry)

    DO iprof = 1, nprofiles
      profiles_ad(iprof)%q = profiles_ad(iprof)%q + profiles_int_ad(iprof)%q
      profiles_ad(iprof)%s2m%q = profiles_ad(iprof)%s2m%q + profiles_int_ad(iprof)%s2m%q
      IF (ASSOCIATED(profiles_int_ad(iprof)%o3) .AND. ASSOCIATED(profiles_ad(iprof)%o3)) THEN
        profiles_ad(iprof)%o3 = profiles_ad(iprof)%o3 + profiles_int_ad(iprof)%o3
        profiles_ad(iprof)%s2m%o = profiles_ad(iprof)%s2m%o + profiles_int_ad(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_int_ad(iprof)%co2) .AND. ASSOCIATED(profiles_ad(iprof)%co2)) &
          profiles_ad(iprof)%co2 = profiles_ad(iprof)%co2 + profiles_int_ad(iprof)%co2
      IF (ASSOCIATED(profiles_int_ad(iprof)%co)  .AND. ASSOCIATED(profiles_ad(iprof)%co))  &
          profiles_ad(iprof)%co  = profiles_ad(iprof)%co  + profiles_int_ad(iprof)%co
      IF (ASSOCIATED(profiles_int_ad(iprof)%ch4) .AND. ASSOCIATED(profiles_ad(iprof)%ch4)) &
          profiles_ad(iprof)%ch4 = profiles_ad(iprof)%ch4 + profiles_int_ad(iprof)%ch4
      IF (ASSOCIATED(profiles_int_ad(iprof)%so2) .AND. ASSOCIATED(profiles_ad(iprof)%so2)) &
          profiles_ad(iprof)%so2 = profiles_ad(iprof)%so2 + profiles_int_ad(iprof)%so2
      IF (ASSOCIATED(profiles_int_ad(iprof)%n2o) .AND. ASSOCIATED(profiles_ad(iprof)%n2o)) &
          profiles_ad(iprof)%n2o = profiles_ad(iprof)%n2o + profiles_int_ad(iprof)%n2o
    ENDDO

  END SELECT

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_profile_units_ad
