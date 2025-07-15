! Description:
!> @file
!!   TL of conversion of units of input gas, cloud and aerosol profiles.
!
!> @brief
!!   TL of conversion of units of input gas, cloud and aerosol profiles.
!!
!! @details
!!   TL of conversion of units of input gas, cloud and aerosol profiles.
!!
!! @param[in]     opts             options to configure the simulations
!! @param[in]     coefs            coefficients structure
!! @param[in]     profiles         profiles structure containing input profiles
!! @param[in]     profiles_tl      profiles structure containing input perturbations
!! @param[in]     profiles_int     profiles structure using internal units
!! @param[in,out] profiles_int_tl  profiles structure containing converted profile perturbations on exit
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
SUBROUTINE rttov_convert_profile_units_tl(opts, coefs, profiles, profiles_tl, &
                                          profiles_int, profiles_int_tl)

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
  TYPE(rttov_profile), INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(rttov_profile), INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile), INTENT(INOUT) :: profiles_int_tl(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles, nlevels, lay, nlayers
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: cldaerfactor(profiles(1)%nlayers)
  REAL(jprb)    :: factor_tl(profiles(1)%nlevels), factor2m_tl
  REAL(jprb)    :: cldaerfactor_tl(profiles(1)%nlayers)
  REAL(jprb)    :: p(profiles(1)%nlayers), t(profiles(1)%nlayers), q(profiles(1)%nlayers)
  REAL(jprb)    :: p_tl(profiles(1)%nlayers), t_tl(profiles(1)%nlayers), q_tl(profiles(1)%nlayers)
  REAL(jprb)    :: qlev(profiles(1)%nlevels)
  REAL(jprb)    :: qlev_tl(profiles(1)%nlevels)
  LOGICAL(jplm) :: docld, doaer
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_TL', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)


  ! Gas units conversion

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!     ppmvdry_gas = kgkgwet_gas / (1._jprb - kgkgwet_h2o) * 1.E+06_jprb * Mair / gas_mass(gas_id)
!     factor = 1.E+06_jprb * Mair / (1._jprb - kgkgwet_h2o)
!     factor_tl = 1.E-06_jprb * kgkgwet_h2o_tl * factor**2 / Mair
!     ppmvdry_gas_tl = (kgkgwet_gas_tl * factor + kgkgwet_gas * factor_tl) / gas_mass(gas_id)

    DO iprof = 1, nprofiles

      factor = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%q)
      factor_tl = 1.E-06_jprb * profiles_tl(iprof)%q * factor**2 / Mair
      IF (opts%rt_all%use_q2m) THEN
        factor2m = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%s2m%q)
        factor2m_tl = 1.E-06_jprb * profiles_tl(iprof)%s2m%q * factor2m**2 / Mair
      ENDIF

      profiles_int_tl(iprof)%q = &
          (profiles_tl(iprof)%q * factor + profiles(iprof)%q * factor_tl) / gas_mass(gas_id_watervapour)
      IF (opts%rt_all%use_q2m) profiles_int_tl(iprof)%s2m%q = &
          (profiles_tl(iprof)%s2m%q * factor2m + profiles(iprof)%s2m%q * factor2m_tl) / gas_mass(gas_id_watervapour)

      IF (ASSOCIATED(profiles_int_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_int_tl(iprof)%o3 = &
            (profiles_tl(iprof)%o3 * factor + profiles(iprof)%o3 * factor_tl) / gas_mass(gas_id_ozone)
        IF (opts%rt_all%use_q2m) profiles_int_tl(iprof)%s2m%o = &
            (profiles_tl(iprof)%s2m%o * factor2m + profiles(iprof)%s2m%o * factor2m_tl) / gas_mass(gas_id_ozone)
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) THEN
        profiles_int_tl(iprof)%co2 = &
            (profiles_tl(iprof)%co2 * factor + profiles(iprof)%co2 * factor_tl) / gas_mass(gas_id_co2)
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  THEN
        profiles_int_tl(iprof)%co = &
            (profiles_tl(iprof)%co * factor + profiles(iprof)%co * factor_tl) / gas_mass(gas_id_co)
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) THEN
        profiles_int_tl(iprof)%ch4 = &
            (profiles_tl(iprof)%ch4 * factor + profiles(iprof)%ch4 * factor_tl) / gas_mass(gas_id_ch4)
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%so2) .AND. ASSOCIATED(profiles_tl(iprof)%so2)) THEN
        profiles_int_tl(iprof)%so2 = &
            (profiles_tl(iprof)%so2 * factor + profiles(iprof)%so2 * factor_tl) / gas_mass(gas_id_so2)
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) THEN
        profiles_int_tl(iprof)%n2o = &
            (profiles_tl(iprof)%n2o * factor + profiles(iprof)%n2o * factor_tl) / gas_mass(gas_id_n2o)
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!     ppmvdry_gas = ppmvwet_gas / (1._jprb - ppmvwet_h2o * 1.E-06_jprb)
!     factor = 1. / (1._jprb - ppmvwet_h2o * 1.E-06_jprb)
!     factor_tl = 1.E-06_jprb * ppmvwet_h2o_tl * factor**2
!     ppmvdry_gas_tl = ppmvwet_gas_tl * factor + ppmvwet_gas * factor_tl

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      factor_tl = 1.E-06_jprb * profiles_tl(iprof)%q * factor**2
      IF (opts%rt_all%use_q2m) THEN
        factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)
        factor2m_tl = 1.E-06_jprb * profiles_tl(iprof)%s2m%q * factor2m**2
      ENDIF

      profiles_int_tl(iprof)%q = &
          profiles_tl(iprof)%q * factor + profiles(iprof)%q * factor_tl
      IF (opts%rt_all%use_q2m) profiles_int_tl(iprof)%s2m%q = &
          profiles_tl(iprof)%s2m%q * factor2m + profiles(iprof)%s2m%q * factor2m_tl

      IF (ASSOCIATED(profiles_int_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_int_tl(iprof)%o3 = &
            profiles_tl(iprof)%o3 * factor + profiles(iprof)%o3 * factor_tl
        IF (opts%rt_all%use_q2m) profiles_int_tl(iprof)%s2m%o = &
            profiles_tl(iprof)%s2m%o * factor2m + profiles(iprof)%s2m%o * factor2m_tl
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) THEN
        profiles_int_tl(iprof)%co2 = &
            profiles_tl(iprof)%co2 * factor + profiles(iprof)%co2 * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  THEN
        profiles_int_tl(iprof)%co = &
            profiles_tl(iprof)%co * factor + profiles(iprof)%co * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) THEN
        profiles_int_tl(iprof)%ch4 = &
            profiles_tl(iprof)%ch4 * factor + profiles(iprof)%ch4 * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%so2) .AND. ASSOCIATED(profiles_tl(iprof)%so2)) THEN
        profiles_int_tl(iprof)%so2 = &
            profiles_tl(iprof)%so2 * factor + profiles(iprof)%so2 * factor_tl
      ENDIF

      IF (ASSOCIATED(profiles_int_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) THEN
        profiles_int_tl(iprof)%n2o = &
            profiles_tl(iprof)%n2o * factor + profiles(iprof)%n2o * factor_tl
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. ppmv dry)

    DO iprof = 1, nprofiles
      profiles_int_tl(iprof)%q     = profiles_tl(iprof)%q
      profiles_int_tl(iprof)%s2m%q = profiles_tl(iprof)%s2m%q
      IF (ASSOCIATED(profiles_int_tl(iprof)%o3) .AND. ASSOCIATED(profiles_tl(iprof)%o3)) THEN
        profiles_int_tl(iprof)%o3    = profiles_tl(iprof)%o3
        profiles_int_tl(iprof)%s2m%o = profiles_tl(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_int_tl(iprof)%co2) .AND. ASSOCIATED(profiles_tl(iprof)%co2)) &
            profiles_int_tl(iprof)%co2 = profiles_tl(iprof)%co2
      IF (ASSOCIATED(profiles_int_tl(iprof)%co)  .AND. ASSOCIATED(profiles_tl(iprof)%co))  &
            profiles_int_tl(iprof)%co  = profiles_tl(iprof)%co
      IF (ASSOCIATED(profiles_int_tl(iprof)%ch4) .AND. ASSOCIATED(profiles_tl(iprof)%ch4)) &
            profiles_int_tl(iprof)%ch4 = profiles_tl(iprof)%ch4
      IF (ASSOCIATED(profiles_int_tl(iprof)%so2) .AND. ASSOCIATED(profiles_tl(iprof)%so2)) &
            profiles_int_tl(iprof)%so2 = profiles_tl(iprof)%so2
      IF (ASSOCIATED(profiles_int_tl(iprof)%n2o) .AND. ASSOCIATED(profiles_tl(iprof)%n2o)) &
            profiles_int_tl(iprof)%n2o = profiles_tl(iprof)%n2o
    ENDDO

  END SELECT


  ! Cloud/aerosol unit conversion

  docld = opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param
  doaer = opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param

  IF (docld .OR. doaer) THEN
    nlevels = profiles(1)%nlevels
    nlayers = profiles(1)%nlayers

    DO iprof = 1, nprofiles
      IF (profiles(iprof)%mmr_cldaer) THEN
        IF (profiles(1)%gas_units == gas_unit_specconc) THEN
          ! Input units were kg/kg moist
          q    = 0.5_jprb * (profiles(iprof)%q(1:nlayers) + profiles(iprof)%q(2:nlevels))
          q_tl = 0.5_jprb * (profiles_tl(iprof)%q(1:nlayers) + profiles_tl(iprof)%q(2:nlevels))
        ELSE
          ! Convert ppmv dry to kg/kg moist
          qlev = profiles_int(iprof)%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(iprof)%q)
          q = 0.5_jprb * (qlev(1:nlayers) + qlev(2:nlevels))

          qlev_tl = profiles_int_tl(iprof)%q * (1._jprb - qlev) / &
                    (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(iprof)%q)
          q_tl = 0.5_jprb * (qlev_tl(1:nlayers) + qlev_tl(2:nlevels))
        ENDIF

        factor(1:nlayers) = (1._jprb + q * (Mair - Mh2o) / Mh2o) * rgc / Mair
        p = 0.5_jprb * 100._jprb * (profiles(iprof)%p(1:nlayers) + profiles(iprof)%p(2:nlevels))
        t = 0.5_jprb * (profiles(iprof)%t(1:nlayers) + profiles(iprof)%t(2:nlevels))
        cldaerfactor(:) = p(:) / (t(:) * factor(1:nlayers))

        factor_tl(1:nlayers) = (q_tl * (Mair - Mh2o) / Mh2o) * rgc / Mair
        t_tl = 0.5_jprb * (profiles_tl(iprof)%t(1:nlayers) + profiles_tl(iprof)%t(2:nlevels))
        cldaerfactor_tl(:) = -t_tl(:) * p(:) / (t(:)**2 * factor(1:nlayers)) - &
                             factor_tl(1:nlayers) * p(:) / (t(:) * factor(1:nlayers)**2)

        IF (opts%interpolation%lgradp) THEN
          p_tl = 0.5_jprb * 100._jprb * (profiles_tl(iprof)%p(1:nlayers) + profiles_tl(iprof)%p(2:nlevels))
          cldaerfactor_tl(:) = cldaerfactor_tl(:) + p_tl(:) / (t(:) * factor(1:nlayers))
        ENDIF

        IF (docld) THEN
          DO lay = 1, nlayers
            profiles_int_tl(iprof)%cloud(:,lay) = &
                profiles_tl(iprof)%cloud(:,lay) * cldaerfactor(lay) + &
                profiles(iprof)%cloud(:,lay) * cldaerfactor_tl(lay)
          ENDDO
        ENDIF

        IF (doaer) THEN
          DO lay = 1, nlayers
            profiles_int_tl(iprof)%aerosols(:,lay) = &
                (profiles_tl(iprof)%aerosols(:,lay) * cldaerfactor(lay) + &
                 profiles(iprof)%aerosols(:,lay) * cldaerfactor_tl(lay)) * coefs%coef_scatt%optp_aer%data(:)%confac
          ENDDO
        ENDIF
      ELSE
        IF (docld) profiles_int_tl(iprof)%cloud    = profiles_tl(iprof)%cloud
        IF (doaer) profiles_int_tl(iprof)%aerosols = profiles_tl(iprof)%aerosols
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_profile_units_tl
