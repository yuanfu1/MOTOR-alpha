! Description:
!> @file
!!   Convert input gas profiles to ppmv wrt dry air and cloud/aerosol
!!   profiles to g/m3 and cm-3 respectively where necessary.
!
!> @brief
!!   Convert input gas profiles to ppmv wrt dry air and cloud/aerosol
!!   profiles to g/m3 and cm-3 respectively where necessary.
!!
!! @details
!!   RTTOV offers a choice of units for input gas, cloud and aerosol
!!   profiles. This subroutine converts the input data to units which
!!   are used consistently within RTTOV.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     coefs          coefficients structure
!! @param[in]     profiles       profiles structure containing input profiles
!! @param[in,out] profiles_int   profiles structure containing converted profiles on exit
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
SUBROUTINE rttov_convert_profile_units(opts, coefs, profiles, profiles_int)

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
  TYPE(rttov_profile), INTENT(INOUT) :: profiles_int(SIZE(profiles))
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles, nlevels, lay, nlayers
  REAL(jprb)    :: factor(profiles(1)%nlevels), factor2m
  REAL(jprb)    :: cldaerfactor(profiles(1)%nlayers)
  REAL(jprb)    :: p(profiles(1)%nlayers), t(profiles(1)%nlayers), q(profiles(1)%nlayers)
  REAL(jprb)    :: qlev(profiles(1)%nlevels)
  LOGICAL(jplm) :: docld, doaer
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS', 0_jpim, ZHOOK_HANDLE)

! These conversion formulae could have been put into separate subroutines
! or functions, but the conversions (and especially the TL/AD/K) are more
! efficiently coded following the pattern below.

  nprofiles = SIZE(profiles)


  ! Gas units conversion

  SELECT CASE (profiles(1)%gas_units)

  CASE (gas_unit_specconc)

!     ppmvdry_gas = kgkgwet_gas / (1._jprb - kgkgwet_h2o) * 1.E+06_jprb * Mair / gas_mass(gas_id)

    DO iprof = 1, nprofiles

      factor = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%q)
      IF (opts%rt_all%use_q2m) factor2m = 1.E+06_jprb * Mair / (1._jprb - profiles(iprof)%s2m%q)

      profiles_int(iprof)%q = factor * profiles(iprof)%q / gas_mass(gas_id_watervapour)
      IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%q = &
                               factor2m * profiles(iprof)%s2m%q / gas_mass(gas_id_watervapour)

      IF (ASSOCIATED(profiles_int(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_int(iprof)%o3 = factor * profiles(iprof)%o3 / gas_mass(gas_id_ozone)
        IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%o = &
                                 factor2m * profiles(iprof)%s2m%o / gas_mass(gas_id_ozone)
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) THEN
        profiles_int(iprof)%co2 = factor * profiles(iprof)%co2 / gas_mass(gas_id_co2)
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co)) THEN
        profiles_int(iprof)%co = factor * profiles(iprof)%co / gas_mass(gas_id_co)
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) THEN
        profiles_int(iprof)%ch4 = factor * profiles(iprof)%ch4 / gas_mass(gas_id_ch4)
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%so2) .AND. ASSOCIATED(profiles(iprof)%so2)) THEN
        profiles_int(iprof)%so2 = factor * profiles(iprof)%so2 / gas_mass(gas_id_so2)
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) THEN
        profiles_int(iprof)%n2o = factor * profiles(iprof)%n2o / gas_mass(gas_id_n2o)
      ENDIF
    ENDDO


  CASE (gas_unit_ppmv)

!     ppmvdry_gas = ppmvwet_gas / (1._jprb - ppmvwet_h2o * 1.E-06_jprb)

    DO iprof = 1, nprofiles

      factor = 1._jprb / (1._jprb - profiles(iprof)%q * 1.E-06_jprb)
      IF (opts%rt_all%use_q2m) factor2m = 1._jprb / (1._jprb - profiles(iprof)%s2m%q * 1.E-06_jprb)

      profiles_int(iprof)%q = factor * profiles(iprof)%q
      IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%q = factor2m * profiles(iprof)%s2m%q

      IF (ASSOCIATED(profiles_int(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_int(iprof)%o3 = factor * profiles(iprof)%o3
        IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%o = factor2m * profiles(iprof)%s2m%o
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) THEN
        profiles_int(iprof)%co2 = factor * profiles(iprof)%co2
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co)) THEN
        profiles_int(iprof)%co = factor * profiles(iprof)%co
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) THEN
        profiles_int(iprof)%ch4 = factor * profiles(iprof)%ch4
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%so2) .AND. ASSOCIATED(profiles(iprof)%so2)) THEN
        profiles_int(iprof)%so2 = factor * profiles(iprof)%so2
      ENDIF

      IF (ASSOCIATED(profiles_int(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) THEN
        profiles_int(iprof)%n2o = factor * profiles(iprof)%n2o
      ENDIF
    ENDDO


  CASE DEFAULT ! Assume no conversion (i.e. ppmv dry)

    DO iprof = 1, nprofiles
      profiles_int(iprof)%q = profiles(iprof)%q
      IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%q = profiles(iprof)%s2m%q
      IF (ASSOCIATED(profiles_int(iprof)%o3) .AND. ASSOCIATED(profiles(iprof)%o3)) THEN
        profiles_int(iprof)%o3 = profiles(iprof)%o3
        IF (opts%rt_all%use_q2m) profiles_int(iprof)%s2m%o = profiles(iprof)%s2m%o
      ENDIF
      IF (ASSOCIATED(profiles_int(iprof)%co2) .AND. ASSOCIATED(profiles(iprof)%co2)) &
            profiles_int(iprof)%co2 = profiles(iprof)%co2
      IF (ASSOCIATED(profiles_int(iprof)%co)  .AND. ASSOCIATED(profiles(iprof)%co))  &
            profiles_int(iprof)%co  = profiles(iprof)%co
      IF (ASSOCIATED(profiles_int(iprof)%ch4) .AND. ASSOCIATED(profiles(iprof)%ch4)) &
            profiles_int(iprof)%ch4 = profiles(iprof)%ch4
      IF (ASSOCIATED(profiles_int(iprof)%so2) .AND. ASSOCIATED(profiles(iprof)%so2)) &
            profiles_int(iprof)%so2 = profiles(iprof)%so2
      IF (ASSOCIATED(profiles_int(iprof)%n2o) .AND. ASSOCIATED(profiles(iprof)%n2o)) &
            profiles_int(iprof)%n2o = profiles(iprof)%n2o
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
          q = 0.5_jprb * (profiles(iprof)%q(1:nlayers) + profiles(iprof)%q(2:nlevels))
        ELSE
          ! Convert ppmv dry to kg/kg moist
          qlev = profiles_int(iprof)%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + profiles_int(iprof)%q)
          q = 0.5_jprb * (qlev(1:nlayers) + qlev(2:nlevels))
        ENDIF

        factor(1:nlayers) = (1._jprb + q * (Mair - Mh2o) / Mh2o) * rgc / Mair
        p = 0.5_jprb * 100._jprb * (profiles(iprof)%p(1:nlayers) + profiles(iprof)%p(2:nlevels))  ! units: Pa
        t = 0.5_jprb * (profiles(iprof)%t(1:nlayers) + profiles(iprof)%t(2:nlevels))
        cldaerfactor(:) = p(:) / (t(:) * factor(1:nlayers))

        IF (docld) THEN
          DO lay = 1, nlayers
            profiles_int(iprof)%cloud(:,lay) = profiles(iprof)%cloud(:,lay) * cldaerfactor(lay)
          ENDDO
        ENDIF

        IF (doaer) THEN
          DO lay = 1, nlayers
            profiles_int(iprof)%aerosols(:,lay) = &
                profiles(iprof)%aerosols(:,lay) * cldaerfactor(lay) * coefs%coef_scatt%optp_aer%data(:)%confac
          ENDDO
        ENDIF
      ELSE
        IF (docld) profiles_int(iprof)%cloud    = profiles(iprof)%cloud
        IF (doaer) profiles_int(iprof)%aerosols = profiles(iprof)%aerosols
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CONVERT_PROFILE_UNITS', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_convert_profile_units
