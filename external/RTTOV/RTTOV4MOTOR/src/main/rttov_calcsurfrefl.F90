! Description:
!> @file
!!   Calculate surface reflectance values for solar-affected channels.
!
!> @brief
!!   Calculate surface reflectance values for solar-affected channels.
!!
!! @details
!!   Two reflectance values are required in the RTE integration:
!!   1) the BRDF for the direct surface-reflected solar beam, only used for
!!      solar radiation
!!   2) the reflectance for downward-emitted or downward-scattered
!!      surface-reflected radiation ("diffuse reflectance"), calculated for
!!      all channels
!!
!!   For solar radiation, the diffuse reflectance represents the specular BRDF
!!   for radiation scattered downwards at the satellite zenith angle which is
!!   then reflected towards the satellite. For thermal radiation, it is the
!!   reflectance used for downwelling surface-reflected atmospheric emission
!!   (which applies under any specularity assumption).
!!
!!   The BRDF values are input/output in the refl_in/refl_out members of the
!!   reflectance structure. In all cases, the refl_out array contains the BRDF
!!   values used by RTTOV for all solar-affected channels.
!!
!!   The diffuse reflectance values are input/output in the diffuse_refl_in/
!!   diffuse_refl_out members of the reflectance structure. The diffuse
!!   reflectance values can only be specified by the user for pure-solar
!!   channels (see below). In all cases, the diffuse_refl_out array contains
!!   reflectance values used by RTTOV for all channels (including thermal-only
!!   channels).
!!
!!   The BRDF and diffuse reflectance values used depend on the channel
!!   wavelength, the surface type, and whether calcrefl is true or false:
!!
!!   -------------------------------------------------------
!!   Pure-solar channels (less than 3 microns):
!!   -------------------------------------------------------
!!   If calcrefl is TRUE:
!!
!!     Sea surface: 
!!       - The BRDF is calculated from the sunglint model (see rttov_refsun)
!!         and the fresnel coefficients (see rttov_fresnel). The BRDF derived
!!         from the USGS water reflectance spectrum (see below) is added to
!!         the calculated sunglint BRDF to avoid unrealistically small BRDFs
!!         away from sunglint.
!!
!!       - The diffuse reflectance is taken from the USGS
!!         water spectra (these are interpolated onto channel wavenumbers in
!!         rttov_init_coef when the coef file is read).
!!
!!     Land/sea-ice surface:
!!       - BRDF computed assuming a fixed albedo: 0.3 for land, 0.8 for
!!         sea-ice (this is not generally recommended, it is better to set
!!         calcrefl to FALSE and supply a BRDF value).
!!
!!       - The diffuse reflectance is the same as the BRDF.
!!
!!   If calcrefl is FALSE:
!!
!!     All surface types:
!!       - The input BRDF value in refl_in is used.
!!
!!       - If the input diffuse_refl_in value is greater than zero, this is
!!         used, otherwise the input BRDF value is used.
!!
!!
!!   -------------------------------------------------------
!!   Mixed solar+thermal channels (between 3 and 5 microns):
!!   -------------------------------------------------------
!!
!!   If calcrefl is TRUE:
!!
!!     Sea surface: 
!!       - The BRDF is calculated from the sunglint model (see rttov_refsun)
!!         and the fresnel coefficients (see rttov_fresnel).
!!
!!       - *The diffuse reflectance is (1-emissivity)/pi.
!!
!!     Land/sea-ice surface:
!!       - *The BRDF is (1-emissivity)/pi.
!!
!!       - *The diffuse reflectance is the same as the BRDF.
!!
!!   If calcrefl is FALSE:
!!
!!     All surface types:
!!       - The input BRDF value in refl_in is used
!!
!!       - *The diffuse reflectance is (1-emissivity)/pi.
!!
!!
!!   -------------------------------------------------------
!!   Thermal channels (greater than 5 microns):
!!   -------------------------------------------------------
!!
!!   - In this case the BRDF is not used.
!!
!!   - *The diffuse reflectance is (1-emissivity)/pi.
!!
!!
!!   *In these cases, the value (1-emissivity)/pi is used for both thermal
!!    emission and solar scattered radiation and has already been calculated by
!!    the surface emissivity calculations so is not modified in this
!!    subroutine. It cannot be modified by the user even if calcrefl is false.
!!
!!
!!   As noted above, for visible/near-IR channels the sea sunglint model gives
!!   extremely small BRDFs away from the sunglint region which results in an
!!   underestimation of the TOA reflectance. Therefore the BRDF derived from the
!!   USGS water reflectance spectrum is added to the calculated sunglint BRDF.
!!   This does not apply to short-wave IR channels.
!!
!!   The BRDF affects only the direct surface-reflected solar term. When this is
!!   applied in rttov_integrate it is multipled by refl_norm which is
!!   COS(sunzenangle). This does not apply where the sunglint model is used
!!   since this explicitly treats the wind-roughened surface: in this case the
!!   BRDF is divided by refl_norm to give a BRDF-like value to be used in the
!!   RTE integration.
!!
!!   For pure-solar channels, if calcrefl is false and the input diffuse_refl_in
!!   value is greater than zero then the input value is used for the diffuse
!!   reflectance. If the input diffuse_refl_in is less than or equal to zero,
!!   then the input surface BRDF is used for the diffuse reflectance. Note
!!   that for mixed thermal+solar channels, the diffuse reflectance cannot be
!!   specified independently (it is set to be consistent with the surface
!!   emissivity).
!!
!!   Unless specified by the user, the diffuse reflectance only differs from the
!!   direct solar BRDF when using the sunglint model and in this case the values
!!   used are fixed (from the USGS reflectance spectra).
!!
!!   The diffuse reflectances are stored in an internal array diffuse_refl which
!!   is used for all channels (solar and thermal). The values in this array
!!   include a factor pi compared to the output diffuse reflectance values:
!!   this factor is divided out where they are used in rttov_integrate.
!!
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in]     fresnrefl      fresnel coefficients
!! @param[in]     solar          flags to indicate channels with solar radiation
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[out]    refl_norm      normalisation factor for direct solar BRDF
!! @param[in]     calcrefl       flags to indicate if BRDF should be provided by RTTOV
!! @param[in]     emissivity     surface emissivity input/output structure
!! @param[in,out] reflectance    surface reflectance input/output structure
!! @param[in,out] diffuse_refl   internal diffuse reflectance array
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
SUBROUTINE rttov_calcsurfrefl( &
              coef,            &
              profiles,        &
              sunglint,        &
              fresnrefl,       &
              solar,           &
              chanprof,        &
              refl_norm,       &
              calcrefl,        &
              emissivity,      &
              reflectance,     &
              diffuse_refl)

  USE rttov_types, ONLY :   &
         rttov_chanprof,    &
         rttov_coef,        &
         rttov_emissivity,  &
         rttov_reflectance, &
         rttov_profile,     &
         rttov_sunglint
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :     &
         min_windsp,          &
         pi, pi_r,            &
         deg2rad,             &
         surftype_sea,        &
         surftype_seaice,     &
         watertype_fresh_water
  USE parkind1, ONLY : jpim
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof),    INTENT(IN)             :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)             :: profiles(:)
  TYPE(rttov_coef),        INTENT(IN)             :: coef
  TYPE(rttov_sunglint),    INTENT(IN)             :: sunglint
  REAL(KIND=jprb),         INTENT(IN)             :: fresnrefl(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)             :: solar(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(OUT)            :: refl_norm(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)             :: calcrefl(SIZE(chanprof))
  TYPE(rttov_emissivity),  INTENT(IN),   OPTIONAL :: emissivity(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT)          :: reflectance(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)          :: diffuse_refl(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: j, prof, chan, nchanprof
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    refl_norm(j) = COS(profiles(prof)%sunzenangle * deg2rad)

    IF (.NOT. calcrefl(j)) THEN
      IF (coef%ss_val_chn(chan) == 2) THEN
        ! Normalisation by 1/pi occurs at point of use in integrate
        IF (reflectance(j)%diffuse_refl_in > 0._jprb) THEN
          diffuse_refl(j) = reflectance(j)%diffuse_refl_in * pi
        ELSE
          diffuse_refl(j) = reflectance(j)%refl_out * pi
        ENDIF
      ENDIF
      CYCLE
    ENDIF

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      reflectance(j)%refl_out = sunglint%s(prof)%glint * fresnrefl(j)

      ! Scale the reflectance by cos(sunzenangle) because the sea surface model
      ! takes the solar zenith angle in account. This ensures the output reflectance
      ! is BRDF-like.
      IF (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2 > min_windsp**2) THEN
        reflectance(j)%refl_out = reflectance(j)%refl_out / refl_norm(j)
      ENDIF

      IF (coef%ss_val_chn(chan) == 2) THEN
        ! For pure solar channels this is the only case where the diffuse and direct solar reflectances differ
        IF (profiles(prof)%skin%watertype == watertype_fresh_water) THEN
          diffuse_refl(j) = coef%refl_visnir_fw(chan)
        ELSE
          diffuse_refl(j) = coef%refl_visnir_ow(chan)
        ENDIF
        reflectance(j)%refl_out = reflectance(j)%refl_out + diffuse_refl(j) * pi_r
      ENDIF

    ELSE

      ! For land and sea-ice, the reflectance values represent directional-hemispherical
      ! albedos. The reflectance(:)%refl_out value is a BRDF, hence the normalisation by 1/pi.

      IF (coef%ss_val_chn(chan) == 1 .AND. PRESENT(emissivity)) THEN
        reflectance(j)%refl_out = (1._jprb - emissivity(j)%emis_out) * pi_r
      ELSE
        IF (profiles(prof)%skin%surftype == surftype_seaice) THEN
          reflectance(j)%refl_out = 0.8_jprb * pi_r  ! Rough value for sea-ice
        ELSE
          reflectance(j)%refl_out = 0.3_jprb * pi_r  ! Rough value for land
        ENDIF
      ENDIF

      ! Normalisation by 1/pi occurs at point of use in integrate
      IF (coef%ss_val_chn(chan) == 2) diffuse_refl(j) = reflectance(j)%refl_out * pi

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl
