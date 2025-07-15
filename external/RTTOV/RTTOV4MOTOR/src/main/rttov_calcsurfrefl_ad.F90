! Description:
!> @file
!!   AD of calculation of surface reflectance values for solar-affected channels.
!
!> @brief
!!   AD of calculation of surface reflectance values for solar-affected channels.
!!
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     sunglint         internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_ad      sunglint model increments
!! @param[in]     fresnrefl        fresnel coefficients
!! @param[in,out] fresnrefl_ad     fresnel coefficient increments
!! @param[in]     solar            flags to indicate channels with solar radiation
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     refl_norm        normalisation factor for direct solar BRDF
!! @param[in]     calcrefl         flags to indicate if BRDF should be provided by RTTOV
!! @param[in,out] emissivity_ad    surface emissivity increments
!! @param[in]     reflectance      surface reflectance input/output structure
!! @param[in,out] reflectance_ad   surface reflectance increments
!! @param[in,out] diffuse_refl_ad  diffuse reflectance increments
!
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
SUBROUTINE rttov_calcsurfrefl_ad( &
              coef,               &
              profiles,           &
              sunglint,           &
              sunglint_ad,        &
              fresnrefl,          &
              fresnrefl_ad,       &
              solar,              &
              chanprof,           &
              refl_norm,          &
              calcrefl,           &
              emissivity_ad,      &
              reflectance,        &
              reflectance_ad,     &
              diffuse_refl_ad)

  USE rttov_types, ONLY :   &
         rttov_chanprof,    &
         rttov_coef,        &
         rttov_emissivity,  &
         rttov_reflectance, &
         rttov_profile,     &
         rttov_sunglint
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : min_windsp, pi, pi_r, surftype_sea
  USE parkind1, ONLY : jpim
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof),    INTENT(IN)              :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)              :: profiles(:)
  TYPE(rttov_coef),        INTENT(IN)              :: coef
  TYPE(rttov_sunglint),    INTENT(IN)              :: sunglint
  TYPE(rttov_sunglint),    INTENT(INOUT)           :: sunglint_ad
  REAL(KIND=jprb),         INTENT(IN)              :: fresnrefl(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)           :: fresnrefl_ad(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)              :: solar(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(IN)              :: refl_norm(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)              :: calcrefl(SIZE(chanprof))
  TYPE(rttov_emissivity),  INTENT(INOUT), OPTIONAL :: emissivity_ad(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(IN)              :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT)           :: reflectance_ad(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)           :: diffuse_refl_ad(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: j, prof, chan, nchanprof
  REAL   (KIND=jprb) :: tmp_refl_ad
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_AD', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    IF (.NOT. calcrefl(j)) THEN
      IF (coef%ss_val_chn(chan) == 2) THEN
        IF (reflectance(j)%diffuse_refl_in > 0._jprb) THEN
          reflectance_ad(j)%diffuse_refl_in = reflectance_ad(j)%diffuse_refl_in + &
                                              diffuse_refl_ad(j) * pi
        ELSE
          reflectance_ad(j)%refl_out = reflectance_ad(j)%refl_out + diffuse_refl_ad(j) * pi
        ENDIF
        diffuse_refl_ad(j) = 0._jprb ! This is important
      ENDIF
      CYCLE
    ENDIF

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      IF (coef%ss_val_chn(chan) == 2) diffuse_refl_ad(j) = 0._jprb

      ! The output reflectance_ad%refl_out should not be scaled by cos(sunzenangle)
      ! Only the refl_ad used in the adjoint sea surface model calculation should be scaled
      IF (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2 > min_windsp**2) THEN
        tmp_refl_ad = reflectance_ad(j)%refl_out / refl_norm(j)
      ELSE
        tmp_refl_ad = reflectance_ad(j)%refl_out
      ENDIF
      reflectance_ad(j)%refl_out = 0._jprb

      sunglint_ad%s(prof)%glint = sunglint_ad%s(prof)%glint + tmp_refl_ad * fresnrefl(j)
      fresnrefl_ad(j)           = fresnrefl_ad(j)           + tmp_refl_ad * sunglint%s(prof)%glint

    ELSE

      IF (coef%ss_val_chn(chan) == 2) THEN
        reflectance_ad(j)%refl_out = reflectance_ad(j)%refl_out + diffuse_refl_ad(j) * pi
        diffuse_refl_ad(j) = 0._jprb ! This is important
      ELSEIF (coef%ss_val_chn(chan) == 1) THEN
        emissivity_ad(j)%emis_out = emissivity_ad(j)%emis_out - reflectance_ad(j)%refl_out * pi_r
        reflectance_ad(j)%refl_out = 0._jprb
      ENDIF

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl_ad
