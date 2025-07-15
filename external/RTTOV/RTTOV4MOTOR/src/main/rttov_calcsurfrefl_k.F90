! Description:
!> @file
!!   Jacobian of calculation of surface reflectance values for solar-affected channels.
!
!> @brief
!!   Jacobian of calculation of surface reflectance values for solar-affected channels.
!!
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     sunglint         internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_k       sunglint model increments
!! @param[in]     fresnrefl        fresnel coefficients
!! @param[in,out] fresnrefl_k      fresnel coefficient increments
!! @param[in]     solar            flags to indicate channels with solar radiation
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     refl_norm        normalisation factor for direct solar BRDF
!! @param[in]     calcrefl         flags to indicate if BRDF should be provided by RTTOV
!! @param[in,out] emissivity_k     surface emissivity increments
!! @param[in]     reflectance      surface reflectance input/output structure
!! @param[in,out] reflectance_k    surface reflectance increments
!! @param[in,out] diffuse_refl_k   diffuse reflectance increments
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
SUBROUTINE rttov_calcsurfrefl_k( &
              coef,              &
              profiles,          &
              sunglint,          &
              sunglint_k,        &
              fresnrefl,         &
              fresnrefl_k,       &
              solar,             &
              chanprof,          &
              refl_norm,         &
              calcrefl,          &
              emissivity_k,      &
              reflectance,       &
              reflectance_k,     &
              diffuse_refl_k)

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
  TYPE(rttov_sunglint),    INTENT(INOUT)           :: sunglint_k
  REAL(KIND=jprb),         INTENT(IN)              :: fresnrefl(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)           :: fresnrefl_k(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)              :: solar(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(IN)              :: refl_norm(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)              :: calcrefl(SIZE(chanprof))
  TYPE(rttov_emissivity),  INTENT(INOUT), OPTIONAL :: emissivity_k(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(IN)              :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT)           :: reflectance_k(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)           :: diffuse_refl_k(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: j, prof, chan, nchanprof
  REAL   (KIND=jprb) :: tmp_refl_k
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_K', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    IF (.NOT. calcrefl(j)) THEN
      IF (coef%ss_val_chn(chan) == 2) THEN
        IF (reflectance(j)%diffuse_refl_in > 0._jprb) THEN
          reflectance_k(j)%diffuse_refl_in = reflectance_k(j)%diffuse_refl_in + &
                                             diffuse_refl_k(j) * pi
        ELSE
          reflectance_k(j)%refl_out = reflectance_k(j)%refl_out + diffuse_refl_k(j) * pi
        ENDIF
        diffuse_refl_k(j) = 0._jprb ! This is important
      ENDIF
      CYCLE
    ENDIF

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      IF (coef%ss_val_chn(chan) == 2) diffuse_refl_k(j) = 0._jprb

      ! The output reflectance_k%refl_out should not be scaled by cos(sunzenangle)
      ! Only the refl_k used in the Jacobian sea surface model calculation should be scaled
      IF (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2 > min_windsp**2) THEN
        tmp_refl_k = reflectance_k(j)%refl_out / refl_norm(j)
      ELSE
        tmp_refl_k = reflectance_k(j)%refl_out
      ENDIF
      reflectance_k(j)%refl_out = 0._jprb

      sunglint_k%s(j)%glint = sunglint_k%s(j)%glint + tmp_refl_k * fresnrefl(j)
      fresnrefl_k(j)        = fresnrefl_k(j)        + tmp_refl_k * sunglint%s(prof)%glint

    ELSE

      IF (coef%ss_val_chn(chan) == 2) THEN
        reflectance_k(j)%refl_out = reflectance_k(j)%refl_out + diffuse_refl_k(j) * pi
        diffuse_refl_k(j) = 0._jprb ! This is important
      ELSEIF (coef%ss_val_chn(chan) == 1) THEN
        emissivity_k(j)%emis_out = emissivity_k(j)%emis_out - reflectance_k(j)%refl_out * pi_r
        reflectance_k(j)%refl_out = 0._jprb
      ENDIF

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl_k
