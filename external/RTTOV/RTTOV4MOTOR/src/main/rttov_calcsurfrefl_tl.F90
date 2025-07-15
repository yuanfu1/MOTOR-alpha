! Description:
!> @file
!!   TL of calculation of surface reflectance values for solar-affected channels.
!
!> @brief
!!   TL of calculation of surface reflectance values for solar-affected channels.
!!
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     sunglint         internal structure for sea surface BRDF model variables
!! @param[in]     sunglint_tl      sunglint model perturbations
!! @param[in]     fresnrefl        fresnel coefficients
!! @param[in]     fresnrefl_tl     fresnel coefficient perturbations
!! @param[in]     solar            flags to indicate channels with solar radiation
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     refl_norm        normalisation factor for direct solar BRDF
!! @param[in]     calcrefl         flags to indicate if BRDF should be provided by RTTOV
!! @param[in]     emissivity_tl    surface emissivity perturbations
!! @param[in]     reflectance      surface reflectance input/output structure
!! @param[in,out] reflectance_tl   surface reflectance perturbations
!! @param[in,out] diffuse_refl_tl  diffuse reflectance perturbations
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
SUBROUTINE rttov_calcsurfrefl_tl( &
              coef,               &
              profiles,           &
              sunglint,           &
              sunglint_tl,        &
              fresnrefl,          &
              fresnrefl_tl,       &
              solar,              &
              chanprof,           &
              refl_norm,          &
              calcrefl,           &
              emissivity_tl,      &
              reflectance,        &
              reflectance_tl,     &
              diffuse_refl_tl)

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

  TYPE(rttov_chanprof),    INTENT(IN)             :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)             :: profiles(:)
  TYPE(rttov_coef),        INTENT(IN)             :: coef
  TYPE(rttov_sunglint),    INTENT(IN)             :: sunglint
  TYPE(rttov_sunglint),    INTENT(IN)             :: sunglint_tl
  REAL(KIND=jprb),         INTENT(IN)             :: fresnrefl(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(IN)             :: fresnrefl_tl(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)             :: solar(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(IN)             :: refl_norm(SIZE(chanprof))
  LOGICAL(KIND=jplm),      INTENT(IN)             :: calcrefl(SIZE(chanprof))
  TYPE(rttov_emissivity),  INTENT(IN),   OPTIONAL :: emissivity_tl(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(IN)             :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(INOUT)          :: reflectance_tl(SIZE(chanprof))
  REAL(KIND=jprb),         INTENT(INOUT)          :: diffuse_refl_tl(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: j, prof, chan, nchanprof
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_TL', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    IF (.NOT. calcrefl(j)) THEN
      IF (coef%ss_val_chn(chan) == 2) THEN
        IF (reflectance(j)%diffuse_refl_in > 0._jprb) THEN
          diffuse_refl_tl(j) = reflectance_tl(j)%diffuse_refl_in * pi
        ELSE
          diffuse_refl_tl(j) = reflectance_tl(j)%refl_out * pi
        ENDIF
      ENDIF
      CYCLE
    ENDIF

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      reflectance_tl(j)%refl_out = sunglint_tl%s(prof)%glint * fresnrefl(j) + &
                                   sunglint%s(prof)%glint * fresnrefl_tl(j)

      IF (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2 > min_windsp**2) THEN
        reflectance_tl(j)%refl_out = reflectance_tl(j)%refl_out / refl_norm(j)
      ENDIF

      IF (coef%ss_val_chn(chan) == 2) THEN
        diffuse_refl_tl(j) = 0._jprb
!         reflectance_tl(j)%refl_out = reflectance_tl(j)%refl_out + diffuse_refl_tl(j) * pi_r  ! Commented to avoid null op
      ENDIF

    ELSE

      IF (coef%ss_val_chn(chan) == 1) THEN
        reflectance_tl(j)%refl_out = - emissivity_tl(j)%emis_out * pi_r
      ELSE
! User should supply TL, same as for ISEM
!         reflectance_tl(j)%refl_out = 0._jprb
      ENDIF

      IF (coef%ss_val_chn(chan) == 2) diffuse_refl_tl(j) = reflectance_tl(j)%refl_out * pi

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl_tl
