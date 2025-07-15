! Description:
!> @file
!!   AD of radiance calculation
!
!> @brief
!!   AD of radiance calculation
!!
!! @details
!!  The derivative of the Planck function with respect to temperature is
!!
!!                                     C2 * Nu
!!              C1 * C2 * Nu**4 * Exp( ------- )
!!                                        T
!! B'(T,Nu) = ------------------------------------- dT
!!                     (      C2 * Nu       )**2
!!               T**2 *( Exp( ------- ) - 1 )
!!                     (         T          )
!!
!!
!! which can be reduced to the following, with
!!  C1 = C1 * Nu**3
!!  C2 = C2 * Nu
!!
!!              C2 * B(T,Nu) * (C1 + B(T,Nu))
!!  B'(T,Nu) =  ----------------------------- dT
!!                        C1 * T**2
!!
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad      profiles increments
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     thermal          flag to indicate channels with thermal emission
!! @param[in]     auxrad           auxiliary radiance structure
!! @param[in]     auxrad_ad        auxiliary radiance increments
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
SUBROUTINE rttov_calcrad_ad( &
              chanprof,     &
              profiles,     &
              profiles_ad,  &
              coef,         &
              thermal,      &
              auxrad,       &
              auxrad_ad)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_profile, rttov_radiance_aux
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_math_mod, ONLY : planck_ad
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof),     INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),      INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),      INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_coef),         INTENT(IN)    :: coef
  LOGICAL(KIND=jplm),       INTENT(IN)    :: thermal(SIZE(chanprof))
  TYPE(rttov_radiance_aux), INTENT(IN)    :: auxrad
  TYPE(rttov_radiance_aux), INTENT(IN)    :: auxrad_ad
!INTF_END

  REAL   (KIND=jprb) :: t_effective_skin_ad, t_effective_s2m_ad
  REAL   (KIND=jprb) :: t_effective_air_ad(SIZE(profiles(1)%t))
  INTEGER(KIND=jpim) :: ichan, chan, prof
  INTEGER(KIND=jpim) :: nchanprof
!- End of header --------------------------------------------------------
  nchanprof = SIZE(chanprof)

  DO ichan = 1, nchanprof
    IF (.NOT. thermal(ichan)) CYCLE
    chan = chanprof(ichan)%chan
    prof = chanprof(ichan)%prof

    CALL planck_ad(coef%planck1(chan), coef%planck2(chan),&
                   auxrad%skin_t_eff(ichan), t_effective_skin_ad, &
                   auxrad%skin(ichan), auxrad_ad%skin(ichan), acc = .FALSE._jplm)

    CALL planck_ad(coef%planck1(chan), coef%planck2(chan), &
                   auxrad%surf_t_eff(ichan), t_effective_s2m_ad, &
                   auxrad%surfair(ichan), auxrad_ad%surfair(ichan), acc = .FALSE._jplm)

    CALL planck_ad(coef%planck1(chan), coef%planck2(chan),&
                   auxrad%air_t_eff(:,ichan), t_effective_air_ad(:), &
                   auxrad%air(:,ichan), auxrad_ad%air(:,ichan), acc = .FALSE._jplm)

    IF (coef%ff_val_bc) THEN
      profiles_ad(prof)%skin%t = profiles_ad(prof)%skin%t + &
        coef%ff_bcs(chan) * t_effective_skin_ad
      profiles_ad(prof)%s2m%t = profiles_ad(prof)%s2m%t + &
        coef%ff_bcs(chan) * t_effective_s2m_ad
      profiles_ad(prof)%t(:) = profiles_ad(prof)%t(:) + &
        coef%ff_bcs(chan) * t_effective_air_ad(:)
    ELSE
      profiles_ad(prof)%skin%t = profiles_ad(prof)%skin%t + &
        t_effective_skin_ad
      profiles_ad(prof)%s2m%t = profiles_ad(prof)%s2m%t + &
        t_effective_s2m_ad
      profiles_ad(prof)%t(:) = profiles_ad(prof)%t(:) + &
        t_effective_air_ad(:)
    ENDIF
  ENDDO

END SUBROUTINE rttov_calcrad_ad
