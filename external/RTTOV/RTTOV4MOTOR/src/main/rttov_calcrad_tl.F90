! Description:
!> @file
!!   TL of radiance calculation
!
!> @brief
!!   TL of radiance calculation
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
!! @param[in]     profiles_tl      profiles perturbations
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     thermal          flag to indicate channels with thermal emission
!! @param[in]     auxrad           auxiliary radiance structure
!! @param[in,out] auxrad_tl        auxiliary radiance perturbations
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
SUBROUTINE rttov_calcrad_tl( &
              chanprof,     &
              profiles,     &
              profiles_tl,  &
              coef,         &
              thermal,      &
              auxrad,       &
              auxrad_tl)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_profile, rttov_radiance_aux
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY : planck_tl
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof),     INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),      INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),      INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(rttov_coef),         INTENT(IN)    :: coef
  LOGICAL(KIND=jplm),       INTENT(IN)    :: thermal(SIZE(chanprof))
  TYPE(rttov_radiance_aux), INTENT(IN)    :: auxrad
  TYPE(rttov_radiance_aux), INTENT(INOUT) :: auxrad_tl
!INTF_END

  REAL   (KIND=jprb) :: t_effective_skin_tl, t_effective_s2m_tl
  REAL   (KIND=jprb) :: t_effective_air_tl(SIZE(profiles(1)%t))
  INTEGER(KIND=jpim) :: chan, prof, ichan
  INTEGER(KIND=jpim) :: nchanprof
!- End of header --------------------------------------------------------
  nchanprof = SIZE(chanprof)

  DO ichan = 1, nchanprof
    IF (.NOT. thermal(ichan)) CYCLE
    chan = chanprof(ichan)%chan
    prof = chanprof(ichan)%prof

    IF (coef%ff_val_bc) THEN
       t_effective_skin_tl = coef%ff_bcs(chan) * profiles_tl(prof)%skin%t
       t_effective_s2m_tl = coef%ff_bcs(chan) * profiles_tl(prof)%s2m%t
       t_effective_air_tl = coef%ff_bcs(chan) * profiles_tl(prof)%t(:)
    ELSE
       t_effective_skin_tl = profiles_tl(prof)%skin%t
       t_effective_s2m_tl = profiles_tl(prof)%s2m%t
       t_effective_air_tl(:) = profiles_tl(prof)%t(:)
    ENDIF

    CALL planck_tl(coef%planck1(chan), coef%planck2(chan),&
                   auxrad%skin_t_eff(ichan), t_effective_skin_tl, &
                   auxrad%skin(ichan), auxrad_tl%skin(ichan))

    CALL planck_tl(coef%planck1(chan), coef%planck2(chan), &
                   auxrad%surf_t_eff(ichan), t_effective_s2m_tl, &
                   auxrad%surfair(ichan), auxrad_tl%surfair(ichan))

    CALL planck_tl(coef%planck1(chan), coef%planck2(chan),&
                   auxrad%air_t_eff(:,ichan), t_effective_air_tl(:), &
                   auxrad%air(:,ichan), auxrad_tl%air(:,ichan))
  ENDDO

END SUBROUTINE rttov_calcrad_tl
