! Description:
!> @file
!!   TL of top-of-atmosphere reflectance calculation
!
!> @brief
!!   TL of top-of-atmosphere reflectance calculation
!!
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     solar_spectrum   top-of-atmosphere solar irradiance values
!! @param[in]     solar            flag to indicate channels with solar radiation
!! @param[in,out] rad_tl           radiance perturbations
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
SUBROUTINE rttov_calcsatrefl_tl(chanprof, profiles, solar_spectrum, solar, rad_tl) 

  USE rttov_types, ONLY : rttov_chanprof, rttov_profile, rttov_radiance
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : pi_r, deg2rad
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  REAL(KIND=jprb),      INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  LOGICAL(KIND=jplm),   INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_radiance), INTENT(INOUT) :: rad_tl
!INTF_END
  REAL   (KIND=jprb) :: toarad_norm
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: i, prof
!- End of header ------------------------------------------------------
  nchanprof = SIZE(chanprof)

  DO i = 1, nchanprof
    IF (solar(i)) THEN
      prof = chanprof(i)%prof

      toarad_norm = solar_spectrum(i) * pi_r * COS(profiles(prof)%sunzenangle * deg2rad)

      rad_tl%refl_clear(i) = rad_tl%clear(i) / toarad_norm
      rad_tl%refl(i)       = rad_tl%total(i) / toarad_norm
    ENDIF
  ENDDO

END SUBROUTINE rttov_calcsatrefl_tl
