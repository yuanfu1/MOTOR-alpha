! Description:
!> @file
!!   Calculate top-of-atmosphere reflectances from radiances
!
!> @brief
!!   Calculate top-of-atmosphere reflectances from radiances
!!
!! @details
!!   The input solar_spectrum contains the top-of-atmosphere solar spectral
!!   irradiance for each channel (corrected for the Earth-sun distance - see
!!   rttov_calc_solar_spec_esd) with units mW/m2/cm-1. This is multiplied by
!!   cos(sol_zen_ang) to obtain the surface irradiance. The factor of 1/pi is
!!   present to convert from a BRDF to the output BRF (bi-directional
!!   reflectance factor).
!!
!!   Note that the normalisation by the ToA solar irradiance means that the
!!   output reflectances are independent of the solar irradiance value (unlike
!!   the simulated radiance).
!!
!!   Reflectances are not calculated for channels at wavelengths above 5um (i.e.
!!   with insignificant solar radiation).
!!
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     solar_spectrum   top-of-atmosphere solar irradiance values
!! @param[in]     solar            flag to indicate channels with solar radiation
!! @param[in,out] rad              radiance structure
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
SUBROUTINE rttov_calcsatrefl(chanprof, profiles, solar_spectrum, solar, rad) 

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
  TYPE(rttov_radiance), INTENT(INOUT) :: rad
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

      rad%refl_clear(i) = rad%clear(i) / toarad_norm
      rad%refl(i)       = rad%total(i) / toarad_norm
    ENDIF
  ENDDO

END SUBROUTINE rttov_calcsatrefl
