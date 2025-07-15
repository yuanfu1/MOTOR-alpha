! Description:
!> @file
!!   Scales RTTOV-SCATT cloud profile increment
!
!> @brief
!!   Scales RTTOV-SCATT cloud profile increment
!!
!! @param[in,out] cld_profiles_inc     array of cloud profiles structures to scale
!! @param[in]     factor               scale factor
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_scale_cld_profile_inc(cld_profiles_inc, factor)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_profile_cloud
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile_cloud), INTENT(INOUT) :: cld_profiles_inc(:)
  REAL(jprb),                INTENT(IN)    :: factor
!INTF_END
  INTEGER(KIND=jpim) :: j, nprofiles

  nprofiles = SIZE(cld_profiles_inc)

  DO j = 1, nprofiles
    cld_profiles_inc(j)%ph = cld_profiles_inc(j)%ph * factor
    cld_profiles_inc(j)%cfrac = cld_profiles_inc(j)%cfrac * factor
    cld_profiles_inc(j)%hydro = cld_profiles_inc(j)%hydro * factor
    cld_profiles_inc(j)%hydro_frac = cld_profiles_inc(j)%hydro_frac * factor
  ENDDO

END SUBROUTINE rttov_scale_cld_profile_inc
