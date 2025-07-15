! Description:
!> @file
!!   Add two RTTOV-SCATT cloud profiles structures.
!
!> @brief
!!   Add two RTTOV-SCATT cloud profiles structures.
!!
!!
!! @param[in,out] cld_profiles       output summed cloud profiles structure
!! @param[in]     cld_profiles1      first input cloud profiles structure
!! @param[in]     cld_profiles2      second input cloud profiles structure
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
SUBROUTINE rttov_add_scatt_prof(cld_profiles, cld_profiles1, cld_profiles2)

  USE rttov_types, ONLY : rttov_profile_cloud
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_profile_cloud), INTENT(INOUT) :: cld_profiles(:)
  TYPE(rttov_profile_cloud), INTENT(IN)    :: cld_profiles1(SIZE(cld_profiles))
  TYPE(rttov_profile_cloud), INTENT(IN)    :: cld_profiles2(SIZE(cld_profiles))
!INTF_END
  INTEGER(KIND=jpim) :: iprof, nprofiles

  nprofiles = SIZE(cld_profiles)

  DO iprof = 1, nprofiles
    cld_profiles(iprof)%cfrac = cld_profiles1(iprof)%cfrac + cld_profiles2(iprof)%cfrac
    cld_profiles(iprof)%ph    = cld_profiles1(iprof)%ph + cld_profiles2(iprof)%ph
    cld_profiles(iprof)%hydro = cld_profiles1(iprof)%hydro + cld_profiles2(iprof)%hydro
    cld_profiles(iprof)%hydro_frac = cld_profiles1(iprof)%hydro_frac + cld_profiles2(iprof)%hydro_frac
  ENDDO
END SUBROUTINE rttov_add_scatt_prof
