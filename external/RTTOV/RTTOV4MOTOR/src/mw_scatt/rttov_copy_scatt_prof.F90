! Description:
!> @file
!!   Copy an RTTOV-SCATT cloud profile structure.
!
!> @brief
!!   Copy an RTTOV-SCATT cloud profile structure.
!!
!!
!! @param[in,out] cld_profiles1      copy of cloud profiles structure
!! @param[in]     cld_profiles2      input cloud profiles structure
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
SUBROUTINE rttov_copy_scatt_prof(cld_profiles1, cld_profiles2)

  USE rttov_types, ONLY : rttov_profile_cloud
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile_cloud), INTENT(INOUT)  :: cld_profiles1(:)
  TYPE(rttov_profile_cloud), INTENT(IN)     :: cld_profiles2(SIZE(cld_profiles1))
!INTF_END
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof

  nprofiles = SIZE(cld_profiles1)
  DO iprof = 1, nprofiles
    cld_profiles1(iprof)%nlevels          = cld_profiles2(iprof)%nlevels
    cld_profiles1(iprof)%nhydro           = cld_profiles2(iprof)%nhydro
    cld_profiles1(iprof)%nhydro_frac      = cld_profiles2(iprof)%nhydro_frac
    cld_profiles1(iprof)%flux_conversion  = cld_profiles2(iprof)%flux_conversion
    cld_profiles1(iprof)%cfrac            = cld_profiles2(iprof)%cfrac

    cld_profiles1(iprof)%ph               = cld_profiles2(iprof)%ph
    cld_profiles1(iprof)%hydro            = cld_profiles2(iprof)%hydro
    cld_profiles1(iprof)%hydro_frac       = cld_profiles2(iprof)%hydro_frac
  ENDDO
END SUBROUTINE rttov_copy_scatt_prof
