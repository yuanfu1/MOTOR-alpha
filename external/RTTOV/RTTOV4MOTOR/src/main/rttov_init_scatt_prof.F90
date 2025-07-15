! Description:
!> @file
!!   Initialise an RTTOV-SCATT cloud profiles structure.
!
!> @brief
!!   Initialise an RTTOV-SCATT cloud profiles structure.
!!
!! @details
!!   The argument is an RTTOV-SCATT cloud profiles array: all
!!   array members will be initialised as well as cfrac.
!!
!! @param[in,out]  cld_profiles   Array of cloud profiles structures
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_init_scatt_prof(cld_profiles)

  USE rttov_types, ONLY : rttov_profile_cloud
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_profile_cloud), INTENT(INOUT) :: cld_profiles(:)
!INTF_END

  INTEGER(jpim) :: iprof, nprofiles

  nprofiles = SIZE(cld_profiles)

  DO iprof = 1, nprofiles
    cld_profiles(iprof)%cfrac      = 0._jprb
    cld_profiles(iprof)%ph         = 0._jprb
    cld_profiles(iprof)%hydro      = 0._jprb
    cld_profiles(iprof)%hydro_frac = 0._jprb
  ENDDO
END SUBROUTINE rttov_init_scatt_prof
