! Description:
!> @file
!!   Initialise an internal profiles structure.
!
!> @brief
!!   Initialise an internal profiles structure.
!!
!! @details
!!   This is used internally by RTTOV, primarily the AD/K,
!!   to initialise the internal profiles structure containing
!!   profile data in the units used internally by RTTOV.
!!
!! @param[in,out]  profiles_int   Array of profiles structures
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
SUBROUTINE rttov_init_prof_internal(profiles_int)

  USE rttov_types, ONLY : rttov_profile
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_profile), INTENT(INOUT) :: profiles_int(:)
!INTF_END

  INTEGER(KIND=jpim) :: j, nprofiles

  nprofiles = SIZE(profiles_int)

  DO j = 1, nprofiles
    profiles_int(j)%s2m%q = 0._jprb
    profiles_int(j)%s2m%o = 0._jprb
    profiles_int(j)%q = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%o3)      ) profiles_int(j)%o3       = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%co2)     ) profiles_int(j)%co2      = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%n2o)     ) profiles_int(j)%n2o      = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%co)      ) profiles_int(j)%co       = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%ch4)     ) profiles_int(j)%ch4      = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%so2)     ) profiles_int(j)%so2      = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%aerosols)) profiles_int(j)%aerosols = 0._jprb
    IF (ASSOCIATED(profiles_int(j)%cloud)   ) profiles_int(j)%cloud    = 0._jprb
  ENDDO
END SUBROUTINE rttov_init_prof_internal
