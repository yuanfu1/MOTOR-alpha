! Description:
!> @file
!!   Initialise a profiles structure.
!
!> @brief
!!   Initialise a profiles structure.
!!
!! @details
!!   The argument is a profiles array: every element of the structure
!!   is initialised except nlevels, nlayers, gas_units and mmr_cldaer
!!   which are initialised by rttov_alloc_prof.
!!
!!   An optional pressure profile may be supplied which will be
!!   copied to the p(:) member of each profile.
!!
!! @param[in,out]  profiles   Array of profiles structures
!! @param[in]      p          Pressure profile to copy to profiles(:)%p, optional
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
SUBROUTINE rttov_init_prof(profiles, p)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_profile
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_profile), INTENT(INOUT)          :: profiles(:)
  REAL(KIND=jprb),     INTENT(IN),   OPTIONAL :: p(:)
!INTF_END

  INTEGER(KIND=jpim) :: j, nprofiles

  nprofiles = SIZE(profiles)

  DO j = 1, nprofiles
    IF (PRESENT(p)) THEN
      profiles(j)%p = p
    ELSE
      profiles(j)%p = 0._jprb
    ENDIF
    profiles(j)%t = 0._jprb
    profiles(j)%q = 0._jprb
    IF (ASSOCIATED(profiles(j)%co2)     ) profiles(j)%co2      = 0._jprb
    IF (ASSOCIATED(profiles(j)%o3)      ) profiles(j)%o3       = 0._jprb
    IF (ASSOCIATED(profiles(j)%n2o)     ) profiles(j)%n2o      = 0._jprb
    IF (ASSOCIATED(profiles(j)%co)      ) profiles(j)%co       = 0._jprb
    IF (ASSOCIATED(profiles(j)%ch4)     ) profiles(j)%ch4      = 0._jprb
    IF (ASSOCIATED(profiles(j)%so2)     ) profiles(j)%so2      = 0._jprb
    IF (ASSOCIATED(profiles(j)%clw)     ) profiles(j)%clw      = 0._jprb
    IF (ASSOCIATED(profiles(j)%aerosols)) profiles(j)%aerosols = 0._jprb
    IF (ASSOCIATED(profiles(j)%cloud)   ) profiles(j)%cloud    = 0._jprb
    IF (ASSOCIATED(profiles(j)%cfrac)   ) profiles(j)%cfrac    = 0._jprb
    IF (ASSOCIATED(profiles(j)%clwde)   ) profiles(j)%clwde    = 0._jprb
    IF (ASSOCIATED(profiles(j)%icede)   ) profiles(j)%icede    = 0._jprb
    profiles(j)%id                 = " "
    profiles(j)%date               = (/ 1950_jpim, 1_jpim, 1_jpim /)
    profiles(j)%time               = (/ 0_jpim, 0_jpim, 0_jpim /)
    profiles(j)%zenangle           = 0._jprb
    profiles(j)%azangle            = 0._jprb
    profiles(j)%sunzenangle        = 0._jprb
    profiles(j)%sunazangle         = 0._jprb
    profiles(j)%elevation          = 0._jprb
    profiles(j)%latitude           = 0._jprb
    profiles(j)%longitude          = 0._jprb
    profiles(j)%ctp                = 0._jprb
    profiles(j)%cfraction          = 0._jprb
    profiles(j)%s2m%t              = 0._jprb
    profiles(j)%s2m%q              = 0._jprb
    profiles(j)%s2m%o              = 0._jprb
    profiles(j)%s2m%p              = 0._jprb
    profiles(j)%s2m%u              = 0._jprb
    profiles(j)%s2m%v              = 0._jprb
    profiles(j)%s2m%wfetc          = 0._jprb
    profiles(j)%skin%t             = 0._jprb
    profiles(j)%skin%fastem        = 0._jprb
    profiles(j)%skin%surftype      = 0_jpim
    profiles(j)%skin%watertype     = 0_jpim
    profiles(j)%skin%salinity      = 0._jprb
    profiles(j)%skin%foam_fraction = 0._jprb
    profiles(j)%skin%snow_fraction = 0._jprb
    profiles(j)%skin%soil_moisture = 0._jprb
    profiles(j)%clwde_param        = 1_jpim ! Currently the only valid option
    profiles(j)%clw_scheme         = -1_jpim
    profiles(j)%icede_param        = -1_jpim
    profiles(j)%ice_scheme         = -1_jpim
    profiles(j)%Be                 = 0._jprb
    profiles(j)%cosbk              = 0._jprb
  ENDDO
END SUBROUTINE
