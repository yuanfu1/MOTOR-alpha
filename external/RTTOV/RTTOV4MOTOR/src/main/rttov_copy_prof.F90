! Description:
!> @file
!!   Copy a profile structure.
!
!> @brief
!!   Copy a profile structure.
!!
!! @details
!!   This subroutine is used in two places: in the direct and TL models
!!   this is used to populate the "profiles_coef" structure which is
!!   on coefficient levels and in RTTOV's internal units for gases.
!!
!!   The larray and lscalar flags allow the array and scalar members
!!   of the structure to be copied independently. If the RTTOV interpolator
!!   is used then this populates the profile arrays and only scalar values
!!   need to be copied into the profiles_coef structure.
!!
!!   If the profiles_gas argument is supplied then the gas arrays and
!!   scalars are copied from this structure instead of the profiles2
!!   argument. This is used to copy the gas data in RTTOV's internal
!!   units to profiles_coef.
!!
!!   The second place this subroutine is used is in the test suite to copy
!!   entire profile structures verbatim. In this case the optional arguments
!!   are not used.
!!
!! @param[in,out] profiles1      copy of profiles structure
!! @param[in]     profiles2      input profiles structure
!! @param[in]     larray         flag to indicate array members should be copied, optional
!! @param[in]     lscalar        flag to indicate scalar members should be copied, optional
!! @param[in]     profiles_gas   input profiles structure from which to take gas profile data, optional
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
SUBROUTINE rttov_copy_prof( &
              profiles1,    &
              profiles2,    &
              larray,       &
              lscalar,      &
              profiles_gas)

  USE rttov_types, ONLY : rttov_profile
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile), INTENT(INOUT)        :: profiles1(:)
  TYPE(rttov_profile), INTENT(IN)           :: profiles2(SIZE(profiles1))
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: larray
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: lscalar
  TYPE(rttov_profile), INTENT(IN), OPTIONAL :: profiles_gas(SIZE(profiles1))
!INTF_END
  LOGICAL(KIND=jplm) :: larray1
  LOGICAL(KIND=jplm) :: lscalar1
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  lscalar1 = .TRUE.
  IF (PRESENT(lscalar)) lscalar1 = lscalar
  larray1 = .TRUE.
  IF (PRESENT(larray)) larray1 = larray
  nprofiles = SIZE(profiles1)
  DO iprof = 1, nprofiles
    IF (PRESENT(profiles_gas)) THEN
      profiles1(iprof)%gas_units     = profiles_gas(iprof)%gas_units
    ELSE
      profiles1(iprof)%gas_units     = profiles2(iprof)%gas_units
    ENDIF
    IF (lscalar1) THEN
      profiles1(iprof)%id            = profiles2(iprof)%id
      profiles1(iprof)%date          = profiles2(iprof)%date
      profiles1(iprof)%time          = profiles2(iprof)%time
      profiles1(iprof)%mmr_cldaer    = profiles2(iprof)%mmr_cldaer
      profiles1(iprof)%zenangle      = profiles2(iprof)%zenangle
      profiles1(iprof)%azangle       = profiles2(iprof)%azangle
      profiles1(iprof)%sunzenangle   = profiles2(iprof)%sunzenangle
      profiles1(iprof)%sunazangle    = profiles2(iprof)%sunazangle
      profiles1(iprof)%latitude      = profiles2(iprof)%latitude
      profiles1(iprof)%longitude     = profiles2(iprof)%longitude
      profiles1(iprof)%ctp           = profiles2(iprof)%ctp
      profiles1(iprof)%cfraction     = profiles2(iprof)%cfraction
      profiles1(iprof)%clw_scheme    = profiles2(iprof)%clw_scheme
      profiles1(iprof)%clwde_param   = profiles2(iprof)%clwde_param
      profiles1(iprof)%ice_scheme    = profiles2(iprof)%ice_scheme
      profiles1(iprof)%icede_param   = profiles2(iprof)%icede_param
      profiles1(iprof)%elevation     = profiles2(iprof)%elevation
      profiles1(iprof)%s2m           = profiles2(iprof)%s2m
      IF (PRESENT(profiles_gas)) THEN
        profiles1(iprof)%s2m%q = profiles_gas(iprof)%s2m%q
        profiles1(iprof)%s2m%o = profiles_gas(iprof)%s2m%o
      ENDIF
      profiles1(iprof)%skin          = profiles2(iprof)%skin
      profiles1(iprof)%be            = profiles2(iprof)%be
      profiles1(iprof)%cosbk         = profiles2(iprof)%cosbk
    ENDIF
    IF (larray1) THEN
      profiles1(iprof)%nlevels = profiles2(iprof)%nlevels
      profiles1(iprof)%nlayers = profiles2(iprof)%nlayers
      profiles1(iprof)%p       = profiles2(iprof)%p
      profiles1(iprof)%t       = profiles2(iprof)%t
      IF (PRESENT(profiles_gas)) THEN
        profiles1(iprof)%gas_units     = profiles_gas(iprof)%gas_units
        profiles1(iprof)%q             = profiles_gas(iprof)%q
        IF (ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles_gas(iprof)%o3))    &
            profiles1(iprof)%o3       = profiles_gas(iprof)%o3
        IF (ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles_gas(iprof)%co2))  &
            profiles1(iprof)%co2      = profiles_gas(iprof)%co2
        IF (ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles_gas(iprof)%n2o))  &
            profiles1(iprof)%n2o      = profiles_gas(iprof)%n2o
        IF (ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles_gas(iprof)%co))    &
            profiles1(iprof)%co       = profiles_gas(iprof)%co
        IF (ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles_gas(iprof)%ch4))  &
            profiles1(iprof)%ch4      = profiles_gas(iprof)%ch4
        IF (ASSOCIATED(profiles1(iprof)%so2) .AND. ASSOCIATED(profiles_gas(iprof)%so2))  &
            profiles1(iprof)%so2      = profiles_gas(iprof)%so2
      ELSE
        profiles1(iprof)%gas_units     = profiles2(iprof)%gas_units
        profiles1(iprof)%q             = profiles2(iprof)%q
        IF (ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles2(iprof)%o3))       &
            profiles1(iprof)%o3       = profiles2(iprof)%o3
        IF (ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles2(iprof)%co2))     &
            profiles1(iprof)%co2      = profiles2(iprof)%co2
        IF (ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles2(iprof)%n2o))     &
            profiles1(iprof)%n2o      = profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles2(iprof)%co))       &
            profiles1(iprof)%co       = profiles2(iprof)%co
        IF (ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles2(iprof)%ch4))     &
            profiles1(iprof)%ch4      = profiles2(iprof)%ch4
        IF (ASSOCIATED(profiles1(iprof)%so2) .AND. ASSOCIATED(profiles2(iprof)%so2))     &
            profiles1(iprof)%so2      = profiles2(iprof)%so2
      ENDIF
      IF (ASSOCIATED(profiles1(iprof)%clw) .AND. ASSOCIATED(profiles2(iprof)%clw))           &
          profiles1(iprof)%clw      = profiles2(iprof)%clw
      IF (ASSOCIATED(profiles1(iprof)%aerosols) .AND. ASSOCIATED(profiles2(iprof)%aerosols)) &
          profiles1(iprof)%aerosols = profiles2(iprof)%aerosols
      IF (ASSOCIATED(profiles1(iprof)%cloud) .AND. ASSOCIATED(profiles2(iprof)%cloud))       &
          profiles1(iprof)%cloud    = profiles2(iprof)%cloud
      IF (ASSOCIATED(profiles1(iprof)%cfrac) .AND. ASSOCIATED(profiles2(iprof)%cfrac))       &
          profiles1(iprof)%cfrac    = profiles2(iprof)%cfrac
      IF (ASSOCIATED(profiles1(iprof)%clwde) .AND. ASSOCIATED(profiles2(iprof)%clwde))       &
          profiles1(iprof)%clwde    = profiles2(iprof)%clwde
      IF (ASSOCIATED(profiles1(iprof)%icede) .AND. ASSOCIATED(profiles2(iprof)%icede))       &
          profiles1(iprof)%icede    = profiles2(iprof)%icede
    ENDIF
  ENDDO
END SUBROUTINE 
