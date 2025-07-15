! Description:
!> @file
!!   Add two profile structures.
!
!> @brief
!!   Add two profile structures.
!!
!! @details
!!   This subroutine is the AD/K counterpart to rttov_copy_prof.
!!   The optional profiles_gas argument fulfils a similar role and, if present,
!!   the gas data in the profiles2 argument are accumulated in profiles_gas
!!   instead of being added to profiles1 in the output profiles argument.
!!
!!   This subroutine is also used in the test suite to add entire profile
!!   structures and in this case the optional arguments are not used.
!!
!! @param[in,out] profiles       output summed profiles structure
!! @param[in]     profiles1      first input profiles structure
!! @param[in]     profiles2      second input profiles structure
!! @param[in]     larray         flag to indicate array members should be added, optional
!! @param[in]     lscalar        flag to indicate scalar members should be added, optional
!! @param[in,out] profiles_gas   profiles structure in which to accumulate gas profile data from profile2, optional
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
SUBROUTINE rttov_add_prof( &
              profiles,     &
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
  TYPE(rttov_profile), INTENT(INOUT)           :: profiles(:)
  TYPE(rttov_profile), INTENT(IN)              :: profiles1(SIZE(profiles))
  TYPE(rttov_profile), INTENT(IN)              :: profiles2(SIZE(profiles))
  LOGICAL(KIND=jplm),  INTENT(IN),    OPTIONAL :: larray
  LOGICAL(KIND=jplm),  INTENT(IN),    OPTIONAL :: lscalar
  TYPE(rttov_profile), INTENT(INOUT), OPTIONAL :: profiles_gas(SIZE(profiles))
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
! Do not add:
!   profiles%id %date %time %gas_units
!
  DO iprof = 1, nprofiles
    IF (lscalar1) THEN
      profiles(iprof)%zenangle            = profiles1(iprof)%zenangle + profiles2(iprof)%zenangle
      profiles(iprof)%azangle             = profiles1(iprof)%azangle + profiles2(iprof)%azangle
      profiles(iprof)%sunzenangle         = profiles1(iprof)%sunzenangle + profiles2(iprof)%sunzenangle
      profiles(iprof)%sunazangle          = profiles1(iprof)%sunazangle + profiles2(iprof)%sunazangle
      profiles(iprof)%latitude            = profiles1(iprof)%latitude + profiles2(iprof)%latitude
      profiles(iprof)%longitude           = profiles1(iprof)%longitude + profiles2(iprof)%longitude
      profiles(iprof)%ctp                 = profiles1(iprof)%ctp + profiles2(iprof)%ctp
      profiles(iprof)%cfraction           = profiles1(iprof)%cfraction + profiles2(iprof)%cfraction
      profiles(iprof)%elevation           = profiles1(iprof)%elevation + profiles2(iprof)%elevation
      profiles(iprof)%s2m%t               = profiles1(iprof)%s2m%t + profiles2(iprof)%s2m%t
      IF (PRESENT(profiles_gas)) THEN
        profiles_gas(iprof)%s2m%q         = profiles_gas(iprof)%s2m%q + profiles2(iprof)%s2m%q
        profiles_gas(iprof)%s2m%o         = profiles_gas(iprof)%s2m%o + profiles2(iprof)%s2m%o
      ELSE
        profiles(iprof)%s2m%q             = profiles1(iprof)%s2m%q + profiles2(iprof)%s2m%q
        profiles(iprof)%s2m%o             = profiles1(iprof)%s2m%o + profiles2(iprof)%s2m%o
      ENDIF
      profiles(iprof)%s2m%p               = profiles1(iprof)%s2m%p + profiles2(iprof)%s2m%p
      profiles(iprof)%s2m%u               = profiles1(iprof)%s2m%u + profiles2(iprof)%s2m%u
      profiles(iprof)%s2m%v               = profiles1(iprof)%s2m%v + profiles2(iprof)%s2m%v
      profiles(iprof)%s2m%wfetc           = profiles1(iprof)%s2m%wfetc + profiles2(iprof)%s2m%wfetc
      profiles(iprof)%skin%t              = profiles1(iprof)%skin%t + profiles2(iprof)%skin%t
      profiles(iprof)%skin%fastem(:)      = profiles1(iprof)%skin%fastem(:) + profiles2(iprof)%skin%fastem(:)
      profiles(iprof)%skin%salinity       = profiles1(iprof)%skin%salinity + profiles2(iprof)%skin%salinity
      profiles(iprof)%skin%foam_fraction  = profiles1(iprof)%skin%foam_fraction + &
                                            profiles2(iprof)%skin%foam_fraction
      profiles(iprof)%skin%snow_fraction  = profiles1(iprof)%skin%snow_fraction + &
                                            profiles2(iprof)%skin%snow_fraction
      profiles(iprof)%skin%soil_moisture  = profiles1(iprof)%skin%soil_moisture + &
                                            profiles2(iprof)%skin%soil_moisture
      profiles(iprof)%be                  = profiles1(iprof)%be + profiles2(iprof)%be
      profiles(iprof)%cosbk               = profiles1(iprof)%cosbk + profiles2(iprof)%cosbk
    ENDIF
    IF (larray1) THEN
      profiles(iprof)%p = profiles1(iprof)%p + profiles2(iprof)%p
      profiles(iprof)%t = profiles1(iprof)%t + profiles2(iprof)%t

      IF (PRESENT(profiles_gas)) THEN
        profiles_gas(iprof)%q            = profiles_gas(iprof)%q + profiles2(iprof)%q
        IF (ASSOCIATED(profiles2(iprof)%o3) .AND. ASSOCIATED(profiles_gas(iprof)%o3))   &
            profiles_gas(iprof)%o3       = profiles_gas(iprof)%o3 + profiles2(iprof)%o3
        IF (ASSOCIATED(profiles2(iprof)%co2) .AND. ASSOCIATED(profiles_gas(iprof)%co2)) &
            profiles_gas(iprof)%co2      = profiles_gas(iprof)%co2 + profiles2(iprof)%co2
        IF (ASSOCIATED(profiles2(iprof)%n2o) .AND. ASSOCIATED(profiles_gas(iprof)%n2o)) &
            profiles_gas(iprof)%n2o      = profiles_gas(iprof)%n2o + profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles2(iprof)%co) .AND. ASSOCIATED(profiles_gas(iprof)%co))   &
            profiles_gas(iprof)%co       = profiles_gas(iprof)%co + profiles2(iprof)%co
        IF (ASSOCIATED(profiles2(iprof)%ch4) .AND. ASSOCIATED(profiles_gas(iprof)%ch4)) &
            profiles_gas(iprof)%ch4      = profiles_gas(iprof)%ch4 + profiles2(iprof)%ch4
        IF (ASSOCIATED(profiles2(iprof)%so2) .AND. ASSOCIATED(profiles_gas(iprof)%so2)) &
            profiles_gas(iprof)%so2      = profiles_gas(iprof)%so2 + profiles2(iprof)%so2
      ELSE
        profiles(iprof)%q            = profiles1(iprof)%q + profiles2(iprof)%q
        IF (ASSOCIATED(profiles(iprof)%o3)  .AND. ASSOCIATED(profiles1(iprof)%o3) .AND. ASSOCIATED(profiles2(iprof)%o3))   &
            profiles(iprof)%o3       = profiles1(iprof)%o3 + profiles2(iprof)%o3
        IF (ASSOCIATED(profiles(iprof)%co2) .AND. ASSOCIATED(profiles1(iprof)%co2) .AND. ASSOCIATED(profiles2(iprof)%co2)) &
            profiles(iprof)%co2      = profiles1(iprof)%co2 + profiles2(iprof)%co2
        IF (ASSOCIATED(profiles(iprof)%n2o) .AND. ASSOCIATED(profiles1(iprof)%n2o) .AND. ASSOCIATED(profiles2(iprof)%n2o)) &
            profiles(iprof)%n2o      = profiles1(iprof)%n2o + profiles2(iprof)%n2o
        IF (ASSOCIATED(profiles(iprof)%co)  .AND. ASSOCIATED(profiles1(iprof)%co) .AND. ASSOCIATED(profiles2(iprof)%co))   &
            profiles(iprof)%co       = profiles1(iprof)%co + profiles2(iprof)%co
        IF (ASSOCIATED(profiles(iprof)%ch4) .AND. ASSOCIATED(profiles1(iprof)%ch4) .AND. ASSOCIATED(profiles2(iprof)%ch4)) &
            profiles(iprof)%ch4      = profiles1(iprof)%ch4 + profiles2(iprof)%ch4
        IF (ASSOCIATED(profiles(iprof)%so2) .AND. ASSOCIATED(profiles1(iprof)%so2) .AND. ASSOCIATED(profiles2(iprof)%so2)) &
            profiles(iprof)%so2      = profiles1(iprof)%so2 + profiles2(iprof)%so2
      ENDIF

      IF (ASSOCIATED(profiles(iprof)%clw) .AND. ASSOCIATED(profiles1(iprof)%clw) .AND. ASSOCIATED(profiles2(iprof)%clw))   &
          profiles(iprof)%clw       = profiles1(iprof)%clw + profiles2(iprof)%clw
      IF (ASSOCIATED(profiles(iprof)%aerosols) .AND. ASSOCIATED(profiles1(iprof)%aerosols) .AND. &
          ASSOCIATED(profiles2(iprof)%aerosols)) &
           profiles(iprof)%aerosols = profiles1(iprof)%aerosols + profiles2(iprof)%aerosols
      IF (ASSOCIATED(profiles(iprof)%cloud) .AND. ASSOCIATED(profiles1(iprof)%cloud) .AND. &
          ASSOCIATED(profiles2(iprof)%cloud)) &
           profiles(iprof)%cloud    = profiles1(iprof)%cloud + profiles2(iprof)%cloud
      IF (ASSOCIATED(profiles(iprof)%cfrac) .AND. ASSOCIATED(profiles1(iprof)%cfrac) .AND. &
          ASSOCIATED(profiles2(iprof)%cfrac)) &
           profiles(iprof)%cfrac    = profiles1(iprof)%cfrac + profiles2(iprof)%cfrac
      IF (ASSOCIATED(profiles(iprof)%clwde) .AND. ASSOCIATED(profiles1(iprof)%clwde) .AND. &
          ASSOCIATED(profiles2(iprof)%clwde)) &
           profiles(iprof)%clwde    = profiles1(iprof)%clwde + profiles2(iprof)%clwde
      IF (ASSOCIATED(profiles(iprof)%icede) .AND. ASSOCIATED(profiles1(iprof)%icede) .AND. &
          ASSOCIATED(profiles2(iprof)%icede)) &
           profiles(iprof)%icede    = profiles1(iprof)%icede + profiles2(iprof)%icede
    ENDIF
  ENDDO
END SUBROUTINE 
