! Description:
!> @file
!!   AD of nearest level calculations.
!
!> @brief
!!   AD of nearest level calculations.
!!
!! @param[in]     do_ctp        switch for CTP nearest lev calculation
!! @param[in]     opts          options to configure the simulations
!! @param[in]     profiles      input profiles
!! @param[in,out] profiles_ad   input profile increments
!! @param[in]     auxs          auxiliary near-level profile data structure
!! @param[in]     auxs_ad       auxiliary near-level profile data increments
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_calc_nearest_lev_ad(do_ctp, opts, profiles, profiles_ad, auxs, auxs_ad)

  USE rttov_types, ONLY : rttov_options, rttov_profile, rttov_profile_aux_s
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(jplm),             INTENT(IN)    :: do_ctp
  TYPE(rttov_options),       INTENT(IN)    :: opts
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_ad(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs_ad(:)
!INTF_END

  INTEGER(jpim) :: nprofiles, iprof
  REAL(jprb)    :: dp, dp_ad
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_AD', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)

  DO iprof = 1, nprofiles

    ! ---------------------------------
    ! Nearest level below cloud top
    ! ---------------------------------
    IF (do_ctp) THEN
      dp = profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - &
           profiles(iprof)%p(auxs(iprof)%nearestlev_ctp - 1)
      dp_ad = 0._jprb
      profiles_ad(iprof)%cfraction = profiles_ad(iprof)%cfraction + auxs_ad(iprof)%cfraction
      profiles_ad(iprof)%ctp       = profiles_ad(iprof)%ctp - auxs_ad(iprof)%pfraction_ctp / dp
      IF (opts%interpolation%lgradp) THEN
        profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp) = profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp) + &
            auxs_ad(iprof)%pfraction_ctp / dp
        dp_ad = dp_ad - (profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - profiles(iprof)%ctp) / dp ** 2 * &
            auxs_ad(iprof)%pfraction_ctp
        profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp) = &
            profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp) + dp_ad
        profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp - 1) = &
            profiles_ad(iprof)%p(auxs(iprof)%nearestlev_ctp - 1) - dp_ad
      ENDIF
    ENDIF

    ! ---------------------------------
    ! Nearest level below surface
    ! ---------------------------------
    dp = profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - &
         profiles(iprof)%p(auxs(iprof)%nearestlev_surf - 1)
    dp_ad = 0._jprb
    profiles_ad(iprof)%s2m%p = profiles_ad(iprof)%s2m%p - auxs_ad(iprof)%pfraction_surf / dp
    IF (opts%interpolation%lgradp) THEN
      profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf) = profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf) + &
          auxs_ad(iprof)%pfraction_surf / dp
      dp_ad = dp_ad - (profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - profiles(iprof)%s2m%p) / dp ** 2 * &
          auxs_ad(iprof)%pfraction_surf
      profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf) = &
          profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf) + dp_ad
      profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf - 1) = &
          profiles_ad(iprof)%p(auxs(iprof)%nearestlev_surf - 1) - dp_ad
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calc_nearest_lev_ad
