! Description:
!> @file
!!   K of nearest level calculations.
!
!> @brief
!!   K of nearest level calculations.
!!
!! @param[in]     do_ctp        switch for CTP nearest lev calculation
!! @param[in]     opts          options to configure the simulations
!! @param[in]     chanprof      channels and profiles to simulate
!! @param[in]     profiles      input profiles
!! @param[in,out] profiles_k    input profile increments
!! @param[in]     auxs          auxiliary near-level profile data structure
!! @param[in]     auxs_k        auxiliary near-level profile data increments
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
SUBROUTINE rttov_calc_nearest_lev_k(do_ctp, opts, chanprof, profiles, profiles_k, auxs, auxs_k)

  USE rttov_types, ONLY : rttov_options, rttov_chanprof, rttov_profile, rttov_profile_aux_s
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(jplm),             INTENT(IN)    :: do_ctp
  TYPE(rttov_options),       INTENT(IN)    :: opts
  TYPE(rttov_chanprof),      INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs_k(:)
!INTF_END

  INTEGER(jpim) :: nchanprof, i, prof
  REAL(jprb)    :: dp, dp_k
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_K', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)

  DO i = 1, nchanprof
    prof = chanprof(i)%prof

    ! ---------------------------------
    ! Nearest level below cloud top
    ! ---------------------------------
    IF (do_ctp) THEN
      dp = profiles(prof)%p(auxs(prof)%nearestlev_ctp) - &
           profiles(prof)%p(auxs(prof)%nearestlev_ctp - 1)
      dp_k = 0._jprb
      profiles_k(i)%cfraction = profiles_k(i)%cfraction + auxs_k(i)%cfraction
      profiles_k(i)%ctp       = profiles_k(i)%ctp - auxs_k(i)%pfraction_ctp / dp
      IF (opts%interpolation%lgradp) THEN
        profiles_k(i)%p(auxs(prof)%nearestlev_ctp) = profiles_k(i)%p(auxs(prof)%nearestlev_ctp) + &
            auxs_k(i)%pfraction_ctp / dp
        dp_k = dp_k - (profiles(prof)%p(auxs(prof)%nearestlev_ctp) - profiles(prof)%ctp) / dp ** 2 * &
            auxs_k(i)%pfraction_ctp
        profiles_k(i)%p(auxs(prof)%nearestlev_ctp) = &
            profiles_k(i)%p(auxs(prof)%nearestlev_ctp) + dp_k
        profiles_k(i)%p(auxs(prof)%nearestlev_ctp - 1) = &
            profiles_k(i)%p(auxs(prof)%nearestlev_ctp - 1) - dp_k
      ENDIF
    ENDIF

    ! ---------------------------------
    ! Nearest level below surface
    ! ---------------------------------
    dp = profiles(prof)%p(auxs(prof)%nearestlev_surf) - &
         profiles(prof)%p(auxs(prof)%nearestlev_surf - 1)
    dp_k = 0._jprb
    profiles_k(i)%s2m%p = profiles_k(i)%s2m%p - auxs_k(i)%pfraction_surf / dp
    IF (opts%interpolation%lgradp) THEN
      profiles_k(i)%p(auxs(prof)%nearestlev_surf) = profiles_k(i)%p(auxs(prof)%nearestlev_surf) + &
          auxs_k(i)%pfraction_surf / dp
      dp_k = dp_k - (profiles(prof)%p(auxs(prof)%nearestlev_surf) - profiles(prof)%s2m%p) / dp ** 2 * &
          auxs_k(i)%pfraction_surf
      profiles_k(i)%p(auxs(prof)%nearestlev_surf) = &
          profiles_k(i)%p(auxs(prof)%nearestlev_surf) + dp_k
      profiles_k(i)%p(auxs(prof)%nearestlev_surf - 1) = &
          profiles_k(i)%p(auxs(prof)%nearestlev_surf - 1) - dp_k
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calc_nearest_lev_k
