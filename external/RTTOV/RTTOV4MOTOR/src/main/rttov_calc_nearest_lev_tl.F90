! Description:
!> @file
!!   TL of nearest level calculations.
!
!> @brief
!!   TL of nearest level calculations.
!!
!! @param[in]     do_ctp        switch for CTP nearest lev calculation
!! @param[in]     opts          options to configure the simulations
!! @param[in]     profiles      input profiles
!! @param[in]     profiles_tl   input profile perturbations
!! @param[in]     auxs          auxiliary near-level profile data structure
!! @param[in,out] auxs_tl       auxiliary near-level profile data perturbations
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
SUBROUTINE rttov_calc_nearest_lev_tl(do_ctp, opts, profiles, profiles_tl, auxs, auxs_tl)

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
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_profile_aux_s), INTENT(INOUT) :: auxs_tl(:)
!INTF_END

  INTEGER(jpim) :: nprofiles, iprof
  REAL(jprb)    :: dp, dp_tl
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_TL', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)

  DO iprof = 1, nprofiles

    ! ---------------------------------
    ! Nearest level below surface
    ! ---------------------------------
    dp = profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - &
         profiles(iprof)%p(auxs(iprof)%nearestlev_surf - 1)
    IF (opts%interpolation%lgradp) THEN
      dp_tl = profiles_tl(iprof)%p(auxs(iprof)%nearestlev_surf) - &
              profiles_tl(iprof)%p(auxs(iprof)%nearestlev_surf - 1)
      auxs_tl(iprof)%pfraction_surf = &
          (profiles_tl(iprof)%p(auxs(iprof)%nearestlev_surf) - profiles_tl(iprof)%s2m%p) / dp -  &
          (profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - profiles(iprof)%s2m%p) * dp_tl / (dp ** 2)
    ELSE
      auxs_tl(iprof)%pfraction_surf =  - profiles_tl(iprof)%s2m%p / dp
    ENDIF

    ! ---------------------------------
    ! Nearest level below cloud top
    ! ---------------------------------
    IF (do_ctp) THEN
      dp = profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - &
           profiles(iprof)%p(auxs(iprof)%nearestlev_ctp - 1)
      IF (opts%interpolation%lgradp) THEN
        dp_tl = profiles_tl(iprof)%p(auxs(iprof)%nearestlev_ctp) - &
                profiles_tl(iprof)%p(auxs(iprof)%nearestlev_ctp - 1)
        auxs_tl(iprof)%pfraction_ctp = &
            (profiles_tl(iprof)%p(auxs(iprof)%nearestlev_ctp) - profiles_tl(iprof)%ctp) / dp -  &
            (profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - profiles(iprof)%ctp) * dp_tl / (dp ** 2)
      ELSE
        auxs_tl(iprof)%pfraction_ctp =  - profiles_tl(iprof)%ctp / dp
      ENDIF
      auxs_tl(iprof)%cfraction = profiles_tl(iprof)%cfraction
    ELSE
      auxs_tl(iprof)%pfraction_ctp = 0._jprb
      auxs_tl(iprof)%cfraction     = 0._jprb
    ENDIF

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calc_nearest_lev_tl
