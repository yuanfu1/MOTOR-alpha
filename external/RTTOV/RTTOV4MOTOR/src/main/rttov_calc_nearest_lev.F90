! Description:
!> @file
!!   Calculate nearest levels to surface and input cloud top pressure.
!
!> @brief
!!   Calculate nearest levels to surface and input cloud top pressure.
!!
!> @details
!!   This is called separately for user levels and coef levels.
!!
!!   For the surface this finds the level at or immediately below
!!   the 2m pressure, unless this lies below the lowest profile level
!!   in which case it returns the bottom level.
!!
!!   For the CTP this finds the level at or immediately below the
!!   CTP. This is not required on coef levels and the simple cloud
!!   scheme is not applicable to MW sensors, so this calculation can
!!   be ommitted by setting the do_ctp argument to false.
!!
!!   For both surface and CTP the fraction of the layer below the
!!   2m/cloud top pressure is also calculated. This is negative when
!!   the surface pressure lies below the bottom profile level.
!!
!! @param[in]     do_ctp        switch for CTP nearest lev calculation
!! @param[in]     profiles      input profiles
!! @param[in,out] auxs          auxiliary near-level profile data structure
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
SUBROUTINE rttov_calc_nearest_lev(do_ctp, profiles, auxs)

  USE rttov_types, ONLY : rttov_profile, rttov_profile_aux_s
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(jplm),             INTENT(IN)    :: do_ctp
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile_aux_s), INTENT(INOUT) :: auxs(:)
!INTF_END

  INTEGER(jpim) :: nprofiles, nlevels, iprof, lev
  REAL(jprb)    :: dp
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlevels = profiles(1)%nlevels

  DO iprof = 1, nprofiles

    ! ---------------------------------
    ! Nearest level below surface
    ! ---------------------------------
    DO lev = nlevels - 1, 1, -1
      IF (profiles(iprof)%s2m%p > profiles(iprof)%p(lev)) EXIT
    ENDDO
    ! case-1: surf lies above lev=nlevels
    !         at exit, lev is first level above surface
    ! case-2: surf lies below lev=nlevels
    !         at exit, lev+1=nlevels, there is no level below surface
    ! case-1: nearestlev_surf is set to first level below surface
    !         auxs(iprof)%pfraction_surf is +ve
    ! case-2: nearestlev_surf is set to first level above surface
    !         auxs(iprof)%pfraction_surf is -ve
    auxs(iprof)%nearestlev_surf = lev + 1
    dp = profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - &
         profiles(iprof)%p(auxs(iprof)%nearestlev_surf - 1)
    auxs(iprof)%pfraction_surf  = &
        (profiles(iprof)%p(auxs(iprof)%nearestlev_surf) - profiles(iprof)%s2m%p) / dp

    ! ---------------------------------
    ! Nearest level below cloud top
    ! ---------------------------------
    auxs(iprof)%nearestlev_ctp = nlevels - 1
    auxs(iprof)%pfraction_ctp  = 0._jprb
    auxs(iprof)%cfraction      = 0._jprb
    IF (do_ctp) THEN
      DO lev = nlevels - 1, 1, -1
        IF (profiles(iprof)%ctp > profiles(iprof)%p(lev)) EXIT
      ENDDO
      IF (lev > 1) THEN
        auxs(iprof)%nearestlev_ctp = lev + 1
        dp = profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - &
             profiles(iprof)%p(auxs(iprof)%nearestlev_ctp - 1)
        auxs(iprof)%pfraction_ctp  = &
            (profiles(iprof)%p(auxs(iprof)%nearestlev_ctp) - profiles(iprof)%ctp) / dp
        auxs(iprof)%cfraction      = profiles(iprof)%cfraction
      ENDIF
    ENDIF

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_CALC_NEAREST_LEV', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calc_nearest_lev
