! Description:
!> @file
!!   Returns BRDF values for the given chanprof array.
!
!> @brief
!!   Returns BRDF values for the given chanprof array.
!!   Optionally, will return bi-hemispherical (back-sky) albedo.
!!
!! @details
!!   The atlas data structure must have been initialised by calling
!!   rttov_setup_brdf_atlas before calling this subroutine.
!!
!!   If you initialised the atlas data with the "single_instrument" option
!!   you must only call this subroutine for that specific instrument. If
!!   "single_instrument" was false then you use the atlas data with any
!!   visible/IR coefficient structure.
!!
!!   Surface types (as determined by profiles(:)\%skin\%surftype):
!!     If the surftype is seaice no BRDFs are returned;
!!     for sea the atlas returns ocean or fresh water reflectance
!!     values but this does not account for sunglint or wind-
!!     roughened water; for land surfaces the atlas returns reflectances
!!     where it has data and this includes climatological snow.
!!
!!   After this subroutine has been called you should check the
!!   returned BRDF values: if the atlas has no value for the given
!!   location it will return values less than or equal to zero.
!!
!!   The important elements of the profiles structure for the BRDF atlas are:
!!   skin\%surftype, skin\%watertype, latitude, longitude,
!!   zenangle, azangle, sunzenangle, sunazangle.
!!
!! @param[out]  err          status on exit
!! @param[in]   opts         options to configure the simulations
!! @param[in]   chanprof     specifies channels and profiles to simulate
!! @param[in]   profiles     input atmospheric profiles and surface variables
!! @param[in]   coefs        coefficients structure for instrument to simulate
!! @param[in]   atlas        atlas data
!! @param[out]  brdf         BRDF values
!! @param[out]  brdf_flag    BRDF atlas flags, optional
!! @param[out]  bh_albedo    Bi-hemispherical (back-sky) albedo, optional
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
SUBROUTINE rttov_get_brdf(  &
                err,        &! out
                opts,       &! in
                chanprof,   &! in
                profiles,   &! in
                coefs,      &! in
                atlas,      &! in
                brdf,       &! out
                brdf_flag,  &! out, optional
                bh_albedo)   ! out, optional
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : &
          rttov_chanprof, &
          rttov_profile,  &
          rttov_coefs,    &
          rttov_options
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data
!INTF_OFF
  USE rttov_const, ONLY :      &
          surftype_seaice,     &
          errorstatus_success, &
          errorstatus_fatal

  USE mod_rttov_brdf_atlas, ONLY : brdf_atlas_id

  USE mod_brdf_atlas, ONLY : rttov_visnirbrdf

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),          INTENT(OUT)           :: err
  TYPE(rttov_options),         INTENT(IN)            :: opts
  TYPE(rttov_chanprof),        INTENT(IN)            :: chanprof(:)
  TYPE(rttov_profile),         INTENT(IN)            :: profiles(:)
  TYPE(rttov_coefs),           INTENT(IN)            :: coefs
  TYPE(rttov_brdf_atlas_data), INTENT(IN)            :: atlas
  REAL(KIND=jprb),             INTENT(OUT)           :: brdf(SIZE(chanprof))
  INTEGER(KIND=jpim),          INTENT(OUT), OPTIONAL :: brdf_flag(SIZE(chanprof))
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: bh_albedo(SIZE(chanprof))
!INTF_END

#include "rttov_errorreport.interface"

  REAL(KIND=jprb)    :: ZHOOK_HANDLE

  INTEGER(KIND=jpim) :: nchanprof                 ! number of channels/profiles
  INTEGER(KIND=jpim) :: prof                      ! index of current profile
  INTEGER(KIND=jpim) :: nchan                     ! number of channels in current profile
  INTEGER(KIND=jpim) :: chans(coefs%coef%fmv_chn) ! channel indexes for current profile

  INTEGER(KIND=jpim) :: k                         ! loop indexes
  INTEGER(KIND=jpim) :: lo, hi                    ! low/high chanprof indexes for current profile

  INTEGER(KIND=jpim) :: instr_brdf_flag

!----------------------------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF', 0_jpim, ZHOOK_HANDLE)

  !-----------------------------
  ! Initialise output arguments
  !-----------------------------

  brdf = -1._jprb
  IF (PRESENT(brdf_flag)) brdf_flag(:) = 0_jpim
  IF (PRESENT(bh_albedo)) bh_albedo(:) = 0._jprb

  IF (.NOT. atlas%init) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Atlas data not initialised')
  ENDIF

  nchanprof = SIZE(chanprof(:))

  ! VIS/NIR sensor

  IF (atlas%atlas_id == brdf_atlas_id) THEN
    ! Nothing to do
  ELSE
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Unknown BRDF atlas ID')
  ENDIF

  !-----------
  ! Main loop
  !-----------

  ! Assume chanprof(:) lists all channels for profile 1, followed by all channels
  ! for profile 2, and so on. Loop through chanprof(:) until the profile number
  ! changes and in this way obtain the list of channels for the current profile.
  k = 1
  DO WHILE (k <= nchanprof)
    prof = chanprof(k)%prof
    lo = k
    k  = lo+1
    DO
      IF (k > nchanprof) EXIT
      IF (prof /= chanprof(k)%prof) EXIT
      k = k+1
    ENDDO
    hi = k-1
    nchan = k-lo
    chans(1:nchan) = chanprof(lo:hi)%chan

    ! Profile index is 'prof'
    ! Number of channels for this profile is 'nchan'
    ! Channel indexes are 'chans(1:nchan)'

    !---------------------------
    ! VIS/NIR atlas
    !---------------------------
    IF (atlas%atlas_id == brdf_atlas_id) THEN

      IF (profiles(prof)%skin%surftype == surftype_seaice) THEN
        brdf(lo:hi) = -999.0_jprb
      ELSE
        IF (PRESENT(bh_albedo)) THEN
          CALL rttov_visnirbrdf(                        &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%watertype,      &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%azangle,             &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%sunazangle,          &! in
                    coefs%coef%ff_cwn(chans(1:nchan)),  &! in
                    chans(1:nchan),                     &! in
                    coefs%coef%id_platform,             &! in
                    coefs%coef%id_sat,                  &! in
                    coefs%coef%id_inst,                 &! in
                    coefs%coef%fmv_chn,                 &! in
                    atlas%brdf_atlas,                   &! in
                    brdf(lo:hi),                        &! out
                    instr_brdf_flag,                    &! out
                    bh_albedo(lo:hi))                    ! out
        ELSE
          CALL rttov_visnirbrdf(                        &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%watertype,      &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%azangle,             &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%sunazangle,          &! in
                    coefs%coef%ff_cwn(chans(1:nchan)),  &! in
                    chans(1:nchan),                     &! in
                    coefs%coef%id_platform,             &! in
                    coefs%coef%id_sat,                  &! in
                    coefs%coef%id_inst,                 &! in
                    coefs%coef%fmv_chn,                 &! in
                    atlas%brdf_atlas,                   &! in
                    brdf(lo:hi),                        &! out
                    instr_brdf_flag)                     ! out
        ENDIF
        THROWM(err .NE. 0, "error in BRDF atlas")

        IF (PRESENT(brdf_flag)) brdf_flag(lo:hi) = instr_brdf_flag

      ENDIF

    ENDIF

  ENDDO ! chanprof

  WHERE (brdf <= 0._jprb)
    brdf = -1._jprb
  ENDWHERE

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_BRDF',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_get_brdf
