! Description:
!> @file
!!   Handle allocation/deallocation of structure(s) for internal RTTOV state.
!
!> @brief
!!   Handle allocation/deallocation of structure(s) for internal RTTOV state.
!!
!! @details
!!   The rttov_traj structure contains internal state for RTTOV.
!!   This can optionally be allocated/deallocated by the user outside of
!!   RTTOV which can sometimes bring performance benefits when making
!!   repeated calls to RTTOV. If it is not passed into RTTOV then it
!!   is allocated/deallocated internally.
!!
!!   This subroutine is called internally by RTTOV. At the end of the
!!   allocation call the traj0 variable (direct, TL, AD, K versions) will
!!   point to the trajectory structure that is to be used.
!!
!!   If the trajectory structure (direct, TL, AD or K version) is allocated
!!   externally then the traj2 argument will be present. In this case traj0
!!   points to traj2 and no allocations are performed. Some checks are carried
!!   out to ensure traj2 is consistent with the current call to RTTOV.
!!
!!   If the trajectory structure is not allocated externally then the traj1
!!   variable is allocated by a call to rttov_alloc_traj and traj0 points
!!   to this.
!!
!!   For deallocation calls nothing is done in the first case (traj2 present).
!!   In the second case traj1 is deallocated.
!!
!! @param[out]    err            status on exit
!! @param[in]     nprofiles      number of profiles being simulated
!! @param[in]     nchanprof      total number of channels being simulated
!! @param[in]     opts           options to configure the simulations
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param         traj0          output pointer to rttov_traj structure for the direct model, optional
!! @param         traj0_tl       output pointer to rttov_traj structure for the TL model, optional
!! @param         traj0_ad       output pointer to rttov_traj structure for the AD model, optional
!! @param         traj0_k        output pointer to rttov_traj structure for the K model, optional
!! @param[in,out] traj1          internally allocated rttov_traj structure for the direct model, optional
!! @param[in,out] traj1_tl       internally allocated rttov_traj structure for the TL model, optional
!! @param[in,out] traj1_ad       internally allocated rttov_traj structure for the AD model, optional
!! @param[in,out] traj1_k        internally allocated rttov_traj structure for the K model, optional
!! @param[in]     traj2          externally allocated rttov_traj structure for the direct model, optional
!! @param[in]     traj2_tl       externally allocated rttov_traj structure for the TL model, optional
!! @param[in]     traj2_ad       externally allocated rttov_traj structure for the AD model, optional
!! @param[in]     traj2_k        externally allocated rttov_traj structure for the K model, optional
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
SUBROUTINE rttov_check_traj( &
       err,             &
       nprofiles,       &
       nchanprof,       &
       opts,            &
       nlevels,         &
       coefs,           &
       asw,             &
       traj0,    traj1,    traj2,    &
       traj0_tl, traj1_tl, traj2_tl, &
       traj0_ad, traj1_ad, traj2_ad, &
       traj0_k,  traj1_k,  traj2_k)

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coefs, rttov_options, rttov_traj

  IMPLICIT NONE

  INTEGER(KIND=jpim),                     INTENT(OUT)   :: err
  INTEGER(KIND=jpim),                     INTENT(IN)    :: nprofiles
  INTEGER(KIND=jpim),                     INTENT(IN)    :: nchanprof
  TYPE(rttov_options),                    INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                     INTENT(IN)    :: nlevels
  TYPE(rttov_coefs),   TARGET,            INTENT(IN)    :: coefs
  INTEGER(KIND=jpim),                     INTENT(IN)    :: asw
  TYPE(rttov_traj),    POINTER, OPTIONAL                :: traj0, traj0_tl, traj0_ad, traj0_k
  TYPE(rttov_traj),    TARGET,  OPTIONAL, INTENT(INOUT) :: traj1, traj1_tl, traj1_ad, traj1_k
  TYPE(rttov_traj),    TARGET,  OPTIONAL, INTENT(IN)    :: traj2, traj2_tl, traj2_ad, traj2_k
!INTF_END

#include "rttov_errorreport.interface"

  REAL(KIND=jprb) :: ZHOOK_HANDLE

  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',0_jpim,ZHOOK_HANDLE)

  IF (PRESENT(traj0)) THEN
    CALL check_traj(        &
        err,                &
        nprofiles,          &
        nchanprof,          &
        opts,               &
        nlevels,            &
        coefs,              &
        asw,                &
        traj0,              &
        traj1,              &
        traj2)
    THROWM(err.NE.0, "check traj0")
  ENDIF

  IF (PRESENT(traj0_tl)) THEN
    CALL check_traj(        &
        err,                &
        nprofiles,          &
        nchanprof,          &
        opts,               &
        nlevels,            &
        coefs,              &
        asw,                &
        traj0_tl,           &
        traj1_tl,           &
        traj2_tl)
    THROWM(err.NE.0, "check traj0_tl")
  ENDIF

  IF (PRESENT(traj0_ad)) THEN
    CALL check_traj(        &
        err,                &
        nprofiles,          &
        nchanprof,          &
        opts,               &
        nlevels,            &
        coefs,              &
        asw,                &
        traj0_ad,           &
        traj1_ad,           &
        traj2_ad)
    THROWM(err.NE.0, "check traj0_ad")
  ENDIF

  IF (PRESENT(traj0_k)) THEN
    CALL check_traj(        &
        err,                &
        nchanprof,          &
        nchanprof,          &
        opts,               &
        nlevels,            &
        coefs,              &
        asw,                &
        traj0_k,            &
        traj1_k,            &
        traj2_k)
     THROWM(err.NE.0, "check traj0_k")
   ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',1_jpim,ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_TEMP',1_jpim,ZHOOK_HANDLE)

  CONTAINS

  SUBROUTINE check_traj(  &
       err,               &
       nprofiles,         &
       nchanprof,         &
       opts,              &
       nlevels,           &
       coefs,             &
       asw,               &
       traj0,             &
       traj1,             &
       traj2)

    INTEGER(KIND=jpim),                     INTENT(OUT)   :: err
    INTEGER(KIND=jpim),                     INTENT(IN)    :: nprofiles
    INTEGER(KIND=jpim),                     INTENT(IN)    :: nchanprof
    TYPE(rttov_options),                    INTENT(IN)    :: opts
    INTEGER(KIND=jpim),                     INTENT(IN)    :: nlevels
    TYPE(rttov_coefs),   TARGET,            INTENT(IN)    :: coefs
    INTEGER(KIND=jpim),                     INTENT(IN)    :: asw
    TYPE(rttov_traj),    POINTER                          :: traj0
    TYPE(rttov_traj),    TARGET,  OPTIONAL, INTENT(INOUT) :: traj1
    TYPE(rttov_traj),    TARGET,  OPTIONAL, INTENT(IN)    :: traj2

#include "rttov_alloc_traj.interface"
#include "rttov_opts_eq.interface"

    REAL(KIND=jprb) :: ZHOOK_HANDLE

    TRY
    IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',0_jpim,ZHOOK_HANDLE)

    IF (asw == 1_jpim) THEN
      IF (PRESENT(traj2)) THEN
        ! If trajectory was externally allocated ensure that dimensions
        ! and other variables are consistent with the current call to RTTOV
        IF ((.NOT. ASSOCIATED(traj2%coefs, coefs))  .OR. &
            (traj2%nchanprof .NE. nchanprof)        .OR. &
            (.NOT. rttov_opts_eq(opts, traj2%opts)) .OR. &
            (traj2%nlevels .NE. nlevels)) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "rttov_check_traj fatal error dimensions mismatch")
        ENDIF
        traj0 => traj2
      ELSE
        ! If trajectory was not externally allocated we allocate it now
        CALL rttov_alloc_traj( &
            err,               &
            nprofiles,         &
            nchanprof,         &
            opts,              &
            nlevels,           &
            coefs,             &
            asw,               &
            traj1)
        THROWM(err.NE.0, "rttov_alloc_traj allocation fatal error")
        traj0 => traj1
      ENDIF
    ELSE
      ! If trajectory was externally allocated we do not deallocate anything
      IF (.NOT. PRESENT(traj2)) THEN
        ! If trajectory was not externally allocated we deallocate it now
        CALL rttov_alloc_traj( &
            err,               &
            nprofiles,         &
            nchanprof,         &
            opts,              &
            nlevels,           &
            coefs,             &
            asw,               &
            traj1)
        THROWM(err.NE.0, "rttov_alloc_traj deallocation fatal error")
      ENDIF
      NULLIFY(traj0)
    ENDIF

    IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',1_jpim,ZHOOK_HANDLE)
    CATCH
    IF (LHOOK) CALL DR_HOOK('CHECK_TEMP',1_jpim,ZHOOK_HANDLE)
  END SUBROUTINE

END SUBROUTINE
