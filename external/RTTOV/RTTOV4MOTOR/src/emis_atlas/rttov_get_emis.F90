! Description:
!> @file
!!   Returns emissivity values for the given chanprof array.
!
!> @brief
!!   Returns emissivity values for the given chanprof array.
!!   Optionally, returns emissivity standard deviations or, for the
!!   TELSEM2 MW atlas, covariances. For the UWIRemis and CAMEL 2007 IR atlases,
!!   it also optionally outputs the Principal Components eigenvalues and
!!   eigenvectors for the output emissivities.
!!
!! @details
!!   The atlas data structure must have been initialised by calling
!!   rttov_setup_emis_atlas on it before calling this subroutine.
!!
!!   For IR instruments where the atlas data were initialised for a
!!   single instrument and for the CNRM MW atlas you must ensure
!!   the coefs and the atlas data correspond to that instrument.
!!   Otherwise any coefficient structure can be used with the
!!   atlas data - but of course an atlas data structure initialised
!!   for an IR atlas can only be used with IR sensors and similarly
!!   for MW atlases and sensors.
!!
!!   Surface types (as determined by profiles(:)\%skin\%surftype):
!!   IR atlases: if the surftype is sea no emissivities are returned;
!!     for sea-ice the atlas returns values interpolated from a
!!     fixed sea-ice emissivity spectrum;
!!     for land surfaces the atlas returns emissivities where it has
!!     data and will linearly combine these with snow emissivity
!!     values if the snow_fraction variable is non-zero.
!!
!!   The CAMEL climatology atlas has a more sophisticated approach to snow:
!!     if snow_fraction = 1: the atlas uses it's snow lab spectra if available
!!       or otherwise uses a fixed snow emissivity spectrum.
!!     if 0 < snow_fraction < 1: if the atlas has any snow-contamination it
!!       simply returns the atlas emissivity, otherwise it linearly combines the
!!       (snow-free) atlas emissivity with a fixed snow spectrum weighted by
!!       snow_fraction.
!!     if snow_fraction = 0 and snow_correction is true: if the atlas emissivity
!!       has some snow-contamination, then it tries to find a snow free lab
!!       spectrum to use (failing that it just returns the snow-contaminated
!!       atlas emissivity)
!!
!!   TELSEM2: this atlas contains climatological emissivity values
!!     for sea-ice. It makes no check on the surftype: if the atlas
!!     contains emissivity data for a given location this is returned
!!     otherwise it returns zeros.
!!
!!   CNRM MW atlas: for sea and sea-ice surftypes zeros are returned.
!!     The atlas returns emissivities for land where it has data.
!!
!!   After this subroutine has been called you should check the
!!   returned emissivity values: if the atlas has no value for the
!!   given location it will return values less than or equal to zero.
!!
!!   Optional arguments provide access to specific features of the atlases.
!!   Note that not all arguments apply to all atlases.
!!
!!   The important elements of the profiles structure for the MW atlases are:
!!   skin\%surftype, latitude, longitude, zenangle.
!!
!!   The important elements of the profiles structure for the IR atlases are:
!!   skin\%surftype, latitude, longitude, skin\%snow_fraction, zenangle, sunzenangle.
!!   The satellite and solar zenith angles are only required if the angle
!!   correction is being applied and sunzenangle is only used to distinguish
!!   between day (< 85 degrees) and night (> 85 degrees).
!!
!!   The output of PCs from the IR atlases is intended for advanced users. Only
!!   the UWIRemis and CAMEL 2007 atlases are currently supported because the
!!   extraction of PCs from the CAMEL climatology atlas is much more
!!   complicated. The PC scores, PC eigenvectors and PC constants can be
!!   returned. The emissivity for chanprof index i and profile index j is
!!   reconstructed as follows:
!!
!!     emis(i) = SUM(pc_eval(1:numpcs,j) * pc_evec(1:numpcs,i)) + pc_const(i)
!!
!!   The first dimension of the pc_eval and pc_evec arrays is the maximum
!!   number of PCs used by the atlas: this is given by the numpcs parameter in
!!   mod_uwiremis_atlas.F90 (6 for UWIRemis) and mod_camel_atlas.F90 (9 for the
!!   CAMEL 2007 atlas). For the CAMEL atlas, not all 9 PCs will be used in
!!   every case: the unused pc_eval and pc_evec values are always set to zero.
!!
!!   Note that pc_eval has dimensions (numpcs,nprofiles): the coefficients
!!   are the same for all channels associated with a given profile.
!!
!!   All pc_* arguments must be supplied together or none of them.
!!
!!   NB The returned PC arrays will be zero if the profile snow_fraction is
!!   greater than zero or if the surface type is not land. Also note that
!!   occasionally the atlas obtains emissivities slightly greater than 1. and
!!   these are capped at 1. In this case the PC outputs are unmodified and so
!!   the emissivity reconstructed from the PC outputs will not match the capped
!!   emissivity returned by the atlas. If the emissivity is exactly 1. and the
!!   PC reconstructed emissivity is greater than 1. then this is the reason.
!!
!! @param[out]    err               status on exit
!! @param[in]     opts              options to configure the simulations
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     profiles          input atmospheric profiles and surface variables
!! @param[in]     coefs             coefficients structure for instrument to simulate
!! @param[in]     atlas             atlas data
!! @param[out]    emissivity        emissivity values
!! @param[out]    emis_std          emissivity errors (standard deviations), IR and TELSEM2 only, optional
!! @param[out]    emis_cov          emissivity covariances, TELSEM2 only, dimensions are (nprof,nchan,nchan)
!!                                    where nchan is the largest number of channels simulated per profile, optional
!! @param[out]    emis_flag         emissivity atlas flags, IR only, optional
!! @param[out]    pc_eval           emissivity atlas PC coefficients (eigenvalues), UWIRemis and CAMEL 2007 only, optional
!! @param[out]    pc_evec           emissivity atlas PC eigenvectors, UWIRemis and CAMEL 2007 only, optional
!! @param[out]    pc_const          emissivity atlas PC constant values, UWIRemis and CAMEL 2007 only, optional
!! @param[in]     resolution        return emissivities at this resolution, TELSEM2 only,
!!                                    units: degrees latitude/longitude, default: 0.25 degrees, optional
!! @param[in]     snow_correction   flag to turn snow correction on/off when snow_frac = 0, but atlas has snow-
!!                                    contamination, CAMEL climatology atlas only (default true)
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
SUBROUTINE rttov_get_emis(      &
                err,            &! out
                opts,           &! in
                chanprof,       &! in
                profiles,       &! in
                coefs,          &! in
                atlas,          &! in
                emissivity,     &! out
                emis_std,       &! out, optional (IR atlases, TELSEM2)
                emis_cov,       &! out, optional (TELSEM2)
                emis_flag,      &! out, optional (IR atlases)
                pc_eval,        &! out, optional (UWIRemis and CAMEL 2007 atlases)
                pc_evec,        &! out, optional (UWIRemis and CAMEL 2007 atlases)
                pc_const,       &! out, optional (UWIRemis and CAMEL 2007 atlases)
                resolution,     &! in,  optional (TELSEM2)
                snow_correction) ! in,  optional (CAMEL climatology atlas only)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_types, ONLY : &
          rttov_chanprof, &
          rttov_profile,  &
          rttov_coefs,    &
          rttov_options
  USE mod_rttov_emis_atlas, ONLY : rttov_emis_atlas_data
!INTF_OFF
  USE rttov_const, ONLY :      &
          deg2rad,             &
          sensor_id_mw,        &
          sensor_id_po,        &
          surftype_land,       &
          surftype_seaice,     &
          pol_v,               &
          pol_h,               &
          errorstatus_success, &
          errorstatus_fatal

  USE mod_rttov_emis_atlas, ONLY : &
    uwiremis_atlas_id, camel_atlas_id, camel_clim_atlas_id, &
    telsem2_atlas_id, cnrm_mw_atlas_id

  USE mod_uwiremis_atlas, ONLY : rttov_uwiremis, numpcs_uwiremis => numpcs
  USE mod_camel_atlas, ONLY : rttov_camel, numpcs_camel_2007 => numpcs
  USE mod_camel_clim_atlas, ONLY : rttov_camel_clim
  USE mod_mwatlas_m2, ONLY : &
          telsem2_emis_interp_ind_mult => emis_interp_ind_mult, &
          telsem2_emis_interp_int_mult => emis_interp_int_mult
  USE mod_cnrm_mw_atlas, ONLY : rttov_cnrmmwemis

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),          INTENT(OUT)           :: err
  TYPE(rttov_options),         INTENT(IN)            :: opts
  TYPE(rttov_chanprof),        INTENT(IN)            :: chanprof(:)
  TYPE(rttov_profile),         INTENT(IN)            :: profiles(:)
  TYPE(rttov_coefs),           INTENT(IN)            :: coefs
  TYPE(rttov_emis_atlas_data), INTENT(IN)            :: atlas
  REAL(KIND=jprb),             INTENT(OUT)           :: emissivity(SIZE(chanprof))
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: emis_std(SIZE(chanprof))
  ! emis_cov dims are (nprofs, nchans, nchans) where nchans is the maximum # channels per profile
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: emis_cov(:,:,:)
  INTEGER(KIND=jpim),          INTENT(OUT), OPTIONAL :: emis_flag(SIZE(chanprof))
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: pc_eval(:,:)
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: pc_evec(:,:)
  REAL(KIND=jprb),             INTENT(OUT), OPTIONAL :: pc_const(SIZE(chanprof))
  REAL(KIND=jprb),             INTENT(IN),  OPTIONAL :: resolution
  LOGICAL(KIND=jplm),          INTENT(IN),  OPTIONAL :: snow_correction
!INTF_END

#include "rttov_errorreport.interface"

  REAL(KIND=jprb)    :: ZHOOK_HANDLE

  INTEGER(KIND=jpim) :: nchanprof                 ! number of channels/profiles
  INTEGER(KIND=jpim) :: prof                      ! index of current profile
  INTEGER(KIND=jpim) :: nchan                     ! number of channels in current profile
  INTEGER(KIND=jpim) :: nchan_coef                ! number of channels in coefficient file
  INTEGER(KIND=jpim), ALLOCATABLE :: chans(:)     ! channel indexes for current profile
  INTEGER(KIND=jpim) :: id_platform, id_sat, id_inst
  REAL(KIND=jprb), POINTER :: chan_wvn(:)
  CHARACTER(LEN=256) :: msg

  LOGICAL(KIND=jplm) :: htfrtc, sensor_mw

  INTEGER(KIND=jpim) :: i, j, k                   ! loop indexes
  INTEGER(KIND=jpim) :: lo, hi                    ! low/high chanprof indexes for current profile
  INTEGER(KIND=jpim) :: chan

  ! Local variables for IR atlas
  INTEGER(KIND=jpim)              :: npcs
  INTEGER(KIND=jpim)              :: instr_emis_flag
  REAL(KIND=jprb), ALLOCATABLE    :: instr_emis_cov(:)
  LOGICAL(KIND=jplm)              :: snow_corr, get_pc

  ! Local variables for TELSEM2
  REAL(KIND=jprb), ALLOCATABLE    :: ev(:), eh(:)
  REAL(KIND=jprb), ALLOCATABLE    :: std(:,:)
  REAL(KIND=jprb)                 :: sinzen, sinview, sinview_sq, cosview_sq
  INTEGER(KIND=jpim), ALLOCATABLE :: pol_id(:)
  REAL(KIND=jprb)   , ALLOCATABLE :: emissfactor_h(:), emissfactor_v(:)

!----------------------------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS', 0_jpim, ZHOOK_HANDLE)

  !-----------------------------
  ! Initialise output arguments
  !-----------------------------

  emissivity(:) = -1._jprb
  IF (PRESENT(emis_std))  emis_std  = 0._jprb
  IF (PRESENT(emis_cov))  emis_cov  = 0._jprb
  IF (PRESENT(emis_flag)) emis_flag = 0_jpim
  IF (PRESENT(pc_eval))   pc_eval   = 0._jprb
  IF (PRESENT(pc_evec))   pc_evec   = 0._jprb
  IF (PRESENT(pc_const))  pc_const  = 0._jprb

  IF (.NOT. atlas%init) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Atlas data not initialised')
  ENDIF

  ! Determine if user is calling HTFRTC
  htfrtc = ASSOCIATED(coefs%coef_htfrtc%freq)

  IF (htfrtc) THEN
    id_platform = 0  ! HTFRTC coefs don't contain these IDs (and in any case
    id_sat      = 0  !   they are irrelevant since the centroid frequencies
    id_inst     = 0  !   are the same regardless of the sensor).
    nchan_coef  = coefs%coef_htfrtc%n_f
    nchanprof   = nchan_coef * SIZE(profiles)
    chan_wvn    => coefs%coef_htfrtc%freq
    sensor_mw   = .FALSE.
  ELSE
    id_platform = coefs%coef%id_platform
    id_sat      = coefs%coef%id_sat
    id_inst     = coefs%coef%id_inst
    nchan_coef  = coefs%coef%fmv_chn
    nchanprof   = SIZE(chanprof)
    chan_wvn    => coefs%coef%ff_cwn
    sensor_mw   = coefs%coef%id_sensor == sensor_id_mw .OR. &
                  coefs%coef%id_sensor == sensor_id_po
  ENDIF
  ALLOCATE(chans(nchan_coef))

  IF (htfrtc) THEN
    nchan = nchan_coef
    DO i = 1, nchan
      chans(i) = i
    ENDDO
  ENDIF

  !-------------------------------------------
  ! Do array allocations before chanprof loop
  !-------------------------------------------

  ! Note that nchan_coef is the largest possible number of channels per profile

  ! Also check optional arguments here and report errors when arguments are supplied
  ! which are not supported by the current atlas.

  IF (sensor_mw) THEN
    ! MW sensor

    IF (.NOT. atlas%is_mw) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Atlas data is initialised for IR sensors')
    ENDIF

    IF (atlas%atlas_id == telsem2_atlas_id) THEN

      ALLOCATE(ev(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of ev')
      ALLOCATE(eh(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of eh')

      IF (PRESENT(emis_std) .OR. PRESENT(emis_cov)) THEN
        ! Emissivity errors or covariances requested
        ALLOCATE(std(2*nchan_coef,2*nchan_coef),stat=err)
        THROWM(err .NE. 0, 'Allocation of std')
      ENDIF

      ALLOCATE(pol_id(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of pol_id')
      ALLOCATE(emissfactor_v(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of emissfactor_v')
      ALLOCATE(emissfactor_h(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of emissfactor_h')

      IF (PRESENT(snow_correction) .AND. opts%config%verbose) THEN
        WARN('snow_correction argument not supported for TELSEM2 atlas')
      ENDIF
      IF (PRESENT(emis_flag) .AND. opts%config%verbose) THEN
        WARN('emis_flag argument not supported for TELSEM2 atlas')
      ENDIF
      IF ((PRESENT(pc_eval) .OR. PRESENT(pc_evec) .OR. PRESENT(pc_const)) &
          .AND. opts%config%verbose) THEN
        WARN('pc_* arguments not supported for TELSEM2 atlas')
      ENDIF

    ELSEIF (atlas%atlas_id == cnrm_mw_atlas_id) THEN

      IF (PRESENT(snow_correction) .AND. opts%config%verbose) THEN
        WARN('snow_correction argument not supported for CNRM MW atlas')
      ENDIF
      IF (PRESENT(resolution) .AND. opts%config%verbose) THEN
        WARN('resolution argument not supported for CNRM MW atlas')
      ENDIF
      IF (PRESENT(emis_std) .AND. opts%config%verbose) THEN
        WARN('emis_std argument not supported for CNRM MW atlas')
      ENDIF
      IF (PRESENT(emis_cov) .AND. opts%config%verbose) THEN
        WARN('emis_cov argument not supported for CNRM MW atlas')
      ENDIF
      IF (PRESENT(emis_flag) .AND. opts%config%verbose) THEN
        WARN('emis_flag argument not supported for CNRM MW atlas')
      ENDIF
      IF ((PRESENT(pc_eval) .OR. PRESENT(pc_evec) .OR. PRESENT(pc_const)) &
          .AND. opts%config%verbose) THEN
        WARN('pc_* arguments not supported for CNRM MW atlas')
      ENDIF

    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Unknown MW atlas ID')
    ENDIF

  ELSE
    ! IR sensor

    IF (atlas%is_mw) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Atlas data is initialised for MW sensors')
    ENDIF

    get_pc = .FALSE.

    IF (atlas%atlas_id == uwiremis_atlas_id .OR. &
        atlas%atlas_id == camel_atlas_id) THEN

      IF (PRESENT(snow_correction) .AND. opts%config%verbose) THEN
        WARN('snow_correction argument only supported by CAMEL climatology atlas')
      ENDIF

      IF (PRESENT(pc_eval) .OR. PRESENT(pc_evec) .OR. PRESENT(pc_const)) THEN

        IF (.NOT. PRESENT(pc_eval) .OR. &
            .NOT. PRESENT(pc_evec) .OR. &
            .NOT. PRESENT(pc_const)) THEN
          err = errorstatus_fatal
          THROWM(err .NE. 0, 'All pc_* output arguments must be supplied or none')
        ENDIF

        ! Check the numpcs dimension of pc_* arguments
        SELECT CASE (atlas%atlas_id)
        CASE (uwiremis_atlas_id)
          npcs = numpcs_uwiremis
        CASE (camel_atlas_id)
          npcs = numpcs_camel_2007
        ENDSELECT
        IF (SIZE(pc_eval, DIM=1) /= npcs) THEN
          err = errorstatus_fatal
          WRITE(msg,'(A,I2)') 'First dimension of pc_eval must be numpcs = ', npcs
          THROWM(err .NE. 0, TRIM(msg))
        ENDIF
        IF (SIZE(pc_evec, DIM=1) /= npcs) THEN
          err = errorstatus_fatal
          WRITE(msg,'(A,I2)') 'First dimension of pc_evec must be numpcs = ', npcs
          THROWM(err .NE. 0, TRIM(msg))
        ENDIF

        get_pc = .TRUE.

      ENDIF

    ELSE

      snow_corr = .TRUE.
      IF (PRESENT(snow_correction)) snow_corr = snow_correction

      IF ((PRESENT(pc_eval) .OR. PRESENT(pc_evec) .OR. PRESENT(pc_const)) &
          .AND. opts%config%verbose) THEN
        WARN('pc_* arguments not supported for CAMEL climatology atlas')
      ENDIF

    ENDIF

    IF (atlas%atlas_id == uwiremis_atlas_id .OR. &
        atlas%atlas_id == camel_atlas_id .OR. &
        atlas%atlas_id == camel_clim_atlas_id) THEN

      ALLOCATE(instr_emis_cov(nchan_coef),stat=err)
      THROWM(err .NE. 0, 'Allocation of instr_emis_cov')

      IF (PRESENT(resolution) .AND. opts%config%verbose) THEN
        WARN('resolution argument not supported for IR atlas')
      ENDIF
      IF (PRESENT(emis_cov) .AND. opts%config%verbose) THEN
        WARN('emis_cov argument not supported for IR atlas')
      ENDIF

    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Unknown IR atlas ID')
    ENDIF

  ENDIF

  !-----------
  ! Main loop
  !-----------

  ! Assume chanprof(:) lists all channels for profile 1, followed by all channels
  ! for profile 2, and so on. Loop through chanprof(:) until the profile number
  ! changes and in this way obtain the list of channels for the current profile.
  k = 1
  prof = 0
  DO WHILE (k <= nchanprof)
    IF (htfrtc) THEN
      prof = prof + 1
      k = k + nchan_coef
      lo = (prof - 1) * nchan_coef + 1
      hi = lo + nchan_coef - 1
    ELSE
      prof = chanprof(k)%prof
      lo = k
      k  = lo + 1
      DO
        IF (k > nchanprof) EXIT
        IF (prof /= chanprof(k)%prof) EXIT
        k = k + 1
      ENDDO
      hi = k - 1
      nchan = k - lo
      IF (nchan > nchan_coef) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, 'Too many (duplicate?) channels specified per profile')
      ENDIF
      chans(1:nchan) = chanprof(lo:hi)%chan
    ENDIF

    ! Profile index is 'prof'
    ! Number of channels for this profile is 'nchan'
    ! Channel indexes are 'chans(1:nchan)'

    IF (sensor_mw) THEN
    !---------------------------
    ! MW atlas
    !---------------------------
      IF (atlas%atlas_id == telsem2_atlas_id) THEN

        ev(:) = 0._jprb
        eh(:) = 0._jprb

        IF (PRESENT(emis_std) .OR. PRESENT(emis_cov)) THEN
          ! Emissivity errors or covariances requested

          std(:,:) = 0._jprb

          IF (PRESENT(resolution)) THEN
            ! User resolution
            CALL telsem2_emis_interp_int_mult(              &
                  profiles(prof)%latitude,                  &! in
                  profiles(prof)%longitude,                 &! in
                  resolution,                               &! in
                  profiles(prof)%zenangle,                  &! in
                  coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  nchan,                                    &! in
                  atlas%telsem2_atlas,                      &! in
                  ev(1:nchan),                              &! out
                  eh(1:nchan),                              &! out
                  std=std(1:(2*nchan),1:(2*nchan)),         &! out
                  verb=0_jpim                               )! in
          ELSE
            ! Atlas resolution
            CALL telsem2_emis_interp_ind_mult(              &
                  profiles(prof)%latitude,                  &! in
                  profiles(prof)%longitude,                 &! in
                  profiles(prof)%zenangle,                  &! in
                  coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  nchan,                                    &! in
                  atlas%telsem2_atlas,                      &! in
                  ev(1:nchan),                              &! out
                  eh(1:nchan),                              &! out
                  std=std(1:(2*nchan),1:(2*nchan)),         &! out
                  verb=0_jpim                               )! in
          ENDIF

        ELSE
          ! No need to supply the std argument

          IF (PRESENT(resolution)) THEN
            ! User resolution
            CALL telsem2_emis_interp_int_mult(              &
                  profiles(prof)%latitude,                  &! in
                  profiles(prof)%longitude,                 &! in
                  resolution,                               &! in
                  profiles(prof)%zenangle,                  &! in
                  coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  nchan,                                    &! in
                  atlas%telsem2_atlas,                      &! in
                  ev(1:nchan),                              &! out
                  eh(1:nchan),                              &! out
                  verb=0_jpim                               )! in
          ELSE
            ! Atlas resolution
            CALL telsem2_emis_interp_ind_mult(              &
                  profiles(prof)%latitude,                  &! in
                  profiles(prof)%longitude,                 &! in
                  profiles(prof)%zenangle,                  &! in
                  coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  nchan,                                    &! in
                  atlas%telsem2_atlas,                      &! in
                  ev(1:nchan),                              &! out
                  eh(1:nchan),                              &! out
                  verb=0_jpim                               )! in
          ENDIF
        ENDIF

        !------------------------------------------------------
        ! Combine H- and V-pol emissivities into output arrays
        !------------------------------------------------------

        sinzen     = SIN(profiles(prof)%zenangle * deg2rad)
        sinview    = sinzen / coefs%coef%ratoe
        sinview_sq = sinview * sinview
        cosview_sq = 1.0_jprb - sinview_sq

        pol_id(1:nchan) = coefs%coef%fastem_polar(chans(1:nchan)) + 1_jpim

        DO j = 1, nchan
          chan = chanprof(j)%chan
          emissfactor_v(j) = pol_v(1,pol_id(j)) * coefs%coef%pol_fac_v(chan) + &
                             pol_v(2,pol_id(j)) * sinview_sq + pol_v(3,pol_id(j)) * cosview_sq
          emissfactor_h(j) = pol_h(1,pol_id(j)) * coefs%coef%pol_fac_h(chan) + &
                             pol_h(2,pol_id(j)) * sinview_sq + pol_h(3,pol_id(j)) * cosview_sq
        ENDDO
        emissivity(lo:hi) = ev(1:nchan) * emissfactor_v(1:nchan) + eh(1:nchan) * emissfactor_h(1:nchan)

        IF (PRESENT(emis_std)) THEN
          ! Generate the stdv

          ! Each dimension of std(:,:) has V-pol values for all channels followed by H-pol values
          ! for all channels.

          ! std(:,:) contains covariances so diagonal elements are variances (not stddev),
          ! but we want output emis_std(:) to be stddev

          DO j = 1, nchan
            emis_std(lo+j-1) = SQRT(emissfactor_v(j) * emissfactor_v(j) * std(j,j)            + &
                                    emissfactor_v(j) * emissfactor_h(j) * 2 * std(j,j+nchan)  + &
                                    emissfactor_h(j) * emissfactor_h(j) * std(j+nchan,j+nchan))
          ENDDO
        ENDIF

        IF (PRESENT(emis_cov)) THEN
          ! Generate the covariances

          ! std(:,:) contains covariances. Combine the H- and V-pol values for each pair of
          ! channels to find the covariance of the combined emissivities in the two channels.

          DO i = 1, nchan
            DO j = 1, nchan
              emis_cov(prof,i,j) = emissfactor_v(i) * emissfactor_v(j) * std(i,       j      ) + &
                                   emissfactor_v(i) * emissfactor_h(j) * std(i,       j+nchan) + &
                                   emissfactor_h(i) * emissfactor_v(j) * std(i+nchan, j      ) + &
                                   emissfactor_h(i) * emissfactor_h(j) * std(i+nchan, j+nchan)
            ENDDO
          ENDDO
        ENDIF

      ENDIF ! TELSEM2

      IF (atlas%atlas_id == cnrm_mw_atlas_id) THEN

        IF (profiles(prof)%skin%surftype == surftype_land) THEN
          CALL rttov_cnrmmwemis(                            &
                  err,                                      &
                  coefs%coef%id_inst,                       &! in
                  nchan,                                    &! in 
                  coefs%coef%frequency_ghz(chans(1:nchan)), &! in
                  coefs%coef%fastem_polar(chans(1:nchan)),  &! in
                  atlas%cnrm_mw_atlas,                      &! in
                  profiles(prof)%latitude,                  &! in
                  profiles(prof)%longitude,                 &! in
                  profiles(prof)%zenangle,                  &! in
                  emissivity(lo:hi))                         ! out
          THROWM(err .NE. 0, "error in CNRM atlas")

        ELSE
          ! Not land
          emissivity(lo:hi) = -1._jprb
        ENDIF

      ENDIF ! CNRM MW atlas

    ELSE
    !---------------------------
    ! IR atlas
    !---------------------------
      IF (atlas%atlas_id == uwiremis_atlas_id) THEN

        IF (profiles(prof)%skin%surftype == surftype_land .OR. &
            profiles(prof)%skin%surftype == surftype_seaice) THEN

          IF (get_pc) THEN
            CALL rttov_uwiremis(                        &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%snow_fraction,  &! in
                    chan_wvn(chans(1:nchan)),           &! in
                    chans(1:nchan),                     &! in
                    id_platform,                        &! in
                    id_sat,                             &! in
                    id_inst,                            &! in
                    nchan_coef,                         &! in
                    atlas%uwiremis_atlas,               &! in
                    emissivity(lo:hi),                  &! out
                    instr_emis_cov(1:nchan),            &! out
                    instr_emis_flag,                    &! out
                    pc_eval(:,prof),                    &! out
                    pc_evec(:,lo:hi),                   &! out
                    pc_const(lo:hi))                     ! out
            THROWM(err .NE. 0, "error in UWIRemis atlas")

          ELSE
            CALL rttov_uwiremis(                        &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%snow_fraction,  &! in
                    chan_wvn(chans(1:nchan)),           &! in
                    chans(1:nchan),                     &! in
                    id_platform,                        &! in
                    id_sat,                             &! in
                    id_inst,                            &! in
                    nchan_coef,                         &! in
                    atlas%uwiremis_atlas,               &! in
                    emissivity(lo:hi),                  &! out
                    instr_emis_cov(1:nchan),            &! out
                    instr_emis_flag)                     ! out
            THROWM(err .NE. 0, "error in UWIRemis atlas")
          ENDIF

          IF (PRESENT(emis_std))  emis_std(lo:hi)  = instr_emis_cov(1:nchan)
          IF (PRESENT(emis_flag)) emis_flag(lo:hi) = instr_emis_flag

        ELSE
          ! Not land or seaice
          emissivity(lo:hi) = -1._jprb
        ENDIF

      ENDIF

      IF (atlas%atlas_id == camel_atlas_id) THEN

        IF (profiles(prof)%skin%surftype == surftype_land .OR. &
            profiles(prof)%skin%surftype == surftype_seaice) THEN

          IF (get_pc) THEN
            CALL rttov_camel(                           &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%snow_fraction,  &! in
                    chan_wvn(chans(1:nchan)),           &! in
                    chans(1:nchan),                     &! in
                    id_platform,                        &! in
                    id_sat,                             &! in
                    id_inst,                            &! in
                    nchan_coef,                         &! in
                    atlas%camel_atlas,                  &! in
                    emissivity(lo:hi),                  &! out
                    instr_emis_cov(1:nchan),            &! out
                    instr_emis_flag,                    &! out
                    pc_eval(:,prof),                    &! out
                    pc_evec(:,lo:hi),                   &! out
                    pc_const(lo:hi))                     ! out
            THROWM(err .NE. 0, "error in CAMEL atlas")

          ELSE
            CALL rttov_camel(                           &
                    err,                                &! out
                    opts%config%verbose,                &! in
                    nchan,                              &! in
                    profiles(prof)%latitude,            &! in
                    profiles(prof)%longitude,           &! in
                    profiles(prof)%zenangle,            &! in
                    profiles(prof)%sunzenangle,         &! in
                    profiles(prof)%skin%surftype,       &! in
                    profiles(prof)%skin%snow_fraction,  &! in
                    chan_wvn(chans(1:nchan)),           &! in
                    chans(1:nchan),                     &! in
                    id_platform,                        &! in
                    id_sat,                             &! in
                    id_inst,                            &! in
                    nchan_coef,                         &! in
                    atlas%camel_atlas,                  &! in
                    emissivity(lo:hi),                  &! out
                    instr_emis_cov(1:nchan),            &! out
                    instr_emis_flag)                     ! out
          ENDIF
          THROWM(err .NE. 0, "error in CAMEL atlas")

          IF (PRESENT(emis_std))  emis_std(lo:hi)  = instr_emis_cov(1:nchan)
          IF (PRESENT(emis_flag)) emis_flag(lo:hi) = instr_emis_flag

        ELSE
          ! Not land or seaice
          emissivity(lo:hi) = -1._jprb
        ENDIF

      ENDIF

      IF (atlas%atlas_id == camel_clim_atlas_id) THEN

        IF (profiles(prof)%skin%surftype == surftype_land .OR. &
            profiles(prof)%skin%surftype == surftype_seaice) THEN

          CALL rttov_camel_clim(                      &
                  err,                                &! out
                  opts%config%verbose,                &! in
                  snow_corr,                          &! in
                  nchan,                              &! in
                  profiles(prof)%latitude,            &! in
                  profiles(prof)%longitude,           &! in
                  profiles(prof)%zenangle,            &! in
                  profiles(prof)%sunzenangle,         &! in
                  profiles(prof)%skin%surftype,       &! in
                  profiles(prof)%skin%snow_fraction,  &! in
                  chan_wvn(chans(1:nchan)),           &! in
                  chans(1:nchan),                     &! in
                  id_platform,                        &! in
                  id_sat,                             &! in
                  id_inst,                            &! in
                  nchan_coef,                         &! in
                  atlas%camel_clim_atlas,             &! in
                  emissivity(lo:hi),                  &! out
                  instr_emis_cov(1:nchan),            &! out
                  instr_emis_flag)                     ! out
          THROWM(err .NE. 0, "error in CAMEL climatology atlas")

          IF (PRESENT(emis_std))  emis_std(lo:hi)  = instr_emis_cov(1:nchan)
          IF (PRESENT(emis_flag)) emis_flag(lo:hi) = instr_emis_flag

        ELSE
          ! Not land or seaice
          emissivity(lo:hi) = -1._jprb
        ENDIF

      ENDIF

    ENDIF ! coefs%coef%id_sensor

  ENDDO ! chanprof

  WHERE (emissivity <= 0._jprb)
    emissivity = -1._jprb
  ENDWHERE


  !------------------
  ! Deallocate arrays
  !------------------
  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS',1_jpim,ZHOOK_HANDLE)

  CATCH

  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_GET_EMIS',1_jpim,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    IF (ALLOCATED(chans)) DEALLOCATE(chans)

    IF (sensor_mw) THEN
      ! MW sensor
      IF (atlas%atlas_id == telsem2_atlas_id) THEN
        IF (ALLOCATED(ev)) DEALLOCATE(ev)
        IF (ALLOCATED(eh)) DEALLOCATE(eh)

        IF (ALLOCATED(std)) DEALLOCATE(std)

        IF (ALLOCATED(pol_id))        DEALLOCATE(pol_id)
        IF (ALLOCATED(emissfactor_v)) DEALLOCATE(emissfactor_v)
        IF (ALLOCATED(emissfactor_h)) DEALLOCATE(emissfactor_h)
      ENDIF

      ! Nothing to do for CNRM atlas

    ELSE
      ! IR sensor
      IF (ALLOCATED(instr_emis_cov)) DEALLOCATE(instr_emis_cov)
    ENDIF
  END SUBROUTINE cleanup
END SUBROUTINE rttov_get_emis
