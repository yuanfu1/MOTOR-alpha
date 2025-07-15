! Description:
!> @file
!!   Main executable for the RTTOV test suite. This should be run via the
!!   rttov_test.pl script.
!
!> @brief
!!   Main executable for the RTTOV test suite. This should be run via the
!!   rttov_test.pl script.
!!
!! @details
!!   When used in conjunction with rttov_test.pl this executable allows
!!   many aspects of RTTOV simulations to be configured at run-time via
!!   command-line options. In addition to simply running the direct, TL,
!!   AD and K models this can also compute Jacobians using the TL and AD,
!!   a brute force Jacobian using the direct model and it can run the
!!   Taylor test to check direct/TL model consistency.
!!
!!   It allows for RTTOV to be run via the standard interface and the
!!   parallel interface.
!!
!!   Input profiles can be multiplied up to test passing many profiles into
!!   RTTOV.
!!
!!   Output is written to ASCII files by default: if selected no output will
!!   be written which is useful for running timing tests.
!!
!!   See the test suite documentation on how to call this program using
!!   rttov_test.pl.
!!
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_test

  USE rttov_const, ONLY :      &
         errorstatus_success,  &
         errorstatus_fatal,    &
         pi_r,                 &
         gas_unit_specconc,    &
         gas_unit_ppmv,        &
         gas_unit_ppmvdry

  USE rttov_types, ONLY :     &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_pccomp,        &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_radiance2,     &
         rttov_reflectivity,  &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance,   &
         rttov_opt_param,     &
         rttov_skin,          &
         rttov_s2m,           &
         rttov_traj,          &
         rttov_scatt_emis_retrieval_type

#include "throw.h"

  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_getoptions, ONLY : getoption
  USE yomhook, ONLY : lhook, dr_hook
  USE rttov_lun, ONLY : rttov_get_lun, rttov_put_lun
  USE rttov_unix_env, ONLY: rttov_exit
  USE rttov_test_mod
  USE rttov_test_k_mod, ONLY : radar_k_lev

#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_tl.interface"
#include "rttov_ad.interface"
#include "rttov_k.interface"
#include "rttov_k_tl.interface"
#include "rttov_k_ad.interface"
#include "rttov_k_bf.interface"
#include "rttov_taylor_test.interface"

#include "rttov_parallel_direct.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_parallel_ad.interface"
#include "rttov_parallel_k.interface"

#include "rttov_scatt.interface"
#include "rttov_scatt_tl.interface"
#include "rttov_scatt_ad.interface"

#include "rttov_parallel_scatt.interface"
#include "rttov_parallel_scatt_tl.interface"
#include "rttov_parallel_scatt_ad.interface"

#include "rttov_make_profile_inc.interface"
#include "rttov_scale_profile_inc.interface"
#include "rttov_make_cld_profile_inc.interface"
#include "rttov_scale_cld_profile_inc.interface"
#include "rttov_make_opt_param_inc.interface"
#include "rttov_scale_opt_param_inc.interface"
#include "rttov_make_radiance_inc.interface"
#include "rttov_scale_radiance_inc.interface"
#include "rttov_make_pccomp_inc.interface"
#include "rttov_scale_pccomp_inc.interface"
#include "rttov_make_emisrefl_inc.interface"

#include "rttov_read_coefs.interface"
#include "rttov_read_coefs_htfrtc.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_alloc_emis_ret_terms.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_scatt_prof.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_scatt_prof.interface"
#include "rttov_copy_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_scale_ref_gas_prof.interface"

#include "rttov_print_radiance_quality.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_cld_profile.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_user_profile_checkinput.interface"

#include "rttov_errorreport.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

! Test definition; we have here all dimensions, flags, etc...
  TYPE rttov_test_defn
! dimensions
    INTEGER(jpim) :: nlevels = -1_jpim, nhydro_frac = -1_jpim, nchannels = -1_jpim, nprofiles = -1_jpim, mult = 1_jpim
    INTEGER(jpim) :: nchannels_rec = -1_jpim, npcscores = -1_jpim, prof_ndigits = -1_jpim
    INTEGER(jpim) :: aer_nphangle = -1_jpim, cld_nphangle = -1_jpim, aer_nmom = -1_jpim, cld_nmom = -1_jpim

! coefficients files
    CHARACTER(LEN=256) :: f_coef = "", f_scaer = "", f_sccld = "", f_pccomp = ""
    CHARACTER(LEN=256) :: f_rttovscatt = "", f_mfasis_cld = "", f_mfasis_aer = ""
    CHARACTER(LEN=256) :: f_htfrtc_coef = "", f_htfrtc_sensor = ""
    CHARACTER(LEN=256) :: coef_prefix

! flags
    TYPE(rttov_options)       :: opts
    TYPE(rttov_options_scatt) :: opts_scatt

! input files
    CHARACTER(LEN=256) :: f_p = "", f_t = "", f_q = "", f_o3 = "", f_co2 = "", f_co = "",  &
                          f_n2o = "", f_ch4  = "", f_so2  = "", f_clw = "", &
                          f_aerosl = "", f_cloud = "", f_cfrac = "", &
                          f_skin = "", f_s2m = "", f_angles = "", f_clw_scheme = "", f_ice_scheme = "", &
                          f_simple_cloud = "", f_flags = "", f_clwde = "", f_icede = "", &
                          f_be = "", f_datetime = "", f_gas_units = "", &
                          f_mmr_cldaer = "", f_prof_id = ""
    CHARACTER(LEN=256) :: f_channels = "", f_pmc = "", &
                          f_calcemis = "", f_emissivity = "", f_specularity = "", &
                          f_calcrefl = "", f_reflectance = "", f_reflectance_diffuse = "", &
                          f_reflectance_cloud_top = "", &
                          f_lprofiles = "", f_pcscores = "", f_channels_rec = ""
    CHARACTER(LEN=256) :: f_aer_opt_param = "", f_cld_opt_param = ""
    CHARACTER(LEN=256) :: f_cld_profiles = ""

! what we will do
    LOGICAL(jplm) :: do_direct = .FALSE., do_tl   = .FALSE., do_ad   = .FALSE., do_k = .FALSE., &
                     do_k_bf   = .FALSE., do_k_tl = .FALSE., do_k_ad = .FALSE.,                 &
                     do_taylor = .FALSE., taylor_by_chan = .FALSE., taylor_on_btrefl = .FALSE.
    LOGICAL(jplm) :: do_rttovscatt = .FALSE., radar = .FALSE., calc_emis_terms = .FALSE., multi_hydro_frac = .FALSE.
    LOGICAL(jplm) :: do_print  = .FALSE., calc_rad2 = .FALSE., ltemp = .FALSE., lallchans = .FALSE.
    LOGICAL(jplm) :: print_opts = .FALSE., print_profiles = .FALSE., print_quality = .FALSE.
    LOGICAL(jplm) :: savehdf5 = .FALSE.
    INTEGER(jpim) :: ntimes = 1, nthreads = 0, parallel_strategy = 0
    LOGICAL(jplm) :: prof_by_prof = .FALSE. ! run RTTOV profile by profile
    LOGICAL(jplm) :: chan_by_chan = .FALSE. ! run RTTOV channel by channel
    LOGICAL(jplm) :: user_check_opts = .FALSE., user_check_prof = .FALSE.
    INTEGER(jpim) :: input_gas_units = -1           ! Gas units in input file
    INTEGER(jpim) :: run_gas_units = gas_unit_ppmv  ! Gas units RTTOV is run with
    REAL(jprb)    :: scale_inc = 1._jprb
    REAL(jprb)    :: scale_out = 1._jprb
    REAL(jprb)    :: fixedemis = -1._jprb
    REAL(jprb)    :: fixedrefl = -1._jprb
    REAL(jprb)    :: fixeddiffuse_refl = -1._jprb
    REAL(jprb)    :: fixedzenang = -1._jprb
    REAL(jprb)    :: fixedsunzenang = -1._jprb
    INTEGER(jpim) :: fixedsurftype = -1
    REAL(jprb)    :: fixedspecularity = -1._jprb
    LOGICAL(jplm) :: fixed2m = .FALSE.
    INTEGER(jpim) :: fixedvisir_clw_scheme = -1
    INTEGER(jpim) :: fixedvisir_ice_scheme = -1
    REAL(jprb)    :: fixedcfraction = -1._jprb
    LOGICAL(jplm) :: zero_cloud = .FALSE.
    LOGICAL(jplm) :: zero_aerosol = .FALSE.
    LOGICAL(jplm) :: no_flux_conv = .FALSE.
    REAL(jprb)    :: o3_col_int = -1._jprb
    REAL(jprb)    :: o3_col_int_du = -1._jprb
    REAL(jprb)    :: o3_max_ppmv = -1._jprb
    REAL(jprb)    :: co2_col_int = -1._jprb
    REAL(jprb)    :: co2_max_ppmv = -1._jprb
    REAL(jprb)    :: n2o_col_int = -1._jprb
    REAL(jprb)    :: n2o_max_ppmv = -1._jprb
    REAL(jprb)    :: co_col_int = -1._jprb
    REAL(jprb)    :: co_max_ppmv = -1._jprb
    REAL(jprb)    :: ch4_col_int = -1._jprb
    REAL(jprb)    :: ch4_max_ppmv = -1._jprb
    REAL(jprb)    :: so2_col_int = -1._jprb
    REAL(jprb)    :: so2_max_ppmv = -1._jprb
    INTEGER(jpim) :: realprec = 6    ! Number of sig figures in real output
  END TYPE

! The basic data structures passed to rttov_direct/rttov_tl/rttov_ad/rttov_k
  TYPE rttov_input_data
    TYPE(rttov_profile),          POINTER :: profiles(:)         => NULL()
    TYPE(rttov_profile_cloud),    POINTER :: cld_profiles(:)     => NULL()
    TYPE(rttov_opt_param)                 :: aer_opt_param
    TYPE(rttov_opt_param)                 :: cld_opt_param
    TYPE(rttov_emissivity),       POINTER :: emissivity(:)       => NULL()
    TYPE(rttov_reflectance),      POINTER :: reflectance(:)      => NULL()
    TYPE(rttov_transmission)              :: transmission
    TYPE(rttov_radiance)                  :: radiance
    TYPE(rttov_radiance)                  :: radiance_saved
    TYPE(rttov_radiance2)                 :: radiance2
    TYPE(rttov_reflectivity)              :: reflectivity
    TYPE(rttov_reflectivity)              :: reflectivity_saved
    TYPE(rttov_pccomp)                    :: pccomp
    TYPE(rttov_pccomp)                    :: pccomp_saved
    REAL(jprb),                   POINTER :: rttovscatt_cfrac(:)  => NULL()
    TYPE(rttov_scatt_emis_retrieval_type) :: emis_terms
  END TYPE

  TYPE rttov_k_input_data
    TYPE(rttov_profile),  POINTER :: profiles_k_pc(:) => NULL()
    TYPE(rttov_profile),  POINTER :: profiles_k_rec(:) => NULL()
  END TYPE

! Test data. All dynamically allocated data goes here.
  TYPE rttov_test_data
    TYPE(rttov_coefs)             :: coefs
    TYPE(rttov_scatt_coef)        :: coefs_scatt
    TYPE(rttov_chanprof), POINTER :: chanprof(:)          => NULL()
    INTEGER(jpim),        POINTER :: frequencies(:)       => NULL()
    INTEGER(jpim)                 :: rttov_errorstatus
    TYPE(rttov_input_data)        :: direct, tl, ad, k, k_bf, k_tl, k_ad
    TYPE(rttov_k_input_data)      :: data_k, data_k_bf, data_k_tl, data_k_ad
    LOGICAL(jplm),        POINTER :: calcemis(:)          => NULL()
    LOGICAL(jplm),        POINTER :: calcrefl(:)          => NULL()
    TYPE(rttov_traj)              :: traj, traj_tl, traj_ad, traj_k
    INTEGER(jpim),        POINTER :: channels_rec(:)      => NULL()
  END TYPE

  TYPE rttov_test_struct
    INTEGER(jpim)         :: rank
    CHARACTER(LEN=256)    :: path
    TYPE(rttov_test_defn) :: defn
    TYPE(rttov_test_data) :: data
    INTEGER(jpim)         :: err = 0_jpim
  END TYPE

  TYPE(rttov_test_struct), POINTER :: ts(:)     => NULL()
  CHARACTER(LEN=256),      POINTER :: paths(:)  => NULL()
  INTEGER(jpim) :: ip, np
  REAL(jprb)    :: zhook_handle
  INTEGER(jpim) :: err, rttov_run_err

  TRY

  IF (lhook) CALL dr_hook('rttov_test',0_jpim,zhook_handle)

  CALL getoption("--test-path-list", paths)

  IF (.NOT. ASSOCIATED(paths)) THEN
    ALLOCATE(paths(1), stat = err)
    THROWM(err.NE.0,"Allocation of paths failed")
    paths(1) = "."
  ENDIF

  np = SIZE(paths)

  ALLOCATE(ts(np), stat = err)
  THROWM(err.NE.0,"Allocation of ts failed")
  DO ip = 1, np
    ts%rank = ip
  ENDDO

  DO ip = 1, np
    CALL setup(ts(ip), paths(ip), err)
    THROW(err.NE.0)
  ENDDO
  DEALLOCATE(paths, stat = err)
  THROWM(err.NE.0,"Deallocation of paths failed")

  rttov_run_err = errorstatus_success

  DO ip = 1, np
    IF (err.NE.0) CYCLE

    IF (ts(ip)%defn%ltemp) THEN
      CALL run(ts(ip), ts(ip)%err, &
        ts(ip)%data%traj,    ts(ip)%data%traj_tl, &
        ts(ip)%data%traj_ad, ts(ip)%data%traj_k)
    ELSE
      CALL run(ts(ip), ts(ip)%err)
    ENDIF
    THROW_L(ts(ip)%err.NE.0,800)
    CYCLE
    CATCH_L(800)

    IF (ts(ip)%err .NE. 0) THEN
!$OMP FLUSH(err)
      rttov_run_err = errorstatus_fatal
    ENDIF
  ENDDO

  DO ip = 1, np
    CALL cleanup(ts(ip), err)
    THROW(err.NE.0)
  ENDDO

  DEALLOCATE(ts, stat = err)
  THROWM(err.NE.0,"Deallocation of ts")

  THROW(rttov_run_err.NE.0)

#ifdef _RTTOV_HDF
  CALL CLOSE_HDF(err)
  THROW(err.NE.0)
#endif

  IF (lhook) CALL dr_hook('rttov_test',1_jpim,zhook_handle)
  PCATCH
  IF (lhook) CALL dr_hook('rttov_test',1_jpim,zhook_handle)
  CALL rttov_exit(1_jpim)

  CONTAINS

  SUBROUTINE run(ts, err, traj, traj_tl, traj_ad, traj_k)
!
  TYPE(rttov_test_struct),    TARGET, INTENT(INOUT) :: ts
  INTEGER(jpim),                      INTENT(OUT)   :: err
  TYPE(rttov_traj), OPTIONAL, TARGET, INTENT(INOUT) :: traj, traj_tl, traj_ad, traj_k
!
  INTEGER(jpim) :: itime, i
  CHARACTER(LEN=80) :: msg
!
  TRY

  ! If we do the user checks here and throw any resulting errors then the cleanup routine
  ! will ensure all memory is deallocated
  IF (ts%defn%user_check_opts) THEN
    CALL rttov_user_options_checkinput(err, ts%defn%opts, ts%data%coefs)
    THROW(err.NE.0)
  ENDIF
  IF (ts%defn%user_check_prof) THEN
    DO i = 1, ts%defn%nprofiles
      CALL rttov_user_profile_checkinput(err, ts%defn%opts, ts%data%coefs, &
            ts%data%direct%profiles(i), ts%data%direct%aer_opt_param, ts%data%direct%cld_opt_param)
      THROW(err.NE.0)
    ENDDO
  ENDIF
  IF (ts%defn%print_opts) THEN
    IF (ts%defn%do_rttovscatt) THEN
      CALL rttov_print_opts_scatt(ts%defn%opts_scatt)
    ELSE
      CALL rttov_print_opts(ts%defn%opts)
    ENDIF
  ENDIF
  IF (ts%defn%print_profiles) THEN
    DO i = 1, ts%defn%nprofiles
      WRITE(msg,'(A,I6)') 'Profile', i
      CALL rttov_print_profile(ts%data%direct%profiles(i), text=TRIM(msg))
      IF (ts%defn%do_rttovscatt) CALL rttov_print_cld_profile(ts%data%direct%cld_profiles(i))
    ENDDO
  ENDIF

  DO itime = 1, ts%defn%ntimes

    IF (ts%defn%do_direct) THEN

      IF (ts%defn%nthreads > 0 .OR. ts%defn%prof_by_prof .OR. ts%defn%chan_by_chan) THEN

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%calc_emis_terms) THEN
            CALL rttov_parallel_scatt(                                      &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,        &
                ts%defn%nlevels,                 ts%data%chanprof,          &
                ts%data%frequencies,             ts%data%direct%profiles,   &
                ts%data%direct%cld_profiles,     ts%data%coefs,             &
                ts%data%coefs_scatt,             ts%data%calcemis,          &
                ts%data%direct%emissivity,       ts%data%direct%radiance,   &
                ts%data%direct%rttovscatt_cfrac, ts%data%direct%emis_terms, &
                nthreads = ts%defn%nthreads,                                &
                strategy = ts%defn%parallel_strategy)
          ELSEIF (ts%defn%radar) THEN
           CALL rttov_parallel_scatt(                                       &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,        &
                ts%defn%nlevels,                 ts%data%chanprof,          &
                ts%data%frequencies,             ts%data%direct%profiles,   &
                ts%data%direct%cld_profiles,     ts%data%coefs,             &
                ts%data%coefs_scatt,             ts%data%calcemis,          &
                ts%data%direct%emissivity,       ts%data%direct%radiance,   &
                reflectivity = ts%data%direct%reflectivity,                 &
                nthreads     = ts%defn%nthreads,                            &
                strategy     = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_scatt(                                      &
                ts%data%rttov_errorstatus,     ts%defn%opts_scatt,          &
                ts%defn%nlevels,               ts%data%chanprof,            &
                ts%data%frequencies,           ts%data%direct%profiles,     &
                ts%data%direct%cld_profiles,   ts%data%coefs,               &
                ts%data%coefs_scatt,           ts%data%calcemis,            &
                ts%data%direct%emissivity,     ts%data%direct%radiance,     &
                ts%data%direct%rttovscatt_cfrac,                            &
                nthreads = ts%defn%nthreads,                                &
                strategy = ts%defn%parallel_strategy)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_parallel_direct(                                        &
                ts%data%rttov_errorstatus,       ts%data%chanprof,             &
                ts%defn%opts,                    ts%data%direct%profiles,      &
                ts%data%coefs,                   ts%data%direct%transmission,  &
                ts%data%direct%radiance,         ts%data%direct%radiance2,     &
                ts%data%calcemis,                ts%data%direct%emissivity,    &
                ts%data%calcrefl,                ts%data%direct%reflectance,   &
                ts%data%direct%aer_opt_param,    ts%data%direct%cld_opt_param, &
                pccomp       = ts%data%direct%pccomp,                          &
                channels_rec = ts%data%channels_rec,                           &
                nthreads     = ts%defn%nthreads,                               &
                strategy     = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_direct(                                        &
                ts%data%rttov_errorstatus,     ts%data%chanprof,               &
                ts%defn%opts,                  ts%data%direct%profiles,        &
                ts%data%coefs,                 ts%data%direct%transmission,    &
                ts%data%direct%radiance,                                       &
                calcemis      = ts%data%calcemis,                              &
                emissivity    = ts%data%direct%emissivity,                     &
                calcrefl      = ts%data%calcrefl,                              &
                reflectance   = ts%data%direct%reflectance,                    &
                aer_opt_param = ts%data%direct%aer_opt_param,                  &
                cld_opt_param = ts%data%direct%cld_opt_param,                  &
                pccomp        = ts%data%direct%pccomp,                         &
                channels_rec  = ts%data%channels_rec,                          &
                nthreads      = ts%defn%nthreads,                              &
                strategy      = ts%defn%parallel_strategy)
          ENDIF
        ENDIF

      ELSE ! parallel or standard interface

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%calc_emis_terms) THEN
            CALL rttov_scatt(                                                 &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,          &
                ts%defn%nlevels,                 ts%data%chanprof,            &
                ts%data%frequencies,             ts%data%direct%profiles,     &
                ts%data%direct%cld_profiles,     ts%data%coefs,               &
                ts%data%coefs_scatt,             ts%data%calcemis,            &
                ts%data%direct%emissivity,       ts%data%direct%radiance,     &
                ts%data%direct%rttovscatt_cfrac, ts%data%direct%emis_terms)
          ELSEIF (ts%defn%radar) THEN
            CALL rttov_scatt(                                                 &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,          &
                ts%defn%nlevels,                 ts%data%chanprof,            &
                ts%data%frequencies,             ts%data%direct%profiles,     &
                ts%data%direct%cld_profiles,     ts%data%coefs,               &
                ts%data%coefs_scatt,             ts%data%calcemis,            &
                ts%data%direct%emissivity,       ts%data%direct%radiance,     &
                reflectivity = ts%data%direct%reflectivity)
          ELSE
            CALL rttov_scatt(                                                 &
                ts%data%rttov_errorstatus,     ts%defn%opts_scatt,            &
                ts%defn%nlevels,               ts%data%chanprof,              &
                ts%data%frequencies,           ts%data%direct%profiles,       &
                ts%data%direct%cld_profiles,   ts%data%coefs,                 &
                ts%data%coefs_scatt,           ts%data%calcemis,              &
                ts%data%direct%emissivity,     ts%data%direct%radiance,       &
                ts%data%direct%rttovscatt_cfrac)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_direct(                                                 &
                ts%data%rttov_errorstatus,       ts%data%chanprof,             &
                ts%defn%opts,                    ts%data%direct%profiles,      &
                ts%data%coefs,                   ts%data%direct%transmission,  &
                ts%data%direct%radiance,         ts%data%direct%radiance2,     &
                ts%data%calcemis,                ts%data%direct%emissivity,    &
                ts%data%calcrefl,                ts%data%direct%reflectance,   &
                ts%data%direct%aer_opt_param,    ts%data%direct%cld_opt_param, &
                traj         = traj,                                           &
                pccomp       = ts%data%direct%pccomp,                          &
                channels_rec = ts%data%channels_rec)
          ELSE
            CALL rttov_direct(                                                 &
                ts%data%rttov_errorstatus,       ts%data%chanprof,             &
                ts%defn%opts,                    ts%data%direct%profiles,      &
                ts%data%coefs,                   ts%data%direct%transmission,  &
                ts%data%direct%radiance,                                       &
                calcemis      = ts%data%calcemis,                              &
                emissivity    = ts%data%direct%emissivity,                     &
                calcrefl      = ts%data%calcrefl,                              &
                reflectance   = ts%data%direct%reflectance,                    &
                aer_opt_param = ts%data%direct%aer_opt_param,                  &
                cld_opt_param = ts%data%direct%cld_opt_param,                  &
                traj          = traj,                                          &
                pccomp        = ts%data%direct%pccomp,                         &
                channels_rec  = ts%data%channels_rec)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ENDIF ! parallel or standard interface

      IF (ts%data%rttov_errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.NE.0)
      ENDIF

      IF (ts%defn%do_print) THEN
        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/direct", "", .FALSE._jplm, &
                         ts%defn%calc_rad2, ts%defn%do_rttovscatt, ts%defn%radar, ts%defn%calc_emis_terms, ts%data%direct)
        THROW(err.NE.0)
      ENDIF

#ifdef _RTTOV_HDF
      IF (ts%defn%savehdf5) THEN
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/options.H5", '/OPTIONS', CREATE=.TRUE., &
          OPTIONS = ts%defn%opts)
       THROW(err.NE.0)
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/chanprof.H5", '/CHANPROF', CREATE=.TRUE., &
          CHANPROF = ts%data%chanprof)
       THROW(err.NE.0)
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/radiance.H5", '/RADIANCE', CREATE=.TRUE., &
          RADIANCE = ts%data%direct%radiance)
       THROW(err.NE.0)
       IF (ts%defn%calc_rad2) THEN
         CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/radiance2.H5", '/RADIANCE2', CREATE=.TRUE., &
            RADIANCE2 = ts%data%direct%radiance2)
         THROW(err.NE.0)
       ENDIF
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/transmission.H5", '/TRANSMISSION', CREATE=.TRUE., &
          TRANSMISSION = ts%data%direct%transmission)
       THROW(err.NE.0)
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/emissivity.H5", '/EMISSIVITY', CREATE=.TRUE., &
          EMISSIVITY = ts%data%direct%emissivity)
       THROW(err.NE.0)
       CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/emissivity.H5", '/EMISSIVITY', CREATE=.FALSE., &
          l1 = ts%data%calcemis, SNAME='CALCEMIS')
       THROW(err.NE.0)
       IF (ts%defn%opts%rt_ir%addsolar) THEN
         CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/reflectance.H5", '/REFLECTANCE', CREATE=.TRUE., &
            REFLECTANCE = ts%data%direct%reflectance)
         THROW(err.NE.0)
         CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/direct/reflectance.H5", '/REFLECTANCE', CREATE=.FALSE., &
            l1 = ts%data%calcrefl, SNAME='CALCREFL')
         THROW(err.NE.0)
       ENDIF
      ENDIF
#endif

      ts%defn%opts%config%verbose = .FALSE.

    ENDIF

    IF (ts%defn%do_tl) THEN

      IF (ts%defn%nthreads > 0 .OR. ts%defn%prof_by_prof .OR. ts%defn%chan_by_chan) THEN

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_parallel_scatt_tl(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%tl%profiles,      &
                ts%data%tl%cld_profiles,         ts%data%tl%emissivity,    &
                ts%data%direct%radiance,         ts%data%tl%radiance,      &
                reflectivity    = ts%data%direct%reflectivity,             &
                reflectivity_tl = ts%data%tl%reflectivity,                 &
                nthreads        = ts%defn%nthreads,                        &
                strategy        = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_scatt_tl(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%tl%profiles,      &
                ts%data%tl%cld_profiles,         ts%data%tl%emissivity,    &
                ts%data%direct%radiance,         ts%data%tl%radiance,      &
                nthreads = ts%defn%nthreads,                               &
                strategy = ts%defn%parallel_strategy)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_parallel_tl(                                        &
                ts%data%rttov_errorstatus,      ts%data%chanprof,          &
                ts%defn%opts,                   ts%data%direct%profiles,   &
                ts%data%tl%profiles,            ts%data%coefs,             &
                ts%data%direct%transmission,    ts%data%tl%transmission,   &
                ts%data%direct%radiance,        ts%data%tl%radiance,       &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_tl    = ts%data%tl%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_tl   = ts%data%tl%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_tl = ts%data%tl%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_tl = ts%data%tl%cld_opt_param,               &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_tl        = ts%data%tl%pccomp,                      &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_tl(                                        &
                ts%data%rttov_errorstatus,      ts%data%chanprof,          &
                ts%defn%opts,                   ts%data%direct%profiles,   &
                ts%data%tl%profiles,            ts%data%coefs,             &
                ts%data%direct%transmission,    ts%data%tl%transmission,   &
                ts%data%direct%radiance,        ts%data%tl%radiance,       &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_tl    = ts%data%tl%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_tl   = ts%data%tl%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_tl = ts%data%tl%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_tl = ts%data%tl%cld_opt_param,               &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_tl        = ts%data%tl%pccomp,                      &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ELSE ! parallel or standard interface

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_scatt_tl(                                           &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%tl%profiles,      &
                ts%data%tl%cld_profiles,         ts%data%tl%emissivity,    &
                ts%data%direct%radiance,         ts%data%tl%radiance,      &
                reflectivity    = ts%data%direct%reflectivity,             &
                reflectivity_tl = ts%data%tl%reflectivity)

          ELSE
            CALL rttov_scatt_tl(                                           &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%tl%profiles,      &
                ts%data%tl%cld_profiles,         ts%data%tl%emissivity,    &
                ts%data%direct%radiance,         ts%data%tl%radiance)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_tl(                                                 &
                ts%data%rttov_errorstatus,      ts%data%chanprof,          &
                ts%defn%opts,                   ts%data%direct%profiles,   &
                ts%data%tl%profiles,            ts%data%coefs,             &
                ts%data%direct%transmission,    ts%data%tl%transmission,   &
                ts%data%direct%radiance,        ts%data%tl%radiance,       &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_tl    = ts%data%tl%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_tl   = ts%data%tl%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_tl = ts%data%tl%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_tl = ts%data%tl%cld_opt_param,               &
                traj             = traj,                                   &
                traj_tl          = traj_tl,                                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_tl        = ts%data%tl%pccomp,                      &
                channels_rec     = ts%data%channels_rec)
          ELSE
            CALL rttov_tl(                                                 &
                ts%data%rttov_errorstatus,      ts%data%chanprof,          &
                ts%defn%opts,                   ts%data%direct%profiles,   &
                ts%data%tl%profiles,            ts%data%coefs,             &
                ts%data%direct%transmission,    ts%data%tl%transmission,   &
                ts%data%direct%radiance,        ts%data%tl%radiance,       &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_tl    = ts%data%tl%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_tl   = ts%data%tl%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_tl = ts%data%tl%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_tl = ts%data%tl%cld_opt_param,               &
                traj             = traj,                                   &
                traj_tl          = traj_tl,                                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_tl        = ts%data%tl%pccomp,                      &
                channels_rec     = ts%data%channels_rec)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ENDIF ! parallel or standard interface

      IF (ts%data%rttov_errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.NE.0)
      ENDIF

      IF (ts%defn%do_print) THEN
        IF (ts%defn%scale_out /= 1.0_jprb) THEN
          ts%data%tl%radiance%total      = ts%data%tl%radiance%total * ts%defn%scale_out
          ts%data%tl%radiance%bt         = ts%data%tl%radiance%bt * ts%defn%scale_out
          ts%data%tl%radiance%refl       = ts%data%tl%radiance%refl * ts%defn%scale_out
          ts%data%tl%radiance%clear      = ts%data%tl%radiance%clear * ts%defn%scale_out
          ts%data%tl%radiance%bt_clear   = ts%data%tl%radiance%bt_clear * ts%defn%scale_out
          ts%data%tl%radiance%refl_clear = ts%data%tl%radiance%refl_clear * ts%defn%scale_out
          ts%data%tl%radiance%overcast   = ts%data%tl%radiance%overcast * ts%defn%scale_out
          ts%data%tl%radiance%cloudy     = ts%data%tl%radiance%cloudy * ts%defn%scale_out

          ts%data%tl%transmission%tau_total           = ts%data%tl%transmission%tau_total * ts%defn%scale_out
          ts%data%tl%transmission%tau_levels          = ts%data%tl%transmission%tau_levels * ts%defn%scale_out
          ts%data%tl%transmission%tausun_total_path1  = ts%data%tl%transmission%tausun_total_path1 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_levels_path1 = ts%data%tl%transmission%tausun_levels_path1 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_total_path2  = ts%data%tl%transmission%tausun_total_path2 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_levels_path2 = ts%data%tl%transmission%tausun_levels_path2 * ts%defn%scale_out
          ts%data%tl%transmission%tau_total_cld       = ts%data%tl%transmission%tau_total_cld * ts%defn%scale_out
          ts%data%tl%transmission%tau_levels_cld      = ts%data%tl%transmission%tau_levels_cld * ts%defn%scale_out

          ts%data%tl%emissivity%emis_out          = ts%data%tl%emissivity%emis_out * ts%defn%scale_out
          ts%data%tl%emissivity%specularity       = ts%data%tl%emissivity%specularity * ts%defn%scale_out
          ts%data%tl%reflectance%refl_out         = ts%data%tl%reflectance%refl_out * ts%defn%scale_out
          ts%data%tl%reflectance%diffuse_refl_out = ts%data%tl%reflectance%diffuse_refl_out * ts%defn%scale_out

          IF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
            ts%data%tl%reflectivity%zef  = ts%data%tl%reflectivity%zef * ts%defn%scale_out
            ts%data%tl%reflectivity%azef = ts%data%tl%reflectivity%azef * ts%defn%scale_out
          ENDIF

          IF (ts%defn%opts%rt_ir%pc%addpc) THEN
            IF (ASSOCIATED(ts%data%tl%pccomp%total_pcscores)) &
               ts%data%tl%pccomp%total_pcscores    = ts%data%tl%pccomp%total_pcscores * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%total_pccomp)) &
               ts%data%tl%pccomp%total_pccomp      = ts%data%tl%pccomp%total_pccomp * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%bt_pccomp)) &
               ts%data%tl%pccomp%bt_pccomp         = ts%data%tl%pccomp%bt_pccomp * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%clear_pcscores)) &
               ts%data%tl%pccomp%clear_pcscores    = ts%data%tl%pccomp%clear_pcscores * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%clear_pccomp)) &
               ts%data%tl%pccomp%clear_pccomp      = ts%data%tl%pccomp%clear_pccomp * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%bt_clear_pccomp)) &
               ts%data%tl%pccomp%bt_clear_pccomp   = ts%data%tl%pccomp%bt_clear_pccomp * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%overcast_pcscores)) &
               ts%data%tl%pccomp%overcast_pcscores = ts%data%tl%pccomp%overcast_pcscores * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%cloudy_pcscores)) &
               ts%data%tl%pccomp%cloudy_pcscores   = ts%data%tl%pccomp%cloudy_pcscores * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%overcast_pccomp)) &
               ts%data%tl%pccomp%overcast_pccomp   = ts%data%tl%pccomp%overcast_pccomp * ts%defn%scale_out
            IF (ASSOCIATED(ts%data%tl%pccomp%cloudy_pccomp)) &
               ts%data%tl%pccomp%cloudy_pccomp     = ts%data%tl%pccomp%cloudy_pccomp * ts%defn%scale_out
          ENDIF
        ENDIF

        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/tl", "", .FALSE._jplm, &
                         ts%defn%calc_rad2, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
        THROW(err.NE.0)
        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/tl", "_tl", .FALSE._jplm, &
                         .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%tl)
        THROW(err.NE.0)
      ENDIF

      ts%defn%opts%config%verbose = .FALSE.

    ENDIF

    IF (ts%defn%do_ad) THEN

      ! Initialise the adjoint variables
      CALL rttov_init_prof(ts%data%ad%profiles)
      IF (ts%defn%opts%rt_ir%user_aer_opt_param) &
        CALL rttov_init_opt_param(err, ts%defn%opts, ts%data%ad%aer_opt_param, zero_only=.TRUE._jplm)
      IF (ts%defn%opts%rt_ir%user_cld_opt_param) &
        CALL rttov_init_opt_param(err, ts%defn%opts, ts%data%ad%cld_opt_param, zero_only=.TRUE._jplm)
      IF (ts%defn%do_rttovscatt) CALL rttov_init_scatt_prof(ts%data%ad%cld_profiles)
      CALL rttov_init_emis_refl(ts%data%ad%emissivity, ts%data%ad%reflectance)
      CALL rttov_init_rad(ts%data%ad%radiance)
      IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        CALL rttov_copy_pccomp(ts%data%ad%pccomp, ts%data%ad%pccomp_saved)
      ELSEIF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
        ts%data%ad%reflectivity%zef(:,:)  = ts%data%ad%reflectivity_saved%zef(:,:)
        ts%data%ad%reflectivity%azef(:,:) = ts%data%ad%reflectivity_saved%azef(:,:)
      ELSE
        ts%data%ad%radiance%total(:) = ts%data%ad%radiance_saved%total(:)
        ts%data%ad%radiance%bt(:)    = ts%data%ad%radiance_saved%bt(:)
      ENDIF

      IF (ts%defn%nthreads > 0 .OR. ts%defn%prof_by_prof .OR. ts%defn%chan_by_chan) THEN

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_parallel_scatt_ad(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%ad%profiles,      &
                ts%data%ad%cld_profiles,         ts%data%ad%emissivity,    &
                ts%data%direct%radiance,         ts%data%ad%radiance,      &
                reflectivity    = ts%data%direct%reflectivity,             &
                reflectivity_ad = ts%data%ad%reflectivity,                 &
                nthreads        = ts%defn%nthreads,                        &
                strategy        = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_scatt_ad(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%ad%profiles,      &
                ts%data%ad%cld_profiles,         ts%data%ad%emissivity,    &
                ts%data%direct%radiance,         ts%data%ad%radiance,      &
                nthreads = ts%defn%nthreads,                               &
                strategy = ts%defn%parallel_strategy)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_parallel_ad(                                        &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%ad%profiles,           ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%ad%transmission,    &
                ts%data%direct%radiance,       ts%data%ad%radiance,        &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_ad    = ts%data%ad%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_ad   = ts%data%ad%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_ad = ts%data%ad%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_ad = ts%data%ad%cld_opt_param,               &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_ad        = ts%data%ad%pccomp,                      &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_ad(                                        &
                ts%data%rttov_errorstatus,       ts%data%chanprof,         &
                ts%defn%opts,                    ts%data%direct%profiles,  &
                ts%data%ad%profiles,             ts%data%coefs,            &
                ts%data%direct%transmission,     ts%data%ad%transmission,  &
                ts%data%direct%radiance,         ts%data%ad%radiance,      &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_ad    = ts%data%ad%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_ad   = ts%data%ad%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_ad = ts%data%ad%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_ad = ts%data%ad%cld_opt_param,               &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_ad        = ts%data%ad%pccomp,                      &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ELSE ! parallel or standard interface

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_scatt_ad(                                           &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%ad%profiles,      &
                ts%data%ad%cld_profiles,         ts%data%ad%emissivity,    &
                ts%data%direct%radiance,         ts%data%ad%radiance,      &
                reflectivity    = ts%data%direct%reflectivity,             &
                reflectivity_ad = ts%data%ad%reflectivity)
          ELSE
            CALL rttov_scatt_ad(                                           &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%ad%profiles,      &
                ts%data%ad%cld_profiles,         ts%data%ad%emissivity,    &
                ts%data%direct%radiance,         ts%data%ad%radiance)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_ad(                                                 &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%ad%profiles,           ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%ad%transmission,    &
                ts%data%direct%radiance,       ts%data%ad%radiance,        &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_ad    = ts%data%ad%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_ad   = ts%data%ad%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_ad = ts%data%ad%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_ad = ts%data%ad%cld_opt_param,               &
                traj             = traj,                                   &
                traj_ad          = traj_ad,                                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_ad        = ts%data%ad%pccomp,                      &
                channels_rec     = ts%data%channels_rec)
          ELSE
            CALL rttov_ad(                                                 &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%ad%profiles,           ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%ad%transmission,    &
                ts%data%direct%radiance,       ts%data%ad%radiance,        &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_ad    = ts%data%ad%emissivity,                  &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_ad   = ts%data%ad%reflectance,                 &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_ad = ts%data%ad%aer_opt_param,               &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_ad = ts%data%ad%cld_opt_param,               &
                traj             = traj,                                   &
                traj_ad          = traj_ad,                                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_ad        = ts%data%ad%pccomp,                      &
                channels_rec     = ts%data%channels_rec)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ENDIF ! parallel or standard interface

      IF (ts%data%rttov_errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.NE.0)
      ENDIF

      IF (ts%defn%do_print) THEN
        IF (ts%defn%scale_out /= 1.0_jprb) THEN
          CALL rttov_scale_profile_inc(ts%data%ad%profiles, ts%defn%scale_out)
          IF (ts%defn%do_rttovscatt) CALL rttov_scale_cld_profile_inc(ts%data%ad%cld_profiles, ts%defn%scale_out)
          IF (ts%defn%opts%rt_ir%user_aer_opt_param) &
            CALL rttov_scale_opt_param_inc(ts%data%ad%aer_opt_param, ts%defn%scale_out)
          IF (ts%defn%opts%rt_ir%user_cld_opt_param) &
            CALL rttov_scale_opt_param_inc(ts%data%ad%cld_opt_param, ts%defn%scale_out)

          ts%data%ad%emissivity%emis_in = ts%data%ad%emissivity%emis_in * ts%defn%scale_out
          ts%data%ad%emissivity%specularity = ts%data%ad%emissivity%specularity * ts%defn%scale_out
          ts%data%ad%reflectance%refl_in = ts%data%ad%reflectance%refl_in * ts%defn%scale_out
          ts%data%ad%reflectance%diffuse_refl_in = ts%data%ad%reflectance%diffuse_refl_in * ts%defn%scale_out
        ENDIF

        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/ad", "", .FALSE._jplm, &
                         ts%defn%calc_rad2, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
        THROW(err.NE.0)
        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/ad", "_ad", .TRUE._jplm, &
                         .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%ad)
        THROW(err.NE.0)
      ENDIF

      ts%defn%opts%config%verbose = .FALSE.

    ENDIF

    IF (ts%defn%do_k) THEN

      CALL rttov_init_prof(ts%data%k%profiles)
      IF (ts%defn%opts%rt_ir%user_aer_opt_param) &
        CALL rttov_init_opt_param(err, ts%defn%opts, ts%data%k%aer_opt_param, zero_only=.TRUE._jplm)
      IF (ts%defn%opts%rt_ir%user_cld_opt_param) &
        CALL rttov_init_opt_param(err, ts%defn%opts, ts%data%k%cld_opt_param, zero_only=.TRUE._jplm)
      IF (ts%defn%do_rttovscatt) CALL rttov_init_scatt_prof(ts%data%k%cld_profiles)
      CALL rttov_init_emis_refl(ts%data%k%emissivity, ts%data%k%reflectance)
      CALL rttov_init_rad(ts%data%k%radiance)
      IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          IF (ts%defn%opts%rt_all%switchrad) THEN
            ts%data%k%pccomp%bt_pccomp = 1._jprb
          ELSE
            ts%data%k%pccomp%total_pccomp = 1._jprb
          ENDIF
        ELSE
          ts%data%k%pccomp%total_pcscores = 1._jprb
        ENDIF
      ELSEIF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
        ts%data%k%reflectivity%zef(:,:) = 0._jprb
        ts%data%k%reflectivity%azef(:,:) = 0._jprb
        ts%data%k%reflectivity%azef(radar_k_lev,:) = 1._jprb
      ELSE
        ts%data%k%radiance%total = 1._jprb ! For any pure-solar channels, increment is always in radiance
        IF (ts%defn%opts%rt_all%switchrad .OR. ts%defn%do_rttovscatt) ts%data%k%radiance%bt = 1._jprb
      ENDIF

      IF (ts%defn%nthreads > 0 .OR. ts%defn%prof_by_prof .OR. ts%defn%chan_by_chan) THEN

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_parallel_scatt_ad(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%k%profiles,       &
                ts%data%k%cld_profiles,          ts%data%k%emissivity,     &
                ts%data%direct%radiance,         ts%data%k%radiance,       &
                reflectivity    = ts%data%direct%reflectivity,             &
                reflectivity_ad = ts%data%k%reflectivity,                  &
                nthreads        = ts%defn%nthreads,                        &
                strategy        = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_scatt_ad(                                  &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,       &
                ts%defn%nlevels,                 ts%data%chanprof,         &
                ts%data%frequencies,             ts%data%direct%profiles,  &
                ts%data%direct%cld_profiles,     ts%data%coefs,            &
                ts%data%coefs_scatt,             ts%data%calcemis,         &
                ts%data%direct%emissivity,       ts%data%k%profiles,       &
                ts%data%k%cld_profiles,          ts%data%k%emissivity,     &
                ts%data%direct%radiance,         ts%data%k%radiance,       &
                nthreads = ts%defn%nthreads,                               &
                strategy = ts%defn%parallel_strategy)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_parallel_k(                                         &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%k%profiles,            ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%k%transmission,     &
                ts%data%direct%radiance,       ts%data%k%radiance,         &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_k     = ts%data%k%emissivity,                   &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_k    = ts%data%k%reflectance,                  &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_k  = ts%data%k%aer_opt_param,                &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_k  = ts%data%k%cld_opt_param,                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_k         = ts%data%k%pccomp,                       &
                profiles_k_rec   = ts%data%data_k%profiles_k_rec,          &
                profiles_k_pc    = ts%data%data_k%profiles_k_pc,           &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ELSE
            CALL rttov_parallel_k(                                         &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%k%profiles,            ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%k%transmission,     &
                ts%data%direct%radiance,       ts%data%k%radiance,         &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_k     = ts%data%k%emissivity,                   &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_k    = ts%data%k%reflectance,                  &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_k  = ts%data%k%aer_opt_param,                &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_k  = ts%data%k%cld_opt_param,                &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_k         = ts%data%k%pccomp,                       &
                profiles_k_rec   = ts%data%data_k%profiles_k_rec,          &
                profiles_k_pc    = ts%data%data_k%profiles_k_pc,           &
                channels_rec     = ts%data%channels_rec,                   &
                nthreads         = ts%defn%nthreads,                       &
                strategy         = ts%defn%parallel_strategy)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ELSE ! parallel or standard interface

        IF (ts%defn%do_rttovscatt) THEN

          IF (ts%defn%radar) THEN
            CALL rttov_scatt_ad(                                          &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,      &
                ts%defn%nlevels,                 ts%data%chanprof,        &
                ts%data%frequencies,             ts%data%direct%profiles, &
                ts%data%direct%cld_profiles,     ts%data%coefs,           &
                ts%data%coefs_scatt,             ts%data%calcemis,        &
                ts%data%direct%emissivity,       ts%data%k%profiles,      &
                ts%data%k%cld_profiles,          ts%data%k%emissivity,    &
                ts%data%direct%radiance,         ts%data%k%radiance,      &
                reflectivity    = ts%data%direct%reflectivity,            &
                reflectivity_ad = ts%data%k%reflectivity)
          ELSE
            CALL rttov_scatt_ad(                                          &
                ts%data%rttov_errorstatus,       ts%defn%opts_scatt,      &
                ts%defn%nlevels,                 ts%data%chanprof,        &
                ts%data%frequencies,             ts%data%direct%profiles, &
                ts%data%direct%cld_profiles,     ts%data%coefs,           &
                ts%data%coefs_scatt,             ts%data%calcemis,        &
                ts%data%direct%emissivity,       ts%data%k%profiles,      &
                ts%data%k%cld_profiles,          ts%data%k%emissivity,    &
                ts%data%direct%radiance,         ts%data%k%radiance)
          ENDIF

        ELSE ! do_rttovscatt or not

          IF (ts%defn%calc_rad2) THEN
            CALL rttov_k(                                                  &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%k%profiles,            ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%k%transmission,     &
                ts%data%direct%radiance,       ts%data%k%radiance,         &
                radiance2        = ts%data%direct%radiance2,               &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_k     = ts%data%k%emissivity,                   &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_k    = ts%data%k%reflectance,                  &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_k  = ts%data%k%aer_opt_param,                &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_k  = ts%data%k%cld_opt_param,                &
                traj             = traj,                                   &
                traj_k           = traj_k,                                 &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_k         = ts%data%k%pccomp,                       &
                profiles_k_rec   = ts%data%data_k%profiles_k_rec,          &
                profiles_k_pc    = ts%data%data_k%profiles_k_pc,           &
                channels_rec     = ts%data%channels_rec)
          ELSE
            CALL rttov_k(                                                  &
                ts%data%rttov_errorstatus,     ts%data%chanprof,           &
                ts%defn%opts,                  ts%data%direct%profiles,    &
                ts%data%k%profiles,            ts%data%coefs,              &
                ts%data%direct%transmission,   ts%data%k%transmission,     &
                ts%data%direct%radiance,       ts%data%k%radiance,         &
                calcemis         = ts%data%calcemis,                       &
                emissivity       = ts%data%direct%emissivity,              &
                emissivity_k     = ts%data%k%emissivity,                   &
                calcrefl         = ts%data%calcrefl,                       &
                reflectance      = ts%data%direct%reflectance,             &
                reflectance_k    = ts%data%k%reflectance,                  &
                aer_opt_param    = ts%data%direct%aer_opt_param,           &
                aer_opt_param_k  = ts%data%k%aer_opt_param,                &
                cld_opt_param    = ts%data%direct%cld_opt_param,           &
                cld_opt_param_k  = ts%data%k%cld_opt_param,                &
                traj             = traj,                                   &
                traj_k           = traj_k,                                 &
                pccomp           = ts%data%direct%pccomp,                  &
                pccomp_k         = ts%data%k%pccomp,                       &
                profiles_k_rec   = ts%data%data_k%profiles_k_rec,          &
                profiles_k_pc    = ts%data%data_k%profiles_k_pc,           &
                channels_rec     = ts%data%channels_rec)
          ENDIF

        ENDIF ! do_rttovscatt or not

      ENDIF ! parallel or standard interface

      IF (ts%data%rttov_errorstatus /= 0) THEN
        err = errorstatus_fatal
        THROW(err.NE.0)
      ENDIF

      IF (ts%defn%do_print) THEN
        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k", "", .FALSE._jplm, &
                         ts%defn%calc_rad2, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
        THROW(err.NE.0)
        CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k", "_k", .TRUE._jplm, &
                         .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%k, ts%data%data_k)
        THROW(err.NE.0)
      ENDIF

#ifdef _RTTOV_HDF
      IF (ts%defn%savehdf5) THEN
        INFO("save kmatrix")
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/KMATRIX', CREATE=.TRUE., &
           KMATRIX = ts%data%k%profiles, OPTS= ts%defn%opts)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/MISC', CREATE=.FALSE., &
           R1 = ts%data%coefs%coef%ff_cwn(ts%data%chanprof(:)%chan), SNAME='WAVENUMBERS', UNITS='cm-1')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/MISC', CREATE=.FALSE., &
           c0 = ts%data%coefs%coef%ID_COMMON_NAME, SNAME='ID_COMMON_NAME')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/MISC', CREATE=.FALSE., &
           c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/PROFILES', CREATE=.FALSE., &
           PROFILES = ts%data%direct%profiles)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/CHANPROF', CREATE=.FALSE., &
           CHANPROF = ts%data%chanprof)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/EMISSIVITY_K', CREATE=.FALSE., &
           EMISSIVITY = ts%data%k%emissivity)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/EMISSIVITY_K', CREATE=.FALSE., &
           l1 = ts%data%calcemis, SNAME='CALCEMIS')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/OPTIONS', CREATE=.FALSE., &
           OPTIONS = ts%defn%opts)
        THROW(err.NE.0)
        IF (ts%defn%opts%rt_ir%addsolar) THEN
          CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/REFLECTANCE_K', CREATE=.FALSE., &
             REFLECTANCE = ts%data%k%reflectance)
          THROW(err.NE.0)
          CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/k/kmatrix.H5", '/REFLECTANCE_K', CREATE=.FALSE., &
           l1 = ts%data%calcrefl, SNAME='CALCREFL')
        THROW(err.NE.0)
        ENDIF
        INFO("fin save")
      ENDIF
#endif

      ts%defn%opts%config%verbose = .FALSE.

    ENDIF

    IF (ts%defn%print_quality) THEN
      DO i = 1, SIZE(ts%data%chanprof)
        WRITE(msg,'(A,I6,A,I6)') 'Profile', ts%data%chanprof(i)%prof, &
                                 ', channel', ts%data%chanprof(i)%chan
        CALL rttov_print_radiance_quality(ts%data%direct%radiance%quality(i), text=TRIM(msg))
      ENDDO
    ENDIF

  ENDDO

  IF (ts%defn%do_k_bf) THEN

    CALL rttov_k_bf(                                       &
        ts%data%rttov_errorstatus,                         &
        ts%data%chanprof,                                  &
        ts%data%frequencies,                               &
        ts%defn%opts,                                      &
        ts%defn%opts_scatt,                                &
        ts%data%direct%profiles,                           &
        ts%data%k_bf%profiles,                             &
        ts%data%direct%cld_profiles,                       &
        ts%data%k_bf%cld_profiles,                         &
        ts%data%coefs,                                     &
        ts%data%coefs_scatt,                               &
        ts%data%direct%transmission,                       &
        ts%data%direct%radiance,                           &
        ts%data%direct%reflectivity,                       &
        ts%data%calcemis,                                  &
        ts%data%direct%emissivity,                         &
        ts%data%k_bf%emissivity,                           &
        ts%data%calcrefl,                                  &
        ts%data%direct%reflectance,                        &
        ts%data%k_bf%reflectance,                          &
        ts%data%direct%aer_opt_param,                      &
        ts%data%k_bf%aer_opt_param,                        &
        ts%data%direct%cld_opt_param,                      &
        ts%data%k_bf%cld_opt_param,                        &
        pccomp         = ts%data%direct%pccomp,            &
        profiles_k_rec = ts%data%data_k_bf%profiles_k_rec, &
        profiles_k_pc  = ts%data%data_k_bf%profiles_k_pc,  &
        channels_rec   = ts%data%channels_rec,             &
        nthreads       = ts%defn%nthreads)

    IF (ts%data%rttov_errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_print) THEN
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_bf", "", .FALSE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
      THROW(err.NE.0)
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_bf", "_k", .TRUE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%k_bf, ts%data%data_k_bf)
      THROW(err.NE.0)
    ENDIF

    ts%defn%opts%config%verbose = .FALSE.

  ENDIF

  IF (ts%defn%do_k_tl) THEN

    CALL rttov_k_tl(                                       &
        ts%data%rttov_errorstatus,                         &
        ts%data%chanprof,                                  &
        ts%data%frequencies,                               &
        ts%defn%opts,                                      &
        ts%defn%opts_scatt,                                &
        ts%data%direct%profiles,                           &
        ts%data%k_tl%profiles,                             &
        ts%data%direct%cld_profiles,                       &
        ts%data%k_tl%cld_profiles,                         &
        ts%data%coefs,                                     &
        ts%data%coefs_scatt,                               &
        ts%data%direct%transmission,                       &
        ts%data%direct%radiance,                           &
        ts%data%direct%reflectivity,                       &
        ts%data%calcemis,                                  &
        ts%data%direct%emissivity,                         &
        ts%data%k_tl%emissivity,                           &
        ts%data%calcrefl,                                  &
        ts%data%direct%reflectance,                        &
        ts%data%k_tl%reflectance,                          &
        ts%data%direct%aer_opt_param,                      &
        ts%data%k_tl%aer_opt_param,                        &
        ts%data%direct%cld_opt_param,                      &
        ts%data%k_tl%cld_opt_param,                        &
        pccomp         = ts%data%direct%pccomp,            &
        profiles_k_rec = ts%data%data_k_tl%profiles_k_rec, &
        profiles_k_pc  = ts%data%data_k_tl%profiles_k_pc,  &
        channels_rec   = ts%data%channels_rec,             &
        nthreads       = ts%defn%nthreads)

    IF (ts%data%rttov_errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_print) THEN
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_tl", "", .FALSE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
      THROW(err.NE.0)
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_tl", "_k", .TRUE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%k_tl, ts%data%data_k_tl)
      THROW(err.NE.0)
    ENDIF

    ts%defn%opts%config%verbose = .FALSE.

  ENDIF

  IF (ts%defn%do_k_ad) THEN

    CALL rttov_k_ad(                                       &
        ts%data%rttov_errorstatus,                         &
        ts%data%chanprof,                                  &
        ts%data%frequencies,                               &
        ts%defn%opts,                                      &
        ts%defn%opts_scatt,                                &
        ts%data%direct%profiles,                           &
        ts%data%k_ad%profiles,                             &
        ts%data%direct%cld_profiles,                       &
        ts%data%k_ad%cld_profiles,                         &
        ts%data%coefs,                                     &
        ts%data%coefs_scatt,                               &
        ts%data%direct%transmission,                       &
        ts%data%direct%radiance,                           &
        ts%data%direct%reflectivity,                       &
        ts%data%calcemis,                                  &
        ts%data%direct%emissivity,                         &
        ts%data%k_ad%emissivity,                           &
        ts%data%calcrefl,                                  &
        ts%data%direct%reflectance,                        &
        ts%data%k_ad%reflectance,                          &
        ts%data%direct%aer_opt_param,                      &
        ts%data%k_ad%aer_opt_param,                        &
        ts%data%direct%cld_opt_param,                      &
        ts%data%k_ad%cld_opt_param,                        &
        pccomp         = ts%data%direct%pccomp,            &
        profiles_k_rec = ts%data%data_k_ad%profiles_k_rec, &
        profiles_k_pc  = ts%data%data_k_ad%profiles_k_pc,  &
        channels_rec   = ts%data%channels_rec,             &
        nthreads       = ts%defn%nthreads)

    IF (ts%data%rttov_errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_print) THEN
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_ad", "", .FALSE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%direct)
      THROW(err.NE.0)
      CALL printoutput(err, ts%defn%opts, TRIM(ts%path)//"/k_ad", "_k", .TRUE._jplm, &
                       .FALSE._jplm, ts%defn%do_rttovscatt, ts%defn%radar, .FALSE._jplm, ts%data%k_ad, ts%data%data_k_ad)
      THROW(err.NE.0)
    ENDIF

    ts%defn%opts%config%verbose = .FALSE.

  ENDIF

  IF (ts%defn%do_taylor) THEN
    CALL rttov_taylor_test(ts%data%rttov_errorstatus,     &
                           ts%path,                       &
                           ts%defn%opts,                  &
                           ts%defn%opts_scatt,            &
                           ts%data%coefs,                 &
                           ts%data%coefs_scatt,           &
                           ts%data%chanprof,              &
                           ts%data%frequencies,           &
                           ts%defn%npcscores,             &
                           ts%data%channels_rec,          &
                           ts%data%direct%profiles,       &
                           ts%data%direct%cld_profiles,   &
                           ts%data%direct%aer_opt_param,  &
                           ts%data%direct%cld_opt_param,  &
                           ts%data%calcemis,              &
                           ts%data%calcrefl,              &
                           ts%data%direct%emissivity,     &
                           ts%data%direct%reflectance,    &
                           ts%defn%radar,                 &
                           ts%defn%taylor_by_chan,        &
                           ts%defn%taylor_on_btrefl)

    IF (ts%data%rttov_errorstatus /= 0) THEN
      err = errorstatus_fatal
      THROW(err.NE.0)
    ENDIF
  ENDIF

  ts%defn%do_print = .FALSE.

  CATCH
  END SUBROUTINE run

  SUBROUTINE alloc_traj(err, ts, asw)

    INTEGER(jpim),           INTENT(OUT)           :: err
    TYPE(rttov_test_struct), INTENT(INOUT), TARGET :: ts
    INTEGER(jpim),           INTENT(IN)            :: asw

    TRY

!     IF (ts%defn%do_direct) THEN
      CALL rttov_alloc_traj(                               &
            err,                  ts%defn%nprofiles,       &
            ts%defn%nchannels,    ts%defn%opts,            &
            ts%defn%nlevels,      ts%data%coefs,           &
            asw,                  traj = ts%data%traj)
      THROW(err.NE.0)
!     ENDIF

    IF (ts%defn%do_tl) THEN
      CALL rttov_alloc_traj(                               &
            err,                  ts%defn%nprofiles,       &
            ts%defn%nchannels,    ts%defn%opts,            &
            ts%defn%nlevels,      ts%data%coefs,           &
            asw,                  traj_tl = ts%data%traj_tl)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_ad) THEN
      CALL rttov_alloc_traj(                               &
            err,                  ts%defn%nprofiles,       &
            ts%defn%nchannels,    ts%defn%opts,            &
            ts%defn%nlevels,      ts%data%coefs,           &
            asw,                  traj_ad = ts%data%traj_ad)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_k) THEN
      CALL rttov_alloc_traj(                               &
            err,                  ts%defn%nprofiles,       &
            ts%defn%nchannels,    ts%defn%opts,            &
            ts%defn%nlevels,      ts%data%coefs,           &
            asw,                  traj_k = ts%data%traj_k)
      THROW(err.NE.0)
    ENDIF

    CATCH
  END SUBROUTINE alloc_traj

  SUBROUTINE setup(ts, path, err)
 !
    TYPE(rttov_test_struct), INTENT(INOUT), TARGET :: ts
    CHARACTER(LEN=*),        INTENT(IN)            :: path
    INTEGER(jpim),           INTENT(OUT)           :: err
 !
    INTEGER(jpim), POINTER :: kchannels(:), &  ! absolute channel number in coef file, full list of channels
                              ichannels(:), &  ! restricted list; each channel number appears once
                              lchannels(:)     ! absolute channel number -> relative channel number
    INTEGER(jpim) :: i, j, lchannelsmax, ichannelsmax
    INTEGER(jpim) :: nprofiles1, nchannels1, nlevels
    INTEGER(jpim) :: imult, k1, k2
    INTEGER(jpim) :: file_id
    INTEGER(jpim) :: file_id_nml
    LOGICAL(jplm) :: file_exists
    TYPE(rttov_test_defn) :: defn
    CHARACTER(LEN=256) :: interp_filename, interp_msg
    LOGICAL(KIND=jplm), ALLOCATABLE :: use_chan(:,:)
    NAMELIST / coef_nml / defn
    NAMELIST / rttov_test_nml / defn

    TRY

    CALL rttov_get_lun(file_id_nml)

    OPEN(file_id_nml, file = TRIM(path)//'/rttov_test.txt', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(path)//'/rttov_test.txt')
    READ(file_id_nml, nml = rttov_test_nml, iostat = err)
    THROWM(err.NE.0,"Cannot read "//TRIM(path)//'/rttov_test.txt')
    CLOSE(file_id_nml, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(path)//'/rttov_test.txt')

    OPEN(file_id_nml, file = TRIM(path)//'/../in/coef.txt', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(path)//'/../in/coef.txt')
    READ(file_id_nml, nml = coef_nml, iostat = err)
    THROWM(err.NE.0,"Cannot read "//TRIM(path)//'/../in/coef.txt')
    CLOSE(file_id_nml, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(path)//'/../in/coef.txt')

    CALL rttov_put_lun (file_id_nml)

    CALL setGformat(defn%realprec)  ! Set format specifier for output

    defn%do_rttovscatt = defn%f_rttovscatt .NE. ""

    defn%opts%rt_all%ozone_data = defn%f_o3     .NE. "" .OR. defn%o3_col_int  > 0._jprb .OR. defn%o3_max_ppmv > 0._jprb .OR. &
                                                             defn%o3_col_int_du > 0._jprb
    defn%opts_scatt%ozone_data = defn%opts%rt_all%ozone_data
    IF (.NOT. defn%do_rttovscatt) THEN
      defn%opts%rt_all%co2_data = defn%f_co2    .NE. "" .OR. defn%co2_col_int > 0._jprb .OR. defn%co2_max_ppmv > 0._jprb
      defn%opts%rt_all%n2o_data = defn%f_n2o    .NE. "" .OR. defn%n2o_col_int > 0._jprb .OR. defn%n2o_max_ppmv > 0._jprb
      defn%opts%rt_all%co_data  = defn%f_co     .NE. "" .OR. defn%co_col_int  > 0._jprb .OR. defn%co_max_ppmv  > 0._jprb
      defn%opts%rt_all%ch4_data = defn%f_ch4    .NE. "" .OR. defn%ch4_col_int > 0._jprb .OR. defn%ch4_max_ppmv > 0._jprb
      defn%opts%rt_all%so2_data = defn%f_so2    .NE. "" .OR. defn%so2_col_int > 0._jprb .OR. defn%so2_max_ppmv > 0._jprb
      defn%opts%rt_mw%clw_data  = defn%f_clw    .NE. ""
      defn%opts%rt_ir%user_aer_opt_param = defn%f_aer_opt_param .NE. ""
      defn%opts%rt_ir%addaerosl = defn%f_aerosl .NE. "" .OR. defn%opts%rt_ir%user_aer_opt_param
      defn%opts%rt_ir%user_cld_opt_param = defn%f_cld_opt_param .NE. ""
      defn%opts%rt_ir%addclouds = (defn%f_cloud .NE. "" .OR. defn%opts%rt_ir%user_cld_opt_param) &
                                    .AND. defn%f_cfrac .NE. ""
      defn%opts%rt_ir%pc%addpc     = defn%f_pcscores     .NE. "" .OR. defn%opts%htfrtc_opts%htfrtc
      defn%opts%rt_ir%pc%addradrec = defn%f_channels_rec .NE. "" .AND. defn%opts%rt_ir%pc%addpc
      defn%opts%htfrtc_opts%reconstruct = defn%opts%rt_ir%pc%addradrec

      ! Switchrad doesn't apply for PC-RTTOV unless addradrec is true: this avoids check_opts reporting an error
      IF (defn%opts%rt_ir%pc%addpc .AND. .NOT. defn%opts%rt_ir%pc%addradrec) defn%opts%rt_all%switchrad = .FALSE.
    ENDIF

    nprofiles1 = defn%nprofiles
    nchannels1 = defn%nchannels

    defn%nprofiles = defn%mult * defn%nprofiles
    defn%nchannels = defn%mult * defn%nchannels

    ts%defn = defn
    ts%path = path

    IF (ts%defn%prof_by_prof) THEN
      ts%defn%parallel_strategy = 1
      ts%defn%nthreads = 1
    ELSEIF (ts%defn%chan_by_chan) THEN
      ts%defn%parallel_strategy = 2
      ts%defn%nthreads = 1
    ELSE
      ts%defn%parallel_strategy = 0
    ENDIF

    IF (ts%defn%opts%htfrtc_opts%htfrtc) THEN
      IF (ts%defn%do_tl .OR. ts%defn%do_k_tl .OR. &
          ts%defn%do_ad .OR. ts%defn%do_k_ad .OR. ts%defn%do_taylor) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0,"HTFRTC: only direct and K models are available")
      ENDIF
    ENDIF

! read input data

    ! This is passed as an argument to various subroutines: NAG v6.1 complains if this is not allocated
    ALLOCATE(ts%data%channels_rec(MAX(ts%defn%nchannels_rec,0)), stat = err)
    THROW(err.NE.0)
    IF (ts%defn%opts%rt_ir%pc%addpc) THEN
      IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
!         ALLOCATE(ts%data%channels_rec(ts%defn%nchannels_rec), stat = err)
!         THROW(err.NE.0)
        CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_channels_rec),  n1 = ts%data%channels_rec, foe = .TRUE._jplm)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (ts%defn%opts%htfrtc_opts%htfrtc) THEN
!
! Read HTFRTC coefficients now because we need n_f (the number of centroids)
!
      IF (ts%defn%opts%htfrtc_opts%reconstruct) THEN
        CALL rttov_read_coefs_htfrtc(err, ts%data%coefs, TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_htfrtc_coef), &
                                     TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_htfrtc_sensor), ts%data%channels_rec)
      ELSE
        CALL rttov_read_coefs_htfrtc(err, ts%data%coefs, TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_htfrtc_coef), &
                                     TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_htfrtc_sensor))
      ENDIF
      THROWM(err.NE.0,"Error reading HTFRTC coefs")

      ! nchannels is number of centroids
      nchannels1 = ts%data%coefs%coef_htfrtc%n_f * nprofiles1
      ts%defn%nchannels = ts%defn%mult * nchannels1
    ENDIF

    ALLOCATE(                                       &
      ts%data%chanprof(ts%defn%nchannels),          &
      ts%data%calcemis(ts%defn%nchannels),          &
      ts%data%calcrefl(ts%defn%nchannels),          &
      kchannels(nchannels1),                        &
      stat = err)
    THROWM(err.NE.0,"Cannot allocate ts%data")

    ! This is passed as an argument to various subroutines: NAG v6.1 complains if this is not allocated
!     IF (ts%defn%do_rttovscatt) THEN
      ALLOCATE(ts%data%frequencies(defn%nchannels), stat = err)
      THROWM(err.NE.0,"Cannot allocate ts%data")
!     ENDIF

    IF (.NOT. ts%defn%opts%htfrtc_opts%htfrtc) THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_channels),  n1 = kchannels, foe = .TRUE._jplm)
      THROW(err.NE.0)
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_lprofiles), n1 = ts%data%chanprof(1:nchannels1)%prof, foe = .TRUE._jplm)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%opts%rt_ir%pc%addpc) THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_pcscores), npcscores = ts%defn%npcscores, &
                 ipcreg = ts%defn%opts%rt_ir%pc%ipcreg, ipcbnd = ts%defn%opts%rt_ir%pc%ipcbnd, foe = .TRUE._jplm)
      THROW(err.NE.0)
      ts%defn%opts%rt_ir%pc%npcscores = ts%defn%npcscores
    ENDIF

    IF (ts%defn%opts%htfrtc_opts%htfrtc) THEN
      ts%defn%opts%htfrtc_opts%n_pc_in = ts%defn%npcscores
      ts%defn%user_check_opts = .FALSE.

      DEALLOCATE(kchannels, stat = err)
      THROWM(err.NE.0,"Deallocation of kchannels failed")
    ELSE
#ifdef _RTTOV_HDF
      CALL OPEN_HDF(.FALSE., err) ! 32 bits
      THROW(err.NE.0)
#endif

!
! Channel indices calculations
!
      ! kchannels(:) holds the full list of channels from channels.txt
      ! lchannels(:) is of size MAXVAL(kchannels): for every channel appearing
      !   in channels.txt, the corresponding entry in lchannels is set non-zero.
      !   This provides a lookup mapping actual channel numbers to the extracted
      !   channel indices.
      ! ichannelsmax is the number of channels to extract
      ! ichannels(:) contains the list of channel numbers to extract
      lchannelsmax = MAXVAL(kchannels)
      ALLOCATE(lchannels(lchannelsmax), stat = err)
      THROWM(err.NE.0,"Cannot allocate lchannels")

      ! First set elements of lchannels to 1 for each channel present
      lchannels = 0
      DO i = 1, nchannels1
        lchannels(kchannels(i)) = 1
      ENDDO

      ! This gives the number of channels to extract
      ichannelsmax = SUM(lchannels)

      ! Now renumber the requested channels from 1 to ichannelsmax
      ! By doing it this way we allow for channels.txt to contain channels
      ! in non-consecutive order (e.g. chan 3 for prof 1, chan 2 for prof 2)
      j = 0
      DO i = 1, lchannelsmax
        IF (lchannels(i) > 0) THEN
          j = j + 1
          lchannels(i) = j
        ENDIF
      ENDDO

      ! Now populate the list of channels to extract by looping over
      ! lchannels and setting ichannel for every non-zero value
      ALLOCATE(ichannels(ichannelsmax), stat = err)
      THROWM(err.NE.0,"Cannot allocate ichannels")
      j = 0
      DO i = 1, lchannelsmax
        IF (lchannels(i) > 0) THEN
          j = j + 1
          ichannels(j) = i
        ENDIF
      ENDDO
!
! read coefficients
!
      IF (ts%defn%lallchans .OR. ts%defn%do_rttovscatt) THEN
        CALL rttov_read_coefs( &
            err, ts%data%coefs, ts%defn%opts, &
            file_coef       = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_coef),       &
            file_scaer      = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_scaer),      &
            file_sccld      = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_sccld),      &
!             file_mfasis_aer = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_mfasis_aer), &
            file_mfasis_cld = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_mfasis_cld), &
            file_pccoef     = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_pccomp))
        THROW(err.NE.0)

        IF (ts%defn%do_rttovscatt) THEN
          CALL rttov_read_scattcoeffs(err, ts%defn%opts_scatt, ts%data%coefs, ts%data%coefs_scatt, &
              file_coef = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_rttovscatt))
          THROW(err.NE.0)
          ts%defn%nhydro_frac = 1
          IF (ts%defn%multi_hydro_frac) ts%defn%nhydro_frac = ts%data%coefs_scatt%nhydro
        ENDIF
      ELSE
        CALL rttov_read_coefs( &
            err, ts%data%coefs, ts%defn%opts, &
            channels        = ichannels(1:ichannelsmax), &
            channels_rec    = ts%data%channels_rec,                                       &
            file_coef       = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_coef),       &
            file_scaer      = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_scaer),      &
            file_sccld      = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_sccld),      &
!             file_mfasis_aer = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_mfasis_aer), &
            file_mfasis_cld = TRIM(ts%defn%coef_prefix)//'/'//TRIM(ts%defn%f_mfasis_cld), &
            file_pccoef     = TRIM(ts%defn%coef_prefix) //'/'//TRIM(ts%defn%f_pccomp))
        THROW(err.NE.0)
      ENDIF

      IF (ts%defn%do_rttovscatt) THEN
        ! use_chan array is dimensioned by the total number of instrument channels
        ALLOCATE(use_chan(nprofiles1,ts%data%coefs%coef%fmv_chn))

        ! Set use_chan to .TRUE. only for required channels
        use_chan(:,:) = .FALSE._jplm
        DO j = 1, nchannels1
          use_chan(ts%data%chanprof(j)%prof,kchannels(j)) = .TRUE._jplm
        ENDDO

        ! Populate chanprof and frequencies arrays
        CALL rttov_scatt_setupindex(err, nprofiles1, ts%data%coefs%coef%fmv_chn, ts%data%coefs, &
              ts%data%coefs_scatt, nchannels1, ts%data%chanprof, ts%data%frequencies, use_chan)
        THROW(err.NE.0)

        DEALLOCATE(use_chan)

      ELSEIF (ts%defn%lallchans) THEN
        ts%data%chanprof(1:nchannels1)%chan = kchannels
      ELSE
!
! Renumber all channels
!
        ts%data%chanprof(1:nchannels1)%chan = lchannels(kchannels)
        IF (ts%defn%opts%rt_ir%pc%addpc) THEN
          IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
            DO i = 1, ts%defn%nchannels_rec
              ts%data%channels_rec(i) = i
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      DEALLOCATE(kchannels, lchannels, ichannels, stat = err)
      THROWM(err.NE.0,"Deallocation of kchannels, lchannels, ... failed")
    ENDIF ! not HTFRTC

    ! Read the aer/cld optical parameter dimensions
    IF (defn%opts%rt_ir%user_aer_opt_param) THEN
      CALL read_opt_param(err, ts, TRIM(ts%path)//'/'//TRIM(defn%f_aer_opt_param), ts%data%direct%aer_opt_param, &
                          ts%defn%aer_nmom, ts%defn%aer_nphangle, .TRUE._jplm)
      THROWM(err.NE.0,"Error initialising aerosol optical parameters")
    ENDIF
    IF (defn%opts%rt_ir%user_cld_opt_param) THEN
      CALL read_opt_param(err, ts, TRIM(ts%path)//'/'//TRIM(defn%f_cld_opt_param), ts%data%direct%cld_opt_param, &
                          ts%defn%cld_nmom, ts%defn%cld_nphangle, .TRUE._jplm)
      THROWM(err.NE.0,"Error initialising cloud optical parameters")
    ENDIF

    CALL alloc_data(err, ts, ts%defn%nprofiles, ts%data%direct)

    IF (ts%defn%do_tl) THEN
      CALL alloc_data(err, ts, ts%defn%nprofiles, ts%data%tl)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_ad) THEN
      CALL alloc_data(err, ts, ts%defn%nprofiles, ts%data%ad, do_ad=.TRUE._jplm)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k) THEN
      CALL alloc_data(err, ts, ts%defn%nchannels, ts%data%k, ts%data%data_k)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_bf) THEN
      CALL alloc_data(err, ts, ts%defn%nchannels, ts%data%k_bf, ts%data%data_k_bf)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_tl) THEN
      CALL alloc_data(err, ts, ts%defn%nchannels, ts%data%k_tl, ts%data%data_k_tl)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_ad) THEN
      CALL alloc_data(err, ts, ts%defn%nchannels, ts%data%k_ad, ts%data%data_k_ad)
      THROW(err.NE.0)
    ENDIF

!
! emissivity
!
    ts%data%calcemis = .TRUE.
    ts%data%direct%emissivity%emis_in = 1.0_jprb
    IF (defn%f_calcemis .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_calcemis), l1 = ts%data%calcemis(1:nchannels1))
      THROW(err.NE.0)
    ENDIF
    IF (defn%f_emissivity .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_emissivity), &
                 r1 = ts%data%direct%emissivity(1:nchannels1)%emis_in)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%fixedemis >= 0._jprb) THEN
      ts%data%calcemis = .FALSE.
      ts%data%direct%emissivity%emis_in = ts%defn%fixedemis
    ENDIF
    ts%data%direct%emissivity%specularity = 0._jprb
    IF (defn%f_specularity .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_specularity), &
                 r1 = ts%data%direct%emissivity(1:nchannels1)%specularity)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%fixedspecularity >= 0._jprb) THEN
      ts%data%direct%emissivity%specularity = ts%defn%fixedspecularity
    ENDIF

!
! reflectance
!
    ts%data%calcrefl = .TRUE.
    ts%data%direct%reflectance%refl_in = 0.3_jprb * pi_r
    ts%data%direct%reflectance%diffuse_refl_in = 0._jprb
    ts%data%direct%reflectance%refl_cloud_top = 0._jprb   ! If zero, RTTOV uses default values
    IF (defn%f_calcrefl .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_calcrefl), l1 = ts%data%calcrefl(1:nchannels1))
      THROW(err.NE.0)
    ENDIF
    IF (defn%f_reflectance .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_reflectance), &
                 r1 = ts%data%direct%reflectance(1:nchannels1)%refl_in)
      THROW(err.NE.0)
    ENDIF
    IF (defn%f_reflectance_diffuse .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_reflectance_diffuse), &
                 r1 = ts%data%direct%reflectance(1:nchannels1)%diffuse_refl_in)
      THROW(err.NE.0)
    ENDIF
    IF (defn%f_reflectance_cloud_top .NE. "") THEN
      CALL reada(err, TRIM(ts%path)//'/'//TRIM(defn%f_reflectance_cloud_top), &
                  r1 = ts%data%direct%reflectance(1:nchannels1)%refl_cloud_top)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%fixedrefl >= 0._jprb) THEN
      ts%data%calcrefl = .FALSE.
      ts%data%direct%reflectance%refl_in = ts%defn%fixedrefl
    ENDIF
    IF (ts%defn%fixeddiffuse_refl >= 0._jprb) THEN
      ts%data%direct%reflectance%diffuse_refl_in = ts%defn%fixeddiffuse_refl
    ENDIF

!
! replicate emissivity and reflectance data
!
    DO imult = 2, ts%defn%mult
      k1 = 1 + (imult - 1) * nchannels1
      k2 = k1 + nchannels1 - 1
      ts%data%direct%emissivity(k1:k2)%emis_in          = ts%data%direct%emissivity(1:nchannels1)%emis_in
      ts%data%direct%emissivity(k1:k2)%specularity      = ts%data%direct%emissivity(1:nchannels1)%specularity
      ts%data%direct%reflectance(k1:k2)%refl_in         = ts%data%direct%reflectance(1:nchannels1)%refl_in
      ts%data%direct%reflectance(k1:k2)%diffuse_refl_in = ts%data%direct%reflectance(1:nchannels1)%diffuse_refl_in
      ts%data%direct%reflectance(k1:k2)%refl_cloud_top  = ts%data%direct%reflectance(1:nchannels1)%refl_cloud_top
      ts%data%chanprof(k1:k2)%chan                      = ts%data%chanprof(1:nchannels1)%chan
      ts%data%chanprof(k1:k2)%prof                      = ts%data%chanprof(1:nchannels1)%prof + (imult - 1) * nprofiles1
      ts%data%calcemis(k1:k2)                           = ts%data%calcemis(1:nchannels1)
      ts%data%calcrefl(k1:k2)                           = ts%data%calcrefl(1:nchannels1)
      IF (ts%defn%do_rttovscatt) THEN
        ts%data%frequencies(k1:k2)                      = ts%data%frequencies(1:nchannels1)
      ENDIF
    ENDDO

    CALL read_profile_data(ts, err)
    THROW(err.NE.0)

    IF (ts%defn%fixed2m) THEN
      nlevels = ts%data%direct%profiles(1)%nlevels
      DO i = 1, defn%nprofiles
        ts%data%direct%profiles(i)%s2m%p = ts%data%direct%profiles(i)%p(nlevels)
        ts%data%direct%profiles(i)%s2m%t = ts%data%direct%profiles(i)%t(nlevels)
        ts%data%direct%profiles(i)%s2m%q = ts%data%direct%profiles(i)%q(nlevels)
        IF (defn%opts%rt_all%ozone_data) ts%data%direct%profiles(i)%s2m%o = ts%data%direct%profiles(i)%o3(nlevels)
      ENDDO
    ENDIF
    IF (ts%defn%fixedzenang >= 0._jprb)     ts%data%direct%profiles(:)%zenangle = ts%defn%fixedzenang
    IF (ts%defn%fixedsunzenang >= 0._jprb)  ts%data%direct%profiles(:)%sunzenangle = ts%defn%fixedsunzenang
    IF (ts%defn%fixedsurftype >= 0)         ts%data%direct%profiles(:)%skin%surftype = ts%defn%fixedsurftype
    IF (ts%defn%fixedvisir_clw_scheme >= 0) ts%data%direct%profiles(:)%clw_scheme = ts%defn%fixedvisir_clw_scheme
    IF (ts%defn%fixedvisir_ice_scheme >= 0) ts%data%direct%profiles(:)%ice_scheme = ts%defn%fixedvisir_ice_scheme
    IF (ts%defn%fixedcfraction >= 0._jprb)  ts%data%direct%profiles(:)%cfraction = ts%defn%fixedcfraction

    IF (ts%defn%do_rttovscatt) THEN
      CALL read_cld_profile_data(ts, err)
      THROW(err.NE.0)
    ENDIF
!
! Check pressure levels to determine whether interpolation is required
!
    ! User can set addinterp on commandline (false by default) and this is read
    ! into ts%defn%opts%interpolation%addinterp, but is overridden if profile
    ! pressure levels differ sufficiently from coef levels (tolerance is more
    ! strict than in rttov_check_reg_limits to ensure we turn interpolator on for
    ! certain tests which require it).

    IF (.NOT. ts%defn%opts%htfrtc_opts%htfrtc) THEN
      CALL rttov_get_lun (file_id)
      interp_filename = TRIM(ts%path)//'/interpolation.log'
      OPEN(file_id, file = interp_filename, form = 'formatted', iostat = err)
      THROWM(err.NE.0,"Cannot open "//TRIM(interp_filename))

      interp_msg = 'Interpolation is OFF'
      IF (ts%defn%opts%interpolation%addinterp) THEN
        interp_msg = 'User has set interpolation ON'
      ELSE
        DO i = 1, nprofiles1
          IF (ts%data%direct%profiles(i)%nlevels /= SIZE(ts%data%coefs%coef%ref_prfl_p(:))) THEN
            IF (.NOT. ts%defn%opts%interpolation%addinterp) &
              interp_msg = 'Number of input levels differs to coef file: switching interpolation ON'
            ts%defn%opts%interpolation%addinterp = .TRUE.
            EXIT
          ELSEIF (ANY(ABS(ts%data%direct%profiles(i)%p(:) - ts%data%coefs%coef%ref_prfl_p(:)) / &
                                     ts%data%coefs%coef%ref_prfl_p(:) > 0.005_jprb)) THEN
            IF (.NOT. ts%defn%opts%interpolation%addinterp) &
              interp_msg = 'Input pressures differ to coef levels: switching interpolation ON'
            ts%defn%opts%interpolation%addinterp = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDIF

      WRITE(file_id, '(a)') TRIM(interp_msg)

      ! lgradp is only valid/useful if interpolation is switched on
      IF (.NOT. ts%defn%opts%interpolation%addinterp .AND. ts%defn%opts%interpolation%lgradp) THEN
        WRITE(file_id, '(a)') 'Interpolation is off: switching lgradp off'
        ts%defn%opts%interpolation%lgradp = .FALSE.
      ENDIF

      CLOSE(file_id, iostat = err)
      THROWM(err.NE.0,"Cannot close "//TRIM(interp_filename))
      CALL rttov_put_lun (file_id)
    ENDIF

    IF (ts%defn%ltemp) THEN
      CALL alloc_traj(err, ts, 1_jpim)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_tl) THEN
      CALL rttov_make_profile_inc(ts%data%tl%profiles, ts%data%direct%profiles, ts%defn%opts)

      IF (ts%defn%opts%rt_ir%user_aer_opt_param) &
        CALL rttov_make_opt_param_inc(ts%data%tl%aer_opt_param, ts%data%direct%aer_opt_param)

      IF (ts%defn%opts%rt_ir%user_cld_opt_param) &
        CALL rttov_make_opt_param_inc(ts%data%tl%cld_opt_param, ts%data%direct%cld_opt_param)

      IF (ts%defn%do_rttovscatt) CALL rttov_make_cld_profile_inc(ts%data%tl%cld_profiles, ts%data%direct%cld_profiles)

      CALL rttov_make_emisrefl_inc(ts%defn%opts, ts%data%coefs, ts%data%direct%profiles, &
            ts%data%chanprof, ts%data%calcemis, &
            ts%data%tl%emissivity%emis_in, ts%data%tl%reflectance%refl_in, &
            ts%data%direct%reflectance%diffuse_refl_in, ts%data%tl%reflectance%diffuse_refl_in, &
            ts%data%direct%emissivity%specularity, ts%data%tl%emissivity%specularity)

      IF (ts%defn%scale_inc /= 1.0_jprb) THEN
        CALL rttov_scale_profile_inc(ts%data%tl%profiles, ts%defn%scale_inc)
        IF (ts%defn%opts%rt_ir%user_aer_opt_param) &
          CALL rttov_scale_opt_param_inc(ts%data%tl%aer_opt_param, ts%defn%scale_inc)
        IF (ts%defn%opts%rt_ir%user_cld_opt_param) &
          CALL rttov_scale_opt_param_inc(ts%data%tl%cld_opt_param, ts%defn%scale_inc)
        IF (ts%defn%do_rttovscatt) CALL rttov_scale_cld_profile_inc(ts%data%tl%cld_profiles, ts%defn%scale_inc)
        ts%data%tl%emissivity(:)%emis_in = ts%data%tl%emissivity(:)%emis_in * ts%defn%scale_inc
        ts%data%tl%emissivity(:)%specularity = ts%data%tl%emissivity(:)%specularity * ts%defn%scale_inc
        ts%data%tl%reflectance(:)%refl_in = ts%data%tl%reflectance(:)%refl_in * ts%defn%scale_inc
        ts%data%tl%reflectance(:)%diffuse_refl_in = ts%data%tl%reflectance(:)%diffuse_refl_in * ts%defn%scale_inc
      ENDIF
    ENDIF

    IF (ts%defn%do_ad) THEN
      IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        CALL rttov_make_pccomp_inc(ts%data%ad%pccomp_saved, ts%defn%opts)
        IF (ts%defn%scale_inc /= 1.0_jprb) &
          CALL rttov_scale_pccomp_inc(ts%data%ad%pccomp_saved, ts%defn%scale_inc, ts%defn%opts)
      ELSEIF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
        CALL rttov_make_reflectivity_inc(ts%data%ad%reflectivity_saved)
        IF (ts%defn%scale_inc /= 1.0_jprb) &
          CALL rttov_scale_reflectivity_inc(ts%data%ad%reflectivity_saved, ts%defn%scale_inc)
      ELSE
        CALL rttov_make_radiance_inc(ts%data%coefs%coef, ts%data%ad%radiance_saved, &
                                     ts%data%chanprof(:)%chan, ts%defn%nchannels)
        IF (ts%defn%scale_inc /= 1.0_jprb) &
          CALL rttov_scale_radiance_inc(ts%data%ad%radiance_saved, ts%defn%scale_inc)
      ENDIF
    ENDIF

    IF (ts%data%coefs%coef%pmc_shift) THEN
      INQUIRE(file=TRIM(ts%path)//'/../in/pmc.txt', exist=file_exists)
      IF (file_exists) THEN
        CALL rttov_get_lun (file_id)
        OPEN(file_id, file = TRIM(path)//'/../in/pmc.txt', iostat = err)
        THROWM(err.NE.0,"Cannot open "//TRIM(path)//'/../in/pmc.txt')
        READ(file_id, *, iostat = err) ts%data%coefs%coef%pmc_ppmc
        THROWM(err.NE.0,"Cannot read "//TRIM(path)//'/../in/pmc.txt')
        CLOSE(file_id, iostat = err)
        THROWM(err.NE.0,"Cannot close "//TRIM(path)//'/../in/pmc.txt')
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,'PMC coef file requires input pmc.txt file')
      ENDIF
    ENDIF

    CATCH
  END SUBROUTINE setup

  SUBROUTINE read_profile_data(ts, err)

    TYPE(rttov_test_struct), INTENT(INOUT), TARGET :: ts
    INTEGER(jpim),           INTENT(OUT)           :: err

    REAL(jprb)  :: becosbk(2)
    INTEGER(jpim) :: datetime(6)
    TYPE(rttov_test_defn), POINTER :: defn
    TYPE(rttov_input_data), POINTER :: direct
    INTEGER(jpim) :: iprof
    INTEGER(jpim) :: kprof
    INTEGER(jpim) :: imult
    CHARACTER(LEN=6)   :: sproffmt
    CHARACTER(LEN=8)   :: sprof
    CHARACTER(LEN=256) :: in_dir
    INTEGER(jpim) :: nprofiles1
    LOGICAL(jplm) :: foe
    LOGICAL(jplm) :: file_exists
    INTEGER(jpim) :: file_id
    LOGICAL(jplm) :: gas_conv_flag
    CHARACTER(LEN=256) :: gas_units_filename

    TYPE(rttov_options) :: opts

#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

    TRY

    defn => ts%defn
    direct => ts%data%direct

    nprofiles1 = ts%defn%nprofiles / ts%defn%mult

    DO iprof = 1, nprofiles1

      foe = (iprof == 1)

      WRITE(sproffmt, "('(i',i1,'.',i1,')')") ts%defn%prof_ndigits, ts%defn%prof_ndigits
      WRITE(sprof, sproffmt) iprof
      in_dir =  TRIM(ts%path)//'/../in/profiles/'//TRIM(sprof)

      CALL reada(err, TRIM(in_dir)//'/ground/skin.txt', k0 = direct%profiles(iprof)%skin, &
                 foe = foe, k0_a = direct%profiles(1)%skin)
      THROW(err.NE.0)
      CALL reada(err, TRIM(in_dir)//'/ground/s2m.txt',  s0 = direct%profiles(iprof)%s2m,  &
                 foe = foe, s0_a = direct%profiles(1)%s2m)
      THROW(err.NE.0)
      CALL reada(err, TRIM(in_dir)//'/angles.txt', &
        zenangle    = direct%profiles(iprof)%zenangle,    azangle    = direct%profiles(iprof)%azangle,    &
        sunzenangle = direct%profiles(iprof)%sunzenangle, sunazangle = direct%profiles(iprof)%sunazangle, &
        latitude    = direct%profiles(iprof)%latitude,    longitude  = direct%profiles(iprof)%longitude,  &
        elevation   = direct%profiles(iprof)%elevation, &
        foe = foe, &
        zenangle_a    = direct%profiles(1)%zenangle,    azangle_a    = direct%profiles(1)%azangle,    &
        sunzenangle_a = direct%profiles(1)%sunzenangle, sunazangle_a = direct%profiles(1)%sunazangle, &
        latitude_a    = direct%profiles(1)%latitude,    longitude_a  = direct%profiles(1)%longitude,  &
        elevation_a   = direct%profiles(1)%elevation)
      THROW(err.NE.0)
      INQUIRE(file=TRIM(in_dir)//'/atm/simple_cloud.txt', exist=file_exists)
      IF (file_exists) THEN
        CALL reada(err, TRIM(in_dir)//'/atm/simple_cloud.txt', &
          ctp = direct%profiles(iprof)%ctp, cfraction = direct%profiles(iprof)%cfraction, &
          ctp_a = direct%profiles(1)%ctp, cfraction_a = direct%profiles(1)%cfraction, &
          foe = foe)
        THROW(err.NE.0)
      ENDIF
      CALL reada(err, TRIM(in_dir)//'/atm/p.txt',  r1 = direct%profiles(iprof)%p, foe = foe)
      THROW(err.NE.0)
      CALL reada(err, TRIM(in_dir)//'/atm/t.txt',  r1 = direct%profiles(iprof)%t, foe = foe)
      THROW(err.NE.0)
      CALL reada(err, TRIM(in_dir)//'/atm/q.txt',  r1 = direct%profiles(iprof)%q, foe = foe)
      THROW(err.NE.0)

      direct%profiles(iprof)%nlevels = ts%defn%nlevels

      IF (defn%f_o3 .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/o3.txt',     &
          r1 = direct%profiles(iprof)%o3,  foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%f_co2 .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/co2.txt',    &
          r1 = direct%profiles(iprof)%co2, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%f_n2o .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/n2o.txt',    &
          r1 = direct%profiles(iprof)%n2o, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%f_co .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/co.txt',     &
          r1 = direct%profiles(iprof)%co,  foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%f_ch4 .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/ch4.txt',    &
          r1 = direct%profiles(iprof)%ch4, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%f_so2 .NE. "") THEN
        CALL reada(err, TRIM(in_dir)//'/atm/so2.txt',    &
          r1 = direct%profiles(iprof)%so2, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%opts%rt_mw%clw_data) THEN
        CALL reada(err, TRIM(in_dir)//'/atm/clw.txt',    &
          r1 = direct%profiles(iprof)%clw, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%opts%rt_ir%addaerosl .AND. .NOT. defn%opts%rt_ir%user_aer_opt_param .AND. .NOT. defn%zero_aerosol) THEN
        CALL reada(err, TRIM(in_dir)//'/atm/mmr_cldaer.txt', &
                    mmr_cldaer = direct%profiles(iprof)%mmr_cldaer, &
                    foe = .TRUE._jplm)
        THROW(err.NE.0)
        CALL reada(err, TRIM(in_dir)//'/atm/aerosl.txt', &
          r2 = direct%profiles(iprof)%aerosols, foe = foe)
        THROW(err.NE.0)
      ENDIF
      IF (defn%opts%rt_ir%addclouds .AND. .NOT. defn%zero_cloud) THEN
        CALL reada(err, TRIM(in_dir)//'/atm/cfrac.txt', &
          r1 = direct%profiles(iprof)%cfrac, foe = foe)
        THROW(err.NE.0)
        IF (.NOT. defn%opts%rt_ir%user_cld_opt_param) THEN
          CALL reada(err, TRIM(in_dir)//'/atm/mmr_cldaer.txt', &
                      mmr_cldaer = direct%profiles(iprof)%mmr_cldaer, &
                      foe = .TRUE._jplm)
          THROW(err.NE.0)
          CALL reada(err, TRIM(in_dir)//'/atm/cloud.txt', &
            r2 = direct%profiles(iprof)%cloud, foe = foe)
          THROW(err.NE.0)
          INQUIRE(file=TRIM(in_dir)//'/atm/clwde.txt', exist=file_exists)
          IF (file_exists) THEN
            CALL reada(err, TRIM(in_dir)//'/atm/clwde.txt', &
              r1 = direct%profiles(iprof)%clwde, foe = foe)
            THROW(err.NE.0)
          ENDIF
          CALL reada(err, TRIM(in_dir)//'/atm/clw_scheme.txt', &
                     clw_scheme = direct%profiles(iprof)%clw_scheme, clwde_param = direct%profiles(iprof)%clwde_param, &
                     clw_scheme_a = direct%profiles(1)%clw_scheme, clwde_param_a = direct%profiles(1)%clwde_param, &
                     foe = foe)
          THROW(err.NE.0)
          INQUIRE(file=TRIM(in_dir)//'/atm/icede.txt', exist=file_exists)
          IF (file_exists) THEN
            CALL reada(err, TRIM(in_dir)//'/atm/icede.txt', &
              r1 = direct%profiles(iprof)%icede, foe = foe)
            THROW(err.NE.0)
          ENDIF
          CALL reada(err, TRIM(in_dir)//'/atm/ice_scheme.txt', &
                     icede_param = direct%profiles(iprof)%icede_param, ice_scheme = direct%profiles(iprof)%ice_scheme, &
                     icede_param_a = direct%profiles(1)%icede_param, ice_scheme_a = direct%profiles(1)%ice_scheme, &
                     foe = foe)
          THROW(err.NE.0)
        ENDIF
      ELSEIF (defn%zero_cloud) THEN
        ! Need to specify valid values for clw/ice schemes even if there is no cloud
        direct%profiles(:)%clw_scheme = 1
        direct%profiles(:)%clwde_param = 1
        direct%profiles(:)%ice_scheme = 1
        direct%profiles(:)%icede_param = 1
      ENDIF

      INQUIRE(file=TRIM(in_dir)//'/datetime.txt', exist=file_exists)
      IF (file_exists) THEN
        CALL reada(err, TRIM(in_dir)//'/datetime.txt', &
          n1 = datetime(:), foe = foe)
        THROW(err.NE.0)
        direct%profiles(iprof)%date(1:3) = datetime(1:3)
        direct%profiles(iprof)%time(1:3) = datetime(4:6)
      ENDIF

      becosbk(1) = direct%profiles(1)%be
      becosbk(2) = direct%profiles(1)%cosbk

      INQUIRE(file=TRIM(in_dir)//'/be.txt', exist=file_exists)
      IF (file_exists) THEN
        CALL reada(err, TRIM(in_dir)//'/be.txt', r1 = becosbk, foe = foe)
        THROW(err.NE.0)
        direct%profiles(iprof)%be    = becosbk(1)
        direct%profiles(iprof)%cosbk = becosbk(2)
      ENDIF

      ! Sort out gas units:
      !   defn%input_gas_units: the units of the input gas profiles in q.txt/etc
      !   defn%gas_units: the gas units with which RTTOV is to be run
      !   Default value for both is ppmv wet

      ! valid values for both input_gas_units and run_gas_units:
      !   0 (ppmv dry), 1 (kg/kg wet), 2 (ppmv wet)

      ! If defn%input_gas_units was specified on commandline that value takes precedence;
      !   otherwise if gas_units.txt present then this is defn%input_gas_units;
      !   otherwise defn%input_gas_units = ppmv wet

      IF (defn%input_gas_units == -1) THEN
        INQUIRE(file=TRIM(in_dir)//'/gas_units.txt', exist=file_exists)
        IF (file_exists) THEN
          CALL reada(err, TRIM(in_dir)//'/gas_units.txt', &
                      gas_units = defn%input_gas_units, &
                      foe = .FALSE._jplm, gas_units_a = gas_unit_ppmv)
          THROW(err.NE.0)
        ELSE
          defn%input_gas_units = gas_unit_ppmv
        ENDIF
      ENDIF

      ! Default for gas_units is ppmv wet: all profiles get the same value
      direct%profiles(iprof)%gas_units = defn%run_gas_units

      ! If input_gas_units and run_gas_units differ do conversion
      gas_conv_flag = .TRUE.
      IF (defn%input_gas_units == gas_unit_ppmv) THEN
        IF (direct%profiles(iprof)%gas_units <= gas_unit_ppmvdry) THEN
          ! ppmv wet -> ppmv dry
          CALL ppmvwet2ppmvdry(direct%profiles(iprof))
        ELSEIF (direct%profiles(iprof)%gas_units == gas_unit_specconc) THEN
          ! ppmv wet -> kgkg wet
          CALL ppmvwet2kgkgwet(direct%profiles(iprof))
        ELSE
          gas_conv_flag = .FALSE.
        ENDIF
      ELSEIF (defn%input_gas_units == gas_unit_specconc) THEN
        IF (direct%profiles(iprof)%gas_units <= gas_unit_ppmvdry) THEN
          ! kgkg wet -> ppmv dry
          CALL kgkgwet2ppmvdry(direct%profiles(iprof))
        ELSEIF (direct%profiles(iprof)%gas_units == gas_unit_ppmv) THEN
          ! kgkg wet -> ppmv wet
          CALL kgkgwet2ppmvwet(direct%profiles(iprof))
        ELSE
          gas_conv_flag = .FALSE.
        ENDIF
      ELSE
        IF (direct%profiles(iprof)%gas_units == gas_unit_ppmv) THEN
          ! ppmv dry -> ppmv wet
          CALL ppmvdry2ppmvwet(direct%profiles(iprof))
        ELSEIF (direct%profiles(iprof)%gas_units == gas_unit_specconc) THEN
          ! ppmv dry -> kgkg wet
          CALL ppmvdry2kgkgwet(direct%profiles(iprof))
        ELSE
          gas_conv_flag = .FALSE.
        ENDIF
      ENDIF

      IF (iprof == 1) THEN
        CALL rttov_get_lun (file_id)
        gas_units_filename = TRIM(ts%path)//'/gas_units.log'
        OPEN(file_id, file = gas_units_filename, form = 'formatted', iostat = err)
        THROWM(err.NE.0,"Cannot open "//TRIM(gas_units_filename))

        IF (defn%input_gas_units == gas_unit_ppmv) THEN
          WRITE(file_id, '(a)') TRIM('Input files have gas units  : ppmv over moist air')
        ELSEIF (defn%input_gas_units == gas_unit_specconc) THEN
          WRITE(file_id, '(a)') TRIM('Input files have gas units  : kg/kg over moist air')
        ELSE
          WRITE(file_id, '(a)') TRIM('Input files have gas units  : ppmv over dry air')
        ENDIF
        IF (defn%run_gas_units == gas_unit_ppmv) THEN
          WRITE(file_id, '(a)') TRIM('Running RTTOV with gas units: ppmv over moist air')
        ELSEIF (defn%run_gas_units == gas_unit_specconc) THEN
          WRITE(file_id, '(a)') TRIM('Running RTTOV with gas units: kg/kg over moist air')
        ELSE
          WRITE(file_id, '(a)') TRIM('Running RTTOV with gas units: ppmv over dry air')
        ENDIF
        IF (gas_conv_flag) WRITE(file_id, '(a)') 'Test suite converted units of input profiles'

        CLOSE(file_id, iostat = err)
        THROWM(err.NE.0,"Cannot close "//TRIM(gas_units_filename))
        CALL rttov_put_lun (file_id)
      ENDIF

      ! Generate scaled reference profiles for optional trace gases
      IF (ts%defn%o3_col_int_du > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), o3_col_int_du=ts%defn%o3_col_int_du)
      ELSE IF (ts%defn%o3_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), o3_col_int=ts%defn%o3_col_int)
      ELSE IF (ts%defn%o3_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), o3_max_ppmv=ts%defn%o3_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference ozone profile')

      IF (ts%defn%co2_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), co2_col_int=ts%defn%co2_col_int)
      ELSE IF (ts%defn%co2_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), co2_max_ppmv=ts%defn%co2_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference CO2 profile')

      IF (ts%defn%n2o_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), n2o_col_int=ts%defn%n2o_col_int)
      ELSE IF (ts%defn%n2o_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), n2o_max_ppmv=ts%defn%n2o_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference N2O profile')

      IF (ts%defn%co_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), co_col_int=ts%defn%co_col_int)
      ELSE IF (ts%defn%co_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), co_max_ppmv=ts%defn%co_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference CO profile')

      IF (ts%defn%ch4_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), ch4_col_int=ts%defn%ch4_col_int)
      ELSE IF (ts%defn%ch4_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), ch4_max_ppmv=ts%defn%ch4_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference CH4 profile')

      IF (ts%defn%so2_col_int > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), so2_col_int=ts%defn%so2_col_int)
      ELSE IF (ts%defn%so2_max_ppmv > 0._jprb) THEN
        CALL rttov_scale_ref_gas_prof(err, ts%data%coefs, direct%profiles(iprof:iprof), so2_max_ppmv=ts%defn%so2_max_ppmv)
      ENDIF
      THROWM(err.NE.0,'Error generating scaled reference SO2 profile')

      !
      ! if multiplicity > 1, then replicate profile data
      !
      DO imult = 2, ts%defn%mult
        kprof = iprof + (imult-1) * nprofiles1
        CALL rttov_copy_prof(direct%profiles(kprof:kprof), direct%profiles(iprof:iprof))
      ENDDO

    ENDDO

    ! read the aer/cld optical parameter files
    IF (defn%opts%rt_ir%user_aer_opt_param) THEN
      CALL read_opt_param(err, ts, TRIM(ts%path)//'/../in/aer_opt_param.txt', ts%data%direct%aer_opt_param, &
                          ts%defn%aer_nmom, ts%defn%aer_nphangle, .FALSE._jplm)
      THROWM(err.NE.0,"Error initialising aerosol optical parameters")
      IF (defn%zero_aerosol) THEN
        ! Need to zero the abs/sca
        ts%data%direct%aer_opt_param%abs = 0.
        ts%data%direct%aer_opt_param%sca = 0.
      ENDIF
    ENDIF
    IF (defn%opts%rt_ir%user_cld_opt_param) THEN
      CALL read_opt_param(err, ts, TRIM(ts%path)//'/../in/cld_opt_param.txt', ts%data%direct%cld_opt_param, &
                          ts%defn%cld_nmom, ts%defn%cld_nphangle, .FALSE._jplm)
      THROWM(err.NE.0,"Error initialising cloud optical parameters")
      IF (defn%zero_cloud) THEN
        ! Need to zero the abs/sca
        ts%data%direct%cld_opt_param%abs = 0.
        ts%data%direct%cld_opt_param%sca = 0.
      ENDIF
    ENDIF

#ifdef _RTTOV_HDF
    ! aer/cld params are saved for all channels and multipicity (shall we remove multiplicity?)
    IF (ts%defn%savehdf5) THEN
      IF (defn%opts%rt_ir%user_aer_opt_param) THEN
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/USER_AEROSOL_OPT_PARAM', &
                CREATE=.TRUE., OPT_PARAM = ts%data%direct%aer_opt_param)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/MISC', CREATE=.FALSE., &
                c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/CHANPROF', CREATE=.FALSE., &
                CHANPROF = ts%data%chanprof)
        THROW(err.NE.0)
      ENDIF

      IF (defn%opts%rt_ir%user_cld_opt_param) THEN
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/USER_CLOUD_OPT_PARAM', &
               CREATE=.TRUE., OPT_PARAM = ts%data%direct%cld_opt_param)
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/MISC', CREATE=.FALSE., &
                c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
        THROW(err.NE.0)
        CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/CHANPROF', CREATE=.FALSE., &
                CHANPROF = ts%data%chanprof)
        THROW(err.NE.0)
      ENDIF

      CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/PROFILES.H5",  '/PROFILES', &
                 CREATE=.TRUE., PROFILES = direct%profiles(1:nprofiles1))
      THROW(err.NE.0)

      ! opts is a local variable, so all fields have the default value
      opts%interpolation%addinterp = .TRUE.
      opts%rt_ir%addaerosl   = defn%opts%rt_ir%addaerosl
      opts%rt_ir%addclouds   = defn%opts%rt_ir%addclouds
      opts%rt_all%ozone_data = defn%opts%rt_all%ozone_data
      opts%rt_all%co2_data   = defn%opts%rt_all%co2_data
      opts%rt_all%n2o_data   = defn%opts%rt_all%n2o_data
      opts%rt_all%co_data    = defn%opts%rt_all%co_data
      opts%rt_all%ch4_data   = defn%opts%rt_all%ch4_data
      opts%rt_all%so2_data   = defn%opts%rt_all%so2_data
      opts%rt_mw%clw_data    = defn%opts%rt_mw%clw_data
      opts%rt_ir%user_aer_opt_param = defn%opts%rt_ir%user_aer_opt_param
      opts%rt_ir%user_cld_opt_param = defn%opts%rt_ir%user_cld_opt_param
      CALL RTTOV_HDF_SAVE(err, TRIM(ts%path)//"/PROFILES.H5",  '/OPTIONS', CREATE=.FALSE., OPTIONS = opts)
      THROW(err.NE.0)
    ENDIF
#endif

    CATCH
  END SUBROUTINE read_profile_data

  SUBROUTINE read_cld_profile_data(ts, err)

    TYPE(rttov_test_struct), INTENT(INOUT), TARGET :: ts
    INTEGER(jpim),           INTENT(OUT)           :: err

    TYPE(rttov_input_data), POINTER :: direct
    INTEGER(jpim) :: iprof, kprof, imult, i, nprofiles1, nhydro, file_id
    CHARACTER(LEN=6)   :: sproffmt
    CHARACTER(LEN=8)   :: sprof
    CHARACTER(LEN=256) :: in_dir, f
    LOGICAL(jplm) :: file_exists

    TRY

    direct => ts%data%direct

    nprofiles1 = ts%defn%nprofiles / ts%defn%mult

    DO iprof = 1, nprofiles1

      WRITE(sproffmt, "('(i',i1,'.',i1,')')") ts%defn%prof_ndigits, ts%defn%prof_ndigits
      WRITE(sprof, sproffmt) iprof
      in_dir =  TRIM(ts%path)//'/../in/profiles/'//TRIM(sprof)

      f = TRIM(in_dir)//'/atm/cld_profiles.txt'
      INQUIRE(file=f, exist=file_exists)
      IF (file_exists) THEN
        CALL rttov_get_lun (file_id)
        OPEN(file_id, file = f, form = 'formatted', status = 'old', iostat = err)
        THROWM(err.NE.0,"Cannot open "//TRIM(f))

        READ(file_id, *) nhydro
        IF (nhydro /= ts%data%coefs_scatt%nhydro) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0,"Inconsistent nhydro in hydrotable and cld_profiles.txt")
        ENDIF
        READ(file_id, *) ts%data%direct%cld_profiles(iprof)%flux_conversion
        READ(file_id, *) ts%data%direct%cld_profiles(iprof)%cfrac
        DO i = 1, ts%defn%nlevels
          READ(file_id, *) ts%data%direct%cld_profiles(iprof)%ph(i), &
                           ts%data%direct%cld_profiles(iprof)%hydro(i,:), &
                           ts%data%direct%cld_profiles(iprof)%hydro_frac(i,:)
        ENDDO
        ts%data%direct%cld_profiles(iprof)%ph(ts%defn%nlevels+1) = ts%data%direct%profiles(iprof)%s2m%p
        IF (ts%defn%no_flux_conv) ts%data%direct%cld_profiles(iprof)%flux_conversion = 0

        CLOSE(file_id, iostat = err)
        THROWM(err.NE.0,"Cannot close "//TRIM(f))

        CALL rttov_put_lun (file_id)
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0, "RTTOV-SCATT simulations require cld_profiles.txt")
      ENDIF

      DO imult = 2, ts%defn%mult
        kprof = iprof + (imult-1) * nprofiles1
        CALL rttov_copy_scatt_prof(direct%cld_profiles(kprof:kprof), direct%cld_profiles(iprof:iprof))
      ENDDO

    ENDDO
    CATCH
  END SUBROUTINE read_cld_profile_data

  SUBROUTINE cleanup(ts, err)

    TYPE(rttov_test_struct), INTENT(INOUT), TARGET :: ts
    INTEGER(jpim),           INTENT(OUT)           :: err

    TRY

    DEALLOCATE(ts%data%chanprof, ts%data%calcemis, ts%data%calcrefl, stat = err)
    THROWM(err.NE.0,"Deallocation of ts%data failed")

    IF (ASSOCIATED(ts%data%frequencies)) THEN
      DEALLOCATE(ts%data%frequencies, stat = err)
      THROWM(err.NE.0,"Deallocation of ts%data failed")
    ENDIF

    IF (ASSOCIATED(ts%data%channels_rec)) THEN
      DEALLOCATE(ts%data%channels_rec, stat = err)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%ltemp) THEN
      CALL alloc_traj(err, ts, 0_jpim)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_tl) THEN
      CALL dealloc_data(err, ts, ts%defn%nprofiles, ts%data%tl)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_ad) THEN
      CALL dealloc_data(err, ts, ts%defn%nprofiles, ts%data%ad,   do_ad=.TRUE._jplm)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k) THEN
      CALL dealloc_data(err, ts, ts%defn%nchannels, ts%data%k,    ts%data%data_k)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_bf) THEN
      CALL dealloc_data(err, ts, ts%defn%nchannels, ts%data%k_bf, ts%data%data_k_bf)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_tl) THEN
      CALL dealloc_data(err, ts, ts%defn%nchannels, ts%data%k_tl, ts%data%data_k_tl)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%do_k_ad) THEN
      CALL dealloc_data(err, ts, ts%defn%nchannels, ts%data%k_ad, ts%data%data_k_ad)
      THROW(err.NE.0)
    ENDIF

    CALL dealloc_data(err, ts, ts%defn%nprofiles, ts%data%direct)
    THROW(err.NE.0)

    IF (ts%defn%do_rttovscatt) THEN
      CALL rttov_dealloc_scattcoeffs(ts%data%coefs_scatt)
    ENDIF

    CALL rttov_dealloc_coefs(err, ts%data%coefs)
    THROW(err.NE.0)

    CATCH
  END SUBROUTINE cleanup

  SUBROUTINE reada(err, f,                                      &
    r1, r1_a, r2, r2_a, n1, n1_a, l1, l1_a, k0, k0_a, s0, s0_a, &
    gas_units, gas_units_a, mmr_cldaer,                         &
    zenangle,   azangle,   sunzenangle,   sunazangle,           &
    latitude,   longitude,   elevation,                         &
    zenangle_a, azangle_a, sunzenangle_a, sunazangle_a,         &
    latitude_a, longitude_a, elevation_a,                       &
    clw_scheme,   clwde_param,                                  &
    clw_scheme_a, clwde_param_a,                                &
    icede_param,   ice_scheme,                                  &
    icede_param_a, ice_scheme_a,                                &
    ctp,   cfraction,                                           &
    ctp_a, cfraction_a,                                         &
    ipcreg, ipcbnd, npcscores,                                  &
    foe)

    INTEGER(jpim),      INTENT(OUT)   :: err
    CHARACTER(LEN=*),   INTENT(IN)    :: f
    REAL(jprb),         INTENT(INOUT), OPTIONAL :: r1(:),   r2(:,:)
    REAL(jprb),         INTENT(IN),    OPTIONAL :: r1_a(:), r2_a(:,:)
    INTEGER(kind=jpim), INTENT(INOUT), OPTIONAL :: n1(:)
    INTEGER(kind=jpim), INTENT(IN),    OPTIONAL :: n1_a(:)
    LOGICAL(jplm),      INTENT(INOUT), OPTIONAL :: l1(:)
    LOGICAL(jplm),      INTENT(IN),    OPTIONAL :: l1_a(:)
    TYPE(rttov_skin),   INTENT(INOUT), OPTIONAL :: k0
    TYPE(rttov_skin),   INTENT(IN),    OPTIONAL :: k0_a
    TYPE(rttov_s2m),    INTENT(INOUT), OPTIONAL :: s0
    TYPE(rttov_s2m),    INTENT(IN),    OPTIONAL :: s0_a
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: gas_units
    INTEGER(jpim),      INTENT(IN),    OPTIONAL :: gas_units_a
    LOGICAL(jplm),      INTENT(INOUT), OPTIONAL :: mmr_cldaer
    REAL(jprb),         INTENT(INOUT), OPTIONAL :: zenangle, azangle, sunzenangle, sunazangle
    REAL(jprb),         INTENT(INOUT), OPTIONAL :: latitude, longitude, elevation
    REAL(jprb),         INTENT(IN),    OPTIONAL :: zenangle_a, azangle_a, sunzenangle_a, sunazangle_a
    REAL(jprb),         INTENT(IN),    OPTIONAL :: latitude_a, longitude_a, elevation_a
    REAL(jprb),         INTENT(INOUT), OPTIONAL :: ctp,   cfraction
    REAL(jprb),         INTENT(IN),    OPTIONAL :: ctp_a, cfraction_a
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: clw_scheme,   clwde_param
    INTEGER(jpim),      INTENT(IN),    OPTIONAL :: clw_scheme_a, clwde_param_a
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: icede_param,   ice_scheme
    INTEGER(jpim),      INTENT(IN),    OPTIONAL :: icede_param_a, ice_scheme_a
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: ipcreg
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: ipcbnd
    INTEGER(jpim),      INTENT(INOUT), OPTIONAL :: npcscores
    LOGICAL(jplm),      INTENT(IN),    OPTIONAL :: foe ! fail on data not present

    LOGICAL(jplm) :: foe1
    INTEGER(jpim) :: file_id, ipcreg1, ipcbnd1, npcscores1
    NAMELIST / pcscores / ipcreg, ipcbnd, npcscores
    NAMELIST / units / gas_units
    NAMELIST / cldaer_units / mmr_cldaer
    NAMELIST / angles / zenangle, azangle, sunzenangle, sunazangle, latitude, longitude, elevation
    NAMELIST / s2m / s0
    NAMELIST / clw_scheme_nml / clw_scheme, clwde_param
    NAMELIST / ice_scheme_nml / icede_param, ice_scheme
    NAMELIST / simple_cloud / ctp, cfraction
    NAMELIST / skin / k0

    TRY

    foe1 = .FALSE.
    IF (PRESENT(foe)) foe1 = foe

    CALL rttov_get_lun (file_id)
    IF (foe1) THEN
      OPEN(file_id, file = TRIM(f), form = 'formatted', status = 'old', iostat = err)
      THROWM(err.NE.0,"Cannot open "//TRIM(f))
    ELSE
      OPEN(file_id, file = TRIM(f), form = 'formatted', status = 'old', iostat = err)
      IF (err .NE. 0) THEN
        IF (PRESENT(r1) .AND. PRESENT(r1_a)) THEN
          r1 = r1_a
        ELSEIF (PRESENT(r2) .AND. PRESENT(r2_a)) THEN
          r2 = r2_a
        ELSEIF (PRESENT(l1) .AND. PRESENT(l1_a)) THEN
          l1 = l1_a
        ELSEIF (PRESENT(n1) .AND. PRESENT(n1_a)) THEN
          n1 = n1_a
        ELSEIF (PRESENT(k0) .AND. PRESENT(k0_a)) THEN
          k0 = k0_a
        ELSEIF (PRESENT(gas_units) .AND. PRESENT(gas_units_a)) THEN
          gas_units = gas_units_a
        ELSEIF (PRESENT(s0) .AND. PRESENT(s0_a)) THEN
          s0 = s0_a
        ELSEIF (PRESENT(zenangle) .AND. PRESENT(zenangle_a)) THEN
          zenangle    = zenangle_a
          azangle     = azangle_a
          sunzenangle = sunzenangle_a
          sunazangle  = sunazangle_a
          latitude    = latitude_a
          longitude   = longitude_a
          elevation   = elevation_a
        ELSEIF (PRESENT(clw_scheme) .AND. PRESENT(clw_scheme_a)) THEN
          clw_scheme = clw_scheme_a
          clwde_param = clwde_param_a
        ELSEIF (PRESENT(ice_scheme) .AND. PRESENT(ice_scheme_a)) THEN
          ice_scheme = ice_scheme_a
          icede_param = icede_param_a
        ELSEIF (PRESENT(ctp) .AND. PRESENT(ctp_a)) THEN
          ctp       = ctp_a
          cfraction = cfraction_a
        ELSE
          THROWM(err.NE.0,"Cannot open "//TRIM(f))
        ENDIF
        err = 0
      ENDIF
    ENDIF

    IF (PRESENT(r1)) THEN
      READ(file_id, *, iostat = err) r1
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(r2)) THEN
      READ(file_id, *, iostat = err) r2
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(l1)) THEN
      READ(file_id, *, iostat = err) l1
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(n1)) THEN
      READ(file_id, *, iostat = err) n1
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(k0)) THEN
      READ(file_id, nml = skin, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(gas_units)) THEN
      READ(file_id, nml = units, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(mmr_cldaer)) THEN
      READ(file_id, nml = cldaer_units, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(s0)) THEN
      READ(file_id, nml = s2m, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(zenangle)) THEN
      READ(file_id, nml = angles, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(clw_scheme)) THEN
      READ(file_id, nml = clw_scheme_nml, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(ice_scheme)) THEN
      READ(file_id, nml = ice_scheme_nml, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(ctp)) THEN
      READ(file_id, nml = simple_cloud, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
    ELSEIF (PRESENT(ipcreg) .AND. PRESENT(ipcbnd) .AND. PRESENT(npcscores)) THEN
      ipcreg1 = ipcreg
      ipcbnd1 = ipcbnd
      npcscores1 = npcscores
      READ(file_id, nml = pcscores, iostat = err)
      THROWM(err.NE.0,"Cannot read from "//TRIM(f))
      IF (ipcreg1 > 0) ipcreg = ipcreg1
      IF (ipcbnd1 > 0) ipcbnd = ipcbnd1
      IF (npcscores1 > 0) npcscores = npcscores1
      IF (ipcbnd < 0) ipcbnd = 1_jpim
    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun (file_id)

    CATCH
  END SUBROUTINE reada

  SUBROUTINE read_opt_param(err, ts, f, opt_param, nmom, nphangle, dims_only)
    INTEGER(jpim),           INTENT(OUT)   :: err
    TYPE(rttov_test_struct), INTENT(IN)    :: ts
    CHARACTER(LEN=*),        INTENT(IN)    :: f
    TYPE(rttov_opt_param),   INTENT(INOUT) :: opt_param
    INTEGER(jpim),           INTENT(INOUT) :: nmom, nphangle
    LOGICAL(jplm),           INTENT(IN)    :: dims_only

    INTEGER(jpim) :: file_id, imult, chanlo, chanhi, nchannels1, lay, i

    TRY

    CALL rttov_get_lun (file_id)
    OPEN(file_id, file = TRIM(f), form = 'formatted', status = 'old', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    READ(file_id, *) nmom, nphangle

    IF (.NOT. dims_only) THEN
      nchannels1 = ts%defn%nchannels / ts%defn%mult

      READ(file_id, *) opt_param%abs(:,1:nchannels1)
      READ(file_id, *) opt_param%sca(:,1:nchannels1)
      READ(file_id, *) opt_param%bpr(:,1:nchannels1)
      opt_param%legcoef = 0.
      opt_param%pha = 0.
      IF (opt_param%nmom > 0) THEN
        DO i = 1, nchannels1
          DO lay = 1, ts%defn%nlevels - 1
            IF (opt_param%sca(lay,i) > 0._jprb) READ(file_id, *) opt_param%legcoef(:,lay,i)
          ENDDO
        ENDDO
      ENDIF
      IF (nphangle > 0) THEN
        READ(file_id, *) opt_param%phangle
        DO i = 1, nchannels1
          IF (ts%data%coefs%coef%ss_val_chn(ts%data%chanprof(i)%chan) > 0) THEN
            DO lay = 1, ts%defn%nlevels - 1
              IF (opt_param%sca(lay,i) > 0._jprb) READ(file_id, *) opt_param%pha(:,lay,i)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun (file_id)

    IF (dims_only) RETURN

    DO imult = 2, ts%defn%mult
       chanlo = (imult - 1) * nchannels1 + 1
       chanhi = imult * nchannels1
       opt_param%abs(:,chanlo:chanhi)  = opt_param%abs(:,1:nchannels1)
       opt_param%sca(:,chanlo:chanhi)  = opt_param%sca(:,1:nchannels1)
       opt_param%bpr(:,chanlo:chanhi)  = opt_param%bpr(:,1:nchannels1)
       opt_param%legcoef(:,:,chanlo:chanhi) = opt_param%legcoef(:,:,1:nchannels1)
       IF (nphangle > 0) THEN
         opt_param%pha(:,:,chanlo:chanhi) = opt_param%pha(:,:,1:nchannels1)
       ENDIF
    ENDDO
    CALL rttov_init_opt_param(err, ts%defn%opts, opt_param)
    THROWM(err.NE.0,"Failure initialising opt_param structure")
    CATCH
  END SUBROUTINE read_opt_param

  SUBROUTINE alloc_data(err, ts, nprofiles, dat, dat_k, do_ad)

    INTEGER(jpim),            INTENT(OUT)             :: err
    TYPE(rttov_test_struct),  INTENT(IN)              :: ts
    INTEGER(jpim),            INTENT(IN)              :: nprofiles
    TYPE(rttov_input_data),   INTENT(INOUT)           :: dat
    TYPE(rttov_k_input_data), INTENT(INOUT), OPTIONAL :: dat_k
    LOGICAL(jplm),            INTENT(IN),    OPTIONAL :: do_ad

    LOGICAL(jplm) :: ldo_ad

    TRY

    ldo_ad = .FALSE._jplm
    IF (PRESENT(do_ad)) ldo_ad = do_ad

    ALLOCATE(                              &
      dat%emissivity(ts%defn%nchannels),   &
      dat%reflectance(ts%defn%nchannels),  &
      dat%profiles(nprofiles), stat = err)
    THROWM(err.NE.0,"Cannot allocate dat")

    IF (ts%defn%do_rttovscatt) THEN
      ALLOCATE(dat%cld_profiles(nprofiles), dat%rttovscatt_cfrac(nprofiles), stat = err)
      THROWM(err.NE.0,"Cannot allocate dat")
      IF (ts%defn%radar) THEN
        CALL rttov_alloc_reflectivity(err, ts%defn%nchannels, dat%reflectivity, ts%defn%nlevels, 1_jpim)
        THROWM(err.NE.0,"Cannot allocate dat")
      ENDIF
      IF (ts%defn%calc_emis_terms) THEN
        CALL rttov_alloc_emis_ret_terms(err, ts%defn%nchannels, dat%emis_terms, 1_jpim)
        THROWM(err.NE.0,"Cannot allocate dat")
      ENDIF
    ENDIF

    CALL rttov_init_emis_refl(dat%emissivity, dat%reflectance)

    CALL rttov_alloc_prof(err, nprofiles, dat%profiles, ts%defn%nlevels, ts%defn%opts, 1_jpim, &
                          init = .TRUE._jplm, coefs = ts%data%coefs)
    THROW(err.NE.0)

    CALL rttov_alloc_rad(err, ts%defn%nchannels, dat%radiance, ts%defn%nlevels, 1_jpim, dat%radiance2, init = .TRUE._jplm)
    THROW(err.NE.0)

    CALL rttov_alloc_transmission(err, dat%transmission, ts%defn%nlevels, ts%defn%nchannels, &
                                  1_jpim, init = .TRUE._jplm)
    THROW(err.NE.0)

    IF (ts%defn%opts%rt_ir%user_aer_opt_param) THEN
      CALL rttov_alloc_opt_param(err, dat%aer_opt_param, ts%defn%nchannels, &
                                 ts%defn%nlevels, ts%defn%aer_nmom, ts%defn%aer_nphangle, 1_jpim)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%opts%rt_ir%user_cld_opt_param) THEN
      CALL rttov_alloc_opt_param(err, dat%cld_opt_param, ts%defn%nchannels, &
                                 ts%defn%nlevels, ts%defn%cld_nmom, ts%defn%cld_nphangle, 1_jpim)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_rttovscatt) THEN
      CALL rttov_alloc_scatt_prof(err, nprofiles, dat%cld_profiles, ts%defn%nlevels, ts%data%coefs_scatt%nhydro, &
            ts%defn%nhydro_frac, 1_jpim, init = .TRUE._jplm)
    ENDIF

    IF (ts%defn%opts%rt_ir%pc%addpc) THEN
      IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
        CALL rttov_alloc_pccomp(err, dat%pccomp, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, init = .TRUE._jplm, &
                                nchannels_rec = ts%defn%nchannels_rec * ts%defn%nprofiles, &
                                opts = ts%defn%opts, nlevels = ts%defn%nlevels)
      ELSE
        CALL rttov_alloc_pccomp(err, dat%pccomp, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, init = .TRUE._jplm, &
                                opts = ts%defn%opts, nlevels = ts%defn%nlevels)
      ENDIF
      THROW(err.NE.0)
    ENDIF

    IF (ldo_ad) THEN
      IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          CALL rttov_alloc_pccomp(err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                  init = .TRUE._jplm, nchannels_rec = ts%defn%nchannels_rec * ts%defn%nprofiles, &
                                  opts = ts%defn%opts, nlevels = ts%defn%nlevels)
        ELSE
          CALL rttov_alloc_pccomp(err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                  init = .TRUE._jplm, opts = ts%defn%opts, nlevels = ts%defn%nlevels)
        ENDIF
      ELSEIF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
        CALL rttov_alloc_reflectivity(err, ts%defn%nchannels, dat%reflectivity_saved, ts%defn%nlevels, 1_jpim, &
                                      init=.TRUE._jplm)
      ELSE
        CALL rttov_alloc_rad(err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels, 1_jpim, &
                             init=.TRUE._jplm)
      ENDIF
      THROW(err.NE.0)
    ENDIF

    IF (PRESENT(dat_k)) THEN
      IF (ts%defn%opts%htfrtc_opts%htfrtc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          ALLOCATE(dat_k%profiles_k_rec(ts%defn%nchannels_rec * ts%defn%nprofiles), stat = err)
          THROW(err.NE.0)
          CALL rttov_alloc_prof(err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, ts%defn%nlevels, &
              ts%defn%opts, 1_jpim, init = .TRUE._jplm)
          THROW(err.NE.0)
        ENDIF
        ALLOCATE(dat_k%profiles_k_pc(ts%defn%npcscores * ts%defn%nprofiles), stat = err)
        THROW(err.NE.0)
        CALL rttov_alloc_prof(err, ts%defn%npcscores * ts%defn%nprofiles, &
            dat_k%profiles_k_pc, &
            ts%defn%nlevels, ts%defn%opts, 1_jpim, init = .TRUE._jplm)
        THROW(err.NE.0)
      ELSE IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          ALLOCATE(dat_k%profiles_k_rec(ts%defn%nchannels_rec * ts%defn%nprofiles), stat = err)
          THROW(err.NE.0)
          CALL rttov_alloc_prof(err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, ts%defn%nlevels, &
              ts%defn%opts, 1_jpim, init = .TRUE._jplm, coefs = ts%data%coefs)
          THROW(err.NE.0)
        ELSE
          ALLOCATE(dat_k%profiles_k_pc(ts%defn%npcscores * ts%defn%nprofiles), stat = err)
          THROW(err.NE.0)
          CALL rttov_alloc_prof(err, ts%defn%npcscores * ts%defn%nprofiles, &
              dat_k%profiles_k_pc, &
              ts%defn%nlevels, ts%defn%opts, 1_jpim, init = .TRUE._jplm, coefs = ts%data%coefs)
          THROW(err.NE.0)
        ENDIF
      ENDIF
    ENDIF

    CATCH
  END SUBROUTINE alloc_data

  SUBROUTINE dealloc_data(err, ts, nprofiles, dat, dat_k, do_ad)

    INTEGER(jpim),            INTENT(OUT)             :: err
    TYPE(rttov_test_struct),  INTENT(IN)              :: ts
    INTEGER(jpim),            INTENT(IN)              :: nprofiles
    TYPE(rttov_input_data),   INTENT(INOUT)           :: dat
    TYPE(rttov_k_input_data), INTENT(INOUT), OPTIONAL :: dat_k
    LOGICAL(jplm),            INTENT(IN),    OPTIONAL :: do_ad

    LOGICAL(jplm) :: ldo_ad

    TRY

    ldo_ad = .FALSE._jplm
    IF (PRESENT(do_ad)) ldo_ad = do_ad

    IF (ts%defn%do_rttovscatt) THEN
      CALL rttov_alloc_scatt_prof(err, nprofiles, dat%cld_profiles, ts%defn%nlevels, ts%data%coefs_scatt%nhydro, &
            ts%defn%nhydro_frac, 0_jpim)
      THROW(err.NE.0)
      IF (ts%defn%radar) THEN
        CALL rttov_alloc_reflectivity(err, ts%defn%nchannels, dat%reflectivity, ts%defn%nlevels, 0_jpim)
        THROWM(err.NE.0,"Cannot allocate dat")
      ENDIF
      IF (ts%defn%calc_emis_terms) THEN
        CALL rttov_alloc_emis_ret_terms(err, ts%defn%nchannels, dat%emis_terms, 0_jpim)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    CALL rttov_alloc_prof(err, nprofiles, dat%profiles, ts%defn%nlevels, ts%defn%opts, 0_jpim)
    THROW(err.NE.0)
    CALL rttov_alloc_rad(err, ts%defn%nchannels, dat%radiance, ts%defn%nlevels, 0_jpim, dat%radiance2)
    THROW(err.NE.0)
    CALL rttov_alloc_transmission(err, dat%transmission, ts%defn%nlevels, ts%defn%nchannels, 0_jpim)
    THROW(err.NE.0)

    IF (ts%defn%opts%rt_ir%user_aer_opt_param) THEN
      CALL rttov_alloc_opt_param(err, dat%aer_opt_param, ts%defn%nchannels, &
                                 ts%defn%nlevels, ts%defn%aer_nmom, 1_jpim, 0_jpim)
      THROW(err.NE.0)
    ENDIF
    IF (ts%defn%opts%rt_ir%user_cld_opt_param) THEN
      CALL rttov_alloc_opt_param(err, dat%cld_opt_param, ts%defn%nchannels, &
                                 ts%defn%nlevels, ts%defn%cld_nmom, 1_jpim, 0_jpim)
      THROW(err.NE.0)
    ENDIF

    IF (ts%defn%do_rttovscatt) THEN
      DEALLOCATE(dat%cld_profiles, dat%rttovscatt_cfrac, stat = err)
      THROWM(err.NE.0,"Cannot allocate dat")
    ENDIF

    DEALLOCATE(dat%emissivity, dat%reflectance, dat%profiles, stat = err)
    THROWM(err.NE.0,"Deallocation of dat failed")

    IF (ts%defn%opts%rt_ir%pc%addpc) THEN
      CALL rttov_alloc_pccomp (err, dat%pccomp, ts%defn%npcscores, 0_jpim)
      THROW(err.NE.0)
    ENDIF

    IF (ldo_ad) THEN
      IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        CALL rttov_alloc_pccomp (err, dat%pccomp_saved, ts%defn%npcscores, 0_jpim)
      ELSEIF (ts%defn%do_rttovscatt .AND. ts%defn%radar) THEN
        CALL rttov_alloc_reflectivity(err, ts%defn%nchannels, dat%reflectivity_saved, ts%defn%nlevels, 0_jpim)
      ELSE
        CALL rttov_alloc_rad(err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels, 0_jpim)
      ENDIF
      THROW(err.NE.0)
    ENDIF

    IF (PRESENT(dat_k)) THEN
      IF (ts%defn%opts%htfrtc_opts%htfrtc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          CALL rttov_alloc_prof(err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, &
                                ts%defn%nlevels, ts%defn%opts, 0_jpim)
          THROW(err.NE.0)

          DEALLOCATE(dat_k%profiles_k_rec, stat = err)
          THROW(err.NE.0)
        ENDIF
        CALL rttov_alloc_prof(err, ts%defn%npcscores * ts%defn%nprofiles, dat_k%profiles_k_pc, &
                              ts%defn%nlevels, ts%defn%opts, 0_jpim)
        THROW(err.NE.0)

        DEALLOCATE(dat_k%profiles_k_pc, stat = err)
        THROW(err.NE.0)
      ELSE IF (ts%defn%opts%rt_ir%pc%addpc) THEN
        IF (ts%defn%opts%rt_ir%pc%addradrec) THEN
          CALL rttov_alloc_prof(err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, &
                                ts%defn%nlevels, ts%defn%opts, 0_jpim)
          THROW(err.NE.0)

          DEALLOCATE(dat_k%profiles_k_rec, stat = err)
          THROW(err.NE.0)
        ELSE
          CALL rttov_alloc_prof(err, ts%defn%npcscores * ts%defn%nprofiles, dat_k%profiles_k_pc, &
                                ts%defn%nlevels, ts%defn%opts, 0_jpim)
          THROW(err.NE.0)

          DEALLOCATE(dat_k%profiles_k_pc, stat = err)
          THROW(err.NE.0)
        ENDIF
      ENDIF
    ENDIF

    CATCH
  END SUBROUTINE dealloc_data

  !> Print all output for a direct, TL, AD or K model test
  !! @param[out]    err             return status
  !! @param[in]     opts            options structure used for test
  !! @param[in]     path            output directory (usually ends "direct", "tl", etc)
  !! @param[in]     suffix          file suffix (empty for direct model, otherwise "_tl", "_ad", etc)
  !! @param[in]     adk             flag to indicate AD or K model data
  !! @param[in]     calcrad2        flag to indicate radiance2 should be output
  !! @param[in]     do_rttovscatt   flag to indicate RTTOV-SCATT simulations
  !! @param[in]     radar           flag to indicate RTTOV-SCATT radar simulations
  !! @param[in]     calcemisterms   flag to indicate RTTOV-SCATT emis terms should be output (direct model only)
  !! @param[in]     data            structure containing data to write out
  !! @param[in]     data_k          additional structure containing PC K output for PC simulations, optional
  SUBROUTINE printoutput(err, opts, path, suffix, adk, calcrad2, do_rttovscatt, radar, calcemisterms, data, data_k)

    USE rttov_unix_env, ONLY: rttov_lower_case, rttov_upper_case

    INTEGER(jpim),              INTENT(OUT)          :: err
    TYPE(rttov_options),        INTENT(IN)           :: opts
    CHARACTER(LEN=*),           INTENT(IN)           :: path
    CHARACTER(LEN=*),           INTENT(IN)           :: suffix
    LOGICAL(jplm),              INTENT(IN)           :: adk
    LOGICAL(jplm),              INTENT(IN)           :: calcrad2
    LOGICAL(jplm),              INTENT(IN)           :: do_rttovscatt
    LOGICAL(jplm),              INTENT(IN)           :: radar
    LOGICAL(jplm),              INTENT(IN)           :: calcemisterms
    TYPE(rttov_input_data),     INTENT(IN)           :: data
    TYPE(rttov_k_input_data),   INTENT(IN), OPTIONAL :: data_k

    CHARACTER(LEN=10) :: suffixuc, suffixlc
    TRY

    CALL rttov_upper_case(suffixuc, suffix)
    CALL rttov_lower_case(suffixlc, suffix)

    IF (calcrad2) THEN
      ! Print radiance and radiance2
      CALL printradiance(err, data%radiance, "RADIANCE"//TRIM(suffixuc), &
              TRIM(path)//"/radiance"//TRIM(suffixlc)//".txt", data%radiance2)
      THROW(err.NE.0)
    ELSE
      IF (.NOT. adk) THEN
        ! Direct/TL - print radiance(_tl) only
        CALL printradiance(err, data%radiance, "RADIANCE"//TRIM(suffixuc), &
                TRIM(path)//"/radiance"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (do_rttovscatt .AND. calcemisterms) THEN
      CALL printemisterms(err, data%emis_terms, "SCATT_EMIS_TERMS"//TRIM(suffixuc), &
              TRIM(path)//"/scatt_emis_terms"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)
    ENDIF

    IF (do_rttovscatt .AND. radar .AND. .NOT. adk) THEN
      CALL printreflectivity(err, data%reflectivity, "REFLECTIVITY"//TRIM(suffixuc), &
              TRIM(path)//"/reflectivity"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)
    ENDIF


    IF (adk) THEN
      ! AD/K - print profiles_ad/k, cld_profiles_ad/k, aer/cld_opt_param_ad/k, emissivity_ad/k, reflectance_ad/k
      CALL printprofiles(err, data%profiles, "PROFILES"//TRIM(suffixuc), &
              TRIM(path)//"/profiles"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)

      IF (do_rttovscatt) THEN
        CALL printcldprofiles(err, data%cld_profiles, "CLD_PROFILES"//TRIM(suffixuc), &
                TRIM(path)//"/cld_profiles"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF

      IF (opts%rt_ir%user_aer_opt_param) THEN
        CALL printoptparam(err, data%aer_opt_param, "AER_OPT_PARAM"//TRIM(suffixuc), &
                TRIM(path)//"/aer_opt_param"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF

      IF (opts%rt_ir%user_cld_opt_param) THEN
        CALL printoptparam(err, data%cld_opt_param, "CLD_OPT_PARAM"//TRIM(suffixuc), &
                TRIM(path)//"/cld_opt_param"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF

      CALL printsurfvar(err, data%emissivity%emis_in, "EMISSIVITY"//TRIM(suffixuc), &
              TRIM(path)//"/emissivity"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)

      IF (opts%rt_ir%addsolar) THEN
        CALL printsurfvar(err, data%reflectance%refl_in, "REFLECTANCE"//TRIM(suffixuc), &
                TRIM(path)//"/reflectance"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)

        CALL printsurfvar(err, data%reflectance%diffuse_refl_in, "REFLECTANCE_DIFFUSE"//TRIM(suffixuc), &
                TRIM(path)//"/reflectance_diffuse"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF
    ELSE
      ! Direct/TL - print transmission(_tl), emissivity_out(_tl), reflectance_out(_tl)
      CALL printtransmission(err, data%transmission, "TRANSMISSION"//TRIM(suffixuc), &
              TRIM(path)//"/transmission"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)

      CALL printsurfvar(err, data%emissivity%emis_out, "EMISSIVITY_OUT"//TRIM(suffixuc), &
              TRIM(path)//"/emissivity_out"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)

      IF (opts%rt_ir%addsolar) THEN
        CALL printsurfvar(err, data%reflectance%refl_out, "REFLECTANCE_OUT"//TRIM(suffixuc), &
                TRIM(path)//"/reflectance_out"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)

        CALL printsurfvar(err, data%reflectance%diffuse_refl_out, "REFLECTANCE_DIFFUSE_OUT"//TRIM(suffixuc), &
                TRIM(path)//"/reflectance_diffuse_out"//TRIM(suffixlc)//".txt")
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (opts%rt_all%do_lambertian) THEN
      CALL printsurfvar(err, data%emissivity%specularity, "SPECULARITY"//TRIM(suffixuc), &
              TRIM(path)//"/specularity"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)
    ENDIF

    IF (opts%rt_ir%pc%addpc) THEN
      ! All models - print pccomp(_tl/ad/k)
      CALL printpcscores(err, data%pccomp, "PCSCORES"//TRIM(suffixuc), &
              TRIM(path)//"/pcscores"//TRIM(suffixlc)//".txt")
      THROW(err.NE.0)

      IF (PRESENT(data_k)) THEN
        ! K only - print PC K or rec rad K
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL printprofiles(err, data_k%profiles_k_rec, "PROFILES_K_REC", &
                  TRIM(path)//"/profiles_k_rec.txt")
          THROW(err.NE.0)
        ELSE
          CALL printprofiles(err, data_k%profiles_k_pc, "PROFILES_K_PC", &
                  TRIM(path)//"/profiles_k_pc.txt")
          THROW(err.NE.0)
        ENDIF
      ENDIF
    ENDIF

    CATCH
  END SUBROUTINE

END PROGRAM
