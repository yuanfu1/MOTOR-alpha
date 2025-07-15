! Description:
!> @file
!!   Creates aer_opt_param.txt and cld_opt_param.txt files for the test suite
!!   from input aerosl.txt and cloud.txt files and the pre-defined aerosol and
!!   cloud optical properties.
!
!> @brief
!!   Creates aer_opt_param.txt and cld_opt_param.txt files for the test suite
!!   from input aerosl.txt and cloud.txt files and the pre-defined aerosol and
!!   cloud optical properties.
!!
!! @details
!!   This is called by rttov_test.pl to automatically generate optical property
!!   input files for the test suite. The output from the direct model for the
!!   original test (using scaer/sccld coefs) and for the test using the
!!   resulting optical property files should be identical.
!!
!!   Usage:
!!   $ rttov_make_opt_param --rtcoef_file  coef_file
!!                          --scaer_file   aer_file
!!                          --sccld_file   cld_file
!!                          --test_dir     test_dir
!!                          --nchanprof    nchanprof
!!                          --nlayers      nlayers
!!   where
!!     coef_file is the instrument rtcoef file
!!     aer_file is the aerosol coef file (mandatory if
!!       an aerosl.txt test file is present)
!!     cld_file is the cloud coef file (mandatory if
!!       a cloud.txt test file is present)
!!     test_dir is the top level test folder containing
!!       cloud and/or aerosol profiles
!!       (e.g. tests.0/seviri/081/in/)
!!     nchanprof is the size of the chanprof array
!!       associated with the test. The channels.txt and
!!       lprofiles.txt files determine the channels/
!!       profiles defined in the output file.
!!     nlayers is the number of layers in the profiles
!!
!!   This outputs an aer_opt_param.txt file if there is
!!   an aerosl.txt present in test_dir/profiles/001/atm/
!!   and it outputs a cld_opt_param.txt file if there is
!!   a cloud.txt present in test_dir/profiles/001/atm/.
!!
!!   NB For clouds it requires that the phase functions in the coef
!!      file are defined on the phangle_hires angle grid defined in
!!      rttov_const: this should always be true since this same grid
!!      is used in the sccld coef generation.
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
PROGRAM rttov_make_opt_param

  USE rttov_getoptions, ONLY : getoption
  USE rttov_unix_env, ONLY: rttov_iargc
  USE rttov_lun, ONLY : rttov_get_lun, rttov_put_lun

  USE parkind1, ONLY : jprb, jpim, jplm

  USE rttov_types, ONLY : &
      rttov_options,      &
      rttov_coefs,        &
      rttov_opt_param,    &
      rttov_chanprof,     &
      rttov_profile,      &
      rttov_profile_aux,  &
      rttov_optp_data,    &
      rttov_skin

  USE rttov_const, ONLY: &
      deg2rad,              &
      nphangle_lores,       &
      phangle_lores,        &
      nphangle_hires,       &
      phangle_hires,        &
      baran_ngauss,         &
      clw_scheme_opac,      &
      ice_scheme_baum,      &
      ice_scheme_baran2018, &
      nwcl_max

  USE rttov_scattering_mod, ONLY: &
      spline_interp,              &
      normalise,                  &
      calc_legendre_coef_gauss

  IMPLICIT NONE

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_aux_prof.interface"
#include "rttov_baran2014_calc_optpar.interface"
#include "rttov_baran2018_calc_optpar.interface"
#include "rttov_baran_calc_phase.interface"
#include "rttov_convert_profile_units.interface"
#include "rttov_profaux_cldaer.interface"

  CHARACTER(256) :: rtcoef_file
  CHARACTER(256) :: scaer_file
  CHARACTER(256) :: sccld_file
  CHARACTER(256) :: test_dir
  INTEGER(jpim)  :: nchanprof
  INTEGER(jpim)  :: nlayers, nlevels
  INTEGER(jpim)  :: nmom, nmomcalc
  INTEGER(jpim)  :: nphangle

  LOGICAL(jplm)  :: exists
  CHARACTER(256) :: pathstr
  INTEGER(jpim)  :: nprof, prof, chan, phchan, lay
  INTEGER(jpim)  :: i, n, iae, iwc, k1, k2
  INTEGER(jpim)  :: naer, ncld
  LOGICAL(jplm)  :: thermal, solar
  INTEGER(jpim)  :: err, file_id
  INTEGER(jpim)  :: gas_units
  LOGICAL(jplm)  :: mmr_cldaer
  INTEGER(jpim)  :: clw_scheme, clwde_param, ice_scheme, icede_param
  REAL(jprb)     :: dgfrac, asym
  INTEGER(jpim)  :: thisnphangle
  REAL(jprb)     :: thisphangle(nphangle_hires), thiscosphangle(nphangle_hires)
  REAL(jprb)     :: baran_phfn_interp(baran_ngauss)
  TYPE(rttov_skin) :: k0

  REAL(jprb), ALLOCATABLE :: absch(:)
  REAL(jprb), ALLOCATABLE :: scach(:)
  REAL(jprb), ALLOCATABLE :: bparh(:)
  REAL(jprb), ALLOCATABLE :: legcoef(:,:)
  REAL(jprb), ALLOCATABLE :: phfn(:,:)

  TYPE(rttov_options)               :: opts
  TYPE(rttov_coefs)                 :: coefs
  TYPE(rttov_opt_param)             :: opt_param
  TYPE(rttov_chanprof), ALLOCATABLE :: chanprof(:)
  TYPE(rttov_profile),  ALLOCATABLE :: profiles(:), profiles_int(:)
  TYPE(rttov_profile_aux)           :: aux
  TYPE(rttov_optp_data), POINTER    :: optp_data

  NAMELIST / clw_scheme_nml / clw_scheme, clwde_param
  NAMELIST / ice_scheme_nml / icede_param, ice_scheme
  NAMELIST / units / gas_units
  NAMELIST / cldaer_units / mmr_cldaer
  NAMELIST / skin / k0

! ----------------------------------------------------------------------------

! ------------------------------------
! Process arguments
! ------------------------------------

  IF (rttov_iargc() == 0) THEN
    PRINT *, 'Usage: --rtcoef_file   input rtcoef file name'
    PRINT *, '       --scaer_file    input aerosol file name (optional)'
    PRINT *, '       --sccld_file    input cloud file name (optional)'
    PRINT *, '       --test_dir      path to test directory'
    PRINT *, '       --nchanprof     integer'
    PRINT *, '       --nlayers       integer'
    STOP
  ENDIF

  CALL getoption('--rtcoef_file', rtcoef_file, mnd=.TRUE._jplm)
  INQUIRE(FILE=rtcoef_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find rtcoef file: '//TRIM(rtcoef_file)
    STOP
  ENDIF

  scaer_file = ''
  CALL getoption('--scaer_file', scaer_file)
  IF (TRIM(scaer_file) /= '') THEN
    INQUIRE(FILE=scaer_file, EXIST=exists)
    IF (.NOT. exists) THEN
      PRINT *, 'Cannot find scaer file: '//TRIM(scaer_file)
      STOP
    ENDIF
  ENDIF

  sccld_file = ''
  CALL getoption('--sccld_file', sccld_file)
  IF (TRIM(sccld_file) /= '') THEN
    INQUIRE(FILE=sccld_file, EXIST=exists)
    IF (.NOT. exists) THEN
      PRINT *, 'Cannot find sccld file: '//TRIM(sccld_file)
      STOP
    ENDIF
  ENDIF

  CALL getoption('--test_dir', test_dir, mnd=.TRUE._jplm)
  INQUIRE(FILE=TRIM(test_dir)//'/in/channels.txt', EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find valid test dir: '//TRIM(test_dir)//'/in/channels.txt'
    STOP
  ENDIF

  CALL getoption('--nchanprof', nchanprof, mnd=.TRUE._jplm)
  IF (nchanprof <= 0) THEN
    nchanprof = 1
  ENDIF

  CALL getoption('--nlayers', nlayers, mnd=.TRUE._jplm)
  IF (nlayers <= 0) THEN
    nlayers = 1
  ENDIF
  nlevels = nlayers + 1


! ------------------------------------
! Aerosols
! ------------------------------------

  IF (TRIM(scaer_file) /= '') THEN

  ! ------------------------------------
  ! Read coefs
  ! ------------------------------------
    opts%rt_ir%addaerosl = .TRUE.
    opts%rt_ir%addclouds = .FALSE.

    CALL rttov_read_coefs(err, coefs, opts, file_coef=rtcoef_file, file_scaer=scaer_file)
    opts%rt_ir%addsolar = (coefs%coef%fmv_model_ver >= 9)

  ! ------------------------------------
  ! Read in data from test directory
  ! ------------------------------------

    CALL rttov_get_lun(file_id)

    ALLOCATE(chanprof(nchanprof))
    OPEN(file_id, file=TRIM(test_dir)//'/in/channels.txt', form='formatted', status='old')
    READ(file_id, *) chanprof(:)%chan
    CLOSE(file_id)
    OPEN(file_id, file=TRIM(test_dir)//'/in/lprofiles.txt', form='formatted', status='old')
    READ(file_id, *) chanprof(:)%prof
    CLOSE(file_id)

    nprof = MAXVAL(chanprof(:)%prof)

    ALLOCATE(profiles(nprof), profiles_int(nprof))
    CALL rttov_alloc_prof(err, nprof, profiles, nlevels, opts, 1_jpim, &
                          coefs, .TRUE._jplm)
    CALL rttov_alloc_prof(err, nprof, profiles_int, nlevels, opts, 1_jpim, &
                          coefs, .TRUE._jplm)
    CALL rttov_alloc_aux_prof(err, nprof, nlevels, aux, opts, coefs%coef, &
                              1_jpim, .TRUE._jplm)

    DO prof = 1, nprof
      WRITE(pathstr, '(A,I3.3)') TRIM(test_dir)//'/in/profiles/', prof

      OPEN(file_id, file=TRIM(pathstr)//'/atm/aerosl.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%aerosols(:,:)
      CLOSE(file_id)

      OPEN(file_id, file=TRIM(pathstr)//'/atm/p.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%p(:)
      CLOSE(file_id)
      profiles(prof)%s2m%p = profiles(prof)%p(nlevels)
      OPEN(file_id, file=TRIM(pathstr)//'/atm/t.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%t(:)
      CLOSE(file_id)
      OPEN(file_id, file=TRIM(pathstr)//'/atm/q.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%q(:)
      CLOSE(file_id)

      OPEN(file_id, file=TRIM(pathstr)//'/atm/mmr_cldaer.txt', form='formatted', status='old')
      READ(file_id, nml=cldaer_units)
      CLOSE(file_id)
      profiles(prof)%mmr_cldaer = mmr_cldaer

      INQUIRE(file=TRIM(pathstr)//'/gas_units.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/gas_units.txt', form='formatted', status='old')
        READ(file_id, nml=units)
        CLOSE(file_id)
        profiles(prof)%gas_units = gas_units
      ENDIF
    ENDDO

    CALL rttov_put_lun(file_id)

    CALL rttov_convert_profile_units(opts, coefs, profiles, profiles_int)

  ! ------------------------------------
  ! Init opt params
  ! ------------------------------------
    naer = SIZE(profiles(1)%aerosols, DIM=1)
    nmom = coefs%coef_scatt%optp_aer%maxnmom
    opts%rt_ir%dom_nstreams = nmom
    nphangle = coefs%coef_scatt%optp_aer%nphangle
    CALL rttov_alloc_opt_param(err, opt_param, nchanprof, nlevels, nmom, &
                               nphangle, 1_jpim)

    ALLOCATE(absch(naer), scach(naer), bparh(naer), &
             legcoef(0:nmom,naer), phfn(nphangle,naer))

  ! ------------------------------------
  ! Compute layer relative humidities (only required for aerosols)
  ! ------------------------------------
    CALL rttov_profaux_cldaer(opts, profiles, profiles_int, aux)

  ! ------------------------------------
  ! Main loop
  ! ------------------------------------

    DO i = 1, nchanprof
      chan = chanprof(i)%chan
      prof = chanprof(i)%prof

      thermal = coefs%coef%ss_val_chn(chan) < 2
      solar   = coefs%coef%ss_val_chn(chan) > 0

      DO lay = 1, nlayers

        ! ------------------------------------
        ! Calculate layer abs, sca and bpr
        ! ------------------------------------

        DO iae = 1, naer
          IF (profiles_int(prof)%aerosols(iae,lay) > 0.) THEN
            CALL aer_interp_relhum(chan, nmom, iae, solar, &
                                   coefs%coef_scatt%optp_aer, aux%relhum(lay,prof), &
                                   absch(iae), scach(iae), bparh(iae), legcoef(:,iae), phfn(:,iae))
          ELSE
            absch(iae) = 0.
            scach(iae) = 0.
            bparh(iae) = 0.
            legcoef(:,iae) = 0.
            phfn(:,iae) = 0.
          ENDIF
        ENDDO

        opt_param%abs(lay,i) = SUM((/ (absch(n) * profiles_int(prof)%aerosols(n,lay), n = 1, naer) /))
        opt_param%sca(lay,i) = SUM((/ (scach(n) * profiles_int(prof)%aerosols(n,lay), n = 1, naer) /))
        IF (thermal .AND. opt_param%sca(lay,i) > 0.) THEN
          opt_param%bpr(lay,i) = SUM((/ (bparh(n) * scach(n) * &
                                         profiles_int(prof)%aerosols(n,lay), n = 1, naer) /)) / opt_param%sca(lay,i)
        ELSE
          opt_param%bpr(lay,i) = 0.
        ENDIF

        opt_param%nmom = nmom
        opt_param%legcoef(:,lay,i) = 0.
        opt_param%pha(:,lay,i) = 0.

        IF (opt_param%sca(lay,i) > 0.) THEN
          DO iae = 1, naer
            opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) + &
                                         legcoef(:,iae) * scach(iae) * profiles_int(prof)%aerosols(iae,lay)
          ENDDO
          opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) / opt_param%sca(lay,i)
        ENDIF

        IF (solar) THEN

          ! ------------------------------------
          ! Calculate layer phase fn
          ! ------------------------------------

          IF (opt_param%sca(lay,i) > 0.) THEN
            DO iae = 1, naer
              opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) + &
                                       phfn(:,iae) * scach(iae) * profiles_int(prof)%aerosols(iae,lay)
            ENDDO
            opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) / opt_param%sca(lay,i)
          ENDIF

        ENDIF

      ENDDO ! layers
    ENDDO ! chanprof

  ! ------------------------------------
  ! Write out opt_param file
  ! ------------------------------------

    CALL rttov_get_lun(file_id)

    OPEN(file_id, file=TRIM(test_dir)//'/in/aer_opt_param.txt', form='formatted')
    WRITE(file_id, *) opt_param%nmom, nphangle
    WRITE(file_id, '(10E20.12)') opt_param%abs(:,:)
    WRITE(file_id, '(10E20.12)') opt_param%sca(:,:)
    WRITE(file_id, '(10E20.12)') opt_param%bpr(:,:)
    IF (opt_param%nmom > 0) THEN
      DO i = 1, nchanprof
        DO lay = 1, nlayers
          IF (opt_param%sca(lay,i) > 0._jprb) WRITE(file_id, '(10E20.12)') opt_param%legcoef(:,lay,i)
        ENDDO
      ENDDO
    ENDIF
    IF (nphangle > 0) THEN
      WRITE(file_id, '(10E20.12)') coefs%coef_scatt%optp_aer%phangle
      DO i = 1, nchanprof
        IF (coefs%coef%ss_val_chn(chanprof(i)%chan) > 0) THEN
          DO lay = 1, nlayers
            IF (opt_param%sca(lay,i) > 0._jprb) WRITE(file_id, '(10E20.12)') opt_param%pha(:,lay,i)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    CLOSE(file_id)

    CALL rttov_put_lun(file_id)

  ! ------------------------------------
  ! Clean up
  ! ------------------------------------

    DEALLOCATE(chanprof, absch, scach, bparh, legcoef, phfn)

    CALL rttov_alloc_opt_param(err, opt_param, nchanprof, nlevels, nmom, &
                               nphangle, 0_jpim)

    CALL rttov_alloc_prof(err, nprof, profiles, nlevels, opts, 0_jpim, coefs)
    DEALLOCATE(profiles)
    CALL rttov_alloc_prof(err, nprof, profiles_int, nlevels, opts, 0_jpim, coefs)
    DEALLOCATE(profiles_int)
    CALL rttov_alloc_aux_prof(err, nprof, nlevels, aux, opts, coefs%coef, 0_jpim)

    CALL rttov_dealloc_coefs(err, coefs)

  ENDIF ! scaer file


! ------------------------------------
! Clouds
! ------------------------------------

  IF (TRIM(sccld_file) /= '') THEN

  ! ------------------------------------
  ! Read coefs
  ! ------------------------------------

    opts%rt_ir%addaerosl = .FALSE.
    opts%rt_ir%addclouds = .TRUE.

    CALL rttov_read_coefs(err, coefs, opts, file_coef=rtcoef_file, file_sccld=sccld_file)
    opts%rt_ir%addsolar = (coefs%coef%fmv_model_ver >= 9)

    IF (coefs%coef_scatt%optp_wcl_opac%nphangle > 0) THEN
      IF (coefs%coef_scatt%optp_wcl_opac%nphangle /= nphangle_hires .OR. &
          coefs%coef_scatt%optp_icl_baum%nphangle /= nphangle_hires) THEN
        PRINT *,'Water and ice cloud phase fn angle grids differ from phangle_hires'
        STOP
      ENDIF
      IF (ANY(ABS(coefs%coef_scatt%optp_wcl_opac%phangle - phangle_hires) > 1.E-8_jprb) .OR. &
          ANY(ABS(coefs%coef_scatt%optp_icl_baum%phangle - phangle_hires) > 1.E-8_jprb)) THEN
        PRINT *,'Water and ice cloud phase fn angle grids differ from phangle_hires'
        STOP
      ENDIF
    ENDIF

    IF (coefs%coef_scatt%optp_wcl_deff%nphangle > 0) THEN
      IF (coefs%coef_scatt%optp_wcl_deff%nphangle /= nphangle_hires) THEN
        PRINT *,'Deff clw scheme phase fn angle grid differs from phangle_hires'
        STOP
      ENDIF
      IF (ANY(ABS(coefs%coef_scatt%optp_wcl_deff%phangle - phangle_hires) > 1.E-8_jprb)) THEN
        PRINT *,'Deff clw scheme phase fn angle grid differs from phangle_hires'
        STOP
      ENDIF
    ENDIF

  ! ------------------------------------
  ! Read in data from test directory
  ! ------------------------------------

    CALL rttov_get_lun(file_id)

    ALLOCATE(chanprof(nchanprof))
    OPEN(file_id, file=TRIM(test_dir)//'/in/channels.txt', form='formatted', status='old')
    READ(file_id, *) chanprof(:)%chan
    CLOSE(file_id)
    OPEN(file_id, file=TRIM(test_dir)//'/in/lprofiles.txt', form='formatted', status='old')
    READ(file_id, *) chanprof(:)%prof
    CLOSE(file_id)

    nprof = MAXVAL(chanprof(:)%prof)

    ALLOCATE(profiles(nprof), profiles_int(nprof))
    CALL rttov_alloc_prof(err, nprof, profiles, nlevels, opts, 1_jpim, &
                          coefs, .TRUE._jplm)
    CALL rttov_alloc_prof(err, nprof, profiles_int, nlevels, opts, 1_jpim, &
                          coefs, .TRUE._jplm)
    CALL rttov_alloc_aux_prof(err, nprof, nlevels, aux, opts, coefs%coef, &
                              1_jpim, .TRUE._jplm)

    DO prof = 1, nprof
      WRITE(pathstr, '(A,I3.3)') TRIM(test_dir)//'/in/profiles/', prof

      OPEN(file_id, file=TRIM(pathstr)//'/atm/cloud.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%cloud(:,:)
      CLOSE(file_id)

      OPEN(file_id, file=TRIM(pathstr)//'/atm/p.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%p(:)
      CLOSE(file_id)
      profiles(prof)%s2m%p = profiles(prof)%p(nlevels)
      OPEN(file_id, file=TRIM(pathstr)//'/atm/t.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%t(:)
      CLOSE(file_id)
      OPEN(file_id, file=TRIM(pathstr)//'/atm/q.txt', form='formatted', status='old')
      READ(file_id, *) profiles(prof)%q(:)
      CLOSE(file_id)

      INQUIRE(file=TRIM(pathstr)//'/atm/clwde.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/atm/clwde.txt', form='formatted', status='old')
        READ(file_id, *) profiles(prof)%clwde(:)
        CLOSE(file_id)
      ENDIF

      INQUIRE(file=TRIM(pathstr)//'/atm/clw_scheme.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/atm/clw_scheme.txt', form='formatted', status='old')
        READ(file_id, nml = clw_scheme_nml)
        CLOSE(file_id)
        profiles(prof)%clw_scheme = clw_scheme
        profiles(prof)%clwde_param = clwde_param
      ENDIF

      INQUIRE(file=TRIM(pathstr)//'/atm/icede.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/atm/icede.txt', form='formatted', status='old')
        READ(file_id, *) profiles(prof)%icede(:)
        CLOSE(file_id)
      ENDIF

      INQUIRE(file=TRIM(pathstr)//'/atm/ice_scheme.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/atm/ice_scheme.txt', form='formatted', status='old')
        READ(file_id, nml = ice_scheme_nml) 
        CLOSE(file_id)
        profiles(prof)%ice_scheme = ice_scheme
        profiles(prof)%icede_param = icede_param
      ENDIF

      OPEN(file_id, file=TRIM(pathstr)//'/atm/mmr_cldaer.txt', form='formatted', status='old')
      READ(file_id, nml=cldaer_units)
      CLOSE(file_id)
      profiles(prof)%mmr_cldaer = mmr_cldaer

      INQUIRE(file=TRIM(pathstr)//'/gas_units.txt', exist=exists)
      IF (exists) THEN
        OPEN(file_id, file=TRIM(pathstr)//'/gas_units.txt', form='formatted', status='old')
        READ(file_id, nml=units)
        CLOSE(file_id)
        profiles(prof)%gas_units = gas_units
      ENDIF

      ! Skin (surftype) required for CLW Deff param
      OPEN(file_id, file=TRIM(pathstr)//'/ground/skin.txt', form='formatted', status='old')
      READ(file_id, nml=skin)
      CLOSE(file_id)
      profiles(prof)%skin = k0
    ENDDO

    CALL rttov_put_lun(file_id)

    CALL rttov_convert_profile_units(opts, coefs, profiles, profiles_int)

  ! ------------------------------------
  ! Init opt params
  ! ------------------------------------
    ncld = SIZE(profiles(1)%cloud, DIM=1)
    nmom = MAX(coefs%coef_scatt%optp_wcl_opac%maxnmom, &
               coefs%coef_scatt%optp_wcl_deff%maxnmom, &
               coefs%coef_scatt%optp_icl_baum%maxnmom)
    opts%rt_ir%dom_nstreams = nmom
    nphangle = coefs%coef_scatt%optp_wcl_opac%nphangle
    CALL rttov_alloc_opt_param(err, opt_param, nchanprof, nlevels, nmom, &
                               nphangle, 1_jpim)

    ALLOCATE(absch(ncld), scach(ncld), bparh(ncld), &
             legcoef(0:nmom,ncld), phfn(nphangle_hires,ncld))

    ! Calculate clw and ice deff
    CALL rttov_profaux_cldaer(opts, profiles, profiles_int, aux)

  ! ------------------------------------
  ! Main loop
  ! ------------------------------------

    DO i = 1, nchanprof
      chan = chanprof(i)%chan
      prof = chanprof(i)%prof

      thermal = coefs%coef%ss_val_chn(chan) < 2
      solar   = coefs%coef%ss_val_chn(chan) > 0

      IF (solar) phchan = coefs%coef_scatt%optp_wcl_opac%chan_pha_index(chan)

      DO lay = 1, nlayers

        ! ------------------------------------
        ! Calculate layer abs, sca and bpr
        ! ------------------------------------
        absch(:) = 0.
        scach(:) = 0.
        bparh(:) = 0.
        legcoef(:,:) = 0.
        phfn(:,:) = 0._jprb

        DO iwc = 1, ncld
          IF (iwc <= nwcl_max) THEN ! Water cloud

            IF (profiles(prof)%clw_scheme == clw_scheme_opac) THEN  ! OPAC
              IF (profiles_int(prof)%cloud(iwc,lay) <= 0._jprb) CYCLE

              optp_data => coefs%coef_scatt%optp_wcl_opac%data(iwc)

              absch(iwc) = optp_data%abs(1,1,chan)
              scach(iwc) = optp_data%sca(1,1,chan)
              bparh(iwc) = optp_data%bpr(1,1,chan)
              legcoef(:,iwc) = optp_data%legcoef(:,1,1,chan)
              IF (solar) THEN
                phfn(:,iwc) = optp_data%pha(:,1,1,phchan)
              ENDIF

            ELSE ! Deff

              IF (iwc > 1) CYCLE
              IF (SUM(profiles_int(prof)%cloud(1:nwcl_max,lay)) <= 0._jprb .OR. aux%clw_dg(lay,prof) <= 0._jprb) CYCLE

              optp_data => coefs%coef_scatt%optp_wcl_deff%data(1)

              IF (aux%clw_dg(lay,prof) >= optp_data%deff(1) .AND. &
                  aux%clw_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                ! Find deff index below this channel
                DO k1 = 1, optp_data%ndeff - 1
                  IF (optp_data%deff(k1+1) > aux%clw_dg(lay,prof)) EXIT
                ENDDO
                k2 = k1 + 1
                dgfrac = (aux%clw_dg(lay,prof) - optp_data%deff(k1)) / &
                          (optp_data%deff(k2) - optp_data%deff(k1))
              ELSE
                ! Take first or last value if deff lies beyond data range
                IF (aux%clw_dg(lay,prof) < optp_data%deff(1)) THEN
                  k1 = 1
                  k2 = 1
                ELSE
                  k1 = optp_data%ndeff
                  k2 = optp_data%ndeff
                ENDIF
                dgfrac = 0._jprb
              ENDIF
              absch(iwc) = SUM(profiles_int(prof)%cloud(1:nwcl_max,lay)) * (optp_data%abs(1,k1,chan) + dgfrac * &
                          (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
              scach(iwc) = SUM(profiles_int(prof)%cloud(1:nwcl_max,lay)) * (optp_data%sca(1,k1,chan) + dgfrac * &
                          (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
              bparh(iwc) = optp_data%bpr(1,k1,chan) + dgfrac * &
                          (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
              legcoef(:,iwc) = optp_data%legcoef(:,1,k1,chan) + dgfrac * &
                          (optp_data%legcoef(:,1,k2,chan) - optp_data%legcoef(:,1,k1,chan))
              IF (solar) THEN
                phfn(:,iwc) = optp_data%pha(:,1,k1,phchan) + dgfrac * &
                          (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
              ENDIF

            ENDIF

          ELSE ! Ice cloud

            IF (profiles_int(prof)%cloud(iwc,lay) <= 0._jprb) CYCLE

            IF (profiles(prof)%ice_scheme == ice_scheme_baum) THEN ! Baum/SSEC

              optp_data => coefs%coef_scatt%optp_icl_baum%data(1)

              IF (aux%ice_dg(lay,prof) >= optp_data%deff(1) .AND. &
                  aux%ice_dg(lay,prof) < optp_data%deff(optp_data%ndeff)) THEN
                ! Find deff index below this channel
                DO k1 = 1, optp_data%ndeff - 1
                  IF (optp_data%deff(k1+1) > aux%ice_dg(lay,prof)) EXIT
                ENDDO
                k2 = k1 + 1
                dgfrac = (aux%ice_dg(lay,prof) - optp_data%deff(k1)) / &
                          (optp_data%deff(k2) - optp_data%deff(k1))
              ELSE
                ! Take first or last value if deff lies beyond data range
                IF (aux%ice_dg(lay,prof) < optp_data%deff(1)) THEN
                  k1 = 1
                  k2 = 1
                ELSE
                  k1 = optp_data%ndeff
                  k2 = optp_data%ndeff
                ENDIF
                dgfrac = 0._jprb
              ENDIF
              absch(iwc) = profiles_int(prof)%cloud(ncld,lay) * (optp_data%abs(1,k1,chan) + dgfrac * &
                          (optp_data%abs(1,k2,chan) - optp_data%abs(1,k1,chan)))
              scach(iwc) = profiles_int(prof)%cloud(ncld,lay) * (optp_data%sca(1,k1,chan) + dgfrac * &
                          (optp_data%sca(1,k2,chan) - optp_data%sca(1,k1,chan)))
              bparh(iwc) = optp_data%bpr(1,k1,chan) + dgfrac * &
                          (optp_data%bpr(1,k2,chan) - optp_data%bpr(1,k1,chan))
              legcoef(:,iwc) = optp_data%legcoef(:,1,k1,chan) + dgfrac * &
                          (optp_data%legcoef(:,1,k2,chan) - optp_data%legcoef(:,1,k1,chan))
              IF (solar) THEN
                phfn(:,iwc) = optp_data%pha(:,1,k1,phchan) + dgfrac * &
                          (optp_data%pha(:,1,k2,phchan) - optp_data%pha(:,1,k1,phchan))
              ENDIF

            ELSE ! Baran

              IF (profiles(prof)%ice_scheme == ice_scheme_baran2018) THEN
                CALL rttov_baran2018_calc_optpar(coefs%coef_scatt%optp_icl_baran2018, chan, &
                    profiles(prof)%t(lay), profiles_int(prof)%cloud(ncld,lay), absch(ncld), &
                    scach(ncld), bparh(ncld), asym)
              ELSE
                CALL rttov_baran2014_calc_optpar(coefs%coef_scatt%optp_icl_baran2014, chan, &
                    profiles(prof)%t(lay), profiles_int(prof)%cloud(ncld,lay), absch(ncld), &
                    scach(ncld), bparh(ncld), asym)
              ENDIF

              IF (solar) THEN
                ! For *all* solar channels use higher resolution angle grid
                ! (for mixed thermal+solar channels with solar scattering we need iphangle and we have this on
                !  the hi-res grid; saves calculating it for both hi-res and lo-res grids - this could be changed)
                thisnphangle = nphangle_hires
                thisphangle = phangle_hires
                thiscosphangle = coefs%coef_scatt%optp_icl_baran2018%phfn_int%cosphangle
              ELSE
                thisnphangle = nphangle_lores
                thisphangle(1:nphangle_lores) = phangle_lores
                thiscosphangle(1:nphangle_lores) = COS(phangle_lores * deg2rad)
              ENDIF

              ! Compute Baran phase fn
              CALL rttov_baran_calc_phase(asym, thisphangle(1:thisnphangle), phfn(1:thisnphangle,ncld))

              ! Compute Legendre coefficients
              CALL spline_interp(thisnphangle, thiscosphangle(thisnphangle:1:-1), &
                                  phfn(thisnphangle:1:-1,ncld), baran_ngauss, &
                                  coefs%coef_scatt%optp_icl_baran2018%q, baran_phfn_interp)
              CALL normalise(baran_ngauss, coefs%coef_scatt%optp_icl_baran2018%w, baran_phfn_interp)
              CALL calc_legendre_coef_gauss(coefs%coef_scatt%optp_icl_baran2018%q, &
                                            coefs%coef_scatt%optp_icl_baran2018%w, &
                                            baran_phfn_interp, nmom, nmom, nmomcalc, legcoef(:,ncld))
            ENDIF
          ENDIF
        ENDDO

        ! Combine liquid and ice properties
        IF (profiles(prof)%clw_scheme == clw_scheme_opac) THEN
          opt_param%abs(lay,i) = &
            SUM((/ (absch(n) * profiles_int(prof)%cloud(n,lay) * &
                    coefs%coef_scatt%optp_wcl_opac%data(n)%confac, n = 1, nwcl_max) /)) + absch(ncld)
          opt_param%sca(lay,i) = &
            SUM((/ (scach(n) * profiles_int(prof)%cloud(n,lay) * &
                    coefs%coef_scatt%optp_wcl_opac%data(n)%confac, n = 1, nwcl_max) /)) + scach(ncld)
          IF (thermal .AND. opt_param%sca(lay,i) > 0.) THEN
            opt_param%bpr(lay,i) = (SUM((/ (bparh(n) * scach(n) * coefs%coef_scatt%optp_wcl_opac%data(n)%confac * &
                                            profiles_int(prof)%cloud(n,lay), n = 1, nwcl_max) /)) + &
                                            bparh(ncld) * scach(ncld)) / opt_param%sca(lay,i)
          ELSE
            opt_param%bpr(lay,i) = 0.
          ENDIF
        ELSE
          opt_param%abs(lay,i) = absch(1) + absch(ncld)
          opt_param%sca(lay,i) = scach(1) + scach(ncld)
          IF (thermal .AND. opt_param%sca(lay,i) > 0.) THEN
            opt_param%bpr(lay,i) = (bparh(1) * scach(1) + bparh(ncld) * scach(ncld)) / opt_param%sca(lay,i)
          ELSE
            opt_param%bpr(lay,i) = 0.
          ENDIF
        ENDIF

        opt_param%nmom = nmom
        opt_param%legcoef(:,lay,i) = 0.
        opt_param%pha(:,lay,i) = 0.

        IF (opt_param%sca(lay,i) > 0.) THEN
          opt_param%legcoef(:,lay,i) = 0.
          IF (profiles(prof)%clw_scheme == clw_scheme_opac) THEN
            DO iwc = 1, nwcl_max
              opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) + &
                  legcoef(:,iwc) * scach(iwc) * profiles_int(prof)%cloud(iwc,lay) * &
                  coefs%coef_scatt%optp_wcl_opac%data(iwc)%confac
            ENDDO
          ELSE
            opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) + legcoef(:,1) * scach(1)
          ENDIF
          opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) + legcoef(:,ncld) * scach(ncld)
          opt_param%legcoef(:,lay,i) = opt_param%legcoef(:,lay,i) / opt_param%sca(lay,i)
        ENDIF

        IF (solar) THEN

          ! ------------------------------------
          ! Calculate layer phase fn
          ! ------------------------------------

          IF (opt_param%sca(lay,i) > 0.) THEN
            IF (profiles(prof)%clw_scheme == clw_scheme_opac) THEN
              DO iwc = 1, nwcl_max
                opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) + &
                    phfn(:,iwc) * scach(iwc) * profiles_int(prof)%cloud(iwc,lay) * &
                    coefs%coef_scatt%optp_wcl_opac%data(iwc)%confac
              ENDDO
            ELSE
              opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) + phfn(:,1) * scach(1)
            ENDIF
            opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) + phfn(:,ncld) * scach(ncld)
            opt_param%pha(:,lay,i) = opt_param%pha(:,lay,i) / opt_param%sca(lay,i)
          ENDIF

        ENDIF

      ENDDO ! layers
    ENDDO ! chanprof

  ! ------------------------------------
  ! Write out opt_param file
  ! ------------------------------------

    CALL rttov_get_lun(file_id)

    OPEN(file_id, file=TRIM(test_dir)//'/in/cld_opt_param.txt', form='formatted')
    WRITE(file_id, *) opt_param%nmom, nphangle
    WRITE(file_id, '(10E20.12)') opt_param%abs(:,:)
    WRITE(file_id, '(10E20.12)') opt_param%sca(:,:)
    WRITE(file_id, '(10E20.12)') opt_param%bpr(:,:)
    IF (opt_param%nmom > 0) THEN
      DO i = 1, nchanprof
        DO lay = 1, nlayers
          IF (opt_param%sca(lay,i) > 0._jprb) WRITE(file_id, '(10E20.12)') opt_param%legcoef(:,lay,i)
        ENDDO
      ENDDO
    ENDIF
    IF (nphangle > 0) THEN
      WRITE(file_id, '(10E20.12)') phangle_hires
      DO i = 1, nchanprof
        IF (coefs%coef%ss_val_chn(chanprof(i)%chan) > 0) THEN
          DO lay = 1, nlayers
            IF (opt_param%sca(lay,i) > 0._jprb) WRITE(file_id, '(10E20.12)') opt_param%pha(:,lay,i)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    CLOSE(file_id)

    CALL rttov_put_lun(file_id)

  ! ------------------------------------
  ! Clean up
  ! ------------------------------------

    DEALLOCATE(chanprof, absch, scach, bparh, legcoef, phfn)

    CALL rttov_alloc_opt_param(err, opt_param, nchanprof, nlevels, nmom, &
                               nphangle, 0_jpim)

    CALL rttov_alloc_prof(err, nprof, profiles, nlevels, opts, 0_jpim, coefs)
    DEALLOCATE(profiles)
    CALL rttov_alloc_prof(err, nprof, profiles_int, nlevels, opts, 0_jpim, coefs)
    DEALLOCATE(profiles_int)
    CALL rttov_alloc_aux_prof(err, nprof, nlevels, aux, opts, coefs%coef, 0_jpim)

    CALL rttov_dealloc_coefs(err, coefs)

  ENDIF ! sccld file

CONTAINS

  ! This subroutine calculates all parameters regardless of options (bpr, legcoefs, phfn)
  ! as opposed to calculating only what is required for the simulation being carried out.
  SUBROUTINE aer_interp_relhum(chan, dom_nstr, iaer, solar, optp, relhum, &
                               absch, scach, bparh, legcoef, phfn)

    USE parkind1, ONLY : jprb, jpim, jplm

    USE rttov_types, ONLY : rttov_optp

    IMPLICIT NONE

    INTEGER(jpim),    INTENT(IN)    :: chan, dom_nstr, iaer
    LOGICAL(jplm),    INTENT(IN)    :: solar
    TYPE(rttov_optp), INTENT(IN)    :: optp
    REAL(jprb),       INTENT(IN)    :: relhum
    REAL(jprb),       INTENT(INOUT) :: absch, scach, bparh
    REAL(jprb),       INTENT(INOUT) :: legcoef(0:)
    REAL(jprb),       INTENT(INOUT) :: phfn(optp%nphangle)

    INTEGER(jpim) :: k, phchan, nmom
    REAL(jprb)    :: frach
    REAL(jprb)    :: afac, sfac, gfac
    REAL(jprb)    :: pfac1(0:dom_nstr)
    REAL(jprb)    :: pfac2(optp%nphangle)
    ! ------------------------------------------------------------------------

    k = optp%data(iaer)%nrelhum
    IF (k /= 1 .AND. aux%relhum(lay,prof) <= optp%data(iaer)%relhum(k)) THEN
      ! Interpolate scattering parameters to actual value of relative humidity
      DO k = 1, optp%data(iaer)%nrelhum - 1
        IF (relhum >= optp%data(iaer)%relhum(k) .AND. &
            relhum <= optp%data(iaer)%relhum(k+1)) THEN

          frach = (relhum - optp%data(iaer)%relhum(k)) / &
                  (optp%data(iaer)%relhum(k+1) - optp%data(iaer)%relhum(k))
          afac  = (optp%data(iaer)%abs(k+1,1,chan) - optp%data(iaer)%abs(k,1,chan))
          absch = optp%data(iaer)%abs(k,1,chan) + afac * frach
          sfac  = (optp%data(iaer)%sca(k+1,1,chan) - optp%data(iaer)%sca(k,1,chan))
          scach = optp%data(iaer)%sca(k,1,chan) + sfac * frach

          gfac  = (optp%data(iaer)%bpr(k+1,1,chan) - optp%data(iaer)%bpr(k,1,chan))
          bparh = optp%data(iaer)%bpr(k,1,chan) + gfac * frach

          nmom = MIN(MAX(optp%data(iaer)%nmom(k,chan), optp%data(iaer)%nmom(k+1,chan)), dom_nstr)
          pfac1(0:nmom) = (optp%data(iaer)%legcoef(1:nmom+1,k+1,1,chan) - &
                           optp%data(iaer)%legcoef(1:nmom+1,k,1,chan))
          legcoef(0:nmom) = optp%data(iaer)%legcoef(1:nmom+1,k,1,chan) + pfac1(0:nmom) * frach

          IF (solar) THEN
            phchan = optp%chan_pha_index(chan)
            pfac2(:) = (optp%data(iaer)%pha(:,k+1,1,phchan) - &
                        optp%data(iaer)%pha(:,k,1,phchan))
            phfn(:) = optp%data(iaer)%pha(:,k,1,phchan) + pfac2(:) * frach
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! Particle doesn't change with rel. hum. (k=1) or rel. hum. exceeds max (k=max_rh_index)
      absch = optp%data(iaer)%abs(k,1,chan)
      scach = optp%data(iaer)%sca(k,1,chan)
      bparh = optp%data(iaer)%bpr(k,1,chan)
      nmom = MIN(optp%data(iaer)%nmom(k,chan), dom_nstr)
      legcoef(0:nmom) = optp%data(iaer)%legcoef(1:nmom+1,k,1,chan)
      IF (solar) THEN
        phchan = optp%chan_pha_index(chan)
        phfn(:) = optp%data(iaer)%pha(:,k,1,phchan)
      ENDIF
    ENDIF
  END SUBROUTINE aer_interp_relhum

END PROGRAM rttov_make_opt_param
