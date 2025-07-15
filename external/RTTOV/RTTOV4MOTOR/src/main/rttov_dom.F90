! Description:
!> @file
!!   RTE integration using Discrete Ordinates Method.
!
!> @brief
!!   RTE integration using Discrete Ordinates Method.
!!
!! @details
!!   The method here closely follows the implementation of DOM in
!!   DISORT. Delta-scaling is applied to the phase function and the
!!   TMS correction is also applied, although here the phase function
!!   is interpolated at the single-scattering angle rather than evaluating
!!   the full Legendre expansion for the TMS correction (as in DISORT).
!!   The same accuracy criterion is applied to terminate the azimuthal
!!   loop (for the solar source term): the loop exits when two radiance
!!   increments have been less than dom_accuracy * current_radiance.
!!
!!   This implementation treats either the thermal emission source term
!!   or the solar source term (according to the dosolar flag). For mixed
!!   thermal+solar channels this subroutine is called twice. This is
!!   done for reasons of efficiency for non-mixed channels.
!!
!!   The surface is assumed to be Lambertian. For the emissive source
!!   term the albedo is (1-emissivity) which is contained in diffuse_refl.
!!   For the solar source term the BRDF is multiplied by pi to obtain an
!!   albedo and this value is capped at 1. to avoid unphysical values.
!!
!!   This implementation assumes the atmosphere is strictly plane-parallel
!!   which is enforced by RTTOV elsewhere.
!!
!!   Currently where the solar zenith angle is close to a quadrature
!!   angle the solar zenith angle is nudged slightly with only very
!!   small impacts on radiances.
!!
!!   The results from intermediate calculations are stored only if the
!!   optional dom_state argument is present: this is only required by the
!!   TL/AD/K models to avoid repeating various direct model calculations.
!!
!!   The code includes a special case for MFASIS training which is switched in
!!   at compile time by supplying the _RTTOV_MFASIS_TRAINING macro. It was done
!!   as a compiler switch to avoid detrimental impact on performance for normal
!!   simulations observed with ifort. This alternative code provides an
!!   efficient way to simulate multiple "geometries" (i.e. sat/sol angles) for
!!   a single profile, which benefits the MFASIS training greatly. It can only
!!   be used with solar simulations. All input profiles must be identical save
!!   for the satellite/solar zenith and azimuth angles. This includes surface
!!   reflectances (they can differ per channel, but for each channel the
!!   reflectances must be identical for all profiles). The same channels must be
!!   simulated for every profile. In this case the DOM code simulates the
!!   multiple geometries in "one go" which is much more efficient. It is
!!   important to note that in order to get strictly identical results to those
!!   from a standard call, you should only vary the relative azimuth among the
!!   input profiles, and/or ensure that the sum of the satellite and solar
!!   zenith angle secants is constant. This is because the zenith angles appear
!!   in the optical depth regression and so otherwise the profiles are not
!!   strictly identical. In this case the dom_accuracy option has no effect.
!!   This code is intended for developers only, primarily for efficient
!!   training of the MFASIS LUTs. It strictly works for the direct model only.
!!
!! @param[out]    err                    status on exit
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     chanprof               specifies channels and profiles to simulate
!! @param[in]     dosolar                flag to indicate if solar source term (true) or emission source term (false)
!!                                       should be calculated
!! @param[in]     chanflag               flags to indicate which channels should be allocated (either channels with
!!                                       significant solar or thermal contributions according to dosolar)
!! @param[in]     nstr                   number of DOM streams; user selected, but RTTOV ensures it is even
!! @param[in]     maxnaz                 max number of azimuthal loops: must be 0 for thermal and nstr-1 for solar
!! @param[in]     profiles               input atmospheric profiles and surface variables
!! @param[in]     profiles_dom           layer optical depth profiles set up for DOM calculations
!! @param[in]     maxnlayers             maximum number of layers in the DOM profiles
!! @param[in]     auxrad                 Planck radiances for emission source term
!! @param[in]     trans_scatt_ir         single-scattering albedos
!! @param[in]     trans_scatt_ir_dyn     phase function Legendre decompositions
!! @param[in]     emissivity             surface emissivities
!! @param[in]     reflectance            direct solar surface BRDFs
!! @param[in]     diffuse_refl           surface albedos
!! @param[in]     solar_spectrum         top of atmosphere solar irradiance
!! @param[in]     raytracing             raytracing structure
!! @param[in]     ircld                  cloud column data
!! @param[in]     rad                    calculated radiances
!! @param         dom_state              pointer to array of rttov_dom_state structures (may be nullified in which case
!!                                       no state is stored)
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
SUBROUTINE rttov_dom(             &
              err,                &
              opts,               &
              chanprof,           &
              dosolar,            &
              chanflag,           &
              nstr,               &
              maxnaz,             &
              profiles,           &
              profiles_dom,       &
              maxnlayers,         &
              auxrad,             &
              trans_scatt_ir,     &
              trans_scatt_ir_dyn, &
              emissivity,         &
              reflectance,        &
              diffuse_refl,       &
              solar_spectrum,     &
              raytracing,         &
              ircld,              &
              rad,                &
              dom_state)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jprb, jpim, jplm

  USE rttov_types, ONLY :          &
      rttov_options,               &
      rttov_chanprof,              &
      rttov_profile,               &
      rttov_profile_dom,           &
      rttov_radiance_aux,          &
      rttov_transmission_scatt_ir, &
      rttov_emissivity,            &
      rttov_reflectance,           &
      rttov_raytracing,            &
      rttov_ircld,                 &
      rttov_radiance,              &
      rttov_dom_state

!INTF_OFF
  USE parkind1, ONLY : jpit

  USE rttov_const, ONLY : &
      pi,                 &
      z4pi_r,             &
      deg2rad

  USE rttov_scattering_mod, ONLY : &
      gauss_quad,                  &
      calc_leg_poly,               &
      asymtx

  USE rttov_lapack_mod, ONLY : &
      dgetrf, dgetrs, dgbtrf, dgbtrs
!INTF_ON

  IMPLICIT NONE

  INTEGER(jpim),                     INTENT(OUT)   :: err
  TYPE(rttov_options),               INTENT(IN)    :: opts
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                     INTENT(IN)    :: dosolar
  LOGICAL(jplm),                     INTENT(IN)    :: chanflag(SIZE(chanprof))
  INTEGER(jpim),                     INTENT(IN)    :: nstr
  INTEGER(jpim),                     INTENT(IN)    :: maxnaz
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile_dom),           INTENT(IN)    :: profiles_dom(0:,:)
  INTEGER(jpim),                     INTENT(IN)    :: maxnlayers
  TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn
  TYPE(rttov_emissivity),            INTENT(IN)    :: emissivity(SIZE(chanprof))
  TYPE(rttov_reflectance),           INTENT(IN)    :: reflectance(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: diffuse_refl(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_radiance),              INTENT(INOUT) :: rad
  TYPE(rttov_dom_state),             POINTER       :: dom_state(:)
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=128) :: msg
  LOGICAL(jplm) :: dothermal, inc_surface !, do_order_n_pi

  INTEGER(jpim) :: n, nmom
#ifdef _RTTOV_MFASIS_TRAINING
  INTEGER(jpim) :: thismaxnaz(SIZE(profiles)), igeom, chan
#else
  INTEGER(jpim) :: thismaxnaz
#endif

  INTEGER(jpim) :: i, j, k, l, m, row, col
  INTEGER(jpim) :: icc, icci
  INTEGER(jpim) :: nlayers, nlevels, lay, lev, laymap
  INTEGER(jpim) :: nchanprof, prof, lastprof
  INTEGER(jpim) :: lo, klu1
  INTEGER(jpim) :: convcount

  INTEGER       :: kl, ku, klku, info, nstri, nstrlay ! LAPACK arguments, no KIND

#ifdef _RTTOV_MFASIS_TRAINING
  INTEGER       :: ngeom                              ! LAPACK argument, no KIND
  REAL(jprb)    :: musat(SIZE(profiles)), musatr(SIZE(profiles))
  REAL(jprb)    :: musun(SIZE(profiles)), musunr(SIZE(profiles)), relazi(SIZE(profiles))
  REAL(jprb)    :: delta_tms(SIZE(profiles)), radinc(SIZE(profiles)), thisrad(SIZE(profiles))
#else
  REAL(jprb)    :: musat, musatr, musun, musunr, relazi
  REAL(jprb)    :: delta_tms, radinc, thisrad
#endif
  REAL(jprb)    :: muscata, taufac, albedo, reflfac, sunsrcfac, radsurfup_sum, radsurfup

  REAL(jprb)              :: muq(nstr/2), muw(nstr/2), muqw(nstr/2), muqr(nstr/2)
  REAL(jprb)              :: legpolymuq(0:nstr-1,nstr,0:maxnaz)
#ifdef _RTTOV_MFASIS_TRAINING
  REAL(jprb)              :: legpolymusat(0:nstr-1,0:0,0:maxnaz,SIZE(profiles))
  REAL(jprb)              :: legpolymusun(0:nstr-1,0:0,0:nstr-1,SIZE(profiles))  ! dims are (0:nmom,0,0:maxnaz,ngeom)
  REAL(jprb)              :: legpolymuscata(0:nstr-1,0:0,0:0,SIZE(profiles))     !(0:nstr-1,0:0,0:0,ngeom)
#else
  REAL(jprb)              :: legpolymusat(0:nstr-1,0:0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymusun(:,:,:) !(0:nstr-1,0:0,0:nstr-1)  ! dims are (0:nmom,0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymuscata(:,:,:) !(0:nstr-1,0:0,0:0)
#endif
  REAL(jprb)              :: legcoef(0:nstr-1,profiles(1)%nlayers,0:1)

  REAL(jprb)              :: f_lay, ftau(maxnlayers), fssa(maxnlayers)
  REAL(jprb)              :: phase_fn(nstr/2,nstr)
  REAL(jprb), ALLOCATABLE :: b0(:), b1(:) !(maxnlayers)
  REAL(jprb)              :: opdep(maxnlayers+1)
#ifdef _RTTOV_MFASIS_TRAINING
  REAL(jprb)              :: tausat(maxnlayers+1,SIZE(profiles))
  REAL(jprb)              :: tausun(maxnlayers+1,SIZE(profiles))
#else
  REAL(jprb)              :: tausat(maxnlayers+1)
  REAL(jprb), ALLOCATABLE :: tausun(:) !(maxnlayers+1)
#endif
  REAL(jprb)              :: tlayer(nstr/2,maxnlayers)
  REAL(jprb)              :: surfterm(nstr)
#ifdef _RTTOV_MFASIS_TRAINING
  REAL(jprb)              :: phasesat(nstr,SIZE(profiles))
#else
  REAL(jprb)              :: phasesat(nstr)
#endif

  REAL(jprb)              :: bp(nstr/2,nstr/2), bm(nstr/2,nstr/2), bpm(nstr/2,nstr/2)
  REAL(jprb)              :: xp(nstr/2,nstr/2), xm(nstr/2,nstr/2)
  REAL(jprb)              :: sp(nstr/2,nstr/2), sm(nstr/2,nstr/2)
!   REAL(jprb)              :: so1, so2, sop(nstr/2), som(nstr/2)

  ! LAPACK arguments, no KIND
  INTEGER                 :: ipiv1(nstr*maxnlayers), ipiv2(nstr)
  DOUBLE PRECISION        :: hh(nstr,nstr), kkpack(9*nstr/2-2,nstr*maxnlayers)
!   DOUBLE PRECISION        :: bpmsave(nstr/2,nstr/2), gp(nstr/2)
#ifdef _RTTOV_MFASIS_TRAINING
  DOUBLE PRECISION        :: x(nstr*maxnlayers,SIZE(profiles))
  DOUBLE PRECISION        :: z(nstr,0:profiles(1)%nlayers,0:1,0:nstr-1,SIZE(profiles))
#else
  DOUBLE PRECISION        :: x(nstr*maxnlayers)
  DOUBLE PRECISION, ALLOCATABLE :: z(:,:,:,:) !(nstr,0:profiles(1)%nlayers,0:1,0:nstr-1)
#endif
  DOUBLE PRECISION, ALLOCATABLE :: y0(:,:) !(nstr,maxnlayers)
  DOUBLE PRECISION, ALLOCATABLE :: y1(:,:) !(nstr,maxnlayers)

  REAL(jprb)              :: eval(nstr/2,0:profiles(1)%nlayers,0:1,0:maxnaz)
#ifdef _RTTOV_MFASIS_TRAINING
  REAL(jprb)              :: fpmuhsup(nstr,profiles(1)%nlayers,0:1,0:maxnaz,SIZE(profiles))
  REAL(jprb)              :: fpmuhsdn(nstr,profiles(1)%nlayers,0:1,0:maxnaz,SIZE(profiles))
  REAL(jprb)              :: fpmupi(profiles(1)%nlayers,0:1,0:nstr-1,SIZE(profiles))
#else
  REAL(jprb)              :: fpmuhsup(nstr,profiles(1)%nlayers,0:1,0:maxnaz)
  REAL(jprb)              :: fpmuhsdn(nstr,profiles(1)%nlayers,0:1,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: fpmupi(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
#endif
  REAL(jprb), ALLOCATABLE :: fpmupi0(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi1(:) !(maxnlayers)

  REAL(jprb), TARGET  :: fp1(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:maxnaz)
  REAL(jprb), TARGET  :: fm1(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:maxnaz)
  REAL(jprb), POINTER :: fp2(:,:,:,:,:), fm2(:,:,:,:,:)
  REAL(jprb), POINTER :: ssa(:,:,:)

  ! This is only an internal logical array so use smallest available KIND
  LOGICAL(jpit) :: done_lay_str_azi(profiles(1)%nlayers,0:1,0:maxnaz)
  LOGICAL(jplm) :: store

  ! --------------------------------------------------------------------------------------

  TRY

!   do_order_n_pi = .FALSE.

  ! --------------------------------------------------------------------------
  ! Initialisation
  ! --------------------------------------------------------------------------
  nchanprof = SIZE(chanprof)
  dothermal = .NOT. dosolar

#ifdef _RTTOV_MFASIS_TRAINING
  IF (dothermal) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "MFASIS training is for solar only")
  ENDIF
  store = .FALSE.
  ngeom = SIZE(profiles)
#else
  store = ASSOCIATED(dom_state)
#endif

  ! nstr is guaranteed to be even so this is safe
  n = nstr / 2_jpim
  nstri = nstr

  ! Number of Legendre moments determined by nstr
  nmom = nstr - 1

  ! Index variables for packing kkpack
  kl = 3 * n - 1
  ku = 3 * n - 1
  klu1 = kl + ku + 1
  klku = 2 * kl + ku + 1

#ifndef _RTTOV_MFASIS_TRAINING
  ! Thermal-/solar-specific allocations to minimise memory usage
  IF (dosolar) THEN
    ALLOCATE(legpolymusun(0:nstr-1,0:0,0:maxnaz), &
             legpolymuscata(0:nstr-1,0:0,0:0), &
             tausun(maxnlayers+1), &
             fpmupi(profiles(1)%nlayers,0:1,0:maxnaz), &
             z(nstr,0:profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
  ELSE
    ALLOCATE(b0(maxnlayers), b1(maxnlayers), &
             fpmupi0(maxnlayers), fpmupi1(maxnlayers), &
             y0(nstr,maxnlayers), y1(nstr,maxnlayers), stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM allocation error')
#endif

  ! Exploit symmetry of eigenvectors using pointers
  fp2 => fm1
  fm2 => fp1

  IF (dosolar) THEN
    ssa => trans_scatt_ir%ssa_solar
  ELSE
    ssa => trans_scatt_ir%ssa_thermal
  ENDIF

  ! --------------------------------------------------------------------------
  ! Calculate Leg. polynomials for Gauss quadrature points
  ! --------------------------------------------------------------------------
  CALL gauss_quad(0._jprb, 1._jprb, muq, muw)
  muqr = 1._jprb / muq
  muqw = muq * muw

  ! Evaluate Legendre polynomials at the quadrature points
  CALL calc_leg_poly(nmom, maxnaz, muq, legpolymuq(:,1:n,:))

  ! Calculate values for -ve muq from those for +ve muq
  DO m = 0, maxnaz
    DO l = m, nstr - 1
      legpolymuq(l,n+1:nstr,m) = (-1) ** (l-m) * legpolymuq(l,1:n,m)
    ENDDO
  ENDDO

  ! --------------------------------------------------------------------------
  ! For clear layers the eigenvalues and eigenvectors are fixed
  ! --------------------------------------------------------------------------
  ! For clear layers the eigenvalue problem is trivial as the matrix is diagonal
  ! so obtain the eigenvalues (eval) and eigenvectors (fp1/fm1) directly
  ! Values are the same for all layers so use layer index zero for this
  DO m = 0, maxnaz
    eval(:,0,0,m) = muqr(:)
    fp1(:,:,0,0,m) = 0._jprb
    fm1(:,:,0,0,m) = 0._jprb
    DO j = 1, n
      fm1(j,j,0,0,m) = -muqr(j)
    ENDDO
    ! See notes below regarding fp2/fm2 (i.e. eigenvectors for -ve eigenvalues)

    ! Particular integral is zero for clear layers for solar channels
    IF (dosolar) THEN
#ifdef _RTTOV_MFASIS_TRAINING
      z(:,0,0,m,:) = 0._jprb
#else
      z(:,0,0,m) = 0._jprb
#endif
      IF (store) THEN
        DO i = 1, nchanprof
          IF (.NOT. chanflag(i)) CYCLE
          dom_state(i)%z(:,0,0,m) = 0._jprb
        ENDDO
      ENDIF
    ENDIF
  ENDDO

  lastprof = -1

  ! --------------------------------------------------------------------------
  ! Channel loop
  ! --------------------------------------------------------------------------
  DO i = 1, nchanprof
    IF (.NOT. chanflag(i)) CYCLE

    prof = chanprof(i)%prof
#ifdef _RTTOV_MFASIS_TRAINING
      IF (prof > 1) EXIT
      chan = chanprof(i)%chan
#endif

    IF (store) THEN
      dom_state(i)%eval = 0._jprb
      dom_state(i)%x = 0._jprb
    ENDIF

    IF (prof /= lastprof) THEN
#ifdef _RTTOV_MFASIS_TRAINING
      musatr = raytracing%pathsat(1,:)
      musat = 1._jprb / musatr

      musunr = raytracing%pathsun(1,:)
      musun = 1._jprb / musunr
      ! Ensure musun doesn't coincide with a quadrature angle
      DO igeom = 1, ngeom
        DO
          IF (.NOT. ANY(ABS(musun(igeom) - muq) < 1.E-10_jprb)) EXIT
          musun(igeom) = musun(igeom) + 2.E-10_jprb
          musunr(igeom) = 1._jprb / musun(igeom)
        ENDDO
      ENDDO

      relazi = (180._jprb - (profiles(:)%azangle - profiles(:)%sunazangle)) * deg2rad

      DO igeom = 1, ngeom
        ! Determine max number of azimuthal components
        IF (musat(igeom) < 1._jprb .AND. musun(igeom) < 1._jprb) THEN
          thismaxnaz(igeom) = maxnaz
        ELSE
          thismaxnaz(igeom) = 0
        ENDIF
      
        ! Evaluate Legendre polynomials at satellite and solar zenith angles and
        ! evaluate Legendre polynomials at scattering angle for TMS correction
        CALL calc_leg_poly(nmom, maxnaz, (/ musat(igeom) /), legpolymusat(:,:,:,igeom))
        CALL calc_leg_poly(nmom, maxnaz, (/ -musun(igeom) /), legpolymusun(:,:,:,igeom))
        muscata = -musat(igeom) * musun(igeom) + &
          SQRT((1._jprb - musun(igeom) ** 2) * (1._jprb - musat(igeom) ** 2)) * COS(relazi(igeom))
        CALL calc_leg_poly(nmom, 0_jpim, (/ muscata /), legpolymuscata(:,:,:,igeom))
      ENDDO
#else
      ! Angles/cosines
      musatr = raytracing%pathsat(1,prof)
      musat = 1._jprb / musatr

      IF (dosolar) THEN
        musunr = raytracing%pathsun(1,prof)
        musun = 1._jprb / musunr
        ! Ensure musun doesn't coincide with a quadrature angle
        DO
          IF (.NOT. ANY(ABS(musun - muq) < 1.E-10_jprb)) EXIT
          musun = musun + 2.E-10_jprb
          musunr = 1._jprb / musun
        ENDDO

        relazi = (180._jprb - (profiles(prof)%azangle - profiles(prof)%sunazangle)) * deg2rad
      ELSE
        musun = 1._jprb
      ENDIF

      ! Determine max number of azimuthal components
      IF (dosolar .AND. musat < 1._jprb .AND. musun < 1._jprb) THEN
        thismaxnaz = maxnaz
      ELSE
        thismaxnaz = 0
      ENDIF

      ! Evaluate Legendre polynomials at satellite and solar zenith angles and
      ! evaluate Legendre polynomials at scattering angle for TMS correction
      CALL calc_leg_poly(nmom, maxnaz, (/ musat /), legpolymusat)
      IF (dosolar) THEN
        CALL calc_leg_poly(nmom, maxnaz, (/ -musun /), legpolymusun)
        muscata = -musat * musun + SQRT((1._jprb - musun ** 2) * (1._jprb - musat ** 2)) * COS(relazi)
        CALL calc_leg_poly(nmom, 0_jpim, (/ muscata /), legpolymuscata)
      ENDIF
#endif
      lastprof = prof
    ENDIF

    ! The top level in profiles_dom is always the space boundary
    opdep(1) = 0._jprb
#ifdef _RTTOV_MFASIS_TRAINING
    tausat(1,:) = 1._jprb
    tausun(1,:) = 1._jprb
#else
    tausat(1) = 1._jprb
    IF (dosolar) tausun(1) = 1._jprb
#endif
    done_lay_str_azi = .FALSE.

    DO icc = 0, ircld%ncolumn(prof)
      thisrad   = 0._jprb
      delta_tms = 0._jprb
      convcount = 0_jpim

      ! Surface reflectance factor
      IF (profiles_dom(icc,i)%surface) THEN
        IF (dosolar) THEN
          ! This prevents unphysical albedos
          albedo = MIN(reflectance(i)%refl_out * pi, 1._jprb)
          reflfac = 2._jprb * albedo
          inc_surface = (reflfac > 0._jprb)
        ELSE
          reflfac = 2._jprb * diffuse_refl(i)
          inc_surface = .TRUE.
        ENDIF
      ELSE
        inc_surface = .FALSE.
      ENDIF

      nlayers = profiles_dom(icc,i)%nlayers
      nlevels = nlayers + 1
      nstrlay = nstr * nlayers

      IF (store .AND. dosolar) dom_state(i)%nazloops(icc) = -1  ! First loop is zero

      ! --------------------------------------------------------------------------
      ! Azimuthal component loop
      ! --------------------------------------------------------------------------
#ifdef _RTTOV_MFASIS_TRAINING
      DO m = 0, MAXVAL(thismaxnaz)
#else
      DO m = 0, thismaxnaz
#endif
        IF (store .AND. dosolar) dom_state(i)%nazloops(icc) = dom_state(i)%nazloops(icc) + 1

        ! Solar source term factor
        IF (dosolar) THEN
          IF (m == 0) THEN
            sunsrcfac = z4pi_r * solar_spectrum(i)
          ELSE
            sunsrcfac = 2._jprb * z4pi_r * solar_spectrum(i)
          ENDIF
        ENDIF

        kkpack = 0._jprb

        ! --------------------------------------------------------------------------
        ! Layer loop
        ! --------------------------------------------------------------------------
        DO lay = 1, nlayers
          lev = lay + 1_jpim

          ! laymap == 0 for clear layers
          laymap = profiles_dom(icc,i)%laymap(lay)
          IF (laymap > 0) THEN
            icci = ircld%icldarr(icc,laymap,prof)
          ELSE
            icci = 0
          ENDIF

          ! --------------------------------------------------------------------------
          ! Set delta-scaling parameters (same for all azimuthal components)
          ! --------------------------------------------------------------------------
          IF (m == 0) THEN
            IF (laymap > 0) THEN
              f_lay = trans_scatt_ir_dyn%phasefn(i)%legcoef(nstr,icci,laymap) / (2._jprb * nstr + 1._jprb)
              ftau(lay) = 1._jprb - ssa(icci,laymap,i) * f_lay
              fssa(lay) = (1._jprb - f_lay) * ssa(icci,laymap,i) / ftau(lay)
            ELSE
              f_lay = 0._jprb
              ftau(lay) = 1._jprb
              fssa(lay) = 0._jprb
            ENDIF
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate optical depths and transmittances
          ! --------------------------------------------------------------------------
          ! Calculate accumulated delta-scaled optical depth and transmittances on
          ! satellite-surface and sun-surface paths.
          IF (m == 0) THEN  ! This calculation applies to all azi components
            opdep(lev) = opdep(lev-1) + ftau(lay) * profiles_dom(icc,i)%layerod(lay)
#ifdef _RTTOV_MFASIS_TRAINING
            tausat(lev,:) = EXP(-opdep(lev) * musatr(:))
            tausun(lev,:) = EXP(-opdep(lev) * musunr(:))
#else
            tausat(lev) = EXP(-opdep(lev) * musatr)
            IF (dosolar) tausun(lev) = EXP(-opdep(lev) * musunr)
#endif
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate emissive source term and interpolated clear-sky emissive source term
          ! --------------------------------------------------------------------------
          IF (dothermal) THEN
            IF (lay == nlayers .AND. profiles_dom(icc,i)%surface) THEN
              b1(lay) = (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) / (opdep(lev-1) - opdep(lev))
            ELSE
              b1(lay) = (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) / (opdep(lev-1) - opdep(lev))
            ENDIF
            b0(lay) = (1._jprb - fssa(lay)) * auxrad%air(lev-1,i) - b1(lay) * opdep(lev-1)

            y1(:,lay) = b1(lay)

            ! For clear layers we have the solution and the interpolated thermal source term is trivial
            IF (.NOT. fssa(lay) > 0._jprb) THEN
              y0(1:n,lay)      = b0(lay) + muq(1:n) * y1(1:n,lay)
              y0(n+1:nstr,lay) = b0(lay) - muq(1:n) * y1(n+1:nstr,lay)
              fpmupi0(lay) = b0(lay)
              fpmupi1(lay) = b1(lay)
            ENDIF
          ENDIF

          ! --------------------------------------------------------------------------
          ! Solve eigenvalue problem and obtain homogenous solution for this layer
          ! --------------------------------------------------------------------------
          ! Note that for clear layers we already wrote the solutions down above

          IF (fssa(lay) > 0._jprb) THEN

            ! --------------------------------------------------------------------------
            ! For scattering layers, only do calculation if we don't already have the results
            ! --------------------------------------------------------------------------
            IF (dothermal .OR. .NOT. done_lay_str_azi(laymap,icci,m)) THEN

              ! --------------------------------------------------------------------------
              ! Calculate the phase function values for each pair of quadrature angles and musun
              ! --------------------------------------------------------------------------
              ! Delta-scale Legendre coefficients for this layer
              ! NB delta-scaled SSA is incorporated into Leg. coefs
              IF (m == 0) THEN  ! This calculation applies to all azi components
                DO l = 0, nmom
                  legcoef(l,laymap,icci) = fssa(lay) * (trans_scatt_ir_dyn%phasefn(i)%legcoef(l,icci,laymap) - &
                                           (2._jprb * l + 1._jprb) * f_lay) / &
                                           (1._jprb - f_lay)
                ENDDO
              ENDIF

              ! Due to symmetry we can calculate just n*nstr values (not nstr*nstr):
              !    phase_fn(i,j) = phase_fn(n+i,n+j)    where n+x => -ve mu
              !    phase_fn(n+i,j) = phase_fn(i,n+j)
              DO j = 1, n
                DO k = j, n
                  phase_fn(k,j)   = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,j,m))
                  phase_fn(k,n+j) = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,n+j,m))

                  ! Exploit symmetry to avoid repeated calculations
                  IF (k > j) THEN
                    phase_fn(j,k)   = phase_fn(k,j)
                    phase_fn(j,n+k) = phase_fn(k,n+j)
                  ENDIF
                ENDDO
              ENDDO

#ifdef _RTTOV_MFASIS_TRAINING
              DO igeom = 1, ngeom
                ! Include quadrature weights in phase fns evaluated for musat and each quadrature angle
                DO j = 1, n
                  phasesat(j,igeom) = muw(j) * SUM(legcoef(m:nmom,laymap,icci) * &
                                                   legpolymusat(m:nmom,0,m,igeom) * legpolymuq(m:nmom,j,m))
                ENDDO
                DO j = n+1, nstr
                  phasesat(j,igeom) = muw(j-n) * SUM(legcoef(m:nmom,laymap,icci) * &
                                                     legpolymusat(m:nmom,0,m,igeom) * legpolymuq(m:nmom,j,m))
                ENDDO
              ENDDO
#else
              ! Include quadrature weights in phase fns evaluated for musat and each quadrature angle
              DO j = 1, n
                phasesat(j) = muw(j) * SUM(legcoef(m:nmom,laymap,icci) * &
                                           legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO
              DO j = n+1, nstr
                phasesat(j) = muw(j-n) * SUM(legcoef(m:nmom,laymap,icci) * &
                                             legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO
#endif
            ENDIF ! dothermal or not done_lay_str_azi

            IF (.NOT. done_lay_str_azi(laymap,icci,m)) THEN

              ! Homogenous solution only done once per layer for scattering layers (applies to all cloud columns)

              ! --------------------------------------------------------------------------
              ! Solve eigenvalue problem and obtain homogenous solution for this layer
              ! --------------------------------------------------------------------------
              DO j = 1, n
                bp(:,j) = 0.5_jprb * muw(j) * phase_fn(:,j) * muqr(:)
                bm(:,j) = 0.5_jprb * muw(j) * phase_fn(:,n+j) * muqr(:)
                bp(j,j) = bp(j,j) - muqr(j)
              ENDDO
              sm = bp - bm
              sp = bp + bm
              bpm = MATMUL(sm, sp)
!               bpmsave = bpm  ! Used for solar order-n PI below

              ! DISORT's ASYMTX routine is faster than LAPACK's DGEEV because it only returns real eigenvalues
              CALL ASYMTX(err, bpm, xp, eval(:,laymap,icci,m), n, n, n)
              IF (err /= 0) THEN
                WRITE(msg,'(a,i5)') 'DOM error: solving eigenvalue problem, ASYMTX error = ', err
                THROWM(err.NE.0,msg)
              ENDIF

              ! ** evals should always be +ve **
              IF (ANY(eval(:,laymap,icci,m) <= 0.)) THEN
                ! Real parts must be positive
                err = errorstatus_fatal
                THROWM(err.NE.0,'DOM error: some eigenvalue(s) <= 0')
              ENDIF

              IF (store) dom_state(i)%xp(:,:,laymap,icci,m) = xp

              ! Calculate the required eigenvalues (+/- each one gives nstr evals)
              eval(:,laymap,icci,m) = SQRT(ABS(eval(:,laymap,icci,m)))

              DO j = 1, n
                xm(:,j) = MATMUL(sp, xp(:,j)) / eval(j,laymap,icci,m)
              ENDDO

              ! Rearrange to obtain F+ and F- (the eigenvectors)
              ! For +ve eigenvalues:
              fp1(:,:,laymap,icci,m) = 0.5_jprb * (xp + xm)
              fm1(:,:,laymap,icci,m) = 0.5_jprb * (xp - xm)

              ! For -ve eigenvalues (eval -> -eval) we could write:
              !   xm = -xm
              !   fp2(:,:,lay) = 0.5_jprb * (xp + xm)   which equals   fm1(:,:,lay)
              !   fm2(:,:,lay) = 0.5_jprb * (xp - xm)   which equals   fp1(:,:,lay)
              ! Instead of this redundancy we have already pointed fp2/fm2 to the relevant
              ! arrays above so we have:
              !   fp2 => fm1
              !   fm2 => fp1


              ! The complete set of eigenvectors can be written as (nstr x nstr) matrix EVECS = /fp1 fp2\
              !                                                                                 \fm1 fm2/
              ! where the columns form the eigenvectors and the rows correspond to the quadrature directions muq.
              !
              ! In other words:
              ! fp1/fp2(i,j) relates to upward radiances (ith component of jth eigenvector, relates to muq(i))
              ! fm1/fm2(i,j) relate to downward radiances (-muq(i))
              ! fp1/fm1 are eigenvectors corresponding to j = 1, ..., n (i.e. the positive eigenvalues)
              ! fp2/fm2 are eigenvectors correspond to j = -1, ..., -n (i.e. the negative eigenvalues)
              ! EVECS(i,:) represents the scattering into direction +/-mu_i from all j=+/1,...,+/-n directions (but
              !   note that the ordering of the j directions depends on the arbitrary order of evals/evecs)


              ! Solar source term particular integral only done once per layer for scattering layers (applies to all cloud columns)
              IF (dosolar) THEN
                ! --------------------------------------------------------------------------
                ! Solve for particular integral for solar source term
                ! --------------------------------------------------------------------------

!                 IF (.NOT. do_order_n_pi) THEN

                  ! Solve PI problem using nstr x nstr system (most straightforward approach)

#ifdef _RTTOV_MFASIS_TRAINING
                  DO igeom = 1, ngeom
#endif
                  ! Construct and solve the linear system for scattering layers
                  DO j = 1, n
                    hh(1:n,j)      = -0.5_jprb * muw(j) * phase_fn(:,j)
                    hh(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry

                    ! gfortran gives wrong answers with this code:
!                     hh(1:n,n+j)      = hh(n+1:nstr,j)   ! Exploit symmetry: matrix structure before   /A B\
!                     hh(n+1:nstr,n+j) = hh(1:n,j)        ! addition of diagonal elements is:           \B A/

                    ! ...so do it in a loop:
                    DO k = 1, n
                      hh(k,n+j)   = hh(n+k,j)
                      hh(n+k,n+j) = hh(k,j)
                    ENDDO

#ifdef _RTTOV_MFASIS_TRAINING
                    hh(j,j)     = hh(j,j) + 1._jprb + muq(j) * musunr(igeom)
                    hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb - muq(j) * musunr(igeom)

                    z(j,laymap,icci,m,igeom)   = sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                                 legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m,igeom))
                    z(n+j,laymap,icci,m,igeom) = sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                                 legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m,igeom))
#else
                    hh(j,j)     = hh(j,j) + 1._jprb + muq(j) * musunr
                    hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb - muq(j) * musunr

                    z(j,laymap,icci,m)   = sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                           legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m))
                    z(n+j,laymap,icci,m) = sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                           legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m))
#endif
                  ENDDO

                  CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
                  IF (info /= 0) THEN
                    WRITE(msg,'(a,3i5)') 'DOM error: solving PI solar, DGETRF error, info = ', info, i, lay
                    err = errorstatus_fatal
                    THROWM(err.NE.0,msg)
                  ENDIF

#ifdef _RTTOV_MFASIS_TRAINING
                  CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, z(:,laymap,icci,m,igeom), nstri, info)
#else
                  CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, z(:,laymap,icci,m), nstri, info)
#endif
                  IF (info /= 0) THEN
                    WRITE(msg,'(a,3i5)') 'DOM error: solving PI solar, DGETRS error, info = ', info, i, lay
                    err = errorstatus_fatal
                    THROWM(err.NE.0,msg)
                  ENDIF
#ifdef _RTTOV_MFASIS_TRAINING
                  ENDDO
#endif

!                 ELSE
!                   ! Solve PI problem using n x n system (instead of nstr x nstr)
!                   ! This offers a small saving in time (e.g. 1% for nstr=16, should get more benefit as nstr increases)
!                   ! By having to store results from homogenous solution it increases memory usage very slightly
!                   ! e.g. by a few times (n x n) (could get rid of hh for solar calcs when using this method)
! 
!                   DO j = 1, n
!                     so1 = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m))
!                     so2 = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m))
!                     sop(j) = sunsrcfac * (so1 + so2) * muqr(j) * musun
!                     som(j) = sunsrcfac * (so1 - so2) * muqr(j) * musun
!                   ENDDO
! 
!                   gp = som + musun * MATMUL(sm, sop)
! 
!                   bpmsave = -musun**2 * bpmsave  ! bpmsave overwritten as it's not required anywhere else
!                   DO j = 1, n
!                     bpmsave(j,j) = 1._jprb + bpmsave(j,j)
!                   ENDDO
! 
!                   CALL DGESV(n, 1, bpmsave, n, ipiv2, gp, n, info)
!                   IF (info /= 0) THEN
!                     PRINT *,'DGESV (solving PI) error, info = ', info, i, lay
!                     err = errorstatus_fatal
!                     THROWM(err.NE.0,'DOM error: solving for particular integral - solar')
!                   ENDIF
! 
!                   ! Re-use som array for calculation of "gm":
!                   som = musun * MATMUL(sp, gp) + sop
!                   z(1:n,laymap,icci,m)      = 0.5_jprb * (gp + som)
!                   z(n+1:nstr,laymap,icci,m) = 0.5_jprb * (gp - som)
! 
!                 ENDIF

#ifndef _RTTOV_MFASIS_TRAINING
                IF (store) dom_state(i)%z(:,laymap,icci,m) = z(:,laymap,icci,m)
#endif
              ENDIF ! dosolar

            ENDIF ! not done_lay_str_azi


            ! Thermal source term done for every layer and *for every cloud column*
            IF (dothermal) THEN
              ! --------------------------------------------------------------------------
              ! Solve for particular integral for thermal source term
              ! --------------------------------------------------------------------------

              ! Construct and solve the linear system for scattering layers
              DO j = 1, n
                hh(1:n,j)      = -0.5_jprb * muw(j) * phase_fn(:,j)
                hh(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry

                ! gfortran gives wrong answers with this code:
!                 hh(1:n,n+j)      = hh(n+1:nstr,j)   ! Exploit symmetry: matrix structure before   /A B\
!                 hh(n+1:nstr,n+j) = hh(1:n,j)        ! addition of diagonal elements is:           \B A/

                ! ...so do it in a loop:
                DO k = 1, n
                  hh(k,n+j)   = hh(n+k,j)
                  hh(n+k,n+j) = hh(k,j)
                ENDDO

                hh(j,j)     = hh(j,j) + 1._jprb
                hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb
              ENDDO

              ! Use DGETRF for LU factorisation (of hh) and then two calls to DGETRS to solve for y1 and y0
              CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving PI thermal, DGETRF error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              ! Solve for y1 (coefficient of tau)
              CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, y1(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving PI thermal y1, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              y0(1:n,lay)      = b0(lay) + muq(1:n) * y1(1:n,lay)
              y0(n+1:nstr,lay) = b0(lay) - muq(1:n) * y1(n+1:nstr,lay)

              ! Solve for y0 (constant coefficient)
              CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, y0(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving PI thermal y0, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

            ENDIF ! dothermal

            ! --------------------------------------------------------------------------
            ! Interpolate HS and PI for final radiance calculation
            ! --------------------------------------------------------------------------

            IF (.NOT. done_lay_str_azi(laymap,icci,m)) THEN
#ifdef _RTTOV_MFASIS_TRAINING
              DO igeom = 1, ngeom
                ! Homogenous part - diffuse radiation, done once per layer (applies to all cloud columns)
                DO j = 1, n
                  fpmuhsup(j,laymap,icci,m,igeom) = 0.5_jprb * (SUM(fp1(:,j,laymap,icci,m) * phasesat(1:n,igeom)) + &
                                                                SUM(fm1(:,j,laymap,icci,m) * phasesat(n+1:nstr,igeom)))

                  fpmuhsdn(j,laymap,icci,m,igeom) = 0.5_jprb * (SUM(fp2(:,j,laymap,icci,m) * phasesat(1:n,igeom)) + &
                                                                SUM(fm2(:,j,laymap,icci,m) * phasesat(n+1:nstr,igeom)))
                ENDDO

                ! Particular integral - solar source term, done once per layer (applies to all cloud columns)
                IF (dosolar) THEN
                  fpmupi(laymap,icci,m,igeom) = 0.5_jprb * SUM(z(:,laymap,icci,m,igeom) * phasesat(:,igeom)) + &
                                                sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                                legpolymusat(m:nmom,0,m,igeom) * legpolymusun(m:nmom,0,m,igeom))
                ENDIF
              ENDDO
#else
              ! Homogenous part - diffuse radiation, done once per layer (applies to all cloud columns)
              DO j = 1, n
                fpmuhsup(j,laymap,icci,m) = 0.5_jprb * (SUM(fp1(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm1(:,j,laymap,icci,m) * phasesat(n+1:nstr)))

                fpmuhsdn(j,laymap,icci,m) = 0.5_jprb * (SUM(fp2(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm2(:,j,laymap,icci,m) * phasesat(n+1:nstr)))
              ENDDO

              ! Particular integral - solar source term, done once per layer (applies to all cloud columns)
              IF (dosolar) THEN
                fpmupi(laymap,icci,m) = 0.5_jprb * SUM(z(:,laymap,icci,m) * phasesat(:)) + &
                                        sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                        legpolymusat(m:nmom,0,m) * legpolymusun(m:nmom,0,m))
              ENDIF
#endif
            ENDIF ! not done_lay_str_azi

#ifndef _RTTOV_MFASIS_TRAINING
            ! Particular integral - thermal source term, done once per layer per cloud column
            IF (dothermal) THEN
              fpmupi0(lay) = 0.5_jprb * SUM(y0(:,lay) * phasesat(:)) + b0(lay)
              fpmupi1(lay) = 0.5_jprb * SUM(y1(:,lay) * phasesat(:)) + b1(lay)
            ENDIF
#endif
            done_lay_str_azi(laymap,icci,m) = .TRUE.

          ENDIF ! fssa > 0 (i.e. scattering layer)

          IF (dothermal .AND. store) THEN
            dom_state(i)%y0(:,lay,icc) = y0(:,lay)
            dom_state(i)%y1(:,lay,icc) = y1(:,lay)
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate layer "transmittances" with opdeps scaled by eigenvalues
          ! --------------------------------------------------------------------------
          tlayer(:,lay) = EXP(-ftau(lay) * profiles_dom(icc,i)%layerod(lay) * eval(:,laymap,icci,m))

          ! --------------------------------------------------------------------------
          ! Construct part of linear system for this layer
          ! --------------------------------------------------------------------------

          ! The matrix is banded and we use LAPACK's optimised solver.
          ! Construct the banded matrix kkpack directly to avoid allocating the full matrix kk:
          !   kkpack(kl+ku+1+i-j,j) = kk(i,j) for max(1,j-ku)<=i<=min(nstr*nlayers,j+kl)
          ! where kl and ku are defined near the top of the subroutine.

          IF (lay == 1) THEN
            row = 0
            col = 0
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j) ! kk(1:n,j)
              kkpack(lo:lo+n-1,j)   = fm1(:,j,laymap,icci,m)
              lo = lo - n                       ! kk(1:n,n+j)
              kkpack(lo:lo+n-1,n+j) = fm2(:,j,laymap,icci,m) * tlayer(j,lay)
            ENDDO
#ifdef _RTTOV_MFASIS_TRAINING
            DO igeom = 1, ngeom
              x(1:n,igeom) = -z(n+1:nstr,laymap,icci,m,igeom)
            ENDDO
#else
            IF (dosolar) THEN
              x(1:n) = -z(n+1:nstr,laymap,icci,m)
            ELSE
              x(1:n) = auxrad%air(1,i) - y0(n+1:nstr,lay) !- y1(n+1:nstr,lay) * (odpep(1) == 0.)
            ENDIF
#endif
          ELSE
            row = (lay-2) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)     ! kk(row+1:row+n,col+j)
              kkpack(lo:lo+n-1,col+j)   = -fp1(:,j,laymap,icci,m)
              lo = lo - n                           ! kk(row+1:row+n,col+n+j)
              kkpack(lo:lo+n-1,col+n+j) = -fp2(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = klu1 + (row + n + 1) - (col + j) ! kk(row+n+1:row+nstr,col+j)
              kkpack(lo:lo+n-1,col+j)   = -fm1(:,j,laymap,icci,m)
              lo = lo - n                           ! kk(row+n+1:row+nstr,col+n+j)
              kkpack(lo:lo+n-1,col+n+j) = -fm2(:,j,laymap,icci,m) * tlayer(j,lay)
            ENDDO
#ifdef _RTTOV_MFASIS_TRAINING
            DO igeom = 1, ngeom
              x(row+1:row+nstr,igeom) = x(row+1:row+nstr,igeom) + z(:,laymap,icci,m,igeom) * tausun(lev-1,igeom)
            ENDDO
#else
            IF (dosolar) THEN
              x(row+1:row+nstr) = x(row+1:row+nstr) + z(:,laymap,icci,m) * tausun(lev-1)
            ELSE
              x(row+1:row+nstr) = x(row+1:row+nstr) + y0(:,lay) + y1(:,lay) * opdep(lev-1)
            ENDIF
#endif
          ENDIF

          IF (lay < nlayers) THEN
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)     ! kk(row+1:row+n,col+j)
              kkpack(lo:lo+n-1,col+j)   = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = lo - n                           ! kk(row+1:row+n,col+n+j)
              kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m)
              lo = klu1 + (row + n + 1) - (col + j) ! kk(row+n+1:row+nstr,col+j)
              kkpack(lo:lo+n-1,col+j)   = fm1(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = lo - n                           ! kk(row+n+1:row+nstr,col+n+j)
              kkpack(lo:lo+n-1,col+n+j) = fm2(:,j,laymap,icci,m)
            ENDDO
#ifdef _RTTOV_MFASIS_TRAINING
            DO igeom = 1, ngeom
              x(row+1:row+nstr,igeom) = -z(:,laymap,icci,m,igeom) * tausun(lev,igeom)
            ENDDO
#else
            IF (dosolar) THEN
              x(row+1:row+nstr) = -z(:,laymap,icci,m) * tausun(lev)
            ELSE
              x(row+1:row+nstr) = -y0(:,lay) - y1(:,lay) * opdep(lev)
            ENDIF
#endif
          ELSE
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            ! Surface reflectance is decomposed azimuthally: for a Lambertian surface, the
            ! azimuthal integration of the downwelling reflected radiance vanishes for m > 0
            IF (m == 0 .AND. inc_surface) THEN
              IF (reflfac > 0._jprb) THEN
                DO j = 1, n
                  surfterm(j)   = reflfac * SUM(muqw(:) * fm1(:,j,laymap,icci,m))
                  surfterm(n+j) = reflfac * SUM(muqw(:) * fm2(:,j,laymap,icci,m))
                ENDDO
              ELSE
                surfterm = 0._jprb
              ENDIF
              DO j = 1, n
                lo = klu1 + (row + 1) - (col + j) ! kk(row+1:row+n,col+j)
                kkpack(lo:lo+n-1,col+j)   = (fp1(:,j,laymap,icci,m) - surfterm(j)) * tlayer(j,lay)
                lo = lo - n                       ! kk(row+1:row+n,col+n+j)
                kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m) - surfterm(n+j)
              ENDDO

#ifdef _RTTOV_MFASIS_TRAINING
              DO igeom = 1, ngeom
                ! The solar reflection here is between musun and muq so use diffuse BRDF
                x(row+1:row+n,igeom) = -tausun(lev,igeom) * (z(1:n,laymap,icci,m,igeom) - &
                                       reflfac * SUM(muqw(:) * z(n+1:nstr,laymap,icci,m,igeom)) - &
                                       (albedo / pi) * solar_spectrum(i) * musun(igeom))
              ENDDO
#else
              IF (dosolar) THEN
                ! The solar reflection here is between musun and muq so use diffuse BRDF
                x(row+1:row+n) = -tausun(lev) * (z(1:n,laymap,icci,m) - &
                                 reflfac * SUM(muqw(:) * z(n+1:nstr,laymap,icci,m)) - &
                                 (albedo / pi) * solar_spectrum(i) * musun)
              ELSE
                x(row+1:row+n) = -(y0(1:n,lay) + y1(1:n,lay) * opdep(lev) - &
                                 reflfac * SUM(muqw(:) * (y0(n+1:nstr,lay) + y1(n+1:nstr,lay) * opdep(lev))) - &
                                 emissivity(i)%emis_out * auxrad%skin(i))
              ENDIF
#endif
            ELSE
              DO j = 1, n
                lo = klu1 + (row + 1) - (col + j) ! kk(row+1:row+n,col+j)
                kkpack(lo:lo+n-1,col+j)   = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
                lo = lo - n                       ! kk(row+1:row+n,col+n+j)
                kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m)
              ENDDO
#ifdef _RTTOV_MFASIS_TRAINING
              DO igeom = 1, ngeom
                x(row+1:row+n,igeom) = -z(1:n,laymap,icci,m,igeom) * tausun(lev,igeom)
              ENDDO
#else
              IF (dosolar) THEN
                x(row+1:row+n) = -z(1:n,laymap,icci,m) * tausun(lev)
              ELSE
                x(row+1:row+n) = -y0(1:n,lay) - y1(1:n,lay) * opdep(lev)
              ENDIF
#endif
            ENDIF
          ENDIF

          ! --------------------------------------------------------------------------
          ! Accumulate TMS correction for this layer
          ! --------------------------------------------------------------------------
          IF (dosolar .AND. m == 0 .AND. fssa(lay) > 0._jprb) THEN
#ifdef _RTTOV_MFASIS_TRAINING
            DO j = 1, nchanprof
              IF (chanprof(j)%chan == chan) THEN
                igeom = chanprof(j)%prof
                taufac = tausun(lev-1,igeom) * tausat(lev-1,igeom) - tausun(lev,igeom) * tausat(lev,igeom)

                ! Single-scattering with unscaled SSA and interpolated phase function minus
                ! single-scattering with scaled SSA and phase function from nmom Legendre moments
                delta_tms(igeom) = delta_tms(igeom) + sunsrcfac / (1._jprb + musat(igeom) * musunr(igeom)) * taufac * &
                  (ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,j) / ftau(lay) - &
                   SUM(legcoef(:,laymap,icci) * legpolymuscata(:,0,0,igeom)))
              ENDIF
            ENDDO
#else
            taufac = tausun(lev-1) * tausat(lev-1) - tausun(lev) * tausat(lev)

            ! Single-scattering with unscaled SSA and interpolated phase function minus
            ! single-scattering with scaled SSA and phase function from nmom Legendre moments
            delta_tms = delta_tms + sunsrcfac / (1._jprb + musat * musunr) * taufac * &
              (ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) / ftau(lay) - &
               SUM(legcoef(:,laymap,icci) * legpolymuscata(:,0,0)))
#endif
          ENDIF

        ENDDO ! layer loop

        ! --------------------------------------------------------------------------
        ! Solve the multi-layer system
        ! --------------------------------------------------------------------------

        ! LAPACK banded matrix solver (much faster than the general solver)
        CALL DGBTRF(nstrlay, nstrlay, kl, ku, kkpack, klku, ipiv1, info)
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving linear system, DGBTRF error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF

#ifdef _RTTOV_MFASIS_TRAINING
        CALL DGBTRS('No transpose', nstrlay, kl, ku, ngeom, kkpack, klku, ipiv1, x(1:nstrlay,:), nstrlay, info)
#else
        CALL DGBTRS('No transpose', nstrlay, kl, ku, 1, kkpack, klku, ipiv1, x, nstrlay, info)
#endif
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving linear system, DGBTRS error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF

#ifndef _RTTOV_MFASIS_TRAINING
        IF (store) dom_state(i)%x(1:nstrlay,m,icc) = x(1:nstrlay)
#endif
        ! --------------------------------------------------------------------------
        ! Calculate TOA radiance
        ! --------------------------------------------------------------------------
        radinc = 0._jprb
        DO lay = 1, nlayers
          lev = lay + 1

          IF (fssa(lay) > 0._jprb) THEN
            laymap = profiles_dom(icc,i)%laymap(lay)
            IF (laymap > 0) THEN
              icci = ircld%icldarr(icc,laymap,prof)
            ELSE
              icci = 0
            ENDIF
            row = (lay-1) * nstr

#ifdef _RTTOV_MFASIS_TRAINING
            DO igeom = 1, ngeom
              ! Homogenous part - diffuse radiation
              DO j = 1, n
                taufac = tausat(lev-1,igeom) - tlayer(j,lay) * tausat(lev,igeom)
                radinc(igeom) = radinc(igeom) + x(row+j,igeom) * fpmuhsup(j,laymap,icci,m,igeom) * &
                                                taufac / (1._jprb + eval(j,laymap,icci,m) * musat(igeom))

                taufac = tlayer(j,lay) * tausat(lev-1,igeom) - tausat(lev,igeom)
                radinc(igeom) = radinc(igeom) + x(row+n+j,igeom) * fpmuhsdn(j,laymap,icci,m,igeom) * &
                                                taufac / (1._jprb - eval(j,laymap,icci,m) * musat(igeom))
              ENDDO

              ! Particular integral - solar source term
              taufac = tausat(lev-1,igeom) * tausun(lev-1,igeom) - tausat(lev,igeom) * tausun(lev,igeom)
              radinc(igeom) = radinc(igeom) + fpmupi(laymap,icci,m,igeom) * &
                                              taufac / (1._jprb + musunr(igeom) * musat(igeom))
            ENDDO
#else
            ! Homogenous part - diffuse radiation
            DO j = 1, n
              taufac = tausat(lev-1) - tlayer(j,lay) * tausat(lev)
              radinc = radinc + x(row+j) * fpmuhsup(j,laymap,icci,m) * &
                                taufac / (1._jprb + eval(j,laymap,icci,m) * musat)

              taufac = tlayer(j,lay) * tausat(lev-1) - tausat(lev)
              radinc = radinc + x(row+n+j) * fpmuhsdn(j,laymap,icci,m) * &
                                taufac / (1._jprb - eval(j,laymap,icci,m) * musat)
            ENDDO

            ! Particular integral - solar source term
            IF (dosolar) THEN
              taufac = tausat(lev-1) * tausun(lev-1) - tausat(lev) * tausun(lev)
              radinc = radinc + fpmupi(laymap,icci,m) * taufac / (1._jprb + musunr * musat)
            ENDIF
#endif
          ENDIF

#ifndef _RTTOV_MFASIS_TRAINING
          ! Particular integral - emissive source term
          IF (dothermal) THEN
            radinc = radinc + fpmupi0(lay) * (tausat(lev-1) - tausat(lev)) + &
                              fpmupi1(lay) * ((opdep(lev-1) + musat) * tausat(lev-1) - &
                                              (opdep(lev) + musat) * tausat(lev))
          ENDIF
#endif
        ENDDO ! nlayers


        ! Surface contribution - assumes surface is Lambertian
        IF (m == 0 .AND. inc_surface) THEN

#ifdef _RTTOV_MFASIS_TRAINING
          DO igeom = 1, ngeom
#endif
          ! Integrate over radiances along downwelling quadrature angles and then add the surface term
          radsurfup_sum = 0._jprb
          row = nstr * (nlayers - 1)
          laymap = profiles_dom(icc,i)%laymap(nlayers)
          IF (laymap > 0) THEN
            icci = ircld%icldarr(icc,laymap,prof)
          ELSE
            icci = 0
          ENDIF

          DO j = 1, n
#ifdef _RTTOV_MFASIS_TRAINING
            radsurfup_sum = radsurfup_sum + muqw(j) * &
              (SUM(x(row+1:row+n,igeom) * fm1(j,:,laymap,icci,m) * tlayer(:,nlayers)) + &
               SUM(x(row+n+1:row+nstr,igeom) * fm2(j,:,laymap,icci,m)))
#else
            radsurfup_sum = radsurfup_sum + muqw(j) * &
              (SUM(x(row+1:row+n) * fm1(j,:,laymap,icci,m) * tlayer(:,nlayers)) + &
               SUM(x(row+n+1:row+nstr) * fm2(j,:,laymap,icci,m)))
#endif
          ENDDO

          IF (dosolar) THEN
#ifdef _RTTOV_MFASIS_TRAINING
            DO j = 1, n
              radsurfup_sum = radsurfup_sum + muqw(j) * z(n+j,laymap,icci,m,igeom) * tausun(nlevels,igeom)
            ENDDO
            radsurfup = reflfac * radsurfup_sum + &
                        (albedo / pi) * solar_spectrum(i) * musun(igeom) * tausun(nlevels,igeom)
#else
            DO j = 1, n
              radsurfup_sum = radsurfup_sum + muqw(j) * z(n+j,laymap,icci,m) * tausun(nlevels)
            ENDDO
            radsurfup = reflfac * radsurfup_sum + &
                        (albedo / pi) * solar_spectrum(i) * musun * tausun(nlevels)
#endif
          ELSE
            DO j = 1, n
              radsurfup_sum = radsurfup_sum + muqw(j) * (y0(n+j,nlayers) + y1(n+j,nlayers) * opdep(nlevels))
            ENDDO
            radsurfup = reflfac * radsurfup_sum + emissivity(i)%emis_out * auxrad%skin(i)
          ENDIF

#ifdef _RTTOV_MFASIS_TRAINING
          radinc(igeom) = radinc(igeom) + radsurfup * tausat(nlevels,igeom)
          ENDDO
#else
          radinc = radinc + radsurfup * tausat(nlevels)
#endif

          IF (store) THEN
            dom_state(i)%radsurfup(icc)     = radsurfup
            dom_state(i)%radsurfup_sum(icc) = radsurfup_sum
          ENDIF
        ENDIF

        ! Multiplier for azimuthal component
        IF (dosolar) radinc = radinc * COS(m * relazi)
#ifdef _RTTOV_MFASIS_TRAINING
        WHERE (thismaxnaz >= m) thisrad = thisrad + radinc
#else
        thisrad = thisrad + radinc
#endif

#ifndef _RTTOV_MFASIS_TRAINING
        ! Quit azimuth loop when two terms have added less than (dom_accuracy * thisrad) to thisrad
        IF (ABS(radinc) < opts%rt_ir%dom_accuracy * thisrad) convcount = convcount + 1
        IF (convcount >= 2) EXIT
#endif
      ENDDO ! azimuth loop

      ! --------------------------------------------------------------------------
      ! Apply TMS correction and accumulate radiance contribution
      ! --------------------------------------------------------------------------
      IF (dosolar) thisrad = thisrad + delta_tms

#ifdef _RTTOV_MFASIS_TRAINING
      DO j = 1, nchanprof
        IF (chanprof(j)%chan == chan) THEN
          igeom = chanprof(j)%prof
          IF (icc == 0) THEN
            rad%clear(j) = rad%clear(j) + thisrad(igeom)
            rad%total(j) = rad%total(j) + thisrad(igeom) * ircld%xcolclr(prof)
          ELSE
            rad%total(j) = rad%total(j) + thisrad(igeom) * (ircld%xcol(icc+1,prof) - ircld%xcol(icc,prof))
          ENDIF
        ENDIF
      ENDDO
#else
      IF (icc == 0) THEN
        rad%clear(i) = rad%clear(i) + thisrad
        rad%total(i) = rad%total(i) + thisrad * ircld%xcolclr(prof)
      ELSE
        rad%total(i) = rad%total(i) + thisrad * (ircld%xcol(icc+1,prof) - ircld%xcol(icc,prof))
      ENDIF

      IF (store) dom_state(i)%thisrad(icc) = thisrad
#endif
    ENDDO ! cloud columns

    IF (store) dom_state(i)%eval = eval

  ENDDO ! channel

  rad%cloudy = rad%total

  ! --------------------------------------------------------------------------
  ! Tidy up
  ! --------------------------------------------------------------------------
  NULLIFY(fp2, fm2)

#ifndef _RTTOV_MFASIS_TRAINING
  IF (dosolar) THEN
    DEALLOCATE(legpolymusun, legpolymuscata, &
               tausun, fpmupi, z, stat=err)
  ELSE
    DEALLOCATE(b0, b1, fpmupi0, fpmupi1, y0, y1, stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM deallocation error')
#endif
CATCH

END SUBROUTINE rttov_dom
