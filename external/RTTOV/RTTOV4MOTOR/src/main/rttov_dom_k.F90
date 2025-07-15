! Description:
!> @file
!!   K of RTE integration using Discrete Ordinates Method.
!
!> @brief
!!   K of RTE integration using Discrete Ordinates Method.
!!
!! @details
!!   See direct model subroutine for details.
!!
!!   The dom_state argument holds intermediate calculation results from
!!   the direct model call which is more efficient than repeating
!!   all the direct model calculations. Nevertheless some repetition is
!!   necessary.
!!
!! @param[out]    err                    status on exit
!! @param[in]     chanprof               specifies channels and profiles to simulate
!! @param[in]     dosolar                flag to indicate if solar source term (true) or emission source term (false)
!!                                       should be calculated
!! @param[in]     chanflag               flags to indicate which channels should be allocated (either channels with
!!                                       significant solar or thermal contributions according to dosolar)
!! @param[in]     nstr                   number of DOM streams; user selected, but RTTOV ensures it is even
!! @param[in]     profiles               input atmospheric profiles and surface variables
!! @param[in]     profiles_dom           layer optical depth profiles set up for DOM calculations
!! @param[in,out] profiles_dom_k         DOM layer optical depth profile increments
!! @param[in]     maxnlayers             maximum number of layers in the DOM profiles
!! @param[in]     auxrad                 Planck radiances for emission source term
!! @param[in,out] auxrad_k               Planck radiance increments
!! @param[in]     trans_scatt_ir         single-scattering albedos
!! @param[in,out] trans_scatt_ir_k       single-scattering albedo increments
!! @param[in]     trans_scatt_ir_dyn     phase function Legendre decompositions
!! @param[in,out] trans_scatt_ir_dyn_k   phase function Legendre decomposition increments
!! @param[in]     emissivity             surface emissivities
!! @param[in,out] emissivity_k           surface emissivity increments
!! @param[in]     reflectance            surface BRDFs
!! @param[in,out] reflectance_k          surface BRDF increments
!! @param[in]     diffuse_refl           surface albedos
!! @param[in,out] diffuse_refl_k         surface albedo increments
!! @param[in]     solar_spectrum         top of atmosphere solar irradiance
!! @param[in]     raytracing             raytracing structure
!! @param[in]     ircld                  cloud column data
!! @param[in,out] ircld_k                cloud column data increments
!! @param[in,out] rad_k                  radiance increments
!! @param         dom_state              pointer to array of rttov_dom_state structures containing intermediate calculation
!!                                       results from the direct model run
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
SUBROUTINE rttov_dom_k(              &
              err,                   &
              chanprof,              &
              dosolar,               &
              chanflag,              &
              nstr,                  &
              profiles,              &
              profiles_dom,          &
              profiles_dom_k,        &
              maxnlayers,            &
              auxrad,                &
              auxrad_k,              &
              trans_scatt_ir,        &
              trans_scatt_ir_k,      &
              trans_scatt_ir_dyn,    &
              trans_scatt_ir_dyn_k,  &
              emissivity,            &
              emissivity_k,          &
              reflectance,           &
              reflectance_k,         &
              diffuse_refl,          &
              diffuse_refl_k,        &
              solar_spectrum,        &
              raytracing,            &
              ircld,                 &
              ircld_k,               &
              rad_k,                 &
              dom_state)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jprb, jpim, jplm

  USE rttov_types, ONLY :          &
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
!       matinv,                      &
      asymtx_ad

  USE rttov_lapack_mod, ONLY : &
      dgetrf, dgetrs, dgbtrf, dgbtrs
!INTF_ON

  IMPLICIT NONE

  INTEGER(jpim),                     INTENT(OUT)   :: err
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                     INTENT(IN)    :: dosolar
  LOGICAL(jplm),                     INTENT(IN)    :: chanflag(SIZE(chanprof))
  INTEGER(jpim),                     INTENT(IN)    :: nstr
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile_dom),           INTENT(IN)    :: profiles_dom(0:,:)
  TYPE(rttov_profile_dom),           INTENT(INOUT) :: profiles_dom_k(0:,:)
  INTEGER(jpim),                     INTENT(IN)    :: maxnlayers
  TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad
  TYPE(rttov_radiance_aux),          INTENT(INOUT) :: auxrad_k
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_k
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_dyn_k
  TYPE(rttov_emissivity),            INTENT(IN)    :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity),            INTENT(INOUT) :: emissivity_k(SIZE(chanprof))
  TYPE(rttov_reflectance),           INTENT(IN)    :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),           INTENT(INOUT) :: reflectance_k(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: diffuse_refl(SIZE(chanprof))
  REAL(jprb),                        INTENT(INOUT) :: diffuse_refl_k(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_ircld),                 INTENT(INOUT) :: ircld_k
  TYPE(rttov_radiance),              INTENT(INOUT) :: rad_k
  TYPE(rttov_dom_state),             INTENT(IN)    :: dom_state(SIZE(chanprof))
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=128) :: msg
  LOGICAL(jplm) :: dothermal, inc_surface

  INTEGER(jpim) :: n, maxnaz, thisnaz, nmom

  INTEGER(jpim) :: i, j, k, l, m, row, col
  INTEGER(jpim) :: icc, icci
  INTEGER(jpim) :: nlayers, nlevels, lay, lev, laymap
  INTEGER(jpim) :: nchanprof, prof, lastprof
  INTEGER(jpim) :: ilo, ihi, lo, hi, klu1

  INTEGER       :: kl, ku, klku, info, nstri, nstrlay ! LAPACK arguments, no KIND

  REAL(jprb)    :: musat, musatr, musun, musunr, relazi, muscata
  REAL(jprb)    :: taufac, albedo, reflfac, sunsrcfac
  REAL(jprb)    :: taufac_k, delta_tms_k, albedo_k, reflfac_k
  REAL(jprb)    :: radsurfup_sum_k, radsurfup_k, radinc_k, thisrad_k
  REAL(jprb)    :: tmpval, tmpval2

  REAL(jprb)              :: muq(nstr/2), muw(nstr/2), muqw(nstr/2), muqr(nstr/2)
  REAL(jprb), ALLOCATABLE :: legpolymuq(:,:,:) !(0:nstr-1,nstr,0:nstr-1)     ! dims are (0:nmom,nstr,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymusat(:,:,:) !(0:nstr-1,0:0,0:nstr-1)    ! dims are (0:nmom,0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymusun(:,:,:) !(0:nstr-1,0:0,0:nstr-1)    ! dims are (0:nmom,0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymuscata(:,:,:) !(0:nstr-1,0:0,0:0)
  REAL(jprb)              :: legcoef(0:nstr-1,profiles(1)%nlayers,0:1)
  REAL(jprb)              :: legcoef_k(0:nstr-1,profiles(1)%nlayers,0:1)

  REAL(jprb)              :: f_lay, ftau(maxnlayers), fssa(maxnlayers)
  REAL(jprb)              :: f_lay_k, ftau_k(maxnlayers), fssa_k(maxnlayers)
  REAL(jprb)              :: phase_fn(nstr/2,nstr)
  REAL(jprb)              :: phase_fn_k(nstr/2,nstr)
  REAL(jprb), ALLOCATABLE :: b0(:), b1(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: b0_k(:), b1_k(:) !(maxnlayers)
  REAL(jprb)              :: opdep(maxnlayers+1)
  REAL(jprb)              :: opdep_k(maxnlayers+1)
  REAL(jprb)              :: tausat(maxnlayers+1)
  REAL(jprb)              :: tausat_k(maxnlayers+1)
  REAL(jprb), ALLOCATABLE :: tausun(:) !(maxnlayers+1)
  REAL(jprb), ALLOCATABLE :: tausun_k(:) !(maxnlayers+1)
  REAL(jprb)              :: tlayer(nstr/2,maxnlayers)
  REAL(jprb)              :: tlayer_k(nstr/2,maxnlayers)
  REAL(jprb)              :: surfterm(nstr)
  REAL(jprb)              :: surfterm_k(nstr)
  REAL(jprb)              :: phasesat(nstr)
  REAL(jprb)              :: phasesat_k(nstr)

  REAL(jprb)              :: bp(nstr/2,nstr/2), bm(nstr/2,nstr/2)!, bpm(nstr/2,nstr/2)
  REAL(jprb)              :: bp_k(nstr/2,nstr/2), bm_k(nstr/2,nstr/2), bpm_k(nstr/2,nstr/2)
  REAL(jprb)              :: xp(nstr/2,nstr/2), xm(nstr/2,nstr/2)
  REAL(jprb)              :: xp_k(nstr/2,nstr/2), xm_k(nstr/2,nstr/2)
  REAL(jprb)              :: sp(nstr/2,nstr/2), sm(nstr/2,nstr/2)
  REAL(jprb)              :: sp_k(nstr/2,nstr/2), sm_k(nstr/2,nstr/2)
  REAL(jprb)              :: kkpack_k(9*nstr/2-2,nstr*maxnlayers)

  ! LAPACK arguments, no KIND
  INTEGER                 :: ipiv1(nstr*maxnlayers), ipiv2(nstr)
  DOUBLE PRECISION        :: hh(nstr,nstr)!, hhinvt(nstr,nstr)
  DOUBLE PRECISION        :: hh_k(nstr,nstr), hhtmp_k(nstr,nstr), x_k(nstr*maxnlayers)
  DOUBLE PRECISION        :: kkpack(9*nstr/2-2,nstr*maxnlayers)
  DOUBLE PRECISION        :: kktmp_k(nstr*maxnlayers,nstr*maxnlayers+1)
  DOUBLE PRECISION, ALLOCATABLE :: z_k(:,:,:,:) !(nstr,0:profiles(1)%nlayers,0:1,0:nstr-1)
  DOUBLE PRECISION, ALLOCATABLE :: y0_k(:,:) !(nstr,maxnlayers)
  DOUBLE PRECISION, ALLOCATABLE :: y1_k(:,:) !(nstr,maxnlayers)

  REAL(jprb)              :: eval(nstr/2)
  REAL(jprb), ALLOCATABLE :: eval_k(:,:,:,:) !(nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsup(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsup_k(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsdn(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsdn_k(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi_k(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi0(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi0_k(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi1(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi1_k(:) !(maxnlayers)

  REAL(jprb), TARGET, ALLOCATABLE :: fp1(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fp1_k(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fm1(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fm1_k(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), POINTER :: fp2(:,:,:,:,:), fm2(:,:,:,:,:)
  REAL(jprb), POINTER :: fp2_k(:,:,:,:,:), fm2_k(:,:,:,:,:)
  REAL(jprb), POINTER :: ssa(:,:,:), ssa_k(:,:,:)

  ! This is only an internal logical array so use smallest available KIND
  LOGICAL(jpit), ALLOCATABLE :: done_lay_str_azi(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
!   LOGICAL(jpit), ALLOCATABLE :: done_lay_str_azi_k(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)

  ! --------------------------------------------------------------------------------------

  TRY

  ! --------------------------------------------------------------------------
  ! Initialisation
  ! --------------------------------------------------------------------------
  nchanprof = SIZE(chanprof)
  dothermal = .NOT. dosolar

  ! nstr is guaranteed to be even so this is safe
  n = nstr / 2_jpim
  nstri = nstr

  ! Number of Legendre moments and azimuthal components are determined by nstr
  nmom = nstr - 1
  IF (dosolar) THEN
    maxnaz = nstr - 1
  ELSE
    maxnaz = 0
  ENDIF

  ! Index variables for packing kkpack
  kl = 3 * n - 1
  ku = 3 * n - 1
  klu1 = kl + ku + 1
  klku = 2 * kl + ku + 1

  ! Allocations to reduce memory usage
  ALLOCATE(legpolymuq(0:nstr-1,nstr,0:maxnaz), &
           legpolymusat(0:nstr-1,0:0,0:maxnaz), &
           eval_k(n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsup(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsup_k(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsdn(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsdn_k(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fp1(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fp1_k(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fm1(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fm1_k(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           done_lay_str_azi(profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
!            done_lay_str_azi_k(profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
  THROWM(err.NE.0, 'DOM allocation error')

  IF (dosolar) THEN
    ALLOCATE(legpolymusun(0:nstr-1,0:0,0:maxnaz), &
             legpolymuscata(0:nstr-1,0:0,0:0), &
             tausun(maxnlayers+1), &
             tausun_k(maxnlayers+1), &
             fpmupi(profiles(1)%nlayers,0:1,0:maxnaz), &
             fpmupi_k(profiles(1)%nlayers,0:1,0:maxnaz), &
             z_k(nstr,0:profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
  ELSE
    ALLOCATE(b0(maxnlayers), b1(maxnlayers), &
             b0_k(maxnlayers), b1_k(maxnlayers), &
             fpmupi0(maxnlayers), fpmupi1(maxnlayers), &
             fpmupi0_k(maxnlayers), fpmupi1_k(maxnlayers), &
             y0_k(nstr,maxnlayers), y1_k(nstr,maxnlayers), stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM allocation error')

  ! Exploit symmetry of eigenvectors using pointers
  fp2 => fm1
  fm2 => fp1
  fp2_k => fm1_k
  fm2_k => fp1_k

  IF (dosolar) THEN
    ssa => trans_scatt_ir%ssa_solar
    ssa_k => trans_scatt_ir_k%ssa_solar
  ELSE
    ssa => trans_scatt_ir%ssa_thermal
    ssa_k => trans_scatt_ir_k%ssa_thermal
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
    fp1(:,:,0,0,m)    = 0._jprb
    fm1(:,:,0,0,m)    = 0._jprb
    DO j = 1, n
      fm1(j,j,0,0,m) = -muqr(j)
    ENDDO
    ! See notes below regarding fp2/fm2 (i.e. eigenvectors for -ve eigenvalues)
  ENDDO

  lastprof = -1

  ! --------------------------------------------------------------------------
  ! Channel loop
  ! --------------------------------------------------------------------------
  DO i = 1, nchanprof
    IF (.NOT. chanflag(i)) CYCLE

    prof = chanprof(i)%prof

    IF (prof /= lastprof) THEN
      ! Angles/cosines - plane parallel means not active in TL/AD/K
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

      ! Evaluate Legendre polynomials at satellite and solar zenith angles and
      ! evaluate Legendre polynomials at scattering angle for TMS correction
      CALL calc_leg_poly(nmom, maxnaz, (/ musat /), legpolymusat)
      IF (dosolar) THEN
        CALL calc_leg_poly(nmom, maxnaz, (/ -musun /), legpolymusun)
        muscata = -musat * musun + SQRT((1._jprb - musun ** 2) * (1._jprb - musat ** 2)) * COS(relazi)
        CALL calc_leg_poly(nmom, 0_jpim, (/ muscata /), legpolymuscata)
      ENDIF

      lastprof = prof
    ENDIF

    ! The top level in profiles_dom is always the space boundary
    opdep(1) = 0._jprb
    tausat(1) = 1._jprb
    IF (dosolar) tausun(1) = 1._jprb


    done_lay_str_azi = .FALSE.
!     done_lay_str_azi_k = .FALSE.

    DO icc = ircld%ncolumn(prof), 0, -1
      profiles_dom_k(icc,i)%layerod(:) = 0._jprb
      thisrad_k   = 0._jprb
      delta_tms_k = 0._jprb
      albedo_k = 0._jprb
      reflfac_k = 0._jprb

      ftau_k = 0._jprb
      fssa_k = 0._jprb
      opdep_k = 0._jprb
      tausat_k = 0._jprb
      IF (dosolar) tausun_k = 0._jprb
      IF (dothermal) THEN
        fpmupi0_k = 0._jprb
        fpmupi1_k = 0._jprb
        b1_k = 0._jprb
        b0_k = 0._jprb
        y0_k = 0._jprb
        y1_k = 0._jprb
      ENDIF

      legcoef_k = 0._jprb
      eval_k = 0._jprb
      fp1_k = 0._jprb
      fm1_k = 0._jprb
      fpmuhsup_k = 0._jprb
      fpmuhsdn_k = 0._jprb
      IF (dosolar) THEN
        fpmupi_k = 0._jprb
        z_k = 0._jprb
      ENDIF

      IF (icc == 0) THEN
        ircld_k%xcolclr(i) = ircld_k%xcolclr(i) + dom_state(i)%thisrad(icc) * rad_k%total(i)
        thisrad_k = thisrad_k + ircld%xcolclr(prof) * rad_k%total(i)

! Do not accumulate from clear - input perturbation is only in total
!         thisrad_k = thisrad_k + rad_k%clear(i)

      ELSE
        ircld_k%xcol(icc+1,i) = ircld_k%xcol(icc+1,i) + dom_state(i)%thisrad(icc) * rad_k%total(i)
        ircld_k%xcol(icc,i) = ircld_k%xcol(icc,i) - dom_state(i)%thisrad(icc) * rad_k%total(i)
        thisrad_k = thisrad_k + rad_k%total(i) * (ircld%xcol(icc+1,prof) - ircld%xcol(icc,prof))
      ENDIF

      IF (dosolar) delta_tms_k = delta_tms_k + thisrad_k


      ! Surface reflectance factor
      IF (profiles_dom(icc,i)%surface) THEN
        IF (dosolar) THEN
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

      IF (dothermal) THEN
        thisnaz = 0
      ELSE
        thisnaz = dom_state(i)%nazloops(icc)
      ENDIF

      ! --------------------------------------------------------------------------
      ! Azimuthal component loop
      ! --------------------------------------------------------------------------
      DO m = thisnaz, 0, -1

        ! Solar source term factor
        IF (dosolar) THEN
          IF (m == 0) THEN
            sunsrcfac = z4pi_r * solar_spectrum(i)
          ELSE
            sunsrcfac = 2._jprb * z4pi_r * solar_spectrum(i)
          ENDIF
        ENDIF

        kkpack = 0._jprb

        ! Direct model layer loop

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
          IF (m == thisnaz) THEN
            IF (laymap > 0) THEN
              f_lay     = trans_scatt_ir_dyn%phasefn(i)%legcoef(nstr,icci,laymap) / (2._jprb * nstr + 1._jprb)
              ftau(lay) = 1._jprb - ssa(icci,laymap,i) * f_lay
              fssa(lay) = (1._jprb - f_lay) * ssa(icci,laymap,i) / ftau(lay)
            ELSE
              f_lay     = 0._jprb
              ftau(lay) = 1._jprb
              fssa(lay) = 0._jprb
            ENDIF
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate optical depths and transmittances
          ! --------------------------------------------------------------------------
          ! Calculate accumulated delta-scaled optical depth and transmittances on
          ! satellite-surface and sun-surface paths.
          IF (m == thisnaz) THEN  ! This calculation applies to all azi components
            opdep(lev) = opdep(lev-1) + ftau(lay) * profiles_dom(icc,i)%layerod(lay)
            tausat(lev) = EXP(-opdep(lev) * musatr)
            IF (dosolar) tausun(lev) = EXP(-opdep(lev) * musunr)
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

            ! For clear layers we have the solution and the interpolated thermal source term is trivial
            IF (.NOT. fssa(lay) > 0._jprb) THEN
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
              IF (m == thisnaz) THEN  ! This calculation applies to all azi components
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

              ! Include quadrature weights in phase fns evaluated for musat and each quadrature angle
              DO j = 1, n
                phasesat(j) = muw(j) * SUM(legcoef(m:nmom,laymap,icci) * &
                                           legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO
              DO j = n+1, nstr
                phasesat(j) = muw(j-n) * SUM(legcoef(m:nmom,laymap,icci) * &
                                             legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO

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
!               bpm = MATMUL(sm, sp)
!               bpmsave    = bpm  ! Used for solar order-n PI below - not implemented in AD (yet)

              xp = dom_state(i)%xp(:,:,laymap,icci,m)

              DO j = 1, n
                xm(:,j) = MATMUL(sp, xp(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
              ENDDO

              ! Rearrange to obtain F+ and F- (the eigenvectors)
              ! For +ve eigenvalues:
              fp1(:,:,laymap,icci,m)    = 0.5_jprb * (xp + xm)
              fm1(:,:,laymap,icci,m)    = 0.5_jprb * (xp - xm)

            ENDIF ! not done_lay_str_azi


            ! --------------------------------------------------------------------------
            ! Interpolate HS and PI for final radiance calculation
            ! --------------------------------------------------------------------------

            IF (.NOT. done_lay_str_azi(laymap,icci,m)) THEN

              ! Homogenous part - diffuse radiation, done once per layer (applies to all cloud columns)
              DO j = 1, n
                fpmuhsup(j,laymap,icci,m) = 0.5_jprb * (SUM(fp1(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm1(:,j,laymap,icci,m) * phasesat(n+1:nstr)))

                fpmuhsdn(j,laymap,icci,m) = 0.5_jprb * (SUM(fp2(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm2(:,j,laymap,icci,m) * phasesat(n+1:nstr)))
              ENDDO

              ! Particular integral - solar source term, done once per layer (applies to all cloud columns)
              IF (dosolar) THEN
                fpmupi(laymap,icci,m) = 0.5_jprb * SUM(dom_state(i)%z(:,laymap,icci,m) * phasesat(:)) + &
                                        sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                        legpolymusat(m:nmom,0,m) * legpolymusun(m:nmom,0,m))
              ENDIF

            ENDIF ! not done_lay_str_azi

            ! Particular integral - thermal source term, done once per layer per cloud column
            IF (dothermal) THEN
              fpmupi0(lay) = 0.5_jprb * SUM(dom_state(i)%y0(:,lay,icc) * phasesat(:)) + b0(lay)
              fpmupi1(lay) = 0.5_jprb * SUM(dom_state(i)%y1(:,lay,icc) * phasesat(:)) + b1(lay)
            ENDIF

            done_lay_str_azi(laymap,icci,m) = .TRUE.

          ENDIF ! fssa > 0 (i.e. scattering layer)

          ! --------------------------------------------------------------------------
          ! Calculate layer "transmittances" with opdeps scaled by eigenvalues
          ! --------------------------------------------------------------------------
          tlayer(:,lay) = EXP(-ftau(lay) * profiles_dom(icc,i)%layerod(lay) * dom_state(i)%eval(:,laymap,icci,m))

          ! --------------------------------------------------------------------------
          ! Construct part of linear system for this layer
          ! --------------------------------------------------------------------------
          IF (lay == 1) THEN
            row = 0
            col = 0
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,j)   = fm1(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,n+j) = fm2(:,j,laymap,icci,m) * tlayer(j,lay)
            ENDDO
          ELSE
            row = (lay-2) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)   = -fp1(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j) = -fp2(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = klu1 + (row + n + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)   = -fm1(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j) = -fm2(:,j,laymap,icci,m) * tlayer(j,lay)
            ENDDO
          ENDIF

          IF (lay < nlayers) THEN
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)   = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m)
              lo = klu1 + (row + n + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)   = fm1(:,j,laymap,icci,m) * tlayer(j,lay)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j) = fm2(:,j,laymap,icci,m)
            ENDDO
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
                lo = klu1 + (row + 1) - (col + j)
                kkpack(lo:lo+n-1,col+j)   = (fp1(:,j,laymap,icci,m) - surfterm(j)) * tlayer(j,lay)
                lo = lo - n
                kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m) - surfterm(n+j)
              ENDDO
            ELSE
              DO j = 1, n
                lo = klu1 + (row + 1) - (col + j)
                kkpack(lo:lo+n-1,col+j)   = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
                lo = lo - n
                kkpack(lo:lo+n-1,col+n+j) = fp2(:,j,laymap,icci,m)
              ENDDO
            ENDIF
          ENDIF

        ENDDO ! layer loop

        ! End of direct model layer loop


        ! Adjoint calculations

        radinc_k = 0._jprb
        tlayer_k = 0._jprb
        x_k = 0._jprb

        ! Multiplier for azimuthal component
        radinc_k = radinc_k + thisrad_k
        IF (dosolar) radinc_k = radinc_k * COS(m * relazi)


        ! Surface contribution - assumes surface is Lambertian
        IF (m == 0 .AND. inc_surface) THEN

          ! Integrate over radiances along downwelling quadrature angles and then add the surface term
          radsurfup_sum_k = 0._jprb
          radsurfup_k = 0._jprb

          row = nstr * (nlayers - 1)
          laymap = profiles_dom(icc,i)%laymap(nlayers)
          IF (laymap > 0) THEN
            icci = ircld%icldarr(icc,laymap,prof)
          ELSE
            icci = 0
          ENDIF

          tausat_k(nlevels) = tausat_k(nlevels) + radinc_k * dom_state(i)%radsurfup(icc)
          radsurfup_k = radsurfup_k + radinc_k * tausat(nlevels)

          IF (dosolar) THEN
            tausun_k(nlevels) = tausun_k(nlevels) + &
                (albedo / pi) * solar_spectrum(i) * musun * radsurfup_k
            albedo_k = albedo_k + &
                solar_spectrum(i) * musun * tausun(nlevels) * radsurfup_k / pi
            radsurfup_sum_k = radsurfup_sum_k + reflfac * radsurfup_k
            reflfac_k = reflfac_k + dom_state(i)%radsurfup_sum(icc) * radsurfup_k

            DO j = n, 1, -1
              tausun_k(nlevels) = tausun_k(nlevels) + muqw(j) * dom_state(i)%z(n+j,laymap,icci,m) * radsurfup_sum_k
              z_k(n+j,laymap,icci,m) = z_k(n+j,laymap,icci,m) + muqw(j) * tausun(nlevels) * radsurfup_sum_k
            ENDDO

          ELSE
            auxrad_k%skin(i) = auxrad_k%skin(i) + emissivity(i)%emis_out * radsurfup_k
            emissivity_k(i)%emis_out = emissivity_k(i)%emis_out + auxrad%skin(i) * radsurfup_k
            radsurfup_sum_k = radsurfup_sum_k + reflfac * radsurfup_k
            reflfac_k = reflfac_k + dom_state(i)%radsurfup_sum(icc) * radsurfup_k

            DO j = n, 1, -1
              opdep_k(nlevels) = opdep_k(nlevels) + muqw(j) * dom_state(i)%y1(n+j,nlayers,icc) * radsurfup_sum_k
              y1_k(n+j,nlayers) = y1_k(n+j,nlayers) + muqw(j) * opdep(nlevels) * radsurfup_sum_k
              y0_k(n+j,nlayers) = y0_k(n+j,nlayers) + muqw(j) * radsurfup_sum_k
            ENDDO

          ENDIF

          DO j = n, 1, -1
            fm2_k(j,:,laymap,icci,m) = fm2_k(j,:,laymap,icci,m) + &
                muqw(j) * dom_state(i)%x(row+n+1:row+nstr,m,icc) * radsurfup_sum_k
            x_k(row+n+1:row+nstr) = x_k(row+n+1:row+nstr) + &
                muqw(j) * fm2(j,:,laymap,icci,m) * radsurfup_sum_k
            tlayer_k(:,nlayers) = tlayer_k(:,nlayers) + &
                muqw(j) * dom_state(i)%x(row+1:row+n,m,icc) * fm1(j,:,laymap,icci,m) * radsurfup_sum_k
            fm1_k(j,:,laymap,icci,m) = fm1_k(j,:,laymap,icci,m) + &
                muqw(j) * dom_state(i)%x(row+1:row+n,m,icc) * tlayer(:,nlayers) * radsurfup_sum_k
            x_k(row+1:row+n) = x_k(row+1:row+n) + &
                muqw(j) * fm1(j,:,laymap,icci,m) * tlayer(:,nlayers) * radsurfup_sum_k
          ENDDO

        ENDIF


        ! --------------------------------------------------------------------------
        ! Calculate TOA radiance
        ! --------------------------------------------------------------------------

        DO lay = nlayers, 1, -1
          lev = lay + 1

          ! Particular integral - emissive source term
          IF (dothermal) THEN
            fpmupi0_k(lay) = fpmupi0_k(lay) + radinc_k * (tausat(lev-1) - tausat(lev))
            tausat_k(lev-1) = tausat_k(lev-1) + radinc_k * fpmupi0(lay)
            tausat_k(lev) = tausat_k(lev) - radinc_k * fpmupi0(lay)

            fpmupi1_k(lay) = fpmupi1_k(lay) + radinc_k * &
                              ((opdep(lev-1) + musat) * tausat(lev-1) - &
                               (opdep(lev) + musat) * tausat(lev))
            opdep_k(lev-1) = opdep_k(lev-1) + radinc_k * fpmupi1(lay) * tausat(lev-1)
            opdep_k(lev) = opdep_k(lev) - radinc_k * fpmupi1(lay) * tausat(lev)
            tausat_k(lev-1) = tausat_k(lev-1) + radinc_k * fpmupi1(lay) * (opdep(lev-1) + musat)
            tausat_k(lev) = tausat_k(lev) - radinc_k * fpmupi1(lay) * (opdep(lev) + musat)
          ENDIF

          IF (fssa(lay) > 0._jprb) THEN
            laymap = profiles_dom(icc,i)%laymap(lay)
            IF (laymap > 0) THEN
              icci = ircld%icldarr(icc,laymap,prof)
            ELSE
              icci = 0
            ENDIF
            row = (lay-1) * nstr

            ! Particular integral - solar source term
            IF (dosolar) THEN
              taufac   = tausat(lev-1) * tausun(lev-1) - tausat(lev) * tausun(lev)

              taufac_k = &!taufac_k +
                         radinc_k * fpmupi(laymap,icci,m) / (1._jprb + musunr * musat)
              fpmupi_k(laymap,icci,m) = fpmupi_k(laymap,icci,m) + radinc_k * &
                                        taufac / (1._jprb + musunr * musat)
              tausat_k(lev-1) = tausat_k(lev-1) + taufac_k * tausun(lev-1)
              tausun_k(lev-1) = tausun_k(lev-1) + taufac_k * tausat(lev-1)
              tausat_k(lev) = tausat_k(lev) - taufac_k * tausun(lev)
              tausun_k(lev) = tausun_k(lev) - taufac_k * tausat(lev)
            ENDIF

            ! Homogenous part - diffuse radiation
            DO j = n, 1, -1
              taufac   = tlayer(j,lay) * tausat(lev-1) - tausat(lev)

              tmpval = 1._jprb / (1._jprb - dom_state(i)%eval(j,laymap,icci,m) * musat)
              eval_k(j,laymap,icci,m) = eval_k(j,laymap,icci,m) + &
                  radinc_k * musat * dom_state(i)%x(row+n+j,m,icc) * &
                  fpmuhsdn(j,laymap,icci,m) * taufac * tmpval**2
              taufac_k = &!taufac_k + &
                  radinc_k * dom_state(i)%x(row+n+j,m,icc) * &
                  fpmuhsdn(j,laymap,icci,m) * tmpval
              fpmuhsdn_k(j,laymap,icci,m) = fpmuhsdn_k(j,laymap,icci,m) + &
                  radinc_k * dom_state(i)%x(row+n+j,m,icc) * taufac * tmpval
              x_k(row+n+j) = x_k(row+n+j) + &
                  radinc_k * fpmuhsdn(j,laymap,icci,m) * taufac * tmpval

              tlayer_k(j,lay) = tlayer_k(j,lay) + taufac_k * tausat(lev-1)
              tausat_k(lev-1) = tausat_k(lev-1) + taufac_k * tlayer(j,lay)
              tausat_k(lev) = tausat_k(lev) - taufac_k


              taufac   = tausat(lev-1) - tlayer(j,lay) * tausat(lev)

              tmpval = 1._jprb / (1._jprb + dom_state(i)%eval(j,laymap,icci,m) * musat)
              eval_k(j,laymap,icci,m) = eval_k(j,laymap,icci,m) - &
                  radinc_k * musat * dom_state(i)%x(row+j,m,icc) * &
                  fpmuhsup(j,laymap,icci,m) * taufac * tmpval**2
              taufac_k = &!taufac_k + &
                  radinc_k * dom_state(i)%x(row+j,m,icc) * &
                  fpmuhsup(j,laymap,icci,m) * tmpval
              fpmuhsup_k(j,laymap,icci,m) = fpmuhsup_k(j,laymap,icci,m) + &
                  radinc_k * dom_state(i)%x(row+j,m,icc) * taufac * tmpval
              x_k(row+j) = x_k(row+j) + &
                  radinc_k * fpmuhsup(j,laymap,icci,m) * taufac * tmpval

              tausat_k(lev-1) = tausat_k(lev-1) + taufac_k
              tlayer_k(j,lay) = tlayer_k(j,lay) - taufac_k * tausat(lev)
              tausat_k(lev) = tausat_k(lev) - taufac_k * tlayer(j,lay)
            ENDDO
          ENDIF
        ENDDO ! nlayers


        ! Direct:  kk * x = b     (in direct model we re-use variable x as b; lapack overwrites the array)
        ! TL: kk * x_tl = b_tl - kk_tl * b

        ! K:  1. kk_k = kk_k - (kk^-1)^T * x_k * x^T
        !     2. let c be solution to kk^T * c = x_k; then b_k = b_k + c   (can overwrite x_k so x_k = c)
        !           (can solve using DGBTRS with "Transpose" option)
        !

        CALL DGBTRF(nstrlay, nstrlay, kl, ku, kkpack, klku, ipiv1, info)
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving K linear system, DGBTRF error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF

        ! 1. Solve each column of kk_k as a linear system. We have:
        !           kk_k = -(kk^-1)^T * x_k * x^T    (no accumulation required)
        ! =>      kk_k^T = -x * x_k^T * kk^-1
        ! => kk_k^T * kk = -x * x_k^T
        ! => kk^T * kk_k = -x_k * x^T
        ! so each column of kk_k(:,j) is the solution c to kk^T * c = -(x_k * x^T)(:,j)
        ! Do this using a temporary array, solving and writing directly to kkpack_k

        ! 2. Solve for x_k

        ! We can do both 1 and 2 in a single call to DGBTRS:
        DO j = 1, nstrlay
          kktmp_k(1:nstrlay,j) = -x_k(1:nstrlay) * dom_state(i)%x(j,m,icc)
        ENDDO
        kktmp_k(1:nstrlay,nstrlay+1) = x_k(1:nstrlay) ! Solve for x_k using the last column of kktmp_k
        CALL DGBTRS('Transpose', nstrlay, kl, ku, nstrlay+1, kkpack, klku, &
                    ipiv1, kktmp_k(1:nstrlay,1:nstrlay+1), nstrlay, info)
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving K linear system, DGBTRS error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF
        DO j = 1, nstrlay
          ilo = MAX(1, j - ku)
          ihi = MIN(nstrlay, j + kl)
          lo = klu1 + ilo - j
          hi = klu1 + ihi - j
          kkpack_k(lo:hi,j) = kktmp_k(ilo:ihi,j)
        ENDDO
        x_k(1:nstrlay) = kktmp_k(1:nstrlay,nstrlay+1) ! Solve for x_k: copy result back


        ! Adjoint layer loop

        ! --------------------------------------------------------------------------
        ! Layer loop
        ! --------------------------------------------------------------------------
        DO lay = nlayers, 1, -1
          lev = lay + 1_jpim

          ! laymap == 0 for clear layers
          laymap = profiles_dom(icc,i)%laymap(lay)
          IF (laymap > 0) THEN
            icci = ircld%icldarr(icc,laymap,prof)
          ELSE
            icci = 0
          ENDIF

          ! f_lay could be stored per layer, but isn't in the direct so recompute here
          IF (m == 0) THEN
            IF (laymap > 0) THEN
              f_lay = trans_scatt_ir_dyn%phasefn(i)%legcoef(nstr,icci,laymap) / (2._jprb * nstr + 1._jprb)
            ELSE
              f_lay = 0._jprb
            ENDIF
          ENDIF
          f_lay_k = 0._jprb

          IF (fssa(lay) > 0._jprb) THEN
            phase_fn_k = 0._jprb
            phasesat_k = 0._jprb
          ENDIF

! Require phase_fn and phasesat on every layer, will have been overwritten in direct loop above **
          IF (fssa(lay) > 0._jprb) THEN
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

            ! Include quadrature weights in phase fns evaluated for musat and each quadrature angle
            DO j = 1, n
              phasesat(j) = muw(j) * SUM(legcoef(m:nmom,laymap,icci) * &
                                         legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
            ENDDO
            DO j = n+1, nstr
              phasesat(j) = muw(j-n) * SUM(legcoef(m:nmom,laymap,icci) * &
                                           legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
            ENDDO
          ENDIF

          ! --------------------------------------------------------------------------
          ! Accumulate TMS correction for this layer
          ! --------------------------------------------------------------------------
          IF (dosolar .AND. m == 0 .AND. fssa(lay) > 0._jprb) THEN

            taufac = tausun(lev-1) * tausat(lev-1) - tausun(lev) * tausat(lev)

            tmpval = sunsrcfac / (1._jprb + musat * musunr)
            tmpval2 = 1._jprb / ftau(lay)

            legcoef_k(:,laymap,icci) = legcoef_k(:,laymap,icci) - &
                delta_tms_k * tmpval * taufac * legpolymuscata(:,0,0)
            ftau_k(lay) = ftau_k(lay) - &
                delta_tms_k * tmpval * taufac * &
                ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) * tmpval2**2
            trans_scatt_ir_k%phup(icci,laymap,i) = trans_scatt_ir_k%phup(icci,laymap,i) + &
                delta_tms_k * tmpval * taufac * &
                ssa(icci,laymap,i) * tmpval2
            ssa_k(icci,laymap,i) = ssa_k(icci,laymap,i) + &
                delta_tms_k * tmpval * taufac * &
                trans_scatt_ir%phup(icci,laymap,i) * tmpval2
            taufac_k = &!taufac_k + &
                delta_tms_k * tmpval * &
                (ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) * tmpval2 - &
                 SUM(legcoef(:,laymap,icci) * legpolymuscata(:,0,0)))

            tausun_k(lev-1) = tausun_k(lev-1) + taufac_k * tausat(lev-1)
            tausat_k(lev-1) = tausat_k(lev-1) + taufac_k * tausun(lev-1)
            tausun_k(lev) = tausun_k(lev) - taufac_k * tausat(lev)
            tausat_k(lev) = tausat_k(lev) - taufac_k * tausun(lev)
          ENDIF


          ! --------------------------------------------------------------------------
          ! Construct part of linear system for this layer
          ! --------------------------------------------------------------------------
          surfterm_k = 0._jprb

          IF (lay < nlayers) THEN
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            IF (dosolar) THEN
              tausun_k(lev) = tausun_k(lev) - SUM(x_k(row+1:row+nstr) * dom_state(i)%z(:,laymap,icci,m))
              z_k(:,laymap,icci,m) = z_k(:,laymap,icci,m) - x_k(row+1:row+nstr) * tausun(lev)
            ELSE
              opdep_k(lev) = opdep_k(lev) - SUM(x_k(row+1:row+nstr) * dom_state(i)%y1(:,lay,icc))
              y1_k(:,lay) = y1_k(:,lay) - x_k(row+1:row+nstr) * opdep(lev)
              y0_k(:,lay) = y0_k(:,lay) - x_k(row+1:row+nstr)
            ENDIF
            DO j = n, 1, -1
              lo = klu1 + (row + 1) - (col + j)
              fp1_k(:,j,laymap,icci,m) = fp1_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+j) * tlayer(j,lay)
              tlayer_k(j,lay) = tlayer_k(j,lay) + SUM(kkpack_k(lo:lo+n-1,col+j) * fp1(:,j,laymap,icci,m))

              lo = lo - n
              fp2_k(:,j,laymap,icci,m) = fp2_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+n+j)

              lo = klu1 + (row + n + 1) - (col + j)
              fm1_k(:,j,laymap,icci,m) = fm1_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+j) * tlayer(j,lay)
              tlayer_k(j,lay) = tlayer_k(j,lay) + SUM(kkpack_k(lo:lo+n-1,col+j) * fm1(:,j,laymap,icci,m))

              lo = lo - n
              fm2_k(:,j,laymap,icci,m) = fm2_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+n+j)
            ENDDO
          ELSE
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            ! Surface reflectance is decomposed azimuthally: for a Lambertian surface, the
            ! azimuthal integration of the downwelling reflected radiance vanishes for m > 0
            IF (m == 0 .AND. inc_surface) THEN
              tmpval = SUM(x_k(row+1:row+n))
              IF (dosolar) THEN
                albedo_k = albedo_k + tausun(lev) * tmpval * solar_spectrum(i) * musun / pi
                z_k(n+1:nstr,laymap,icci,m) = z_k(n+1:nstr,laymap,icci,m) + &
                    tausun(lev) * tmpval * reflfac * muqw(:)
                reflfac_k = reflfac_k + &
                    tausun(lev) * tmpval * SUM(muqw(:) * dom_state(i)%z(n+1:nstr,laymap,icci,m))
                z_k(1:n,laymap,icci,m) = z_k(1:n,laymap,icci,m) - tausun(lev) * x_k(row+1:row+n)
                tausun_k(lev) = tausun_k(lev) - SUM(x_k(row+1:row+n) * (dom_state(i)%z(1:n,laymap,icci,m) - &
                                    reflfac * SUM(muqw(:) * dom_state(i)%z(n+1:nstr,laymap,icci,m)) - &
                                    (albedo / pi) * solar_spectrum(i) * musun))
              ELSE
                auxrad_k%skin(i) = auxrad_k%skin(i) + tmpval * emissivity(i)%emis_out
                emissivity_k(i)%emis_out = emissivity_k(i)%emis_out + tmpval * auxrad%skin(i)

                opdep_k(lev) = opdep_k(lev) + &
                    reflfac * SUM(x_k(row+1:row+n)) * SUM(muqw(:) * dom_state(i)%y1(n+1:nstr,lay,icc))
                y1_k(n+1:nstr,lay) = y1_k(n+1:nstr,lay) + reflfac * muqw(:) * tmpval * opdep(lev)
                y0_k(n+1:nstr,lay) = y0_k(n+1:nstr,lay) + reflfac * muqw(:) * tmpval

                reflfac_k = reflfac_k + tmpval * SUM(muqw(:) * (dom_state(i)%y0(n+1:nstr,lay,icc) + &
                                                                dom_state(i)%y1(n+1:nstr,lay,icc) * opdep(lev)))

                opdep_k(lev) = opdep_k(lev) - SUM(x_k(row+1:row+n) * dom_state(i)%y1(1:n,lay,icc))
                y1_k(1:n,lay) = y1_k(1:n,lay) - x_k(row+1:row+n) * opdep(lev)
                y0_k(1:n,lay) = y0_k(1:n,lay) - x_k(row+1:row+n)
              ENDIF

              DO j = n, 1, -1
                lo = klu1 + (row + 1) - (col + j)
                tlayer_k(j,lay) = tlayer_k(j,lay) + &
                    SUM(kkpack_k(lo:lo+n-1,col+j) * (fp1(:,j,laymap,icci,m) - surfterm(j)))
                surfterm_k(j) = surfterm_k(j) - SUM(kkpack_k(lo:lo+n-1,col+j)) * tlayer(j,lay)
                fp1_k(:,j,laymap,icci,m) = fp1_k(:,j,laymap,icci,m) + &
                    kkpack_k(lo:lo+n-1,col+j) * tlayer(j,lay)

                lo = lo - n
                fp2_k(:,j,laymap,icci,m) = fp2_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+n+j)
                surfterm_k(n+j) = surfterm_k(n+j) - SUM(kkpack_k(lo:lo+n-1,col+n+j))
              ENDDO

              IF (reflfac > 0._jprb) THEN
                DO j = n, 1, -1
                  fm2_k(:,j,laymap,icci,m) = fm2_k(:,j,laymap,icci,m) + surfterm_k(n+j) * reflfac * muqw(:)
                  fm1_k(:,j,laymap,icci,m) = fm1_k(:,j,laymap,icci,m) + surfterm_k(j) * reflfac * muqw(:)

                  reflfac_k = reflfac_k + surfterm_k(n+j) * surfterm(n+j) / reflfac
                  reflfac_k = reflfac_k + surfterm_k(j) * surfterm(j) / reflfac
                ENDDO
              ELSE
                surfterm_k = 0._jprb
              ENDIF
            ELSE
              IF (dosolar) THEN
                tausun_k(lev) = tausun_k(lev) - SUM(x_k(row+1:row+n) * dom_state(i)%z(1:n,laymap,icci,m))
                z_k(1:n,laymap,icci,m) = z_k(1:n,laymap,icci,m) - x_k(row+1:row+n) * tausun(lev)
              ELSE
                opdep_k(lev) = opdep_k(lev) - SUM(x_k(row+1:row+n) * dom_state(i)%y1(1:n,lay,icc))
                y1_k(1:n,lay) = y1_k(1:n,lay) - x_k(row+1:row+n) * opdep(lev)
                y0_k(1:n,lay) = y0_k(1:n,lay) - x_k(row+1:row+n)
              ENDIF
              DO j = n, 1, -1
                lo = klu1 + (row + 1) - (col + j)
                tlayer_k(j,lay) = tlayer_k(j,lay) + SUM(kkpack_k(lo:lo+n-1,col+j) * fp1(:,j,laymap,icci,m))
                fp1_k(:,j,laymap,icci,m) = fp1_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+j) * tlayer(j,lay)

                lo = lo - n
                fp2_k(:,j,laymap,icci,m) = fp2_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,col+n+j)
              ENDDO
            ENDIF
          ENDIF

          IF (lay == 1) THEN
            row = 0
            col = 0
            IF (dosolar) THEN
              z_k(n+1:nstr,laymap,icci,m) = z_k(n+1:nstr,laymap,icci,m) - x_k(1:n)
            ELSE
              y0_k(n+1:nstr,lay) = y0_k(n+1:nstr,lay) - x_k(1:n)
              auxrad_k%air(1,i) = auxrad_k%air(1,i) + SUM(x_k(1:n))
            ENDIF
            DO j = n, 1, -1
              lo = klu1 + (row + 1) - (col + j)
              fm1_k(:,j,laymap,icci,m) = fm1_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,j)

              lo = lo - n
              tlayer_k(j,lay) = tlayer_k(j,lay) + SUM(kkpack_k(lo:lo+n-1,n+j) * fm2(:,j,laymap,icci,m))
              fm2_k(:,j,laymap,icci,m) = fm2_k(:,j,laymap,icci,m) + kkpack_k(lo:lo+n-1,n+j) * tlayer(j,lay)
            ENDDO
          ELSE
            row = (lay-2) * nstr + n
            col = (lay-1) * nstr
            IF (dosolar) THEN
              tausun_k(lev-1) = tausun_k(lev-1) + SUM(x_k(row+1:row+nstr) * dom_state(i)%z(:,laymap,icci,m))
              z_k(:,laymap,icci,m) = z_k(:,laymap,icci,m) + x_k(row+1:row+nstr) * tausun(lev-1)
            ELSE
              opdep_k(lev-1) = opdep_k(lev-1) + SUM(x_k(row+1:row+nstr) * dom_state(i)%y1(:,lay,icc))
              y1_k(:,lay) = y1_k(:,lay) + x_k(row+1:row+nstr) * opdep(lev-1)
              y0_k(:,lay) = y0_k(:,lay) + x_k(row+1:row+nstr)
            ENDIF
            DO j = n, 1, -1
              lo = klu1 + (row + 1) - (col + j)
              fp1_k(:,j,laymap,icci,m) = fp1_k(:,j,laymap,icci,m) - kkpack_k(lo:lo+n-1,col+j)

              lo = lo - n
              tlayer_k(j,lay) = tlayer_k(j,lay) - &
                  SUM(kkpack_k(lo:lo+n-1,col+n+j) * fp2(:,j,laymap,icci,m))
              fp2_k(:,j,laymap,icci,m) = fp2_k(:,j,laymap,icci,m) - &
                  kkpack_k(lo:lo+n-1,col+n+j) * tlayer(j,lay)

              lo = klu1 + (row + n + 1) - (col + j)
              fm1_k(:,j,laymap,icci,m) = fm1_k(:,j,laymap,icci,m) - kkpack_k(lo:lo+n-1,col+j)

              lo = lo - n
              tlayer_k(j,lay) = tlayer_k(j,lay) - &
                  SUM(kkpack_k(lo:lo+n-1,col+n+j) * fm2(:,j,laymap,icci,m))
              fm2_k(:,j,laymap,icci,m) = fm2_k(:,j,laymap,icci,m) - &
                  kkpack_k(lo:lo+n-1,col+n+j) * tlayer(j,lay)
            ENDDO
          ENDIF


          ! --------------------------------------------------------------------------
          ! Calculate layer "transmittances" with opdeps scaled by eigenvalues
          ! --------------------------------------------------------------------------
          eval_k(:,laymap,icci,m) = eval_k(:,laymap,icci,m) - &
              tlayer_k(:,lay) * tlayer(:,lay) * ftau(lay) * profiles_dom(icc,i)%layerod(lay)
          tmpval = SUM(tlayer_k(:,lay) * tlayer(:,lay) * dom_state(i)%eval(:,laymap,icci,m))
          profiles_dom_k(icc,i)%layerod(lay) = profiles_dom_k(icc,i)%layerod(lay) - tmpval * ftau(lay)
          ftau_k(lay) = ftau_k(lay) - tmpval * profiles_dom(icc,i)%layerod(lay)


          ! --------------------------------------------------------------------------
          ! Solve eigenvalue problem and obtain homogenous solution for this layer
          ! --------------------------------------------------------------------------
          ! Note that for clear layers we already wrote the solutions down above

          IF (fssa(lay) > 0._jprb) THEN

            ! Particular integral - thermal source term, done once per layer per cloud column
            IF (dothermal) THEN
              b1_k(lay) = b1_k(lay) + fpmupi1_k(lay)
              phasesat_k(:) = phasesat_k(:) + 0.5_jprb * fpmupi1_k(lay) * dom_state(i)%y1(:,lay,icc)
              y1_k(:,lay) = y1_k(:,lay) + 0.5_jprb * fpmupi1_k(lay) * phasesat(:)

              b0_k(lay) = b0_k(lay) + fpmupi0_k(lay)
              phasesat_k(:) = phasesat_k(:) + 0.5_jprb * fpmupi0_k(lay) * dom_state(i)%y0(:,lay,icc)
              y0_k(:,lay) = y0_k(:,lay) + 0.5_jprb * fpmupi0_k(lay) * phasesat(:)
            ENDIF


            ! --------------------------------------------------------------------------
            ! Interpolate HS and PI for final radiance calculation
            ! --------------------------------------------------------------------------
!             IF (.NOT. done_lay_str_azi_k(laymap,icci,m)) THEN

              ! Particular integral - solar source term, done once per layer (applies to all cloud columns)
              IF (dosolar) THEN
                legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                    sunsrcfac * fpmupi_k(laymap,icci,m) * legpolymusat(m:nmom,0,m) * legpolymusun(m:nmom,0,m)
                phasesat_k(:) = phasesat_k(:) + &
                    0.5_jprb * fpmupi_k(laymap,icci,m) * dom_state(i)%z(:,laymap,icci,m)
                z_k(:,laymap,icci,m) = z_k(:,laymap,icci,m) + &
                    0.5_jprb * fpmupi_k(laymap,icci,m) * phasesat(:)
              ENDIF

              ! Homogenous part - diffuse radiation, done once per layer (applies to all cloud columns)
              DO j = n, 1, -1
                phasesat_k(n+1:nstr) = phasesat_k(n+1:nstr) + &
                    0.5_jprb * fpmuhsdn_k(j,laymap,icci,m) * fm2(:,j,laymap,icci,m)
                fm2_k(:,j,laymap,icci,m) = fm2_k(:,j,laymap,icci,m) + &
                    0.5_jprb * fpmuhsdn_k(j,laymap,icci,m) * phasesat(n+1:nstr)
                phasesat_k(1:n) = phasesat_k(1:n) + &
                    0.5_jprb * fpmuhsdn_k(j,laymap,icci,m) * fp2(:,j,laymap,icci,m)
                fp2_k(:,j,laymap,icci,m) = fp2_k(:,j,laymap,icci,m) + &
                    0.5_jprb * fpmuhsdn_k(j,laymap,icci,m) * phasesat(1:n)
                phasesat_k(n+1:nstr) = phasesat_k(n+1:nstr) + &
                    0.5_jprb * fpmuhsup_k(j,laymap,icci,m) * fm1(:,j,laymap,icci,m)
                fm1_k(:,j,laymap,icci,m) = fm1_k(:,j,laymap,icci,m) + &
                    0.5_jprb * fpmuhsup_k(j,laymap,icci,m) * phasesat(n+1:nstr)
                phasesat_k(1:n) = phasesat_k(1:n) + &
                    0.5_jprb * fpmuhsup_k(j,laymap,icci,m) * fp1(:,j,laymap,icci,m)
                fp1_k(:,j,laymap,icci,m) = fp1_k(:,j,laymap,icci,m) + &
                    0.5_jprb * fpmuhsup_k(j,laymap,icci,m) * phasesat(1:n)
              ENDDO

!             ENDIF ! not done_lay_str_azi_k


            ! Thermal source term done for every layer and *for every cloud column*
            IF (dothermal) THEN
              ! --------------------------------------------------------------------------
              ! Solve for particular integral for thermal source term
              ! --------------------------------------------------------------------------

              ! Construct and solve the linear system for scattering layers
              DO j = 1, n
                hh(1:n,j)      = -0.5_jprb * muw(j) * phase_fn(:,j)
                hh(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry

                DO k = 1, n
                  hh(k,n+j)   = hh(n+k,j)
                  hh(n+k,n+j) = hh(k,j)
                ENDDO

                hh(j,j)     = hh(j,j) + 1._jprb
                hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb
              ENDDO

!               hhinvt = TRANSPOSE(matinv(hh, err))
!               IF (err /= 0) THEN
!                 PRINT *,'Matrix inversion error ', err
!                 THROWM(err.NE.0,'DOM error: solving for K particular integral - thermal')
!               ENDIF
!               ! No accumulation required
!               hh_k = - MATMUL(hhinvt, &
!                   MATMUL(y0_k(:,lay:lay), TRANSPOSE(dom_state(i)%y0(:,lay:lay,icc))))

!               y0_k(:,lay) = MATMUL(hhinvt, y0_k(:,lay))

              ! Use DGETRF for LU factorisation (of hh) and then two calls to DGETRS to solve for y1 and y0
              CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving K PI thermal, DGETRF error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              hh_k = -MATMUL(y0_k(:,lay:lay), TRANSPOSE(dom_state(i)%y0(:,lay:lay,icc)))
              CALL DGETRS('Transpose', nstri, nstri, hh, nstri, ipiv2, hh_k, nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving K PI thermal y1, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              CALL DGETRS('Transpose', nstri, 1, hh, nstri, ipiv2, y0_k(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving K PI thermal y0, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              y1_k(n+1:nstr,lay) = y1_k(n+1:nstr,lay) - muq(1:n) * y0_k(n+1:nstr,lay)
              y1_k(1:n,lay) = y1_k(1:n,lay) + muq(1:n) * y0_k(1:n,lay)
              b0_k(lay) = b0_k(lay) + SUM(y0_k(1:nstr,lay))

!               hh_k = hh_k - MATMUL(hhinvt, &
!                   MATMUL(y1_k(:,lay:lay), TRANSPOSE(dom_state(i)%y1(:,lay:lay,icc))))

!               y1_k(:,lay) = MATMUL(hhinvt, y1_k(:,lay))

              hhtmp_k = -MATMUL(y1_k(:,lay:lay), TRANSPOSE(dom_state(i)%y1(:,lay:lay,icc)))
              CALL DGETRS('Transpose', nstri, nstri, hh, nstri, ipiv2, hhtmp_k, nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving K PI thermal y1, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF
              hh_k = hh_k + hhtmp_k

              CALL DGETRS('Transpose', nstri, 1, hh, nstri, ipiv2, y1_k(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving K PI thermal y0, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              ! Construct and solve the linear system for scattering layers
              DO j = n, 1, -1
                DO k = n, 1, -1
                  hh_k(k,j) = hh_k(k,j) + hh_k(n+k,n+j)
                  hh_k(n+k,j) = hh_k(n+k,j) + hh_k(k,n+j)
                ENDDO

                phase_fn_k(:,n+j) = phase_fn_k(:,n+j) - 0.5_jprb * muw(j) * hh_k(n+1:nstr,j)
                phase_fn_k(:,j) = phase_fn_k(:,j) - 0.5_jprb * muw(j) * hh_k(1:n,j)
              ENDDO

            ENDIF ! dothermal


!             IF (.NOT. done_lay_str_azi_k(laymap,icci,m)) THEN

              ! Homogenous solution only done once per layer for scattering layers (applies to all cloud columns)

              ! --------------------------------------------------------------------------
              ! Solve eigenvalue problem and obtain homogenous solution for this layer
              ! --------------------------------------------------------------------------

              ! Solar source term particular integral only done once per layer for scattering layers (applies to all cloud columns)
              IF (dosolar) THEN
                ! --------------------------------------------------------------------------
                ! Solve for particular integral for solar source term
                ! --------------------------------------------------------------------------

                ! Solve PI problem using nstr x nstr system (most straightforward approach)

                ! Construct and solve the linear system for scattering layers
                DO j = 1, n
                  hh(1:n,j)      = -0.5_jprb * muw(j) * phase_fn(:,j)
                  hh(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry

                  DO k = 1, n
                    hh(k,n+j)   = hh(n+k,j)
                    hh(n+k,n+j) = hh(k,j)
                  ENDDO

                  hh(j,j)     = hh(j,j) + 1._jprb + muq(j) * musunr
                  hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb - muq(j) * musunr
                ENDDO

!                 ! The first step is tricky to do efficiently, but we can make an initial stab at solving directly:
!                 hhinvt = TRANSPOSE(matinv(hh, err))
!                 IF (err /= 0) THEN
!                   PRINT *,'Matrix inversion error ', err
!                   THROWM(err.NE.0,'DOM error: solving for K particular integral - solar')
!                 ENDIF
!                 ! No accumulation required
!                 hh_k = - MATMUL(hhinvt, &
!                     MATMUL(z_k(:,laymap:laymap,icci,m), TRANSPOSE(dom_state(i)%z(:,laymap:laymap,icci,m))))
! 
!                 z_k(:,laymap,icci,m) = MATMUL(hhinvt, z_k(:,laymap,icci,m))

                CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
                IF (info /= 0) THEN
                  WRITE(msg,'(a,3i5)') 'DOM error: solving K PI solar, DGETRF error, info = ', info, i, lay
                  err = errorstatus_fatal
                  THROWM(err.NE.0,msg)
                ENDIF

                hh_k = -MATMUL(z_k(:,laymap:laymap,icci,m), TRANSPOSE(dom_state(i)%z(:,laymap:laymap,icci,m)))
                CALL DGETRS('Transpose', nstri, nstri, hh, nstri, ipiv2, hh_k, nstri, info)
                IF (info /= 0) THEN
                  WRITE(msg,'(a,3i5)') 'DOM error: solving K PI solar, DGETRS error, info = ', info, i, lay
                  err = errorstatus_fatal
                  THROWM(err.NE.0,msg)
                ENDIF

                CALL DGETRS('Transpose', nstri, 1, hh, nstri, ipiv2, z_k(:,laymap,icci,m), nstri, info)
                IF (info /= 0) THEN
                  WRITE(msg,'(a,3i5)') 'DOM error: solving K PI solar, DGETRS error, info = ', info, i, lay
                  err = errorstatus_fatal
                  THROWM(err.NE.0,msg)
                ENDIF

                DO j = n, 1, -1
                  legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                      sunsrcfac * z_k(n+j,laymap,icci,m) * legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m)
                  legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                      sunsrcfac * z_k(j,laymap,icci,m) * legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m)

                  DO k = n, 1, -1
                    hh_k(k,j) = hh_k(k,j) + hh_k(n+k,n+j)
                    hh_k(n+k,j) = hh_k(n+k,j) + hh_k(k,n+j)
                  ENDDO

                  phase_fn_k(:,n+j) = phase_fn_k(:,n+j) - 0.5_jprb * muw(j) * hh_k(n+1:nstr,j)
                  phase_fn_k(:,j) = phase_fn_k(:,j) - 0.5_jprb * muw(j) * hh_k(1:n,j)
                ENDDO

              ENDIF ! dosolar

              DO j = 1, n
                bp(:,j) = 0.5_jprb * muw(j) * phase_fn(:,j) * muqr(:)
                bm(:,j) = 0.5_jprb * muw(j) * phase_fn(:,n+j) * muqr(:)
                bp(j,j) = bp(j,j) - muqr(j)
              ENDDO
              sm = bp - bm
              sp = bp + bm
!               bpm        = MATMUL(sm, sp)
!               bpmsave    = bpm  ! Used for solar order-n PI below

              xp = dom_state(i)%xp(:,:,laymap,icci,m)
! ** Recompute xm or store above? **
              DO j = 1, n
                xm(:,j) = MATMUL(sp, xp(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
              ENDDO

              ! First use of xp_k, xm_k:
              xm_k = 0.5_jprb * (fp1_k(:,:,laymap,icci,m) - fm1_k(:,:,laymap,icci,m))
              xp_k = 0.5_jprb * (fp1_k(:,:,laymap,icci,m) + fm1_k(:,:,laymap,icci,m))

              sp_k = 0._jprb
              DO j = n, 1, -1
                eval_k(j,laymap,icci,m) = eval_k(j,laymap,icci,m) - &
                    SUM(xm_k(:,j) * xm(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
                sp_k = sp_k + MATMUL(xm_k(:,j:j), TRANSPOSE(xp(:,j:j))) / dom_state(i)%eval(j,laymap,icci,m)
                xp_k(:,j) = xp_k(:,j) + MATMUL(TRANSPOSE(sp), xm_k(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
              ENDDO

              ! Calculate the required eigenvalues (+/- each one gives nstr evals)
              eval_k(:,laymap,icci,m) = 0.5_jprb * eval_k(:,laymap,icci,m) / dom_state(i)%eval(:,laymap,icci,m)

              ! The eigenvalues output by ASYMTX are required here (before we take the square root)
!               bpm_k = 0._jprb ! Not accumulated in ASYMTX_AD
              eval = dom_state(i)%eval(:,laymap,icci,m)**2
              CALL ASYMTX_AD(err, n, xp, eval, xp_k, eval_k(:,laymap,icci,m), bpm_k)
              IF (err /= 0) THEN
                WRITE(msg,'(a,i5)') 'DOM error: solving K eigenvalue problem, ASYMTX_AD error = ', err
                THROWM(err.NE.0,msg)
              ENDIF

!               bpm_k = bpm_k + bpmsave_k ! AD of order n PI not implemented (yet)
              ! First use of sm_k, bp_k, bm_k
              sm_k = MATMUL(bpm_k, TRANSPOSE(sp))
              sp_k = sp_k + MATMUL(TRANSPOSE(sm), bpm_k)
              bp_k = sp_k + sm_k
              bm_k = sp_k - sm_k
              DO j = n, 1, -1
                phase_fn_k(:,n+j) = phase_fn_k(:,n+j) + 0.5_jprb * muw(j) * bm_k(:,j) * muqr(:)
                phase_fn_k(:,j) = phase_fn_k(:,j) + 0.5_jprb * muw(j) * bp_k(:,j) * muqr(:)
              ENDDO

!             ENDIF ! not done_lay_str_azi_k


            ! --------------------------------------------------------------------------
            ! For scattering layers, only do calculation if we don't already have the results
            ! --------------------------------------------------------------------------
!             IF (dothermal .OR. .NOT. done_lay_str_azi_k(laymap,icci,m)) THEN

              DO j = nstr, n+1, -1
                legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                    phasesat_k(j) * muw(j-n) * legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m)
              ENDDO
              DO j = n, 1, -1
                legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                    phasesat_k(j) * muw(j) * legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m)
              ENDDO

              ! Due to symmetry we can calculate just n*nstr values (not nstr*nstr):
              !    phase_fn(i,j) = phase_fn(n+i,n+j)    where n+x => -ve mu
              !    phase_fn(n+i,j) = phase_fn(i,n+j)
              DO j = n, 1, -1
                DO k = n, j, -1

                  ! Exploit symmetry to avoid repeated calculations
                  IF (k > j) THEN
                    phase_fn_k(k,n+j) = phase_fn_k(k,n+j) + phase_fn_k(j,n+k)
                    phase_fn_k(k,j) = phase_fn_k(k,j) + phase_fn_k(j,k)
                  ENDIF

                  legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                      phase_fn_k(k,n+j) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,n+j,m)
                  legcoef_k(m:nmom,laymap,icci) = legcoef_k(m:nmom,laymap,icci) + &
                      phase_fn_k(k,j) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,j,m)
                ENDDO
              ENDDO

              ! --------------------------------------------------------------------------
              ! Calculate the phase function values for each pair of quadrature angles and musun
              ! --------------------------------------------------------------------------
              ! Delta-scale Legendre coefficients for this layer
              ! NB delta-scaled SSA is incorporated into Leg. coefs
              IF (m == 0) THEN  ! This calculation applies to all azi components
                DO l = nmom, 0, -1
                  f_lay_k = f_lay_k + legcoef_k(l,laymap,icci) * &
                      fssa(lay) * (trans_scatt_ir_dyn%phasefn(i)%legcoef(l,icci,laymap) - &
                      (2._jprb * l + 1._jprb) * f_lay) / (1._jprb - f_lay)**2
                  f_lay_k = f_lay_k - legcoef_k(l,laymap,icci) * &
                      fssa(lay) * (2._jprb * l + 1._jprb) / (1._jprb - f_lay)
                  trans_scatt_ir_dyn_k%phasefn(i)%legcoef(l,icci,laymap) = &
                      trans_scatt_ir_dyn_k%phasefn(i)%legcoef(l,icci,laymap) + &
                      legcoef_k(l,laymap,icci) * fssa(lay) / (1._jprb - f_lay)
                  fssa_k(lay) = fssa_k(lay) + legcoef_k(l,laymap,icci) * &
                      (trans_scatt_ir_dyn%phasefn(i)%legcoef(l,icci,laymap) - &
                      (2._jprb * l + 1._jprb) * f_lay) / (1._jprb - f_lay)
                ENDDO
              ENDIF

!             ENDIF ! dothermal or not done_lay_str_azi_k

!             done_lay_str_azi_k(laymap,icci,m) = .TRUE.

          ENDIF ! fssa > 0 (i.e. scattering layer)


          ! --------------------------------------------------------------------------
          ! Calculate emissive source term and interpolated clear-sky emissive source term
          ! --------------------------------------------------------------------------
          IF (dothermal) THEN
            ! For clear layers we have the solution and the interpolated thermal source term is trivial
            IF (.NOT. fssa(lay) > 0._jprb) THEN
              b1_k(lay) = b1_k(lay) + fpmupi1_k(lay)
              b0_k(lay) = b0_k(lay) + fpmupi0_k(lay)
              y1_k(n+1:nstr,lay) = y1_k(n+1:nstr,lay) - y0_k(n+1:nstr,lay) * muq(1:n)
              y1_k(1:n,lay) = y1_k(1:n,lay) + y0_k(1:n,lay) * muq(1:n)
              b0_k(lay) = b0_k(lay) + SUM(y0_k(1:nstr,lay))
            ENDIF

            b1_k(lay) = b1_k(lay) + SUM(y1_k(:,lay))

            opdep_k(lev-1) = opdep_k(lev-1) - b0_k(lay) * b1(lay)
            b1_k(lay) = b1_k(lay) - b0_k(lay) * opdep(lev-1)
            auxrad_k%air(lev-1,i) = auxrad_k%air(lev-1,i) + b0_k(lay) * (1._jprb - fssa(lay))
            fssa_k(lay) = fssa_k(lay) - b0_k(lay) * auxrad%air(lev-1,i)

            tmpval = 1._jprb / (opdep(lev-1) - opdep(lev))
            IF (lay == nlayers .AND. profiles_dom(icc,i)%surface) THEN
              opdep_k(lev) = opdep_k(lev) + b1_k(lay) * &
                  (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) * tmpval**2
              opdep_k(lev-1) = opdep_k(lev-1) - b1_k(lay) * &
                  (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) * tmpval**2
              auxrad_k%surfair(i) = auxrad_k%surfair(i) - b1_k(lay) * &
                  (1._jprb - fssa(lay)) * tmpval
              auxrad_k%air(lev-1,i) = auxrad_k%air(lev-1,i) + b1_k(lay) * &
                  (1._jprb - fssa(lay)) * tmpval
              fssa_k(lay) = fssa_k(lay) - b1_k(lay) * &
                  (auxrad%air(lev-1,i) - auxrad%surfair(i)) * tmpval
            ELSE
              opdep_k(lev) = opdep_k(lev) + b1_k(lay) * &
                  (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) * tmpval**2
              opdep_k(lev-1) = opdep_k(lev-1) - b1_k(lay) * &
                  (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) * tmpval**2
              auxrad_k%air(lev,i) = auxrad_k%air(lev,i) - b1_k(lay) * &
                  (1._jprb - fssa(lay)) * tmpval
              auxrad_k%air(lev-1,i) = auxrad_k%air(lev-1,i) + b1_k(lay) * &
                  (1._jprb - fssa(lay)) * tmpval
              fssa_k(lay) = fssa_k(lay) - b1_k(lay) * &
                  (auxrad%air(lev-1,i) - auxrad%air(lev,i)) * tmpval
            ENDIF

          ENDIF


          ! --------------------------------------------------------------------------
          ! Calculate optical depths and transmittances
          ! --------------------------------------------------------------------------
          ! Calculate accumulated delta-scaled optical depth and transmittances on
          ! satellite-surface and sun-surface paths.
          IF (m == 0) THEN  ! This calculation applies to all azi components
            IF (dosolar) opdep_k(lev) = opdep_k(lev) - tausun_k(lev) * musunr * tausun(lev)
            opdep_k(lev) = opdep_k(lev) - tausat_k(lev) * musatr * tausat(lev)
            profiles_dom_k(icc,i)%layerod(lay) = profiles_dom_k(icc,i)%layerod(lay) + &
                opdep_k(lev) * ftau(lay)
            ftau_k(lay) = ftau_k(lay) + opdep_k(lev) * profiles_dom(icc,i)%layerod(lay)
            opdep_k(lev-1) = opdep_k(lev-1) + opdep_k(lev)
          ENDIF


          ! --------------------------------------------------------------------------
          ! Set delta-scaling parameters (same for all azimuthal components)
          ! --------------------------------------------------------------------------
          IF (m == 0) THEN
            IF (laymap > 0) THEN
              ftau_k(lay) = ftau_k(lay) - fssa_k(lay) * (1._jprb - f_lay) * ssa(icci,laymap,i) / ftau(lay)**2
              ssa_k(icci,laymap,i) = ssa_k(icci,laymap,i) + fssa_k(lay) * (1._jprb - f_lay) / ftau(lay)
              f_lay_k = f_lay_k - fssa_k(lay) * ssa(icci,laymap,i) / ftau(lay)

              f_lay_k = f_lay_k - ftau_k(lay) * ssa(icci,laymap,i)
              ssa_k(icci,laymap,i) = ssa_k(icci,laymap,i) - ftau_k(lay) * f_lay
              trans_scatt_ir_dyn_k%phasefn(i)%legcoef(nstr,icci,laymap) = &
                  trans_scatt_ir_dyn_k%phasefn(i)%legcoef(nstr,icci,laymap) + f_lay_k / (2._jprb * nstr + 1._jprb)
!             ELSE
!               f_lay_k     = 0._jprb
!               ftau_k(lay) = 0._jprb
!               fssa_k(lay) = 0._jprb
            ENDIF
          ENDIF

        ENDDO ! layer loop

      ENDDO ! azimuth loop

      ! Surface reflectance factor
      IF (profiles_dom(icc,i)%surface) THEN
        IF (dosolar) THEN
          albedo_k = albedo_k + 2._jprb * reflfac_k
          IF (diffuse_refl(i) > 1._jprb) THEN
            albedo_k = 0._jprb
          ELSE
            reflectance_k(i)%refl_out = reflectance_k(i)%refl_out + albedo_k * pi
          ENDIF
        ELSE
          diffuse_refl_k(i) = diffuse_refl_k(i) + 2._jprb * reflfac_k
        ENDIF
      ENDIF

    ENDDO ! cloud columns

  ENDDO ! channel


  ! --------------------------------------------------------------------------
  ! Tidy up
  ! --------------------------------------------------------------------------
  NULLIFY(fp2, fm2, fp2_k, fm2_k)

  IF (dosolar) THEN
    DEALLOCATE(legpolymusun, legpolymuscata, &
               tausun, fpmupi, &
               tausun_k, fpmupi_k, z_k, stat=err)
  ELSE
    DEALLOCATE(b0, b1, fpmupi0, fpmupi1, &
               b0_k, b1_k, fpmupi0_k, fpmupi1_k, y0_k, y1_k, stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM deallocation error')

  DEALLOCATE(legpolymuq, legpolymusat, fpmuhsup, fpmuhsdn, fp1, fm1, &
             eval_k, fpmuhsup_k, fpmuhsdn_k, fp1_k, fm1_k, &
             done_lay_str_azi, stat=err) !done_lay_str_azi_k
  THROWM(err.NE.0, 'DOM deallocation error')

CATCH

END SUBROUTINE rttov_dom_k
