! Description:
!> @file
!!   TL of RTE integration using Discrete Ordinates Method.
!
!> @brief
!!   TL of RTE integration using Discrete Ordinates Method.
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
!! @param[in]     profiles_dom_tl        DOM layer optical depth profile perturbations
!! @param[in]     maxnlayers             maximum number of layers in the DOM profiles
!! @param[in]     auxrad                 Planck radiances for emission source term
!! @param[in]     auxrad_tl              Planck radiance perturbations
!! @param[in]     trans_scatt_ir         single-scattering albedos
!! @param[in]     trans_scatt_ir_tl      single-scattering albedo perturbations
!! @param[in]     trans_scatt_ir_dyn     phase function Legendre decompositions
!! @param[in]     trans_scatt_ir_dyn_tl  phase function Legendre decomposition perturbations
!! @param[in]     emissivity             surface emissivities
!! @param[in]     emissivity_tl          surface emissivity perturbations
!! @param[in]     reflectance            surface BRDFs
!! @param[in]     reflectance_tl         surface BRDF perturbations
!! @param[in]     diffuse_refl           surface albedos
!! @param[in]     diffuse_refl_tl        surface albedo perturbations
!! @param[in]     solar_spectrum         top of atmosphere solar irradiance
!! @param[in]     raytracing             raytracing structure
!! @param[in]     ircld                  cloud column data
!! @param[in]     ircld_tl               cloud column data perturbations
!! @param[in,out] rad_tl                 calculated radiance perturbations
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
SUBROUTINE rttov_dom_tl(             &
              err,                   &
              chanprof,              &
              dosolar,               &
              chanflag,              &
              nstr,                  &
              profiles,              &
              profiles_dom,          &
              profiles_dom_tl,       &
              maxnlayers,            &
              auxrad,                &
              auxrad_tl,             &
              trans_scatt_ir,        &
              trans_scatt_ir_tl,     &
              trans_scatt_ir_dyn,    &
              trans_scatt_ir_dyn_tl, &
              emissivity,            &
              emissivity_tl,         &
              reflectance,           &
              reflectance_tl,        &
              diffuse_refl,          &
              diffuse_refl_tl,       &
              solar_spectrum,        &
              raytracing,            &
              ircld,                 &
              ircld_tl,              &
              rad_tl,                &
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
      asymtx_tl

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
  TYPE(rttov_profile_dom),           INTENT(IN)    :: profiles_dom_tl(0:,:)
  INTEGER(jpim),                     INTENT(IN)    :: maxnlayers
  TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad
  TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_dyn_tl
  TYPE(rttov_emissivity),            INTENT(IN)    :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity),            INTENT(IN)    :: emissivity_tl(SIZE(chanprof))
  TYPE(rttov_reflectance),           INTENT(IN)    :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),           INTENT(IN)    :: reflectance_tl(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: diffuse_refl(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: diffuse_refl_tl(SIZE(chanprof))
  REAL(jprb),                        INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld_tl
  TYPE(rttov_radiance),              INTENT(INOUT) :: rad_tl
  TYPE(rttov_dom_state),             INTENT(IN)    :: dom_state(SIZE(chanprof))
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=128) :: msg
  LOGICAL(jplm) :: dothermal, inc_surface!, do_order_n_pi

  INTEGER(jpim) :: n, maxnaz, thisnaz, nmom

  INTEGER(jpim) :: i, j, k, l, m, row, col
  INTEGER(jpim) :: icc, icci
  INTEGER(jpim) :: nlayers, nlevels, lay, lev, laymap
  INTEGER(jpim) :: nchanprof, prof, lastprof
  INTEGER(jpim) :: ilo, ihi, lo, hi, klu1

  INTEGER       :: kl, ku, klku, info, nstri, nstrlay ! LAPACK arguments, no KIND

  REAL(jprb)    :: musat, musatr, musun, musunr, relazi, muscata
  REAL(jprb)    :: taufac, albedo, reflfac, sunsrcfac
  REAL(jprb)    :: taufac_tl, delta_tms_tl, albedo_tl, reflfac_tl
  REAL(jprb)    :: radsurfup_sum_tl, radsurfup_tl, radinc_tl, thisrad_tl
  REAL(jprb)    :: tmpval, tmpval2

  REAL(jprb)              :: muq(nstr/2), muw(nstr/2), muqw(nstr/2), muqr(nstr/2)
  REAL(jprb), ALLOCATABLE :: legpolymuq(:,:,:) !(0:nstr-1,nstr,0:nstr-1)     ! dims are (0:nmom,nstr,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymusat(:,:,:) !(0:nstr-1,0:0,0:nstr-1)    ! dims are (0:nmom,0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymusun(:,:,:) !(0:nstr-1,0:0,0:nstr-1)    ! dims are (0:nmom,0,0:maxnaz)
  REAL(jprb), ALLOCATABLE :: legpolymuscata(:,:,:) !(0:nstr-1,0:0,0:0)
  REAL(jprb)              :: legcoef(0:nstr-1,profiles(1)%nlayers,0:1)
  REAL(jprb)              :: legcoef_tl(0:nstr-1,profiles(1)%nlayers,0:1)

  REAL(jprb)              :: f_lay, ftau(maxnlayers), fssa(maxnlayers)
  REAL(jprb)              :: f_lay_tl, ftau_tl(maxnlayers), fssa_tl(maxnlayers)
  REAL(jprb)              :: phase_fn(nstr/2,nstr)
  REAL(jprb)              :: phase_fn_tl(nstr/2,nstr)
  REAL(jprb), ALLOCATABLE :: b0(:), b1(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: b0_tl(:), b1_tl(:) !(maxnlayers)
  REAL(jprb)              :: opdep(maxnlayers+1)
  REAL(jprb)              :: opdep_tl(maxnlayers+1)
  REAL(jprb)              :: tausat(maxnlayers+1)
  REAL(jprb)              :: tausat_tl(maxnlayers+1)
  REAL(jprb), ALLOCATABLE :: tausun(:) !(maxnlayers+1)
  REAL(jprb), ALLOCATABLE :: tausun_tl(:) !(maxnlayers+1)
  REAL(jprb)              :: tlayer(nstr/2,maxnlayers)
  REAL(jprb)              :: tlayer_tl(nstr/2,maxnlayers)
  REAL(jprb)              :: surfterm(nstr)
  REAL(jprb)              :: surfterm_tl(nstr)
  REAL(jprb)              :: phasesat(nstr)
  REAL(jprb)              :: phasesat_tl(nstr)

  REAL(jprb)              :: bp(nstr/2,nstr/2), bm(nstr/2,nstr/2)!, bpm(nstr/2,nstr/2)
  REAL(jprb)              :: bp_tl(nstr/2,nstr/2), bm_tl(nstr/2,nstr/2), bpm_tl(nstr/2,nstr/2)
  REAL(jprb)              :: xp(nstr/2,nstr/2), xm(nstr/2,nstr/2)
  REAL(jprb)              :: xp_tl(nstr/2,nstr/2), xm_tl(nstr/2,nstr/2)
  REAL(jprb)              :: sp(nstr/2,nstr/2), sm(nstr/2,nstr/2)
  REAL(jprb)              :: sp_tl(nstr/2,nstr/2), sm_tl(nstr/2,nstr/2)!, bpmsave_tl(nstr/2,nstr/2)
!   REAL(jprb)              :: so1, so2, sop(nstr/2), som(nstr/2)
!   REAL(jprb)              :: so1_tl, so2_tl, sop_tl(nstr/2), som_tl(nstr/2)
  REAL(jprb)              :: kkpack_tl(9*nstr/2-2,nstr*maxnlayers)

  ! LAPACK arguments, no KIND
  INTEGER                 :: ipiv1(nstr*maxnlayers), ipiv2(nstr)
  DOUBLE PRECISION        :: hh(nstr,nstr)
  DOUBLE PRECISION        :: hh_tl(nstr,nstr), x_tl(nstr*maxnlayers)
  DOUBLE PRECISION        :: kkpack(9*nstr/2-2,nstr*maxnlayers)
!   DOUBLE PRECISION        :: bpmsave(nstr/2,nstr/2), gp(nstr/2), gp_tl(nstr/2)
  DOUBLE PRECISION, ALLOCATABLE :: z_tl(:,:,:,:) !(nstr,0:profiles(1)%nlayers,0:1,0:nstr-1)
  DOUBLE PRECISION, ALLOCATABLE :: y0_tl(:,:) !(nstr,maxnlayers)
  DOUBLE PRECISION, ALLOCATABLE :: y1_tl(:,:) !(nstr,maxnlayers)

  REAL(jprb)              :: eval(nstr/2)
  REAL(jprb), ALLOCATABLE :: eval_tl(:,:,:,:) !(nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsup(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsup_tl(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsdn(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmuhsdn_tl(:,:,:,:) !(nstr,profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi_tl(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), ALLOCATABLE :: fpmupi0(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi0_tl(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi1(:) !(maxnlayers)
  REAL(jprb), ALLOCATABLE :: fpmupi1_tl(:) !(maxnlayers)

  REAL(jprb), TARGET, ALLOCATABLE :: fp1(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fp1_tl(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fm1(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), TARGET, ALLOCATABLE :: fm1_tl(:,:,:,:,:) !(nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1)
  REAL(jprb), POINTER :: fp2(:,:,:,:,:), fm2(:,:,:,:,:)
  REAL(jprb), POINTER :: fp2_tl(:,:,:,:,:), fm2_tl(:,:,:,:,:)
  REAL(jprb), POINTER :: ssa(:,:,:), ssa_tl(:,:,:)

  ! This is only an internal logical array so use smallest available KIND
  LOGICAL(jpit), ALLOCATABLE :: done_lay_str_azi(:,:,:) !(profiles(1)%nlayers,0:1,0:nstr-1)

  ! --------------------------------------------------------------------------------------

  TRY

!   do_order_n_pi = .FALSE.

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
           eval_tl(n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsup(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsup_tl(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsdn(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fpmuhsdn_tl(nstr,profiles(1)%nlayers,0:1,0:maxnaz), &
           fp1(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fp1_tl(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fm1(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           fm1_tl(n,n,0:profiles(1)%nlayers,0:1,0:maxnaz), &
           done_lay_str_azi(profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
  THROWM(err.NE.0, 'DOM allocation error')

  IF (dosolar) THEN
    ALLOCATE(legpolymusun(0:nstr-1,0:0,0:maxnaz), &
             legpolymuscata(0:nstr-1,0:0,0:0), &
             tausun(maxnlayers+1), &
             tausun_tl(maxnlayers+1), &
             fpmupi(profiles(1)%nlayers,0:1,0:maxnaz), &
             fpmupi_tl(profiles(1)%nlayers,0:1,0:maxnaz), &
             z_tl(nstr,0:profiles(1)%nlayers,0:1,0:maxnaz), stat=err)
  ELSE
    ALLOCATE(b0(maxnlayers), b1(maxnlayers), &
             b0_tl(maxnlayers), b1_tl(maxnlayers), &
             fpmupi0(maxnlayers), fpmupi1(maxnlayers), &
             fpmupi0_tl(maxnlayers), fpmupi1_tl(maxnlayers), &
             y0_tl(nstr,maxnlayers), y1_tl(nstr,maxnlayers), stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM allocation error')

  ! Exploit symmetry of eigenvectors using pointers
  fp2 => fm1
  fm2 => fp1
  fp2_tl => fm1_tl
  fm2_tl => fp1_tl

  IF (dosolar) THEN
    ssa => trans_scatt_ir%ssa_solar
    ssa_tl => trans_scatt_ir_tl%ssa_solar
  ELSE
    ssa => trans_scatt_ir%ssa_thermal
    ssa_tl => trans_scatt_ir_tl%ssa_thermal
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

    eval_tl(:,0,0,m)  = 0._jprb
    fp1_tl(:,:,0,0,m) = 0._jprb
    fm1_tl(:,:,0,0,m) = 0._jprb
    DO j = 1, n
      fm1_tl(j,j,0,0,m) = 0._jprb
    ENDDO
    ! See notes below regarding fp2/fm2 (i.e. eigenvectors for -ve eigenvalues)

    ! Particular integral is zero for clear layers for solar channels
    IF (dosolar) z_tl(:,0,0,m) = 0._jprb
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

    opdep_tl(1) = 0._jprb
    tausat_tl(1) = 0._jprb
    IF (dosolar) tausun_tl(1) = 0._jprb

    done_lay_str_azi = .FALSE.

    DO icc = 0, ircld%ncolumn(prof)
      thisrad_tl   = 0._jprb
      delta_tms_tl = 0._jprb

      ! Surface reflectance factor
      IF (profiles_dom(icc,i)%surface) THEN
        IF (dosolar) THEN
          ! This prevents unphysical albedos
          albedo = reflectance(i)%refl_out * pi
          IF (albedo > 1._jprb) THEN
            albedo = 1._jprb
            albedo_tl = 0._jprb
          ELSE
            albedo_tl = reflectance_tl(i)%refl_out * pi
          ENDIF
          reflfac = 2._jprb * albedo
          reflfac_tl = 2._jprb * albedo_tl
          inc_surface = (reflfac > 0._jprb)
        ELSE
          reflfac = 2._jprb * diffuse_refl(i)
          reflfac_tl = 2._jprb * diffuse_refl_tl(i)
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
      DO m = 0, thisnaz

        ! Solar source term factor
        IF (dosolar) THEN
          IF (m == 0) THEN
            sunsrcfac = z4pi_r * solar_spectrum(i)
          ELSE
            sunsrcfac = 2._jprb * z4pi_r * solar_spectrum(i)
          ENDIF
        ENDIF

        kkpack = 0._jprb
        kkpack_tl = 0._jprb

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
              f_lay        = trans_scatt_ir_dyn%phasefn(i)%legcoef(nstr,icci,laymap) / (2._jprb * nstr + 1._jprb)
              ftau(lay)    = 1._jprb - ssa(icci,laymap,i) * f_lay
              fssa(lay)    = (1._jprb - f_lay) * ssa(icci,laymap,i) / ftau(lay)

              f_lay_tl     = trans_scatt_ir_dyn_tl%phasefn(i)%legcoef(nstr,icci,laymap) / (2._jprb * nstr + 1._jprb)
              ftau_tl(lay) = - ssa_tl(icci,laymap,i) * f_lay - &
                               ssa(icci,laymap,i) * f_lay_tl
              fssa_tl(lay) = - f_lay_tl * ssa(icci,laymap,i) / ftau(lay) + &
                             (1._jprb - f_lay) * ssa_tl(icci,laymap,i) / ftau(lay) - &
                             ftau_tl(lay) * (1._jprb - f_lay) * ssa(icci,laymap,i) / ftau(lay)**2
            ELSE
              f_lay        = 0._jprb
              ftau(lay)    = 1._jprb
              fssa(lay)    = 0._jprb

              f_lay_tl     = 0._jprb
              ftau_tl(lay) = 0._jprb
              fssa_tl(lay) = 0._jprb
            ENDIF
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate optical depths and transmittances
          ! --------------------------------------------------------------------------
          ! Calculate accumulated delta-scaled optical depth and transmittances on
          ! satellite-surface and sun-surface paths.
          IF (m == 0) THEN  ! This calculation applies to all azi components
            opdep(lev) = opdep(lev-1) + ftau(lay) * profiles_dom(icc,i)%layerod(lay)
            tausat(lev) = EXP(-opdep(lev) * musatr)
            IF (dosolar) tausun(lev) = EXP(-opdep(lev) * musunr)

            opdep_tl(lev) = opdep_tl(lev-1) + ftau_tl(lay) * profiles_dom(icc,i)%layerod(lay) + &
                                              ftau(lay) * profiles_dom_tl(icc,i)%layerod(lay)
            tausat_tl(lev) = -opdep_tl(lev) * musatr * tausat(lev)
            IF (dosolar) tausun_tl(lev) = -opdep_tl(lev) * musunr * tausun(lev)
          ENDIF

          ! --------------------------------------------------------------------------
          ! Calculate emissive source term and interpolated clear-sky emissive source term
          ! --------------------------------------------------------------------------
          IF (dothermal) THEN
            IF (lay == nlayers .AND. profiles_dom(icc,i)%surface) THEN
              b1(lay) = (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) / (opdep(lev-1) - opdep(lev))

              b1_tl(lay) = - fssa_tl(lay) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) / (opdep(lev-1) - opdep(lev)) + &
                (1._jprb - fssa(lay)) * (auxrad_tl%air(lev-1,i) - auxrad_tl%surfair(i)) / (opdep(lev-1) - opdep(lev)) - &
                (opdep_tl(lev-1) - opdep_tl(lev)) * &
                (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%surfair(i)) / (opdep(lev-1) - opdep(lev))**2
            ELSE
              b1(lay) = (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) / (opdep(lev-1) - opdep(lev))

              b1_tl(lay) = - fssa_tl(lay) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) / (opdep(lev-1) - opdep(lev)) + &
                (1._jprb - fssa(lay)) * (auxrad_tl%air(lev-1,i) - auxrad_tl%air(lev,i)) / (opdep(lev-1) - opdep(lev)) - &
                (opdep_tl(lev-1) - opdep_tl(lev)) * &
                (1._jprb - fssa(lay)) * (auxrad%air(lev-1,i) - auxrad%air(lev,i)) / (opdep(lev-1) - opdep(lev))**2
            ENDIF
            b0(lay) = (1._jprb - fssa(lay)) * auxrad%air(lev-1,i) - b1(lay) * opdep(lev-1)
            b0_tl(lay) = - fssa_tl(lay) * auxrad%air(lev-1,i) + (1._jprb - fssa(lay)) * auxrad_tl%air(lev-1,i) - &
                         b1_tl(lay) * opdep(lev-1) - b1(lay) * opdep_tl(lev-1)

            y1_tl(:,lay) = b1_tl(lay)

            ! For clear layers we have the solution and the interpolated thermal source term is trivial
            IF (.NOT. fssa(lay) > 0._jprb) THEN
              y0_tl(1:n,lay)      = b0_tl(lay) + muq(1:n) * y1_tl(1:n,lay)
              y0_tl(n+1:nstr,lay) = b0_tl(lay) - muq(1:n) * y1_tl(n+1:nstr,lay)
              fpmupi0(lay)        = b0(lay)
              fpmupi1(lay)        = b1(lay)
              fpmupi0_tl(lay)     = b0_tl(lay)
              fpmupi1_tl(lay)     = b1_tl(lay)
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
                  legcoef_tl(l,laymap,icci) = fssa_tl(lay) * (trans_scatt_ir_dyn%phasefn(i)%legcoef(l,icci,laymap) - &
                      (2._jprb * l + 1._jprb) * f_lay) / &
                      (1._jprb - f_lay) + &
                      fssa(lay) * (trans_scatt_ir_dyn_tl%phasefn(i)%legcoef(l,icci,laymap) - &
                      (2._jprb * l + 1._jprb) * f_lay_tl) / &
                      (1._jprb - f_lay) + &
                      f_lay_tl * fssa(lay) * (trans_scatt_ir_dyn%phasefn(i)%legcoef(l,icci,laymap) - &
                      (2._jprb * l + 1._jprb) * f_lay) / &
                      (1._jprb - f_lay)**2
                ENDDO
              ENDIF

              ! Due to symmetry we can calculate just n*nstr values (not nstr*nstr):
              !    phase_fn(i,j) = phase_fn(n+i,n+j)    where n+x => -ve mu
              !    phase_fn(n+i,j) = phase_fn(i,n+j)
              DO j = 1, n
                DO k = j, n
                  phase_fn(k,j)      = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,j,m))
                  phase_fn(k,n+j)    = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,n+j,m))
                  phase_fn_tl(k,j)   = SUM(legcoef_tl(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,j,m))
                  phase_fn_tl(k,n+j) = SUM(legcoef_tl(m:nmom,laymap,icci) * legpolymuq(m:nmom,k,m) * legpolymuq(m:nmom,n+j,m))

                  ! Exploit symmetry to avoid repeated calculations
                  IF (k > j) THEN
                    phase_fn(j,k)      = phase_fn(k,j)
                    phase_fn(j,n+k)    = phase_fn(k,n+j)
                    phase_fn_tl(j,k)   = phase_fn_tl(k,j)
                    phase_fn_tl(j,n+k) = phase_fn_tl(k,n+j)
                  ENDIF
                ENDDO
              ENDDO

              ! Include quadrature weights in phase fns evaluated for musat and each quadrature angle
              DO j = 1, n
                phasesat(j) = muw(j) * SUM(legcoef(m:nmom,laymap,icci) * &
                                           legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
                phasesat_tl(j) = muw(j) * SUM(legcoef_tl(m:nmom,laymap,icci) * &
                                              legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO
              DO j = n+1, nstr
                phasesat(j) = muw(j-n) * SUM(legcoef(m:nmom,laymap,icci) * &
                                             legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
                phasesat_tl(j) = muw(j-n) * SUM(legcoef_tl(m:nmom,laymap,icci) * &
                                                legpolymusat(m:nmom,0,m) * legpolymuq(m:nmom,j,m))
              ENDDO

            ENDIF ! dothermal or not done_lay_str_azi

            IF (.NOT. done_lay_str_azi(laymap,icci,m)) THEN

              ! Homogenous solution only done once per layer for scattering layers (applies to all cloud columns)

              ! --------------------------------------------------------------------------
              ! Solve eigenvalue problem and obtain homogenous solution for this layer
              ! --------------------------------------------------------------------------
              DO j = 1, n
                bp(:,j)    = 0.5_jprb * muw(j) * phase_fn(:,j) * muqr(:)
                bm(:,j)    = 0.5_jprb * muw(j) * phase_fn(:,n+j) * muqr(:)
                bp(j,j)    = bp(j,j) - muqr(j)

                bp_tl(:,j) = 0.5_jprb * muw(j) * phase_fn_tl(:,j) * muqr(:)
                bm_tl(:,j) = 0.5_jprb * muw(j) * phase_fn_tl(:,n+j) * muqr(:)
              ENDDO
              sm    = bp - bm
              sp    = bp + bm
!               bpm        = MATMUL(sm, sp)
!               bpmsave    = bpm  ! Used for solar order-n PI below

              sm_tl = bp_tl - bm_tl
              sp_tl = bp_tl + bm_tl
              bpm_tl     = MATMUL(sm_tl, sp) + MATMUL(sm, sp_tl)
!               bpmsave_tl = bpm_tl

              xp = dom_state(i)%xp(:,:,laymap,icci,m)

              ! The eigenvalues output by ASYMTX are required here (before we take the square root)
              eval = dom_state(i)%eval(:,laymap,icci,m)**2
              CALL ASYMTX_TL(err, n, xp, eval, bpm_tl, xp_tl, eval_tl(:,laymap,icci,m))
              IF (err /= 0) THEN
                WRITE(msg,'(a,i5)') 'DOM error: solving TL eigenvalue problem, ASYMTX_TL error = ', err
                THROWM(err.NE.0,msg)
              ENDIF

              ! Calculate the required eigenvalues (+/- each one gives nstr evals)
              eval_tl(:,laymap,icci,m) = 0.5_jprb * eval_tl(:,laymap,icci,m) / dom_state(i)%eval(:,laymap,icci,m)

              DO j = 1, n
                xm(:,j) = MATMUL(sp, xp(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
                xm_tl(:,j) = (MATMUL(sp_tl, xp(:,j)) + MATMUL(sp, xp_tl(:,j)) - &
                             eval_tl(j,laymap,icci,m) * xm(:,j)) / dom_state(i)%eval(j,laymap,icci,m)
              ENDDO

              ! Rearrange to obtain F+ and F- (the eigenvectors)
              ! For +ve eigenvalues:
              fp1(:,:,laymap,icci,m)    = 0.5_jprb * (xp + xm)
              fm1(:,:,laymap,icci,m)    = 0.5_jprb * (xp - xm)

              fp1_tl(:,:,laymap,icci,m) = 0.5_jprb * (xp_tl + xm_tl)
              fm1_tl(:,:,laymap,icci,m) = 0.5_jprb * (xp_tl - xm_tl)

              ! Solar source term particular integral only done once per layer for scattering layers (applies to all cloud columns)
              IF (dosolar) THEN
                ! --------------------------------------------------------------------------
                ! Solve for particular integral for solar source term
                ! --------------------------------------------------------------------------

!                 IF (.NOT. do_order_n_pi) THEN

                  ! Solve PI problem using nstr x nstr system (most straightforward approach)

                  ! Construct and solve the linear system for scattering layers
                  DO j = 1, n
                    hh(1:n,j)         = -0.5_jprb * muw(j) * phase_fn(:,j)
                    hh(n+1:nstr,j)    = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry
                    hh_tl(1:n,j)      = -0.5_jprb * muw(j) * phase_fn_tl(:,j)
                    hh_tl(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn_tl(:,n+j)   ! Use phase fn symmetry

                    DO k = 1, n
                      hh(k,n+j)      = hh(n+k,j)
                      hh(n+k,n+j)    = hh(k,j)
                      hh_tl(k,n+j)   = hh_tl(n+k,j)
                      hh_tl(n+k,n+j) = hh_tl(k,j)
                    ENDDO

                    hh(j,j)     = hh(j,j) + 1._jprb + muq(j) * musunr
                    hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb - muq(j) * musunr

                    z_tl(j,laymap,icci,m)   = sunsrcfac * SUM(legcoef_tl(m:nmom,laymap,icci) * &
                                                              legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m))
                    z_tl(n+j,laymap,icci,m) = sunsrcfac * SUM(legcoef_tl(m:nmom,laymap,icci) * &
                                                              legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m))
                  ENDDO

                  ! Compute LU factorisation of hh
                  CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
                  IF (info /= 0) THEN
                    WRITE(msg,'(a,3i5)') 'DOM error: solving TL PI solar, DGETRF error, info = ', info, i, lay
                    err = errorstatus_fatal
                    THROWM(err.NE.0,msg)
                  ENDIF

                  ! Solve TL linear system
                  z_tl(:,laymap,icci,m) = z_tl(:,laymap,icci,m) - MATMUL(hh_tl, dom_state(i)%z(:,laymap,icci,m))
                  CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, z_tl(:,laymap,icci,m), nstri, info)
                  IF (info /= 0) THEN
                    WRITE(msg,'(a,3i5)') 'DOM error: solving TL PI solar, DGETRS error, info = ', info, i, lay
                    err = errorstatus_fatal
                    THROWM(err.NE.0,msg)
                  ENDIF


!                 ELSE
!                   ! Solve PI problem using n x n system (instead of nstr x nstr)
!                   ! This offers a small saving in time (e.g. 1% for nstr=16, should get more benefit as nstr increases)
!                   ! By having to store results from homogenous solution it increases memory usage very slightly
!                   ! e.g. by a few times (n x n) (could get rid of hh for solar calcs when using this method)
! 
!                   DO j = 1, n
!                     so1 = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m))
!                     so1_tl = SUM(legcoef_tl(m:nmom,laymap,icci) * legpolymuq(m:nmom,j,m) * legpolymusun(m:nmom,0,m))
!                     so2 = SUM(legcoef(m:nmom,laymap,icci) * legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m))
!                     so2_tl = SUM(legcoef_tl(m:nmom,laymap,icci) * legpolymuq(m:nmom,n+j,m) * legpolymusun(m:nmom,0,m))
!                     sop(j) = sunsrcfac * (so1 + so2) * muqr(j) * musun
!                     sop_tl(j) = sunsrcfac * (so1_tl + so2_tl) * muqr(j) * musun
!                     som(j) = sunsrcfac * (so1 - so2) * muqr(j) * musun
!                     som_tl(j) = sunsrcfac * (so1_tl - so2_tl) * muqr(j) * musun
!                   ENDDO
! 
!                   gp = som + musun * MATMUL(sm, sop)
!                   gp_tl = som_tl + musun * (MATMUL(sm_tl, sop) + MATMUL(sm, sop_tl))
! 
!                   bpmsave = -musun**2 * bpmsave  ! bpmsave overwritten as it's not required anywhere else
!                   bpmsave_tl = -musun**2 * bpmsave_tl
!                   DO j = 1, n
!                     bpmsave(j,j) = 1._jprb + bpmsave(j,j)
!                   ENDDO
! 
!                   ! For a linear system of equations  A.x = y
!                   ! So for the TL we have  A.x_tl = y_tl - A_tl.x
!                   ! i.e. we solve a linear system with the same coefficient matrix A
! 
!                   ! In this case A == bpmsave  and  y == gp, so the RHS becomes gp_tl - MATMUL(bpmsave_tl, gp)
! 
!                   ! Compute LU factorisation of bpmsave
!                   CALL DGETRF(n, n, bpmsave, n, ipiv2, info)
!                   IF (info /= 0) THEN
!                     PRINT *,'DGETRF (solving PI) error, info = ', info, i, lay
!                     err = errorstatus_fatal
!                     THROWM(err.NE.0,'DOM error: solving for particular integral - solar')
!                   ENDIF
! 
!                   ! Solve direct linear system
!                   CALL DGETRS('No transpose', n, 1, bpmsave, n, ipiv2, gp, n, info)
!                   IF (info /= 0) THEN
!                     PRINT *,'DGETRS (solving PI) error, info = ', info, i, lay
!                     err = errorstatus_fatal
!                     THROWM(err.NE.0,'DOM error: solving for direct particular integral - solar')
!                   ENDIF
! 
!                   ! Solve TL linear system
!                   gp_tl = gp_tl - MATMUL(bpmsave_tl, gp)
!                   CALL DGETRS('No transpose', n, 1, bpmsave, n, ipiv2, gp_tl, n, info)
!                   IF (info /= 0) THEN
!                     PRINT *,'DGETRS (solving PI) error, info = ', info, i, lay
!                     err = errorstatus_fatal
!                     THROWM(err.NE.0,'DOM error: solving for TL particular integral - solar')
!                   ENDIF
! 
!                   ! Re-use som array for calculation of "gm":
!                   som = musun * MATMUL(sp, gp) + sop
!                   som_tl = musun * (MATMUL(sp_tl, gp) + MATMUL(sp, gp_tl)) + sop_tl
!                   z_tl(1:n,laymap,icci,m)      = 0.5_jprb * (gp_tl + som_tl)
!                   z_tl(n+1:nstr,laymap,icci,m) = 0.5_jprb * (gp_tl - som_tl)
! 
!                 ENDIF

              ENDIF ! dosolar

            ENDIF ! not done_lay_str_azi


            ! Thermal source term done for every layer and *for every cloud column*
            IF (dothermal) THEN
              ! --------------------------------------------------------------------------
              ! Solve for particular integral for thermal source term
              ! --------------------------------------------------------------------------

              ! Construct and solve the linear system for scattering layers
              DO j = 1, n
                hh(1:n,j)         = -0.5_jprb * muw(j) * phase_fn(:,j)
                hh(n+1:nstr,j)    = -0.5_jprb * muw(j) * phase_fn(:,n+j)   ! Use phase fn symmetry
                hh_tl(1:n,j)      = -0.5_jprb * muw(j) * phase_fn_tl(:,j)
                hh_tl(n+1:nstr,j) = -0.5_jprb * muw(j) * phase_fn_tl(:,n+j)   ! Use phase fn symmetry

                DO k = 1, n
                  hh(k,n+j)      = hh(n+k,j)
                  hh(n+k,n+j)    = hh(k,j)
                  hh_tl(k,n+j)   = hh_tl(n+k,j)
                  hh_tl(n+k,n+j) = hh_tl(k,j)
                ENDDO

                hh(j,j)     = hh(j,j) + 1._jprb
                hh(n+j,n+j) = hh(n+j,n+j) + 1._jprb
              ENDDO

              CALL DGETRF(nstri, nstri, hh, nstri, ipiv2, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving TL PI thermal, DGETRF error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              y1_tl(:,lay) = y1_tl(:,lay) - MATMUL(hh_tl, dom_state(i)%y1(:,lay,icc))
              CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, y1_tl(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving TL PI thermal y1, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

              y0_tl(1:n,lay)      = b0_tl(lay) + muq(1:n) * y1_tl(1:n,lay)
              y0_tl(n+1:nstr,lay) = b0_tl(lay) - muq(1:n) * y1_tl(n+1:nstr,lay)

              y0_tl(:,lay) = y0_tl(:,lay) - MATMUL(hh_tl, dom_state(i)%y0(:,lay,icc))
              CALL DGETRS('No transpose', nstri, 1, hh, nstri, ipiv2, y0_tl(:,lay), nstri, info)
              IF (info /= 0) THEN
                WRITE(msg,'(a,3i5)') 'DOM error: solving TL PI thermal y0, DGETRS error, info = ', info, i, lay
                err = errorstatus_fatal
                THROWM(err.NE.0,msg)
              ENDIF

            ENDIF ! dothermal

            ! --------------------------------------------------------------------------
            ! Interpolate HS and PI for final radiance calculation
            ! --------------------------------------------------------------------------

            IF (.NOT. done_lay_str_azi(laymap,icci,m)) THEN

              ! Homogenous part - diffuse radiation, done once per layer (applies to all cloud columns)
              DO j = 1, n
                fpmuhsup(j,laymap,icci,m) = 0.5_jprb * (SUM(fp1(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm1(:,j,laymap,icci,m) * phasesat(n+1:nstr)))

                fpmuhsup_tl(j,laymap,icci,m) = 0.5_jprb * (SUM(fp1_tl(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                           SUM(fp1(:,j,laymap,icci,m) * phasesat_tl(1:n)) + &
                                                           SUM(fm1_tl(:,j,laymap,icci,m) * phasesat(n+1:nstr)) + &
                                                           SUM(fm1(:,j,laymap,icci,m) * phasesat_tl(n+1:nstr)))

                fpmuhsdn(j,laymap,icci,m) = 0.5_jprb * (SUM(fp2(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                        SUM(fm2(:,j,laymap,icci,m) * phasesat(n+1:nstr)))

                fpmuhsdn_tl(j,laymap,icci,m) = 0.5_jprb * (SUM(fp2_tl(:,j,laymap,icci,m) * phasesat(1:n)) + &
                                                           SUM(fp2(:,j,laymap,icci,m) * phasesat_tl(1:n)) + &
                                                           SUM(fm2_tl(:,j,laymap,icci,m) * phasesat(n+1:nstr)) + &
                                                           SUM(fm2(:,j,laymap,icci,m) * phasesat_tl(n+1:nstr)))
              ENDDO

              ! Particular integral - solar source term, done once per layer (applies to all cloud columns)
              IF (dosolar) THEN
                fpmupi(laymap,icci,m) = 0.5_jprb * SUM(dom_state(i)%z(:,laymap,icci,m) * phasesat(:)) + &
                                        sunsrcfac * SUM(legcoef(m:nmom,laymap,icci) * &
                                                        legpolymusat(m:nmom,0,m) * legpolymusun(m:nmom,0,m))
                fpmupi_tl(laymap,icci,m) = 0.5_jprb * SUM(z_tl(:,laymap,icci,m) * phasesat(:)) + &
                                           0.5_jprb * SUM(dom_state(i)%z(:,laymap,icci,m) * phasesat_tl(:)) + &
                                           sunsrcfac * SUM(legcoef_tl(m:nmom,laymap,icci) * &
                                                           legpolymusat(m:nmom,0,m) * legpolymusun(m:nmom,0,m))
              ENDIF

            ENDIF ! not done_lay_str_azi

            ! Particular integral - thermal source term, done once per layer per cloud column
            IF (dothermal) THEN
              fpmupi0(lay)    = 0.5_jprb * SUM(dom_state(i)%y0(:,lay,icc) * phasesat(:)) + b0(lay)
              fpmupi0_tl(lay) = 0.5_jprb * (SUM(y0_tl(:,lay) * phasesat(:)) + &
                                            SUM(dom_state(i)%y0(:,lay,icc) * phasesat_tl(:))) + b0_tl(lay)
              fpmupi1(lay)    = 0.5_jprb * SUM(dom_state(i)%y1(:,lay,icc) * phasesat(:)) + b1(lay)
              fpmupi1_tl(lay) = 0.5_jprb * (SUM(y1_tl(:,lay) * phasesat(:)) + &
                                            SUM(dom_state(i)%y1(:,lay,icc) * phasesat_tl(:))) + b1_tl(lay)
            ENDIF

            done_lay_str_azi(laymap,icci,m) = .TRUE.

          ENDIF ! fssa > 0 (i.e. scattering layer)

          ! --------------------------------------------------------------------------
          ! Calculate layer "transmittances" with opdeps scaled by eigenvalues
          ! --------------------------------------------------------------------------
          tlayer(:,lay) = EXP(-ftau(lay) * profiles_dom(icc,i)%layerod(lay) * dom_state(i)%eval(:,laymap,icci,m))
          tlayer_tl(:,lay) = -tlayer(:,lay) * &
                              (ftau_tl(lay) * profiles_dom(icc,i)%layerod(lay) * dom_state(i)%eval(:,laymap,icci,m) + &
                               ftau(lay) * profiles_dom_tl(icc,i)%layerod(lay) * dom_state(i)%eval(:,laymap,icci,m) + &
                               ftau(lay) * profiles_dom(icc,i)%layerod(lay) * eval_tl(:,laymap,icci,m))

          ! --------------------------------------------------------------------------
          ! Construct part of linear system for this layer
          ! --------------------------------------------------------------------------
          IF (lay == 1) THEN
            row = 0
            col = 0
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,j)      = fm1(:,j,laymap,icci,m)
              kkpack_tl(lo:lo+n-1,j)   = fm1_tl(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,n+j)    = fm2(:,j,laymap,icci,m) * tlayer(j,lay)
              kkpack_tl(lo:lo+n-1,n+j) = fm2_tl(:,j,laymap,icci,m) * tlayer(j,lay) + &
                                         fm2(:,j,laymap,icci,m) * tlayer_tl(j,lay)
            ENDDO
            IF (dosolar) THEN
              x_tl(1:n) = -z_tl(n+1:nstr,laymap,icci,m)
            ELSE
              x_tl(1:n) = auxrad_tl%air(1,i) - y0_tl(n+1:nstr,lay)
            ENDIF
          ELSE
            row = (lay-2) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)      = -fp1(:,j,laymap,icci,m)
              kkpack_tl(lo:lo+n-1,col+j)   = -fp1_tl(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j)    = -fp2(:,j,laymap,icci,m) * tlayer(j,lay)
              kkpack_tl(lo:lo+n-1,col+n+j) = -fp2_tl(:,j,laymap,icci,m) * tlayer(j,lay) - &
                                              fp2(:,j,laymap,icci,m) * tlayer_tl(j,lay)
              lo = klu1 + (row + n + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)      = -fm1(:,j,laymap,icci,m)
              kkpack_tl(lo:lo+n-1,col+j)   = -fm1_tl(:,j,laymap,icci,m)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j)    = -fm2(:,j,laymap,icci,m) * tlayer(j,lay)
              kkpack_tl(lo:lo+n-1,col+n+j) = -fm2_tl(:,j,laymap,icci,m) * tlayer(j,lay) - &
                                              fm2(:,j,laymap,icci,m) * tlayer_tl(j,lay)
            ENDDO
            IF (dosolar) THEN
              x_tl(row+1:row+nstr) = x_tl(row+1:row+nstr) + z_tl(:,laymap,icci,m) * tausun(lev-1) + &
                                                            dom_state(i)%z(:,laymap,icci,m) * tausun_tl(lev-1)
            ELSE
              x_tl(row+1:row+nstr) = x_tl(row+1:row+nstr) + y0_tl(:,lay) + &
                                     y1_tl(:,lay) * opdep(lev-1) + dom_state(i)%y1(:,lay,icc) * opdep_tl(lev-1)
            ENDIF
          ENDIF

          IF (lay < nlayers) THEN
            row = (lay-1) * nstr + n
            col = (lay-1) * nstr
            DO j = 1, n
              lo = klu1 + (row + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)      = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
              kkpack_tl(lo:lo+n-1,col+j)   = fp1_tl(:,j,laymap,icci,m) * tlayer(j,lay) + &
                                             fp1(:,j,laymap,icci,m) * tlayer_tl(j,lay)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j)    = fp2(:,j,laymap,icci,m)
              kkpack_tl(lo:lo+n-1,col+n+j) = fp2_tl(:,j,laymap,icci,m)
              lo = klu1 + (row + n + 1) - (col + j)
              kkpack(lo:lo+n-1,col+j)      = fm1(:,j,laymap,icci,m) * tlayer(j,lay)
              kkpack_tl(lo:lo+n-1,col+j)   = fm1_tl(:,j,laymap,icci,m) * tlayer(j,lay) + &
                                             fm1(:,j,laymap,icci,m) * tlayer_tl(j,lay)
              lo = lo - n
              kkpack(lo:lo+n-1,col+n+j)    = fm2(:,j,laymap,icci,m)
              kkpack_tl(lo:lo+n-1,col+n+j) = fm2_tl(:,j,laymap,icci,m)
            ENDDO
            IF (dosolar) THEN
              x_tl(row+1:row+nstr) = -z_tl(:,laymap,icci,m) * tausun(lev) - dom_state(i)%z(:,laymap,icci,m) * tausun_tl(lev)
            ELSE
              x_tl(row+1:row+nstr) = -y0_tl(:,lay) - y1_tl(:,lay) * opdep(lev) - dom_state(i)%y1(:,lay,icc) * opdep_tl(lev)
            ENDIF
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

                  ! Using surfterm/reflfac is more efficient than recomputing the SUMs above
                  surfterm_tl(j)   = reflfac_tl * surfterm(j) / reflfac + &
                                     reflfac * SUM(muqw(:) * fm1_tl(:,j,laymap,icci,m))
                  surfterm_tl(n+j) = reflfac_tl * surfterm(n+j) / reflfac + &
                                     reflfac * SUM(muqw(:) * fm2_tl(:,j,laymap,icci,m))
                ENDDO
              ELSE
                surfterm    = 0._jprb
                surfterm_tl = 0._jprb
              ENDIF
              DO j = 1, n
                lo = klu1 + (row + 1) - (col + j)
                kkpack(lo:lo+n-1,col+j)      = (fp1(:,j,laymap,icci,m) - surfterm(j)) * tlayer(j,lay)
                kkpack_tl(lo:lo+n-1,col+j)   = (fp1_tl(:,j,laymap,icci,m) - surfterm_tl(j)) * tlayer(j,lay) + &
                                               (fp1(:,j,laymap,icci,m) - surfterm(j)) * tlayer_tl(j,lay)
                lo = lo - n
                kkpack(lo:lo+n-1,col+n+j)    = fp2(:,j,laymap,icci,m) - surfterm(n+j)
                kkpack_tl(lo:lo+n-1,col+n+j) = fp2_tl(:,j,laymap,icci,m) - surfterm_tl(n+j)
              ENDDO

              IF (dosolar) THEN
                tmpval = SUM(muqw(:) * dom_state(i)%z(n+1:nstr,laymap,icci,m))
                x_tl(row+1:row+n) = -tausun_tl(lev) * (dom_state(i)%z(1:n,laymap,icci,m) - &
                                    reflfac * tmpval - (albedo / pi) * solar_spectrum(i) * musun) - &
                                    tausun(lev) * (z_tl(1:n,laymap,icci,m) - reflfac_tl * tmpval - &
                                    reflfac * SUM(muqw(:) * z_tl(n+1:nstr,laymap,icci,m)) - &
                                    (albedo_tl / pi) * solar_spectrum(i) * musun)
              ELSE
                x_tl(row+1:row+n) = -(y0_tl(1:n,lay) + y1_tl(1:n,lay) * opdep(lev) + &
                                      dom_state(i)%y1(1:n,lay,icc) * opdep_tl(lev) - &
                                    reflfac_tl * SUM(muqw(:) * (dom_state(i)%y0(n+1:nstr,lay,icc) + &
                                                     dom_state(i)%y1(n+1:nstr,lay,icc) * opdep(lev))) - &
                                    reflfac * SUM(muqw(:) * (y0_tl(n+1:nstr,lay) + &
                                                             y1_tl(n+1:nstr,lay) * opdep(lev) + &
                                                             dom_state(i)%y1(n+1:nstr,lay,icc) * opdep_tl(lev))) - &
                                    emissivity_tl(i)%emis_out * auxrad%skin(i) - emissivity(i)%emis_out * auxrad_tl%skin(i))
              ENDIF
            ELSE
              DO j = 1, n
                lo = klu1 + (row + 1) - (col + j)
                kkpack(lo:lo+n-1,col+j)      = fp1(:,j,laymap,icci,m) * tlayer(j,lay)
                kkpack_tl(lo:lo+n-1,col+j)   = fp1_tl(:,j,laymap,icci,m) * tlayer(j,lay) + &
                                               fp1(:,j,laymap,icci,m) * tlayer_tl(j,lay)
                lo = lo - n
                kkpack(lo:lo+n-1,col+n+j)    = fp2(:,j,laymap,icci,m)
                kkpack_tl(lo:lo+n-1,col+n+j) = fp2_tl(:,j,laymap,icci,m)
              ENDDO
              IF (dosolar) THEN
                x_tl(row+1:row+n) = -z_tl(1:n,laymap,icci,m) * tausun(lev) - &
                                     dom_state(i)%z(1:n,laymap,icci,m) * tausun_tl(lev)
              ELSE
                x_tl(row+1:row+n) = -y0_tl(1:n,lay) - y1_tl(1:n,lay) * opdep(lev) - dom_state(i)%y1(1:n,lay,icc) * opdep_tl(lev)
              ENDIF
            ENDIF
          ENDIF


          ! --------------------------------------------------------------------------
          ! Accumulate TMS correction for this layer
          ! --------------------------------------------------------------------------
          IF (dosolar .AND. m == 0 .AND. fssa(lay) > 0._jprb) THEN

            taufac = tausun(lev-1) * tausat(lev-1) - tausun(lev) * tausat(lev)
            taufac_tl = tausun_tl(lev-1) * tausat(lev-1) + tausun(lev-1) * tausat_tl(lev-1) - &
                        tausun_tl(lev) * tausat(lev) - tausun(lev) * tausat_tl(lev)

            tmpval = sunsrcfac / (1._jprb + musat * musunr)
            tmpval2 = 1._jprb / ftau(lay)

            delta_tms_tl = delta_tms_tl + tmpval * taufac_tl * &
              (ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) * tmpval2 - &
               SUM(legcoef(:,laymap,icci) * legpolymuscata(:,0,0))) + &
              tmpval * taufac * &
              (ssa_tl(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) * tmpval2 + &
               ssa(icci,laymap,i) * trans_scatt_ir_tl%phup(icci,laymap,i) * tmpval2 - &
               ssa(icci,laymap,i) * trans_scatt_ir%phup(icci,laymap,i) * ftau_tl(lay) * tmpval2**2 - &
               SUM(legcoef_tl(:,laymap,icci) * legpolymuscata(:,0,0)))
          ENDIF

        ENDDO ! layer loop


        ! --------------------------------------------------------------------------
        ! Solve the multi-layer system
        ! --------------------------------------------------------------------------

        CALL DGBTRF(nstrlay, nstrlay, kl, ku, kkpack, klku, ipiv1, info)
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving TL linear system, DGBTRF error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF

        ! x_tl = x_tl - kk_tl * x
        DO j = 1, nstrlay
          ilo = MAX(1, j - ku)
          ihi = MIN(nstr * nlayers, j + kl)
          lo = klu1 + ilo - j
          hi = klu1 + ihi - j
          x_tl(ilo:ihi) = x_tl(ilo:ihi) - kkpack_tl(lo:hi,j) * dom_state(i)%x(j,m,icc)
        ENDDO
        CALL DGBTRS('No transpose', nstrlay, kl, ku, 1, kkpack, klku, ipiv1, x_tl, nstrlay, info)
        IF (info /= 0) THEN
          WRITE(msg,'(a,i5)') 'DOM error: solving TL linear system, DGBTRS error, info = ', info
          err = errorstatus_fatal
          THROWM(err.NE.0,msg)
        ENDIF


        ! --------------------------------------------------------------------------
        ! Calculate TOA radiance
        ! --------------------------------------------------------------------------
        radinc_tl = 0._jprb
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

            ! Homogenous part - diffuse radiation
            DO j = 1, n
              taufac    = tausat(lev-1) - tlayer(j,lay) * tausat(lev)
              taufac_tl = tausat_tl(lev-1) - tlayer_tl(j,lay) * tausat(lev) - tlayer(j,lay) * tausat_tl(lev)

              tmpval = 1._jprb / (1._jprb + dom_state(i)%eval(j,laymap,icci,m) * musat)
              radinc_tl = radinc_tl + x_tl(row+j) * fpmuhsup(j,laymap,icci,m) * taufac * tmpval + &
                                      dom_state(i)%x(row+j,m,icc) * fpmuhsup_tl(j,laymap,icci,m) * taufac * tmpval + &
                                      dom_state(i)%x(row+j,m,icc) * fpmuhsup(j,laymap,icci,m) * taufac_tl * tmpval - &
                                      eval_tl(j,laymap,icci,m) * musat * &
                                      dom_state(i)%x(row+j,m,icc) * fpmuhsup(j,laymap,icci,m) * taufac * tmpval**2

              taufac    = tlayer(j,lay) * tausat(lev-1) - tausat(lev)
              taufac_tl = tlayer_tl(j,lay) * tausat(lev-1) + tlayer(j,lay) * tausat_tl(lev-1) - tausat_tl(lev)

              tmpval = 1._jprb / (1._jprb - dom_state(i)%eval(j,laymap,icci,m) * musat)
              radinc_tl = radinc_tl + x_tl(row+n+j) * fpmuhsdn(j,laymap,icci,m) * taufac * tmpval + &
                                      dom_state(i)%x(row+n+j,m,icc) * fpmuhsdn_tl(j,laymap,icci,m) * taufac * tmpval + &
                                      dom_state(i)%x(row+n+j,m,icc) * fpmuhsdn(j,laymap,icci,m) * taufac_tl * tmpval + &
                                      eval_tl(j,laymap,icci,m) * musat * &
                                      dom_state(i)%x(row+n+j,m,icc) * fpmuhsdn(j,laymap,icci,m) * taufac * tmpval**2
            ENDDO

            ! Particular integral - solar source term
            IF (dosolar) THEN
              taufac    = tausat(lev-1) * tausun(lev-1) - tausat(lev) * tausun(lev)
              taufac_tl = tausat_tl(lev-1) * tausun(lev-1) + tausat(lev-1) * tausun_tl(lev-1) - &
                          tausat_tl(lev) * tausun(lev) - tausat(lev) * tausun_tl(lev)
              radinc_tl = radinc_tl + fpmupi_tl(laymap,icci,m) * taufac / (1._jprb + musunr * musat) + &
                                      fpmupi(laymap,icci,m) * taufac_tl / (1._jprb + musunr * musat)
            ENDIF
          ENDIF

          ! Particular integral - emissive source term
          IF (dothermal) THEN
            radinc_tl = radinc_tl + fpmupi0_tl(lay) * (tausat(lev-1) - tausat(lev)) + &
                                    fpmupi0(lay) * (tausat_tl(lev-1) - tausat_tl(lev)) + &
                                    fpmupi1_tl(lay) * ((opdep(lev-1) + musat) * tausat(lev-1) - &
                                                       (opdep(lev) + musat) * tausat(lev)) + &
                                    fpmupi1(lay) * (opdep_tl(lev-1) * tausat(lev-1) - &
                                                    opdep_tl(lev) * tausat(lev) + &
                                                    (opdep(lev-1) + musat) * tausat_tl(lev-1) - &
                                                    (opdep(lev) + musat) * tausat_tl(lev))
          ENDIF

        ENDDO ! nlayers


        ! Surface contribution - assumes surface is Lambertian
        IF (m == 0 .AND. inc_surface) THEN

          ! Integrate over radiances along downwelling quadrature angles and then add the surface term
          radsurfup_sum_tl = 0._jprb
          row = nstr * (nlayers - 1)
          laymap = profiles_dom(icc,i)%laymap(nlayers)
          IF (laymap > 0) THEN
            icci = ircld%icldarr(icc,laymap,prof)
          ELSE
            icci = 0
          ENDIF

          DO j = 1, n
            radsurfup_sum_tl = radsurfup_sum_tl + muqw(j) * &
              (SUM(x_tl(row+1:row+n) * fm1(j,:,laymap,icci,m) * tlayer(:,nlayers) + &
                   dom_state(i)%x(row+1:row+n,m,icc) * fm1_tl(j,:,laymap,icci,m) * tlayer(:,nlayers) + &
                   dom_state(i)%x(row+1:row+n,m,icc) * fm1(j,:,laymap,icci,m) * tlayer_tl(:,nlayers)) + &
               SUM(x_tl(row+n+1:row+nstr) * fm2(j,:,laymap,icci,m) + &
                   dom_state(i)%x(row+n+1:row+nstr,m,icc) * fm2_tl(j,:,laymap,icci,m)))
          ENDDO

          IF (dosolar) THEN
            DO j = 1, n
              radsurfup_sum_tl = radsurfup_sum_tl + muqw(j) * &
                (z_tl(n+j,laymap,icci,m) * tausun(nlevels) + dom_state(i)%z(n+j,laymap,icci,m) * tausun_tl(nlevels))
            ENDDO

            radsurfup_tl = reflfac_tl * dom_state(i)%radsurfup_sum(icc) + reflfac * radsurfup_sum_tl + &
                           (albedo_tl / pi) * solar_spectrum(i) * musun * tausun(nlevels) + &
                           (albedo / pi) * solar_spectrum(i) * musun * tausun_tl(nlevels)
          ELSE
            DO j = 1, n
              radsurfup_sum_tl = radsurfup_sum_tl + muqw(j) * &
                (y0_tl(n+j,nlayers) + y1_tl(n+j,nlayers) * opdep(nlevels) + dom_state(i)%y1(n+j,nlayers,icc) * opdep_tl(nlevels))
            ENDDO
            radsurfup_tl = reflfac_tl * dom_state(i)%radsurfup_sum(icc) + reflfac * radsurfup_sum_tl + &
                           emissivity_tl(i)%emis_out * auxrad%skin(i) + emissivity(i)%emis_out * auxrad_tl%skin(i)
          ENDIF

          radinc_tl = radinc_tl + radsurfup_tl * tausat(nlevels) + &
                      dom_state(i)%radsurfup(icc) * tausat_tl(nlevels)
        ENDIF

        ! Multiplier for azimuthal component
        IF (dosolar) radinc_tl = radinc_tl * COS(m * relazi)
        thisrad_tl = thisrad_tl + radinc_tl

      ENDDO ! azimuth loop

      ! --------------------------------------------------------------------------
      ! Apply TMS correction and accumulate radiance contribution
      ! --------------------------------------------------------------------------
      IF (dosolar) thisrad_tl = thisrad_tl + delta_tms_tl

      IF (icc == 0) THEN
        rad_tl%clear(i) = rad_tl%clear(i) + thisrad_tl
        rad_tl%total(i) = rad_tl%total(i) + thisrad_tl * ircld%xcolclr(prof) + &
                          dom_state(i)%thisrad(icc) * ircld_tl%xcolclr(prof)
      ELSE
        rad_tl%total(i) = rad_tl%total(i) + thisrad_tl * (ircld%xcol(icc+1,prof) - ircld%xcol(icc,prof)) + &
                          dom_state(i)%thisrad(icc) * (ircld_tl%xcol(icc+1,prof) - ircld_tl%xcol(icc,prof))
      ENDIF

    ENDDO ! cloud columns

    rad_tl%cloudy(i) = rad_tl%total(i)

  ENDDO ! channel

  ! --------------------------------------------------------------------------
  ! Tidy up
  ! --------------------------------------------------------------------------
  NULLIFY(fp2, fm2, fp2_tl, fm2_tl)

  IF (dosolar) THEN
    DEALLOCATE(legpolymusun, legpolymuscata, &
               tausun, fpmupi, &
               tausun_tl, fpmupi_tl, z_tl, stat=err)
  ELSE
    DEALLOCATE(b0, b1, fpmupi0, fpmupi1, &
               b0_tl, b1_tl, fpmupi0_tl, fpmupi1_tl, y0_tl, y1_tl, stat=err)
  ENDIF
  THROWM(err.NE.0, 'DOM deallocation error')

  DEALLOCATE(legpolymuq, legpolymusat, fpmuhsup, fpmuhsdn, fp1, fm1, &
             eval_tl, fpmuhsup_tl, fpmuhsdn_tl, fp1_tl, fm1_tl, &
             done_lay_str_azi, stat=err)
  THROWM(err.NE.0, 'DOM deallocation error')

CATCH

END SUBROUTINE rttov_dom_tl
