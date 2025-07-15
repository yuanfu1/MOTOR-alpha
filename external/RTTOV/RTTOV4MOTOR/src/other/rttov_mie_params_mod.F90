! Description:
!> @file
!!   Various routines for visible/IR aerosol/cloud optical property
!!   calculations
!
!> @brief
!!   Various routines for visible/IR aerosol/cloud optical property
!!   calculations
!!
!! @details
!!   This module contains subroutines that are used for calculating cloud and
!!   aerosol optical properties for visible/IR sensors. These include:
!!   - read_opac_ref_index: read OPAC refractive index file
!!   - read_segelstein_ref_index: read Segelstein (1981) refractive index data
!!   - interpolate_ref_index: interpolate refractive index data onto given wavelengths
!!   - make_radius_grid: creates a grid of radius values
!!   - mie_sphere: calculate Mie optical properties
!!   - calc_mie_ext_sca: calculates Mie extinction and scattering coefficients
!!   - calc_mie_phasefn: calculates Mie phase function
!!   - calc_mass: calculates mass integral over size distribution
!!   - lognorm: log-normal size distribution function
!!   - gammadist: modified gamma size distribution function
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
MODULE rttov_mie_params_mod

#include "throw.h"

  USE rttov_scattering_mod, ONLY : integrate, inter

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_opac_ref_index, read_segelstein_ref_index, &
            interpolate_ref_index, &
            lognorm, gammadist, make_radius_grid, calc_mass, &
            mie_sphere, calc_mie_ext_sca, calc_mie_phasefn, &
            init_new_optp

#include "rttov_errorreport.interface"
#include "rttov_nullify_optp_data.interface"

CONTAINS

!> Initialise a new aerosol/cloud optical property structure
!! @param[in]      coef        optical depth coefficients structure
!! @param[in,out]  optp        optical property structure to initialise
!! @param[in]      ntypes      number of particle types
!! @param[in]      chou_only   only compute Chou-scaling properties if true
SUBROUTINE init_new_optp(coef, optp, ntypes, chou_only)
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_coef, rttov_optp

  TYPE(rttov_coef), INTENT(IN)    :: coef
  TYPE(rttov_optp), INTENT(INOUT) :: optp
  INTEGER(jpim),    INTENT(IN)    :: ntypes
  LOGICAL(jplm),    INTENT(IN)    :: chou_only

  INTEGER(jpim) :: i
  INTEGER(jpim), ALLOCATABLE :: solar_chanlist(:)

  ! Number of channels
  optp%nchan = coef%fmv_chn

  optp%nchan_pha = 0

  IF (.NOT. chou_only) THEN
    ! Create the index lookup for solar channels
    ALLOCATE(optp%chan_pha_index(coef%fmv_chn))
    ALLOCATE(solar_chanlist(coef%fmv_chn))
    optp%chan_pha_index(:) = 0
    solar_chanlist(:) = 0

    DO i = 1, coef%fmv_chn
      IF (coef%ss_val_chn(i) > 0) THEN
        optp%nchan_pha = optp%nchan_pha + 1
        solar_chanlist(optp%nchan_pha) = i
        optp%chan_pha_index(i) = optp%nchan_pha
      ENDIF
    ENDDO
    
    IF (optp%nchan_pha > 0) THEN
      ALLOCATE(optp%chan_pha(optp%nchan_pha))
      optp%chan_pha(:) = solar_chanlist(1:optp%nchan_pha)
    ENDIF
    DEALLOCATE(solar_chanlist)
  ENDIF

  ! We now have:
  ! optp%nchan_pha = number of solar channels (i.e. with phase functions)
  ! optp%chan_pha_index(:) = index into phase array for each solar channel
  ! optp%chan_pha(:) = list of solar channel indexes

  optp%ntypes = ntypes
  ALLOCATE(optp%data(ntypes))

  DO i = 1, ntypes
    CALL rttov_nullify_optp_data(optp%data(i))
  ENDDO

END SUBROUTINE init_new_optp

!> Interpolate complex refractive index dataset onto a different set of wavelengths
!! @param[in]      n           number of input wavelengths
!! @param[in]      winp        input wavelengths (microns)
!! @param[in]      minp        input complex refractive indices
!! @param[in]      wout        output wavelengths (microns)
!! @param[in,out]  mout        interpolated complex refractive indices
SUBROUTINE interpolate_ref_index(n, winp, minp, wout, mout)

  USE parkind1, ONLY : jpim, jprb

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: winp(:)
  COMPLEX(jprb), INTENT(IN)    :: minp(:)
  REAL(jprb),    INTENT(IN)    :: wout(:)
  COMPLEX(jprb), INTENT(INOUT) :: mout(:)

  INTEGER(jpim) :: iwav
  REAL(jprb)    :: moutreal, moutimag, hr

  DO iwav = 1, SIZE(wout)
    IF (wout(iwav) <= winp(1)) THEN
      mout(iwav) = minp(1)
    ELSE IF (wout(iwav) >= winp(n)) THEN
      mout(iwav) = minp(n)
    ELSE
      CALL inter(n, n, 2, wout(iwav), winp(1:n), REAL(minp(1:n)),  moutreal, hr)
      CALL inter(n, n, 2, wout(iwav), winp(1:n), AIMAG(minp(1:n)), moutimag, hr)
      mout(iwav) = CMPLX(moutreal, moutimag, KIND=jprb)
    ENDIF
  ENDDO

END SUBROUTINE interpolate_ref_index

!> Read OPAC refractive index dataset and interpolate onto given wavelengths
!! @param[out]     err         error status
!! @param[in]      fname       full path to OPAC file
!! @param[in]      nopac       number of OPAC wavelengths
!! @param[in]      nskip       number of lines of header to skip
!! @param[in]      wout        wavelengths on which to output refractive indices (microns)
!! @param[in,out]  mout        complex refractive indices (imaginary parts are positive)
SUBROUTINE read_opac_ref_index(err, fname, nopac, nskip, wout, mout)

  USE parkind1, ONLY : jpim, jprb

  INTEGER(jpim),    INTENT(OUT)   :: err
  CHARACTER(LEN=*), INTENT(IN)    :: fname
  INTEGER(jpim),    INTENT(IN)    :: nopac
  INTEGER(jpim),    INTENT(IN)    :: nskip
  REAL(jprb),       INTENT(IN)    :: wout(:)
  COMPLEX(jprb),    INTENT(INOUT) :: mout(:)

  INTEGER(jpim) :: i, lun
  REAL(jprb)    :: extc, scac, absc, sisca, asym, extnor
  REAL(jprb)    :: wopac(nopac)
  COMPLEX(jprb) :: mopac(nopac)
  LOGICAL       :: lopen

  TRY

  ! Look for a free logical unit (not thread-safe)
  DO lun = 9, 99
    INQUIRE(lun, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO
  IF (lopen) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'Cannot find a free lun')
  ENDIF

  OPEN(lun, file=TRIM(fname), status='old', action='read', iostat=err)
  THROWM(err.NE.0, 'Error opening OPAC file '//TRIM(fname))
  DO i = 1, nskip
    READ(lun, *)
  ENDDO
  DO i = 1, nopac
    READ (lun, '(2x,7e10.3,2e11.3)') wopac(i), extc, scac, absc, &
                                     sisca, asym, extnor, mopac(i)
    mopac(i) = CONJG(mopac(i))
  ENDDO
  CLOSE(lun)

  CALL interpolate_ref_index(nopac, wopac, mopac, wout, mout)

  CATCH
END SUBROUTINE read_opac_ref_index

!> Read Segelstein (1981) refractive index dataset and interpolate onto given wavelengths
!! @param[out]     err         error status
!! @param[in]      fname       full path to file containing segelstein refractive index data
!! @param[in]      wout        wavelengths on which to output refractive indices (microns)
!! @param[in,out]  mout        complex refractive indices (imaginary parts are positive)
SUBROUTINE read_segelstein_ref_index(err, fname, wout, mout)

  USE parkind1, ONLY : jpim, jprb

  INTEGER(jpim),    INTENT(OUT)   :: err
  CHARACTER(LEN=*), INTENT(IN)    :: fname
  REAL(jprb),       INTENT(IN)    :: wout(:)
  COMPLEX(jprb),    INTENT(INOUT) :: mout(:)

  INTEGER(jpim)              :: i, lun, n
  REAL(jprb),    ALLOCATABLE :: ssw(:), ssre(:), ssim(:)
  COMPLEX(jprb), ALLOCATABLE :: ssm(:)
  LOGICAL                    :: lopen

  TRY

  ! Look for a free logical unit (not thread-safe)
  DO lun = 9, 99
    INQUIRE(lun, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO
  IF (lopen) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'Cannot find a free lun')
  ENDIF

  OPEN(lun, file=TRIM(fname), status='old', action='read', iostat=err)
  THROWM(err.NE.0, 'Error opening Segelstein refractive index file '//TRIM(fname))
  READ(lun, *) n

  ALLOCATE(ssw(n), ssre(n), ssim(n), ssm(n))
  DO i = 1, n
    READ (lun, *) ssw(i), ssre(i), ssim(i)
  ENDDO
  CLOSE(lun)

  ssm = CMPLX(ssre, ssim, KIND=jprb)
  ssw = ssw * 1E-3_jprb
  CALL interpolate_ref_index(n, ssw, ssm, wout, mout)
  DEALLOCATE(ssw, ssre, ssim, ssm)

  CATCH
END SUBROUTINE read_segelstein_ref_index

SUBROUTINE make_radius_grid(err, rmin, rmax, rfac, ntot, r)

  ! Compute the radius grid between rmin and rmax

  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : errorstatus_fatal

  INTEGER(jpim), INTENT(OUT) :: err
  REAL(jprb),    INTENT(IN)  :: rmin, rmax, rfac
  INTEGER(jpim), INTENT(OUT) :: ntot
  REAL(jprb),    INTENT(OUT) :: r(:)

  INTEGER(jpim) :: i, maxnradii

  TRY
  maxnradii = SIZE(r)
  r(1) = rmin
  DO i = 2, maxnradii
    r(i) = r(i - 1) * 10**rfac
    IF (r(i) > rmax) THEN
      r(i) = rmax
      EXIT
    ENDIF
  ENDDO
  IF (r(i) < rmax) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Increase maxnradii')
  ENDIF
  ntot = i
  CATCH
END SUBROUTINE make_radius_grid

SUBROUTINE mie_sphere(x, kl, r_angle, m, q_sct, q_ext,f, g, q_bsct, phasef)

!     Based on original by Marco Matricardi, ECMWF, 17-5-2012

  USE parkind1, ONLY : jprb, jpim
  IMPLICIT NONE

!-----Common variables
  REAL(KIND=jprb),    INTENT (IN)  :: x, kl, r_angle             ! The size parameter
  COMPLEX(KIND=jprb), INTENT (IN)  :: m                          ! The refractive index
  REAL(KIND=jprb),    INTENT (OUT) :: q_sct, q_ext, f, g, q_bsct, phasef

!-----Local variables
  REAL(KIND=jprb)    :: n1, n2, d0_re, d0_im, f_n, g_n
  COMPLEX(KIND=jprb) :: mc, mx, c_bsct, ss1, ss2
  INTEGER(KIND=jpim) :: n, n_end

  COMPLEX(KIND=jprb), ALLOCATABLE :: dn(:), wn(:), an(:), bn(:), s1n(:), s2n(:)
  REAL(KIND=jprb),    ALLOCATABLE :: pin(:), taun(:)

!-----How many iterations
  n_end = ANINT(x + 4.05_jprb * x ** (1.0_jprb / 3.0_jprb) + 2.0_jprb)

  IF (n_end < 5) n_end = 5
!      IF (n_end > 1000) n_end = 1000

  ALLOCATE(dn(-1:n_end+1))
  ALLOCATE(wn(-1:n_end+1))
  ALLOCATE(an(-1:n_end+1))
  ALLOCATE(bn(-1:n_end+1))
  ALLOCATE(s1n(-1:n_end+1))
  ALLOCATE(s2n(-1:n_end+1))
  ALLOCATE(pin(-1:n_end+1))
  ALLOCATE(taun(-1:n_end+1))

  mc = CONJG(m)
  mx = mc * x
  n1 = REAL(m) * x
  n2 = AIMAG(m) * x

  d0_re = SIN(n1) * COS(n1) / (SIN(n1) * SIN(n1) + SINH(n2) * SINH(n2))
  d0_im = SINH(n2) * COSH(n2) / (SIN(n1) * SIN(n1) + SINH(n2) * SINH(n2))

!-----Initialize
  pin(0)  = 0._jprb
  pin(1)  = 1._jprb
  pin(2)  = 3._jprb * COS(r_angle)
  taun(0) = 0._jprb
  taun(1) = COS(r_angle)
  taun(2) = 3._jprb * COS(2._jprb * r_angle)

  dn( 0) = CMPLX(d0_re, d0_im, KIND=jprb)
  wn(-1) = CMPLX(COS(x), -1._jprb * SIN(x), KIND=jprb)
  wn( 0) = CMPLX(SIN(x), COS(x), KIND=jprb)

  q_ext  = 0._jprb
  q_sct  = 0._jprb
  q_bsct = 0._jprb
  f      = 0._jprb
  g      = 0._jprb
  c_bsct = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  ss1    = (0._jprb, 0._jprb)
  ss2    = (0._jprb, 0._jprb)

  DO n = 1, n_end + 1

    IF (n > 2) THEN
      pin(n) = ((2_jprb * n - 1._jprb) * COS(r_angle) * pin(n - 1) / (n - 1._jprb)) - &
               n * pin(n - 2) / (n - 1._jprb)
      taun(n) = n * COS(r_angle) * pin(n) - (n + 1._jprb) * pin(n - 1)
    ENDIF

    f_n = 2._jprb * n - 1._jprb
    g_n = 2._jprb * n + 1._jprb

    wn (n) = f_n / x * wn(n - 1) - wn(n - 2)
    dn (n) = 1._jprb / (n / mx - dn(n - 1)) - n / mx

    an (n) = ((dn(n) / mc + n / x) * REAL(wn(n)) - REAL(wn(n - 1))) / &
             ((dn(n) / mc + n / x) * wn(n) - wn(n - 1))
    bn (n) = ((dn(n) * mc + n / x) * REAL(wn(n)) - REAL(wn(n - 1))) / &
             ((dn(n) * mc + n / x) * wn(n) - wn(n - 1))

    q_ext  = q_ext  + g_n * (REAL(an(n)) + REAL(bn(n)))
    q_sct  = q_sct  + g_n * (ABS(an(n)) ** 2 + ABS(bn(n)) ** 2)
    c_bsct = c_bsct + g_n * (-1) ** n * (an(n) - bn(n))

    s1n(n) = (2._jprb * n + 1._jprb) * (an(n) * pin(n) + bn(n) * taun(n)) / &
             (n * (n + 1._jprb))
    s2n(n) = (2._jprb * n + 1._jprb) * (bn(n) * pin(n) + an(n) * taun(n)) / &
             (n * (n + 1._jprb))

    ss1 = ss1 + s1n(n)
    ss2 = ss2 + s2n(n)
  END DO

  phasef = (ABS(ss1) ** 2 + ABS(ss2) ** 2)
  phasef = phasef / (2._jprb * kl ** 2)

  q_ext  = q_ext * 2._jprb / (x * x)
  q_sct  = q_sct * 2._jprb / (x * x)
  IF (q_sct > q_ext) q_sct = q_ext

  DEALLOCATE(an, bn, dn, s1n, s2n, wn)
  DEALLOCATE(pin, taun)

END SUBROUTINE mie_sphere

SUBROUTINE calc_mie_ext_sca(ntot, r, n, wavl, refindex, ecoef, scoef, nthreads)

  ! Compute Mie extinction and scattering coefficients

  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : pi

  INTEGER(jpim), INTENT(IN)   :: ntot
  REAL(jprb),    INTENT(IN)   :: r(ntot), n(ntot)
  REAL(jprb),    INTENT(IN)   :: wavl
  COMPLEX(jprb), INTENT(IN)   :: refindex
  REAL(jprb),    INTENT(OUT)  :: ecoef, scoef
  INTEGER(jpim), INTENT(IN)   :: nthreads

  INTEGER(jpim) :: i, ifail
  REAL(jprb)    :: x, kl, r_angle, q_sct, q_ext, f, g, q_bsct, phasef, error
  REAL(jprb)    :: yarr(ntot), zarr(ntot)

#ifndef RTTOV_NAG53
!$OMP PARALLEL DO NUM_THREADS(nthreads) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC) &
!$OMP             SHARED(ntot, r, wavl, refindex, n, yarr, zarr)
#endif
  DO i = 1, ntot
    x  = r(i) * 2._jprb * pi / wavl
    kl = 2._jprb * pi / wavl
    r_angle = 0._jprb

    CALL mie_sphere(x, kl, r_angle, refindex, q_sct, q_ext, f, g, q_bsct, phasef)

    yarr(i) = q_ext * n(i) * r(i)**2 * pi
    zarr(i) = q_sct * n(i) * r(i)**2 * pi
  ENDDO
#ifndef RTTOV_NAG53
!$OMP END PARALLEL DO
#endif

  CALL integrate(ntot, r(1:ntot), yarr(1:ntot), ecoef, error, ifail)
  CALL integrate(ntot, r(1:ntot), zarr(1:ntot), scoef, error, ifail)

END SUBROUTINE calc_mie_ext_sca

SUBROUTINE calc_mie_phasefn(nphangle, phangle, ntot, r, n, wavl, refindex, scoef, phasefn, nthreads)

  ! Compute Mie phase function

  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : pi, deg2rad

  INTEGER(jpim), INTENT(IN)   :: nphangle
  REAL(jprb),    INTENT(IN)   :: phangle(nphangle)
  INTEGER(jpim), INTENT(IN)   :: ntot
  REAL(jprb),    INTENT(IN)   :: r(ntot), n(ntot)
  REAL(jprb),    INTENT(IN)   :: wavl
  COMPLEX(jprb), INTENT(IN)   :: refindex
  REAL(jprb),    INTENT(IN)   :: scoef
  REAL(jprb),    INTENT(OUT)  :: phasefn(nphangle)
  INTEGER(jpim), INTENT(IN)   :: nthreads

  INTEGER(jpim) :: iang, i, ifail
  REAL(jprb)    :: x, kl, r_angle, q_sct, q_ext, f, g, q_bsct, phasef, phfnc, error
  REAL(jprb)    :: parr(ntot)

#ifndef RTTOV_NAG53
!$OMP PARALLEL DO NUM_THREADS(nthreads) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC) &
!$OMP             SHARED(nphangle, phangle, ntot, r, wavl, refindex, n, phasefn, scoef)
#endif
  DO iang = 1, nphangle

    r_angle = phangle(iang) * deg2rad
    DO i = 1, ntot
      x  = r(i) * 2._jprb * pi / wavl
      kl = 2._jprb * pi / wavl
      CALL mie_sphere(x, kl, r_angle, refindex, q_sct, q_ext, f, g, q_bsct, phasef)
      parr(i) = phasef * n(i)
    ENDDO

    CALL integrate(ntot, r, parr, phfnc, error, ifail)
    phasefn(iang) = 4._jprb * pi * phfnc / scoef

  ENDDO
#ifndef RTTOV_NAG53
!$OMP END PARALLEL DO
#endif
END SUBROUTINE calc_mie_phasefn

SUBROUTINE calc_mass(ntot, r, n, rho, mass)

  ! Compute integrated mass

  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : pi

  INTEGER(jpim), INTENT(IN)   :: ntot
  REAL(jprb),    INTENT(IN)   :: r(ntot), n(ntot)
  REAL(jprb),    INTENT(IN)   :: rho
  REAL(jprb),    INTENT(OUT)  :: mass

  INTEGER(jpim) :: ifail
  REAL(jprb)    :: error, varr(ntot)

  varr = n * r**3
  CALL integrate(ntot, r, varr, mass, error, ifail)
  mass = mass * rho * 4._jprb * pi / 3._jprb

END SUBROUTINE calc_mass


SUBROUTINE lognorm(ntot, rarr, asigma, armod, sqrt2pi, n)

  USE parkind1, ONLY : jprb, jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)  :: ntot
  REAL(KIND=jprb),    INTENT(IN)  :: rarr(1:ntot)
  REAL(KIND=jprb),    INTENT(IN)  :: asigma
  REAL(KIND=jprb),    INTENT(IN)  :: armod
  REAL(KIND=jprb),    INTENT(IN)  :: sqrt2pi
  REAL(KIND=jprb),    INTENT(OUT) :: n(1:ntot)

  INTEGER(KIND=jpim) :: i
  REAL(KIND=jprb)    :: r

  DO i = 1, ntot
    r = rarr(i)
    n(i) = (1._jprb / (sqrt2pi * r * LOG10(asigma) * LOG(10._jprb))) * &
           EXP(-0.5_jprb * ((LOG10(r) - LOG10(armod)) / LOG10(asigma))**2)
  ENDDO

  RETURN
END SUBROUTINE lognorm

SUBROUTINE gammadist(ntot, rarr, aacoef, aalpha, abcoef, agamma, n)

  USE parkind1, ONLY : jprb, jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)  :: ntot
  REAL(KIND=jprb),    INTENT(IN)  :: rarr(1:ntot)
  REAL(KIND=jprb),    INTENT(IN)  :: aacoef
  REAL(KIND=jprb),    INTENT(IN)  :: aalpha
  REAL(KIND=jprb),    INTENT(IN)  :: abcoef
  REAL(KIND=jprb),    INTENT(IN)  :: agamma
  REAL(KIND=jprb),    INTENT(OUT) :: n(1:ntot)

  INTEGER(KIND=jpim) :: i
  REAL(KIND=jprb)    :: r

  DO i = 1, ntot
    r = rarr(i)
    n(i) = aacoef * (r ** aalpha) * EXP(-abcoef * (r ** agamma))
  ENDDO

  RETURN
END SUBROUTINE gammadist

END MODULE rttov_mie_params_mod
