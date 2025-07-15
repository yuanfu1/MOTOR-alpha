! Description:
!> @file
!!   Various routines related to visible/IR scattering calculations
!
!> @brief
!!   Various routines related to visible/IR scattering calculations
!!
!! @details
!!   This module contains subroutines that are used internally by RTTOV
!!   for visible and IR scattering calculations. These include:
!!   - calc_rel_hum(+TL/AD/K): calculation of relative humidity used
!!        to interpolate certain aerosol particle optical properties
!!   - calc_leg_poly: evaluate Legendre polynomials given expansion coefs
!!   - gauss_quad: calculate Gaussian quadrature
!!   - calc_legendre_coef_gauss+tl/ad: calculate Legendre coefficients using
!!        Gaussian quadrature
!!   - normalise+tl/ad: normalise a phase function calculated by Gaussian quadrature
!!   - spline_interp/spline/splint+tl/ad: spline interpolation
!!
!!   The following mathematical subroutines/functions are used in aerosol/cloud
!!   optical property calculations, and also elsewhere in RTTOV:
!!   - integrate: implementation of NAG library subroutine d01GAF
!!   - inter/xl: subroutines for interpolation
!!
!!   The following subroutines/functions were taken from external sources and
!!   are used in the DOM algorithm:
!!   - asymtx: from DISORT, for a real asymmetric matrix this calculates
!!        eigenvalues and eigenvectors which are known a priori to be real
!!   - asymtx_tl/ad: from CRTM, TL/AD of ASYMTX
!!   - matinv: from CRTM, matrix inversion
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE rttov_scattering_mod

#include "throw.h"

  USE rttov_const, ONLY : pi
  USE parkind1, ONLY : jprb

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_rel_hum, calc_rel_hum_tl, calc_rel_hum_ad, calc_rel_hum_k, &
            calc_leg_poly, gauss_quad, calc_legendre_coef_gauss, &
            calc_legendre_coef_gauss_tl, calc_legendre_coef_gauss_ad, &
            normalise, normalise_tl, normalise_ad, &
            spline_interp, spline_interp_tl, spline_interp_ad, &
            integrate, inter, asymtx, asymtx_tl, asymtx_ad

  ! This comes from the LAPACK D1MACH subroutine and is used by ASYMTX
  ! Need to convert RADIX to REAL hence "1.0 *"
  REAL(jprb), PARAMETER :: d1mach4 = (1.0_jprb * RADIX(1._jprb)) ** (1 - DIGITS(1._jprb))

  ! Used by asymtx_tl/ad
  REAL(jprb), PARAMETER :: eigen_threshold = 1.E-20_jprb

#include "rttov_errorreport.interface"

CONTAINS

!> Calculate relative humidity profiles - used to interpolate aerosol
!! optical properties to appropriate relative humidity value.
!! @param[in]     profiles        profiles structure
!! @param[in]     profiles_dry    profiles structure containing gas profiles in units of ppmv dry
!! @param[in,out] aux             auxiliary profile data structure
SUBROUTINE calc_rel_hum(profiles, profiles_dry, aux)

  USE rttov_types, ONLY :  rttov_profile, rttov_profile_aux
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : e00, t00, ti

  IMPLICIT NONE

  TYPE(rttov_profile),     INTENT(IN)              :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)              :: profiles_dry(:)
  TYPE(rttov_profile_aux), INTENT(INOUT)           :: aux

  INTEGER(jpim) :: nprofiles, nlayers, j, lay, lev
  ! ------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

  DO j = 1, nprofiles
    DO lay = 1, nlayers
      lev = lay + 1
      aux%tave(lay,j)     = 0.5_jprb * (profiles(j)%t(lev-1) + profiles(j)%t(lev))
      aux%wmixave(lay,j)  = 0.5_jprb * (profiles_dry(j)%q(lev-1) + profiles_dry(j)%q(lev))
      aux%xpresave(lay,j) = 0.5_jprb * (profiles(j)%p(lev-1) + profiles(j)%p(lev))

      ! Saturated vapour pressure
      aux%esw(lay,j) = e00 * EXP(17.502_jprb * (aux%tave(lay,j) - t00) / &
                                   (aux%tave(lay,j) - 32.19_jprb))
      aux%esi(lay,j) = e00 * EXP(22.587_jprb * (aux%tave(lay,j) - t00) / &
                                   (aux%tave(lay,j) + 0.7_jprb))

      IF (aux%tave(lay,j) > t00) THEN
        ! Water phase
        aux%ppv(lay,j) = aux%esw(lay,j)
      ELSE IF (aux%tave(lay,j) > ti .AND. aux%tave(lay,j) <= t00) THEN
        ! Mixed phase
        aux%ppv(lay,j) = aux%esi(lay,j) +     &
          (aux%esw(lay,j) - aux%esi(lay,j)) * &
          ((aux%tave(lay,j) - ti) / (t00 - ti)) ** 2
      ELSE !IF (aux%tave(lay,j) <= ti) THEN
        ! Ice phase
        aux%ppv(lay,j) = aux%esi(lay,j)
      ENDIF
      aux%ppv(lay,j) = aux%ppv(lay,j) / 100._jprb

      ! Layer average relative humidity
      aux%relhum(lay,j) = 100._jprb * aux%wmixave(lay,j) *     &
          1.E-6_jprb * 0.622_jprb * aux%xpresave(lay,j) /  &
          (aux%ppv(lay,j) * 0.622_jprb * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb))
    ENDDO ! layers
  ENDDO ! profiles
END SUBROUTINE calc_rel_hum

!> Calculate relative humidity profiles - TL.
!! @param[in]     opts            options structure
!! @param[in]     profiles        profiles structure
!! @param[in]     profiles_tl     profiles structure perturbations
!! @param[in]     profiles_dry_tl profiles structure containing gas perturbations in units of ppmv dry
!! @param[in]     aux             auxiliary profile data structure
!! @param[in,out] aux_tl          perturbations in auxiliary profile data
SUBROUTINE calc_rel_hum_tl(opts, profiles, profiles_tl, profiles_dry_tl, aux, aux_tl)

  USE rttov_types, ONLY :  rttov_options, rttov_profile, rttov_profile_aux
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : t00, ti

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_dry_tl(:)
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_tl

  INTEGER(jpim) :: nprofiles, nlayers, j, lay, lev
  ! ------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

  DO j = 1, nprofiles
    DO lay = 1, nlayers
      lev = lay + 1
      aux_tl%tave(lay,j)     = 0.5_jprb * (profiles_tl(j)%t(lev-1) + profiles_tl(j)%t(lev))
      aux_tl%wmixave(lay,j)  = 0.5_jprb * (profiles_dry_tl(j)%q(lev-1) + profiles_dry_tl(j)%q(lev))
      IF (opts%interpolation%lgradp) THEN
        aux_tl%xpresave(lay,j) = 0.5_jprb * (profiles_tl(j)%p(lev-1) + profiles_tl(j)%p(lev))
      ELSE
        aux_tl%xpresave(lay,j) = 0._jprb
      ENDIF

      ! Saturated vapour pressure
      aux_tl%esw(lay,j) = aux_tl%tave(lay,j) * aux%esw(lay,j) * &
                            17.502_jprb * (t00 - 32.19_jprb) / (aux%tave(lay,j) - 32.19_jprb) ** 2_jpim
      aux_tl%esi(lay,j) = aux_tl%tave(lay,j) * aux%esi(lay,j) * &
                            22.587_jprb * (t00 + 0.7_jprb) / (aux%tave(lay,j) + 0.7_jprb) ** 2_jpim

      IF (aux%tave(lay,j) > t00) THEN
        ! Water phase
        aux_tl%ppv(lay,j) = aux_tl%esw(lay,j)
      ELSE IF (aux%tave(lay,j) > ti .AND. aux%tave(lay,j) <= t00) THEN
        ! Mixed phase
        aux_tl%ppv(lay,j) = aux_tl%esi(lay,j) + &
          ((aux_tl%esw(lay,j) - aux_tl%esi(lay,j)) * (aux%tave(lay,j) - ti) ** 2_jpim + &
           (aux%esw(lay,j) - aux%esi(lay,j)) * &
           2_jpim * aux_tl%tave(lay,j) * (aux%tave(lay,j) - ti)) / &
             ((t00 - ti) ** 2_jpim)
      ELSE !IF (aux%tave(lay,j) <= ti) THEN
        ! Ice phase
        aux_tl%ppv(lay,j) = aux_tl%esi(lay,j)
      ENDIF
      aux_tl%ppv(lay,j) = aux_tl%ppv(lay,j) / 100._jprb

      ! Layer average relative humidity
      ! NB aux_tl%xpresave is set to zero above if .not. lgradp
      aux_tl%relhum(lay,j) = 100._jprb * 1.E-6_jprb * 0.622_jprb * &
          ((aux_tl%wmixave(lay,j) * aux%xpresave(lay,j) + aux%wmixave(lay,j) * aux_tl%xpresave(lay,j)) /  &
           (aux%ppv(lay,j) * 0.622_jprb * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb)) - &
           (aux%wmixave(lay,j) * aux%xpresave(lay,j)) * &
           0.622_jprb * (aux_tl%ppv(lay,j) * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb) + &
                         aux%ppv(lay,j) * aux_tl%wmixave(lay,j) * 1.E-6_jprb) / &
           (aux%ppv(lay,j) * 0.622_jprb * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb)) ** 2_jpim)
    ENDDO ! layers
  ENDDO ! profiles

END SUBROUTINE calc_rel_hum_tl

!> Calculate relative humidity profiles - AD.
!! @param[in]     opts            options structure
!! @param[in]     profiles        profiles structure
!! @param[in,out] profiles_ad     profiles structure increments
!! @param[in,out] profiles_dry_ad profiles structure containing gas increments in units of ppmv dry
!! @param[in]     aux             auxiliary profile data structure
!! @param[in,out] aux_ad          increments of auxiliary profile data
SUBROUTINE calc_rel_hum_ad(opts, profiles, profiles_ad, profiles_dry_ad, aux, aux_ad)

  USE rttov_types, ONLY : rttov_options, rttov_profile, rttov_profile_aux
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : t00, ti

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_ad(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_dry_ad(:)
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_ad

  INTEGER(jpim) :: nprofiles, nlayers, j, lay, lev
  ! ------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

  DO j = 1, nprofiles
    DO lay = nlayers, 1, -1
      lev = lay + 1
      aux_ad%wmixave(lay,j) = aux_ad%wmixave(lay,j) +  &
          aux_ad%relhum(lay,j) * 100._jprb * 1.E-6_jprb * aux%xpresave(lay,j) / &
          (aux%ppv(lay,j) * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb))
      aux_ad%wmixave(lay,j) = aux_ad%wmixave(lay,j) - &
          aux_ad%relhum(lay,j) * 100._jprb * 1.E-6_jprb ** 2 * aux%xpresave(lay,j) * &
          aux%wmixave(lay,j) * aux%ppv(lay,j) / &
          (aux%ppv(lay,j) * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb)) ** 2
      aux_ad%ppv(lay,j) = aux_ad%ppv(lay,j) - &
          aux_ad%relhum(lay,j) * 100._jprb * aux%wmixave(lay,j) * &
          1.E-6_jprb * aux%xpresave(lay,j) * &
          (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb) / &
          (aux%ppv(lay,j) * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb)) ** 2
      IF (opts%interpolation%lgradp) aux_ad%xpresave(lay,j) = aux_ad%xpresave(lay,j) + &
            100._jprb * aux%wmixave(lay,j) * 1.e-6_jprb * aux_ad%relhum(lay,j) / &
            (aux%ppv(lay,j) * (1._jprb + aux%wmixave(lay,j) * 1.E-6_jprb))

      aux_ad%ppv(lay,j) = aux_ad%ppv(lay,j) / 100._jprb
      IF (aux%tave(lay,j) > t00) THEN
        aux_ad%esw(lay,j) = aux_ad%esw(lay,j) + aux_ad%ppv(lay,j)
      ELSE IF (aux%tave(lay,j) > ti .AND. aux%tave(lay,j) <= t00) THEN
        aux_ad%esi(lay,j) = aux_ad%esi(lay,j) + aux_ad%ppv(lay,j)
        aux_ad%esw(lay,j) = &
          & aux_ad%esw(lay,j) + aux_ad%ppv(lay,j) * ((aux%tave(lay,j) - ti) / (t00 - ti)) ** 2
        aux_ad%esi(lay,j) = &
          & aux_ad%esi(lay,j) - aux_ad%ppv(lay,j) * ((aux%tave(lay,j) - ti) / (t00 - ti)) ** 2
        aux_ad%tave(lay,j) = aux_ad%tave(lay,j) + &
          & aux_ad%ppv(lay,j) * (aux%esw(lay,j) - aux%esi(lay,j)) * 2 * &
          & ((aux%tave(lay,j) - ti) / (t00 - ti) ** 2)
      ELSE IF (aux%tave(lay,j) <= ti) THEN
        aux_ad%esi(lay,j) = aux_ad%esi(lay,j) + aux_ad%ppv(lay,j)
      ENDIF
      aux_ad%tave(lay,j) = aux_ad%tave(lay,j) + &
        & aux_ad%esw(lay,j) * aux%esw(lay,j) * 17.502_jprb * (t00 - 32.19_jprb) / &
        & (aux%tave(lay,j) - 32.19_jprb) ** 2
      aux_ad%tave(lay,j) = aux_ad%tave(lay,j) +      &
        & aux_ad%esi(lay,j) * aux%esi(lay,j) * 22.587_jprb * (t00 + 0.7_jprb) / &
        & (aux%tave(lay,j) + 0.7_jprb) ** 2

      profiles_dry_ad(j)%q(lev-1) = profiles_dry_ad(j)%q(lev-1) + 0.5_jprb * aux_ad%wmixave(lay,j)
      profiles_dry_ad(j)%q(lev)   = profiles_dry_ad(j)%q(lev) + 0.5_jprb * aux_ad%wmixave(lay,j)
      profiles_ad(j)%t(lev-1)     = profiles_ad(j)%t(lev-1) + 0.5_jprb * aux_ad%tave(lay,j)
      profiles_ad(j)%t(lev)       = profiles_ad(j)%t(lev) + 0.5_jprb * aux_ad%tave(lay,j)
      IF (opts%interpolation%lgradp) THEN
        profiles_ad(j)%p(lev-1) = profiles_ad(j)%p(lev-1) + 0.5_jprb * aux_ad%xpresave(lay, j)
        profiles_ad(j)%p(lev)   = profiles_ad(j)%p(lev) + 0.5_jprb * aux_ad%xpresave(lay, j)
      ENDIF
    ENDDO ! layers
  ENDDO ! profiles

END SUBROUTINE calc_rel_hum_ad

!> Calculate relative humidity profiles - K.
!! @param[in]     opts            options structure
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     profiles        profiles structure
!! @param[in,out] profiles_k      profiles structure increments
!! @param[in,out] profiles_dry_k  profiles structure containing gas increments in units of ppmv dry
!! @param[in]     aux             auxiliary profile data structure
!! @param[in,out] aux_k           increments of auxiliary profile data
SUBROUTINE calc_rel_hum_k(opts, chanprof, profiles, profiles_k, profiles_dry_k, aux, aux_k)

  USE rttov_types, ONLY : rttov_options, rttov_chanprof, rttov_profile, rttov_profile_aux
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : t00, ti

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_dry_k(:)
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_k

  INTEGER(jpim) :: nchanprof, nlayers, j, lay, lev, prof
  ! ------------------------------------------------------------------------

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    DO lay = nlayers, 1, -1
      lev = lay + 1
      aux_k%wmixave(lay,j) = aux_k%wmixave(lay,j) +  &
          aux_k%relhum(lay,j) * 100._jprb * 1.E-6_jprb * aux%xpresave(lay,prof) / &
          (aux%ppv(lay,prof) * (1._jprb + aux%wmixave(lay,prof) * 1.E-6_jprb))
      aux_k%wmixave(lay,j) = aux_k%wmixave(lay,j) - &
          aux_k%relhum(lay,j) * 100._jprb * 1.E-6_jprb ** 2 * aux%xpresave(lay,prof) * &
          aux%wmixave(lay,prof) * aux%ppv(lay,prof) / &
          (aux%ppv(lay,prof) * (1._jprb + aux%wmixave(lay,prof) * 1.E-6_jprb)) ** 2
      aux_k%ppv(lay,j) = aux_k%ppv(lay,j) - &
          aux_k%relhum(lay,j) * 100._jprb * aux%wmixave(lay,prof) * &
          1.E-6_jprb * aux%xpresave(lay,prof) * &
          (1._jprb + aux%wmixave(lay,prof) * 1.E-6_jprb) / &
          (aux%ppv(lay,prof) * (1._jprb + aux%wmixave(lay,prof) * 1.E-6_jprb)) ** 2
      IF (opts%interpolation%lgradp) aux_k%xpresave(lay,j) = aux_k%xpresave(lay,j) + &
            100._jprb * aux%wmixave(lay,prof) * 1.e-6_jprb * aux_k%relhum(lay,j) / &
            (aux%ppv(lay,prof) * (1._jprb + aux%wmixave(lay,prof) * 1.E-6_jprb))

      aux_k%ppv(lay,j) = aux_k%ppv(lay,j) / 100._jprb
      IF (aux%tave(lay,prof) > t00) THEN
        aux_k%esw(lay,j) = aux_k%esw(lay,j) + aux_k%ppv(lay,j)
      ELSE IF (aux%tave(lay,prof) > ti .AND. aux%tave(lay,prof) <= t00) THEN
        aux_k%esi(lay,j) = aux_k%esi(lay,j) + aux_k%ppv(lay,j)
        aux_k%esw(lay,j) = &
          & aux_k%esw(lay,j) + aux_k%ppv(lay,j) * ((aux%tave(lay,prof) - ti) / (t00 - ti)) ** 2
        aux_k%esi(lay,j) = &
          & aux_k%esi(lay,j) - aux_k%ppv(lay,j) * ((aux%tave(lay,prof) - ti) / (t00 - ti)) ** 2
        aux_k%tave(lay,j) = aux_k%tave(lay,j) + &
          & aux_k%ppv(lay,j) * (aux%esw(lay,prof) - aux%esi(lay,prof)) * 2 * &
          & ((aux%tave(lay,prof) - ti) / (t00 - ti) ** 2)
      ELSE IF (aux%tave(lay,prof) <= ti) THEN
        aux_k%esi(lay,j) = aux_k%esi(lay,j) + aux_k%ppv(lay,j)
      ENDIF
      aux_k%tave(lay,j) = aux_k%tave(lay,j) + &
        & aux_k%esw(lay,j) * aux%esw(lay,prof) * 17.502_jprb * (t00 - 32.19_jprb) / &
        & (aux%tave(lay,prof) - 32.19_jprb) ** 2
      aux_k%tave(lay,j) = aux_k%tave(lay,j) +      &
        & aux_k%esi(lay,j) * aux%esi(lay,prof) * 22.587_jprb * (t00 + 0.7_jprb) / &
        & (aux%tave(lay,prof) + 0.7_jprb) ** 2

      profiles_dry_k(j)%q(lev-1) = profiles_dry_k(j)%q(lev-1) + 0.5_jprb * aux_k%wmixave(lay,j)
      profiles_dry_k(j)%q(lev)   = profiles_dry_k(j)%q(lev) + 0.5_jprb * aux_k%wmixave(lay,j)
      profiles_k(j)%t(lev-1)     = profiles_k(j)%t(lev-1) + 0.5_jprb * aux_k%tave(lay,j)
      profiles_k(j)%t(lev)       = profiles_k(j)%t(lev) + 0.5_jprb * aux_k%tave(lay,j)
      IF (opts%interpolation%lgradp) THEN
        profiles_k(j)%p(lev-1) = profiles_k(j)%p(lev-1) + 0.5_jprb * aux_k%xpresave(lay, j)
        profiles_k(j)%p(lev)   = profiles_k(j)%p(lev) + 0.5_jprb * aux_k%xpresave(lay, j)
      ENDIF
    ENDDO ! layers
  ENDDO ! profiles

END SUBROUTINE calc_rel_hum_k


!> Evaluate associated Legendre polynomials of order m, l (l = m...nmom) for
!! m = 0...naz at values x(:) using equations 58a, 58d and 58h from
!! Stamnes et al 2000 DISORT, a General-Purpose Fortran Program for Discrete-
!! Ordinate-Method Radiative Transfer in Scattering and Emitting Layered
!! Media: Documentation of Methodology (version 1.1). The outputs are
!! normalised Legendre polynomials.
!! @param[in]       nmom      number of Legendre moments (l = m...nmom)
!! @param[in]       naz       number of azimuthal terms  (m = 0...naz)
!! @param[in]       x         array of values at which to evaluate Legendre polynomials
!! @param[in,out]   legpolyx  output normalised Legendre polynomials at x(:) P_l^m
SUBROUTINE calc_leg_poly(nmom, naz, x, legpolyx)

  USE parkind1, ONLY : jprb, jpim

  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: nmom
  INTEGER(jpim), INTENT(IN)    :: naz
  REAL(jprb),    INTENT(IN)    :: x(:)
  REAL(jprb),    INTENT(INOUT) :: legpolyx(0:nmom,SIZE(x),0:naz)

  INTEGER(jpim) :: l, m
  REAL(jprb)    :: legpoly(SIZE(x),3)

  DO m = 0, naz
    legpolyx(0:m-1,:,m) = 0._jprb
    DO l = m, nmom
      IF (l == m) THEN
        IF (l == 0) THEN
          legpoly(:,1) = 1._jprb
        ELSE
          legpoly(:,1) = -SQRT((2 * m - 1) * (1 - x(:)**2) / (2 * m)) * legpolyx(m-1,:,m-1)
        ENDIF
      ELSE IF (l == m + 1) THEN
        legpoly(:,2) = legpoly(:,1)
        legpoly(:,1) = x(:) * SQRT(2._jprb * m + 1._jprb) * legpoly(:,2)
      ELSE
        legpoly(:,3) = legpoly(:,2)
        legpoly(:,2) = legpoly(:,1)
        legpoly(:,1) = ((2._jprb * l - 1._jprb) * x(:) * legpoly(:,2) - &
                        SQRT((l + m - 1._jprb) * (l - m - 1)) * legpoly(:,3)) / &
                        SQRT(REAL((l - m) * (l + m), jprb))
      ENDIF
      legpolyx(l,:,m) = legpoly(:,1)
    ENDDO
  ENDDO

END SUBROUTINE calc_leg_poly

!> Calculate a Gaussian quadrature on the interval [x1, x2]
!! @param[in]     x1      interval lower bound
!! @param[in]     x2      interval upper bound
!! @param[in,out] x       calculated quadrature ordinates
!! @param[in,out] w       calculated quadrature weights
SUBROUTINE gauss_quad(x1, x2, x, w)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  REAL(jprb),    INTENT(IN)    :: x1, x2
  REAL(jprb),    INTENT(INOUT) :: x(:), w(:)

  INTEGER(jpim)             :: i, j, m, n

  REAL(jprb), PARAMETER     :: EPS = 3.e-14
  REAL(jprb)                :: p1, p2, p3, pp, xl, xm, z, z1

  n = SIZE(x)
  m = (n + 1) / 2
  xm = 0.5_jprb * (x2 + x1)
  xl = 0.5_jprb * (x2 - x1)
  DO i = 1, m
    z = COS(pi * (i - .25_jprb) / (n + .5_jprb))
    z1 = -999._jprb
    DO WHILE (ABS(z - z1) > EPS)
      p1 = 1._jprb
      p2 = 0._jprb
      DO j = 1, n
        p3 = p2
        p2 = p1
        p1 = ((2._jprb * j - 1._jprb) * z * p2 - (j - 1._jprb) * p3) / j
      ENDDO
      pp = n * (z * p1 - p2) / (z * z - 1._jprb)
      z1 = z
      z = z1 - p1 / pp
    ENDDO
!     IF (ABS(z - z1) > EPS) GOTO 1
    x(i) = xm - xl * z
    x(n+1-i) = xm + xl * z
    w(i) = 2._jprb * xl / ((1._jprb - z * z) * pp * pp)
    w(n+1-i) = w(i)
  ENDDO

END SUBROUTINE gauss_quad

!> Calculates Legendre coefficients for given phase function f using a
!! precomputed Gaussian quadrature u with weights w. At least minnmom
!! and at most maxnmom coefficients will be calculated. This calculates
!! fixednmom coefficients if supplied, otherwise the expansion stops when
!! the magnitude of the coefficients drops below 1E-6. In any case
!! minnmom/maxnmom are always respected.
!! @param[in]     u           Gaussian quadrature ordinates
!! @param[in]     w           Gaussian quadrature weights
!! @param[in]     f           phase function on quadrature ordinates
!! @param[in]     minnmom     minimum number of coefficients to compute (0...minnmom)
!! @param[in]     maxnmom     maximum number of coefficients to compute (0...maxnmom)
!! @param[in,out] nmom        number of coefficients computed (0...nmom)
!! @param[in,out] lcoef       calculated Legendre coefficients
!! @param[in]     fixednmom   number of coefficients to compute (0...fixednmom), optional
!! @param[in,out] f_recomp    phase function recomputed from Legendre expansion, optional
SUBROUTINE calc_legendre_coef_gauss(u, w, f, minnmom, maxnmom, nmom, lcoef, fixednmom, f_recomp)

  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  REAL(jprb),    INTENT(IN)              :: u(:)
  REAL(jprb),    INTENT(IN)              :: w(SIZE(u))
  REAL(jprb),    INTENT(IN)              :: f(SIZE(u))
  INTEGER(jpim), INTENT(IN)              :: minnmom
  INTEGER(jpim), INTENT(IN)              :: maxnmom
  INTEGER(jpim), INTENT(INOUT)           :: nmom
  REAL(jprb),    INTENT(INOUT)           :: lcoef(0:)
  INTEGER(jpim), INTENT(IN),    OPTIONAL :: fixednmom
  REAL(jprb),    INTENT(INOUT), OPTIONAL :: f_recomp(SIZE(u))

  INTEGER(jpim) :: ngauss, j, k, thismaxnmom
  REAL(jprb)    :: y1(SIZE(u))
  REAL(jprb)    :: p00(SIZE(u),-1:maxnmom+1)

  lcoef(:) = 0._jprb
  thismaxnmom = MIN(SIZE(lcoef)-1, maxnmom)
  IF (thismaxnmom < minnmom) THEN
    PRINT *, 'lcoef is too small or minnmom is too large'
    nmom = 0
    IF (PRESENT(f_recomp)) f_recomp = 1.E38_jprb
    RETURN
  ENDIF

  !**************************************************************************
  !         Compute lcoef for F11 and F44
  !  Rq: P00(l) stand for the generalized spherical fonction
  !       l
  !      P (u) see de Haan et al., 1987 (Astronomy and astrophysics)
  !       m,m
  !  From B.B. verifier avec l'article ci-dessus (2/2/2005)
  !**************************************************************************
  ngauss = SIZE(u)
  p00(:,-1) = 0._jprb
  p00(:,0) = 1._jprb
  y1 = f * w
  IF (PRESENT(f_recomp)) f_recomp = 0._jprb

  DO k = 0, thismaxnmom
    DO j = 1, ngauss
      p00(j,k+1) = ((2._jprb * k + 1._jprb) * u(j) * p00(j,k) - k * p00(j,k-1)) / (k + 1._jprb)
      lcoef(k) = lcoef(k) + y1(j) * p00(j,k)
    ENDDO
    lcoef(k) = (2._jprb * k + 1._jprb) * lcoef(k) / 2._jprb

    IF (PRESENT(f_recomp)) f_recomp = f_recomp + lcoef(k) * p00(:,k)

    nmom = k

    IF (PRESENT(fixednmom)) THEN
      IF (k == MAX(minnmom, fixednmom)) EXIT
    ELSE
      ! Stop if the coefficients have grown sufficiently small
      IF (lcoef(k) < 1.E-6_jprb .AND. k >= minnmom) EXIT
    ENDIF
  ENDDO

END SUBROUTINE calc_legendre_coef_gauss

!> Calculate Legendre coefficients for given phase function - TL
!! This call carries out both the forward and TL model calculations.
!! @param[in]     u           Gaussian quadrature ordinates
!! @param[in]     w           Gaussian quadrature weights
!! @param[in]     f           phase function on quadrature ordinates
!! @param[in]     f_tl        phase function perturbations on quadrature ordinates
!! @param[in]     minnmom     minimum number of coefficients to compute (0...minnmom)
!! @param[in]     maxnmom     maximum number of coefficients to compute (0...maxnmom)
!! @param[in,out] nmom        number of coefficients computed (0...nmom)
!! @param[in,out] lcoef       calculated Legendre coefficients
!! @param[in,out] lcoef_tl    calculated Legendre coefficient perturbations
!! @param[in]     fixednmom   number of coefficients to compute (0...fixednmom), optional
SUBROUTINE calc_legendre_coef_gauss_tl(u, w, f, f_tl, minnmom, maxnmom, nmom, lcoef, lcoef_tl, fixednmom)

  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  REAL(jprb),    INTENT(IN)              :: u(:)
  REAL(jprb),    INTENT(IN)              :: w(SIZE(u))
  REAL(jprb),    INTENT(IN)              :: f(SIZE(u))
  REAL(jprb),    INTENT(IN)              :: f_tl(SIZE(u))
  INTEGER(jpim), INTENT(IN)              :: minnmom
  INTEGER(jpim), INTENT(IN)              :: maxnmom
  INTEGER(jpim), INTENT(INOUT)           :: nmom
  REAL(jprb),    INTENT(INOUT)           :: lcoef(0:)
  REAL(jprb),    INTENT(INOUT)           :: lcoef_tl(0:)
  INTEGER(jpim), INTENT(IN),    OPTIONAL :: fixednmom

  INTEGER(jpim) :: ngauss, j, k, thismaxnmom
  REAL(jprb)    :: y1(SIZE(u)), y1_tl(SIZE(u))
  REAL(jprb)    :: p00(SIZE(u),-1:maxnmom+1)
  REAL(jprb)    :: p00_tl(SIZE(u),-1:maxnmom+1)

  lcoef(:) = 0._jprb
  lcoef_tl(:) = 0._jprb
  thismaxnmom = MIN(SIZE(lcoef)-1, maxnmom)
  IF (thismaxnmom < minnmom) THEN
    PRINT *, 'lcoef is too small or minnmom is too large'
    nmom = 0
    RETURN
  ENDIF

  ngauss = SIZE(u)
  p00(:,-1) = 0._jprb
  p00(:,0) = 1._jprb
  p00_tl(:,-1:0) = 0._jprb
  y1 = f * w
  y1_tl = f_tl * w

  DO k = 0, thismaxnmom
    DO j = 1, ngauss
      p00(j,k+1) = ((2._jprb * k + 1._jprb) * u(j) * p00(j,k) - k * p00(j,k-1)) / (k + 1._jprb)
      p00_tl(j,k+1) = ((2._jprb * k + 1._jprb) * u(j) * p00_tl(j,k) - k * p00_tl(j,k-1)) / (k + 1._jprb)
      lcoef(k) = lcoef(k) + y1(j) * p00(j,k)
      lcoef_tl(k) = lcoef_tl(k) + y1_tl(j) * p00(j,k) + y1(j) * p00_tl(j,k)
    ENDDO
    lcoef(k) = (2._jprb * k + 1._jprb) * lcoef(k) / 2._jprb
    lcoef_tl(k) = (2._jprb * k + 1._jprb) * lcoef_tl(k) / 2._jprb

    nmom = k

    IF (PRESENT(fixednmom)) THEN
      IF (k == MAX(minnmom, fixednmom)) EXIT
    ELSE
      ! Stop if the coefficients have grown sufficiently small
      IF (lcoef(k) < 1.E-6_jprb .AND. k >= minnmom) EXIT
    ENDIF
  ENDDO

END SUBROUTINE calc_legendre_coef_gauss_tl

!> Calculate Legendre coefficients for given phase function - AD
!! This call does NOT carry out the forward model calculations (cf TL).
!! @param[in]     u           Gaussian quadrature ordinates
!! @param[in]     w           Gaussian quadrature weights
!! @param[in]     f           phase function on quadrature ordinates
!! @param[in,out] f_ad        phase function increments on quadrature ordinates
!! @param[in]     nmom        number of coefficients (0...nmom)
!! @param[in,out] lcoef_ad    calculated Legendre coefficient increments
SUBROUTINE calc_legendre_coef_gauss_ad(u, w, f, f_ad, nmom, lcoef_ad)

  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  REAL(jprb),    INTENT(IN)              :: u(:)
  REAL(jprb),    INTENT(IN)              :: w(SIZE(u))
  REAL(jprb),    INTENT(IN)              :: f(SIZE(u))
  REAL(jprb),    INTENT(INOUT)           :: f_ad(SIZE(u))
  INTEGER(jpim), INTENT(IN)              :: nmom
  REAL(jprb),    INTENT(INOUT)           :: lcoef_ad(0:)

  INTEGER(jpim) :: ngauss, j, k
!   REAL(jprb)    :: y1(SIZE(u))
  REAL(jprb)    :: y1_ad(SIZE(u))
  REAL(jprb)    :: p00(SIZE(u),-1:nmom+1)
!   REAL(jprb)    :: p00_ad(SIZE(u),-1:nmom+1)

  ! direct calculations

  ngauss = SIZE(u)
  p00(:,-1) = 0._jprb
  p00(:,0) = 1._jprb
!   y1 = f * w

  DO k = 0, nmom
    DO j = 1, ngauss
      p00(j,k+1) = ((2._jprb * k + 1._jprb) * u(j) * p00(j,k) - k * p00(j,k-1)) / (k + 1._jprb)
    ENDDO
  ENDDO

  ! adjoint calculation

  f_ad = 0._jprb
  y1_ad = 0._jprb
!   p00_ad(:,:) = 0._jprb

  DO k = nmom, 0, -1
    lcoef_ad(k) = (2._jprb * k + 1._jprb) * lcoef_ad(k) / 2._jprb
    DO j = ngauss, 1, -1
      y1_ad(j) = y1_ad(j) + lcoef_ad(k) * p00(j,k)
!       p00_ad(j,k) = p00_ad(j,k) + lcoef_ad(k) * y1(j)
!       p00_ad(j,k) = p00_ad(j,k) + (2._jprb * k + 1._jprb) * u(j) * p00_ad(j,k+1) / (k + 1._jprb)
!       p00_ad(j,k-1) = p00_ad(j,k-1) - k * p00_ad(j,k+1) / (k + 1._jprb)
    ENDDO
  ENDDO

  f_ad = f_ad + y1_ad * w

END SUBROUTINE calc_legendre_coef_gauss_ad

!> Normalise a phase function which has been interpolated onto a Gaussian
!! quadrature
!! @param[in]     n     size of quadrature
!! @param[in]     w     quadrature weights
!! @param[in,out] f     phase function, normalised on output
SUBROUTINE normalise(n, w, f)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: w(n)
  REAL(jprb),    INTENT(INOUT) :: f(n)

  f = 2._jprb * f / SUM(f * w)

END SUBROUTINE normalise

!> Normalise a phase function interpolated onto a Gaussian quadrature - TL
!! This call carries out both the forward and TL model calculations.
!! @param[in]     n     size of quadrature
!! @param[in]     w     quadrature weights
!! @param[in,out] f     UN-NORMALISED phase function, normalised on output
!! @param[in,out] f_tl  phase function perturbations
SUBROUTINE normalise_tl(n, w, f, f_tl)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: w(n)
  REAL(jprb),    INTENT(INOUT) :: f(n), f_tl(n)

  REAL(jprb) :: tmp, tmp_tl

  tmp = SUM(f * w)
  tmp_tl = SUM(f_tl * w)

  ! NB f and f_tl overwritten here (cf direct): the order of the lines is important
  f_tl = 2._jprb * f_tl / tmp - 2._jprb * f * tmp_tl / tmp**2
  f = 2._jprb * f / tmp

END SUBROUTINE normalise_tl

!> Calculate Legendre coefficients for given phase function - AD
!! This call does NOT carry out the forward model calculations (cf TL).
!! @param[in]     n     size of quadrature
!! @param[in]     w     quadrature weights
!! @param[in]     f     UN-NORMALISED phase function (i.e. before call to normalise)
!! @param[in,out] f_ad  phase function increments
SUBROUTINE normalise_ad(n, w, f, f_ad)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: w(n), f(n)
  REAL(jprb),    INTENT(INOUT) :: f_ad(n)

  REAL(jprb) :: tmp, tmp_ad

  tmp = SUM(f * w)
  ! NB the order of the following lines is important as f_ad is overwritten (cf the TL above)
  tmp_ad = - 2._jprb * SUM(f * f_ad) / tmp**2
  f_ad = 2._jprb * f_ad / tmp + tmp_ad * w

END SUBROUTINE normalise_ad


SUBROUTINE spline_interp(n, x, y, nn, xn, yn)
! ---------spline fit to derive the yn value at point xn
!   Inputs:
!      n:    the length of x and y
!      x(n): the x values which x(1) < x(2) ... < x(n)
!      y(n): the y value which correspondent to x(n)
!      nn:  the length of vector xn and yn
!      xn:  the x value at which y value is wanted
!
!   Outputs:
!      yn: the wanted y value from the fitting
!
!   Internal variables:
!      yp1: the derivative of y over x at x(1), for natural bc, yp1=1.e31
!      ypn: the derivative of y over x at x(n), for natural bc, ypn=1.e31
!      y2(n): the second derivatives
!
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n, nn
  REAL(jprb),    INTENT(IN)    :: x(n), y(n), xn(nn)
  REAL(jprb),    INTENT(INOUT) :: yn(nn)

  INTEGER(jpim) :: i
  REAL(jprb)    :: y2(n), xx, yy, yp1, ypn

  yp1 = 1.e+31_jprb
  ypn = 1.e+31_jprb
  CALL spline(x, y, n, yp1, ypn, y2)

  DO i = 1, nn
    xx = xn(i)
    CALL splint(x, y, y2, n, xx, yy)
    yn(i) = yy
  ENDDO

END SUBROUTINE spline_interp

SUBROUTINE spline(x, y, n, yp1, ypn, y2)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: yp1, ypn, x(n), y(n)
  REAL(jprb),    INTENT(INOUT) :: y2(n)

  INTEGER(jpim) :: i, k
  REAL(jprb)    :: p, qn, sig, un, u(n)

  IF (yp1 > .99e+30_jprb) THEN
    y2(1) = 0._jprb
    u(1) = 0._jprb
  ELSE
    y2(1) = -0.5_jprb
    u(1) = (3._jprb / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
  ENDIF

  DO i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig * y2(i-1) + 2._jprb
    y2(i) = (sig - 1._jprb) / p
    u(i) = (6._jprb * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                       (y(i) - y(i-1)) / (x(i) - x(i-1))) / &
                       (x(i+1) - x(i-1)) - sig * u(i-1)) / p
  ENDDO

  IF (ypn > .99e+30_jprb) THEN
    qn = 0._jprb
    un = 0._jprb
  ELSE
    qn = 0.5_jprb
    un = (3._jprb / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
  ENDIF

  y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1._jprb)
  DO k = n-1, 1, -1
    y2(k) = y2(k) * y2(k+1) + u(k)
  ENDDO

END SUBROUTINE spline

SUBROUTINE splint(xa, ya, y2a, n, x, y)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: x, xa(n), y2a(n), ya(n)
  REAL(jprb),    INTENT(INOUT) :: y

  INTEGER(jpim) :: k, khi, klo
  REAL(jprb)    :: a, b, h

  klo = 1
  khi = n
  DO
    IF (khi-klo <= 1) EXIT
    k = (khi + klo) / 2
    IF (xa(k) > x) THEN
      khi = k
    ELSE
      klo = k
    ENDIF
  ENDDO

  h = xa(khi) - xa(klo)

  IF (h == 0._jprb) PRINT *, 'bad xa input in splint'

  a = (xa(khi) - x) / h
  b = (x - xa(klo)) / h
  y = a * ya(klo) + b * ya(khi) + &
      ((a ** 3 - a) * y2a(klo) + (b ** 3 - b) * y2a(khi)) * (h * h) / 6._jprb

END SUBROUTINE splint


SUBROUTINE spline_interp_tl(n, x, y, y_tl, nn, xn, yn, yn_tl)
  ! Direct and TL of spline_interp; x and xn are fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n, nn
  REAL(jprb),    INTENT(IN)    :: x(n), y(n), y_tl(n), xn(nn)
  REAL(jprb),    INTENT(INOUT) :: yn(nn), yn_tl(nn)

  INTEGER(jpim) :: i
  REAL(jprb)    :: y2(n), y2_tl(n), xx, yy, yy_tl, yp1, ypn

  yp1 = 1.e+31_jprb
  ypn = 1.e+31_jprb
  CALL spline_tl(x, y, y_tl, n, yp1, ypn, y2, y2_tl)

  DO i = 1, nn
    xx = xn(i)
    CALL splint_tl(x, y, y_tl, y2, y2_tl, n, xx, yy, yy_tl)
    yn(i) = yy
    yn_tl(i) = yy_tl
  ENDDO

END SUBROUTINE spline_interp_tl

SUBROUTINE spline_tl(x, y, y_tl, n, yp1, ypn, y2, y2_tl)
  ! Direct and TL of spline; x is fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: yp1, ypn, x(n), y(n), y_tl(n)
  REAL(jprb),    INTENT(INOUT) :: y2(n), y2_tl(n)

  INTEGER(jpim) :: i, k
  REAL(jprb)    :: p, p_tl, qn, sig, un, un_tl, u(n), u_tl(n), y3(n), y3_tl(n)

  IF (yp1 > .99e+30_jprb) THEN
    y2(1) = 0._jprb
    y2_tl(1) = 0._jprb
    u(1) = 0._jprb
    u_tl(1) = 0._jprb
  ELSE
    y2(1) = -0.5_jprb
    y2_tl(1) = 0._jprb
    u(1) = (3._jprb / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
    u_tl(1) = 3._jprb * (y_tl(2) - y_tl(1)) / (x(2) - x(1))**2
  ENDIF

  DO i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig * y2(i-1) + 2._jprb
    p_tl = sig * y2_tl(i-1)
    y2(i) = (sig - 1._jprb) / p
    y2_tl(i) = -p_tl * y2(i) / p
    u(i) = (6._jprb * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                       (y(i) - y(i-1)) / (x(i) - x(i-1))) / &
                       (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    u_tl(i) = (6._jprb * ((y_tl(i+1) - y_tl(i)) / (x(i+1) - x(i)) - &
                          (y_tl(i) - y_tl(i-1)) / (x(i) - x(i-1))) / &
                          (x(i+1) - x(i-1)) - sig * u_tl(i-1)) / p - p_tl * u(i) / p
  ENDDO

  IF (ypn > .99e+30_jprb) THEN
    qn = 0._jprb
    un = 0._jprb
    un_tl = 0._jprb
  ELSE
    qn = 0.5_jprb
    un = (3._jprb / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
    un_tl = 3._jprb * ( -(y_tl(n) - y_tl(n-1)) / (x(n) - x(n-1))**2)
  ENDIF

  y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1._jprb)
  y2_tl(n) = (un_tl - qn * u_tl(n-1)) / (qn * y2(n-1) + 1._jprb) - qn * y2_tl(n-1) * y2(n) / (qn * y2(n-1) + 1._jprb)
  y3 = y2
  y3_tl = y2_tl
  DO k = n-1, 1, -1
    y2(k) = y3(k) * y2(k+1) + u(k)
    y2_tl(k) = y3_tl(k) * y2(k+1) + y3(k) * y2_tl(k+1) + u_tl(k)
  ENDDO

END SUBROUTINE spline_tl

SUBROUTINE splint_tl(xa, ya, ya_tl, y2a, y2a_tl, n, x, y, y_tl)
  ! Direct and TL of splint; x and xa are fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: x, xa(n), y2a(n), ya(n), y2a_tl(n), ya_tl(n)
  REAL(jprb),    INTENT(INOUT) :: y, y_tl

  INTEGER(jpim) :: k, khi, klo
  REAL(jprb)    :: a, b, h

  klo = 1
  khi = n
  DO
    IF (khi-klo <= 1) EXIT
    k = (khi + klo) / 2
    IF (xa(k) > x) THEN
      khi = k
    ELSE
      klo = k
    ENDIF
  ENDDO

  h = xa(khi) - xa(klo)

  IF (h == 0._jprb) PRINT *, 'bad xa input in splint_tl'

  a = (xa(khi) - x) / h
  b = (x - xa(klo)) / h
  y = a * ya(klo) + b * ya(khi) + &
      ((a ** 3 - a) * y2a(klo) + (b ** 3 - b) * y2a(khi)) * (h * h) / 6._jprb
  y_tl = a * ya_tl(klo) + b * ya_tl(khi) + &
         ((a ** 3 - a) * y2a_tl(klo) + (b ** 3 - b) * y2a_tl(khi)) * (h * h) / 6._jprb

END SUBROUTINE splint_tl


SUBROUTINE spline_interp_ad(n, x, y, y_ad, nn, xn, yn_ad)
  ! AD of spline_interp; x and xn are fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n, nn
  REAL(jprb),    INTENT(IN)    :: x(n), y(n), xn(nn), yn_ad(nn)
  REAL(jprb),    INTENT(INOUT) :: y_ad(n)

  INTEGER(jpim) :: i
  REAL(jprb)    :: y2_ad(n), xx, yy_ad, yp1, ypn

  yp1 = 1.e+31_jprb
  ypn = 1.e+31_jprb

  y2_ad = 0._jprb
  DO i = nn, 1, -1
    xx = xn(i)
    yy_ad = yn_ad(i)
    CALL splint_ad(x, y_ad, y2_ad, n, xx, yy_ad)
  ENDDO

  CALL spline_ad(x, y, y_ad, n, yp1, ypn, y2_ad)

END SUBROUTINE spline_interp_ad

SUBROUTINE spline_ad(x, y, y_ad, n, yp1, ypn, y2_ad)
  ! AD of spline; x is fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: yp1, ypn, x(n), y(n)
  REAL(jprb),    INTENT(INOUT) :: y_ad(n), y2_ad(n)

  INTEGER(jpim) :: i, k
  REAL(jprb)    :: p, p_ad, qn, sig, un, un_ad, u(n), u_ad(n), y2(n), y3(n), y3_ad(n)

  ! direct calculations

  IF (yp1 > .99e+30_jprb) THEN
    y2(1) = 0._jprb
    u(1) = 0._jprb
  ELSE
    y2(1) = -0.5_jprb
    u(1) = (3._jprb / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
  ENDIF

  DO i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig * y2(i-1) + 2._jprb
    y2(i) = (sig - 1._jprb) / p
    u(i) = (6._jprb * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                       (y(i) - y(i-1)) / (x(i) - x(i-1))) / &
                       (x(i+1) - x(i-1)) - sig * u(i-1)) / p
  ENDDO

  IF (ypn > .99e+30_jprb) THEN
    qn = 0._jprb
    un = 0._jprb
  ELSE
    qn = 0.5_jprb
    un = (3._jprb / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
  ENDIF

  y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1._jprb)
  y3 = y2
  DO k = n-1, 1, -1
    y2(k) = y3(k) * y2(k+1) + u(k)
  ENDDO

  ! adjoint

  y3_ad = 0._jprb
  u_ad = 0._jprb
  un_ad = 0._jprb
  DO k = 1, n-1
    u_ad(k) = u_ad(k) + y2_ad(k)
    y2_ad(k+1) = y2_ad(k+1) + y3(k) * y2_ad(k)
    y3_ad(k) = y3_ad(k) + y2(k+1) * y2_ad(k)
  ENDDO

  y2_ad = y2_ad + y3_ad

  un_ad = un_ad + y2_ad(n) / (qn * y3(n-1) + 1._jprb)
  u_ad(n-1) = u_ad(n-1) - y2_ad(n) * qn / (qn * y3(n-1) + 1._jprb)
  y2_ad(n-1) = y2_ad(n-1) - y2_ad(n) * qn * y3(n) / (qn * y3(n-1) + 1._jprb)

  IF (ypn > .99e+30_jprb) THEN
    un_ad = 0._jprb
  ELSE
    y_ad(n) = y_ad(n) - un_ad * 3._jprb / (x(n) - x(n-1))**2
    y_ad(n-1) = y_ad(n-1) + un_ad * 3._jprb / (x(n) - x(n-1))**2
  ENDIF

  DO i = n-1, 2, -1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig * y3(i-1) + 2._jprb

    y_ad(i+1) = y_ad(i+1) + 6._jprb * u_ad(i) / (p * (x(i+1) - x(i)) * (x(i+1) - x(i-1)))
    y_ad(i)   = y_ad(i)   - 6._jprb * u_ad(i) / (p * (x(i+1) - x(i)) * (x(i+1) - x(i-1)))
    y_ad(i)   = y_ad(i)   - 6._jprb * u_ad(i) / (p * (x(i) - x(i-1)) * (x(i+1) - x(i-1)))
    y_ad(i-1) = y_ad(i-1) + 6._jprb * u_ad(i) / (p * (x(i) - x(i-1)) * (x(i+1) - x(i-1)))
    u_ad(i-1) = u_ad(i-1) - u_ad(i) * sig / p
    p_ad = -u_ad(i) * u(i) / p

    p_ad = p_ad - y2_ad(i) * y3(i) / p

    y2_ad(i-1) = y2_ad(i-1) + sig * p_ad
  ENDDO

  IF (yp1 > .99e+30_jprb) THEN
    y2_ad(1) = 0._jprb
    u_ad(1) = 0._jprb
  ELSE
    y2_ad(1) = 0._jprb

    y_ad(2) = y_ad(2) + u_ad(1) * 3._jprb / (x(2) - x(1))**2
    y_ad(1) = y_ad(1) - u_ad(1) * 3._jprb / (x(2) - x(1))**2
  ENDIF
END SUBROUTINE spline_ad

SUBROUTINE splint_ad(xa, ya_ad, y2a_ad, n, x, y_ad)
  ! AD of splint; x and xa are fixed

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: x, xa(n), y_ad
  REAL(jprb),    INTENT(INOUT) :: ya_ad(n), y2a_ad(n)

  INTEGER(jpim) :: k, khi, klo
  REAL(jprb)    :: a, b, h

  klo = 1
  khi = n
  DO
    IF (khi-klo <= 1) EXIT
    k = (khi + klo) / 2
    IF (xa(k) > x) THEN
      khi = k
    ELSE
      klo = k
    ENDIF
  ENDDO

  h = xa(khi) - xa(klo)

  IF (h == 0._jprb) PRINT *, 'bad xa input in splint_ad'

  a = (xa(khi) - x) / h
  b = (x - xa(klo)) / h

  ya_ad(klo) = ya_ad(klo) + a * y_ad
  ya_ad(khi) = ya_ad(khi) + b * y_ad
  y2a_ad(klo) = y2a_ad(klo) + (a ** 3 - a) * y_ad * (h * h) / 6._jprb
  y2a_ad(khi) = y2a_ad(khi) + (b ** 3 - b) * y_ad * (h * h) / 6._jprb
END SUBROUTINE splint_ad


SUBROUTINE integrate(n, x, f, intg, error, fail)

  ! This replaces the NAG library d01GAF subroutine to avoid dependency on an external
  ! library. The differences in the calculated quantities between this and d01GAF are
  ! extremely small.

  ! Uses the method of Gill and Miller 1972 (I copied the code from the paper) and adds
  ! the estimated error to the output value for consistency with d01GAF (the NAG folks
  ! found this generally results in better agreement with the true value)

  USE parkind1, ONLY : jprb, jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)  :: n
  REAL(KIND=jprb),    INTENT(IN)  :: x(1:n)
  REAL(KIND=jprb),    INTENT(IN)  :: f(1:n)
  REAL(KIND=jprb),    INTENT(OUT) :: intg
  REAL(KIND=jprb),    INTENT(OUT) :: error
  INTEGER(KIND=jpim), INTENT(OUT) :: fail

  INTEGER(KIND=jpim) :: i, j, k
  REAL(KIND=jprb)    :: h1, h2, h3, h4, r1, r2, r3, r4, d1, d2, d3, c, s

  fail = 0

  IF (n < 4) THEN
    fail = 1
    intg = 0._jprb
    error = 0._jprb
    PRINT *,'Requires at least 4 data points'
    RETURN
  ENDIF

  intg  = 0._jprb
  error = 0._jprb
  s     = 0._jprb
  c     = 0._jprb
  r4    = 0._jprb

  j = 3
  k = n - 1
  DO i = j, k
    IF (i == j) THEN
      h2 = x(j-1) - x(j-2)
      d3 = (f(j-1) - f(j-2))/h2
      h3 = x(j) - x(j-1)
      d1 = (f(j) - f(j-1))/h3
      h1 = h2 + h3
      d2 = (d1 - d3)/h1
      h4 = x(j+1) - x(j)
      r1 = (f(j+1) - f(j))/h4
      r2 = (r1 - d1)/(h4 + h3)
      h1 = h1 + h4
      r3 = (r2 - d2)/h1
      intg = h2 * (f(1) + h2 * (d3/2._jprb - h2 * (d2/6._jprb - (h2 + 2._jprb * h3) * r3/12._jprb)))
      s    = -h2**3 * (h2 * (3._jprb * h2 + 5._jprb * h4) + 10._jprb * h3 * h1) / 60._jprb
    ELSE
      h4 = x(i+1) - x(i)
      r1 = (f(i+1) - f(i))/h4
      r4 = h4 + h3
      r2 = (r1 - d1)/r4
      r4 = r4 + h2
      r3 = (r2 - d2)/r4
      r4 = r4 + h1
      r4 = (r3 - d3)/r4
    ENDIF

    intg = intg + h3 * ((f(i) + f(i-1))/2._jprb - h3 * h3 * (d2 + r2 + (h2 - h4) * r3)/12._jprb)
    c    = h3**3 * (2._jprb * h3 * h3 +  5._jprb * (h3 * (h4 + h2) + 2._jprb * h4 * h2))/120._jprb
    error = error + (c + s) * r4
    IF (i == j) THEN
      s = 2._jprb*c + s
    ELSE
      s = c
    ENDIF

    IF (i == k) THEN
      intg = intg + h4 * (f(n) - h4 * (r1/2._jprb + h4 * (r2/6._jprb + (2._jprb * h3 + h4) * r3/12._jprb)))
      error = error - h4**3 * r4 * (h4 * (3._jprb * h4 + 5._jprb * h2) + 10._jprb * h3 * (h2 + h3 + h4))/60._jprb
      error = error + s * r4
    ELSE
      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    ENDIF
  ENDDO
  intg = intg + error  ! The error is added to the answer in the NAG lib function
END SUBROUTINE integrate


!***********************************************************************
!  PROGRAM        INTER   SUBROUTINE
!-----------------------------------------------------------------------
!  PURPOSE        TO INTERPOLATE AT THE POINT ARG BETWEEN THE ARRAY ARX
!                 AND THE CORRESPONDING FUNCTION VALUES ARY.
!                 EXPONENTIAL INTERPOLATION IS USED FOR PRESSURE AND
!                 NUMBER DENSITY, LINEAR INTERPOLATION FOR TEMPERATURE.
!                 IF THE MODE NO IS <0 NEGATIVE NUMBERS ARE SET = 0.
!-----------------------------------------------------------------------
!  VERSION        4.0   D.P. EDWARDS   15/08/94
!                 Last changed 96/11/05
!-----------------------------------------------------------------------
!  ARGUMENTS      MXDIM   I*4  I/P  ARRAY DIMENSION
!                 NREC    I*4  I/P  NO OF ELEMENTS IN ARRAYS ARX AND ARY
!                 MODE    I*4  I/P  INTERPOLATION MODE
!                 ARG     R*4  I/P  INTERPOLATION ARGUMENT
!                 ARX     R*4  I/P  X VALUE ARRAY
!                 ARY     R*4  I/P  Y VALUE FUNCTION ARRAY
!                 SS      R*4  O/P  INTERPOLATED FUNCTION VALUE AT ARG
!                 H       R*4  O/P  GRADIENT OR SCALE HEIGHT VALUE
!***********************************************************************
       SUBROUTINE INTER(MXDIM,NREC,MODE,ARG,ARX,ARY,SS,H)
!-----------------------------------------------------------------------
       USE parkind1, ONLY : jprb, jpim

       INTEGER(KIND=JPIM)   MXDIM, NREC, MODE
       REAL(KIND=JPRB)      ARG
       REAL(KIND=JPRB)      ARX(NREC),ARY(NREC)
       REAL(KIND=JPRB)      SS, H

       INTEGER(KIND=JPIM)   IR, KL, NO, NOO, KS, KF, K
       REAL(KIND=JPRB)      AA, BB
!***********************************************************************
!
       IF (ARG .GE. ARX(1) .AND. ARG .LE. ARX(NREC)) THEN
         DO 10 IR=1,NREC-1
           IF (ARG .GE. ARX(IR) .AND. ARG .LT. ARX(IR+1)) KL = IR
   10    CONTINUE
         IF (ARG .EQ. ARX(NREC)) KL = NREC - 1
       ELSEIF (ARG .LT. ARX(1)) THEN
         KL = 1
       ELSEIF (ARG .GT. ARX(NREC)) THEN
         KL = NREC - 1
       ENDIF
!
!  INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KL TO KL+1
!
!  EXPONENTIAL INTERPOLATION
!
       IF (ABS(MODE) .EQ. 1 .OR. ABS(MODE) .EQ. 4) THEN
         H = -(ARX(KL+1) - ARX(KL))/LOG(ARY(KL+1)/ARY(KL))
         SS = ARY(KL)*EXP(-(ARG-ARX(KL))/H)
         IF (ABS(MODE) .EQ. 4)  SS = 1.0/SS
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LINEAR INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 2 .OR. ABS(MODE) .EQ. 5) THEN
         H = (ARY(KL+1) - ARY(KL))/(ARX(KL+1) - ARX(KL))
         SS = ARY(KL) + H*(ARG - ARX(KL))
         IF (ABS(MODE) .EQ. 5) SS = 1.0/SS
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LINEAR-LOG INTERPOLATION
       ELSEIF (ABS(MODE) .EQ. 7) THEN
         H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(ARX(KL+1)-ARX(KL))
         SS = EXP(LOG(ARY(KL)) + H*(ARG - ARX(KL)))
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LOG-LOG INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 8) THEN
         H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(LOG(ARX(KL+1))-LOG(ARX(KL)))
         SS = EXP(LOG(ARY(KL)) + H*(LOG(ARG) - LOG(ARX(KL))))
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LOGARITHMIC INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 3) THEN
         AA = ARX(KL+1)/ARX(KL)
         BB = ARG/ARX(KL)
         IF (AA .EQ. BB) THEN
           SS = ARY(KL+1)
         ELSE
           H = (ARY(KL+1) - ARY(KL))/LOG(AA)
           SS = ARY(KL) + H*LOG(BB)
         ENDIF
         IF (MODE .LT. 0 .AND. SS .LT. 0.0) SS = 0.0
!
!  LAGRANGIAN INTERPOLATION
!
       ELSEIF (ABS(MODE) .EQ. 6) THEN
!
!  NUMBER OF DATA POINT POINTS TO INTERPOLATE BETWEEN
!
         NO = 4
!
!  FIND DATA POINTS BETWEEN WHICH TO INTERPOLATE
!
         NOO = NO
   20    IF (ARG .LT. ARX(1)) THEN
           NOO = 2
           KS = 1
           KF = 1 + NOO - 1
         ELSEIF (ARG .GT. ARX(NREC)) THEN
           NOO = 2
           KS = NREC - NOO + 1
           KF = NREC
         ELSE
           IF (MOD(NOO,2) .EQ. 0) THEN
             KS = KL - 0.5*NOO + 1
             KF = KL + 0.5*NOO
           ELSE
             KS = KL - 0.5*(NOO - 1) + 1
             KF = KL + 0.5*(NOO + 1)
           ENDIF
           IF (KS .LT. 1) KS = 1
           IF (KF .GT. NREC) KF = NREC
         ENDIF
!
!  INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KS TO KF
!
         SS = 0.0
         DO 30 K=KS,KF
           SS = SS + XL(MXDIM,KS,KF,K,ARG,ARX)*ARY(K)     &
          /XL(MXDIM,KS,KF,K,ARX(K),ARX)
   30    CONTINUE
         H = NOO
!
!  IF INTERPOLATION HAS OVERSHOT, REDUCE ORDER AND RETRY
!
         IF (ARG .LT. ARX(1) .OR. ARG .GT. ARX(NREC)) THEN
           IF (SS .LT. 0.0) THEN
             IF (NOO .EQ. 2) THEN
               SS = 0.0
             ELSE
               NOO = NOO - 1
               GOTO 20
             ENDIF
           ENDIF
         ELSE
           IF (((ARY(KL) .LE. ARY(KL+1)) .AND.                       &
          (SS .LT. ARY(KL) .OR. SS .GT. ARY(KL+1))) .OR.            &
          ((ARY(KL) .GT. ARY(KL+1)) .AND.                           &
          (SS .GT. ARY(KL) .OR. SS .LT. ARY(KL+1)))) THEN
             NOO = NOO-1
             GOTO 20
           ENDIF
         ENDIF
!
       ENDIF
!
!-----------------------------------------------------------------------
       RETURN
       END SUBROUTINE
!***********************************************************************
!
!  PROGRAM        XL  FUNCTION
!
!  PURPOSE        TO COMPUTE LAGRANGE INTERPOLATION COEFFICIENTS
!
!  VERSION        3.0   D.P. EDWARDS   01/01/89
!
!  ARGUMENTS      MXDIM   I*4  I/P  ARRAY DIMENSION
!                 KS      I*4  I/P  LOWER LIMIT OF LAGRANGE SUM
!                 KF      I*4  I/P  UPPER LIMIT OF LAGRANGE SUM
!                 K       I*4  I/P  CURRENT INDEX OF LAGRANGE SUM
!                 ARG     R*4  I/P  INTERPOLATION ARGUMENT
!                 ARR     R*4  I/P  ARRAY TO INTERPOLATE BETWEEN
!
!***********************************************************************
!
       FUNCTION XL(MXDIM,KS,KF,K,ARG,ARR)
!-----------------------------------------------------------------------
       USE parkind1, ONLY : jprb, jpim

       INTEGER(KIND=JPIM) MXDIM
       INTEGER(KIND=JPIM) KS, KF, K
       REAL(KIND=JPRB)    ARG
       REAL(KIND=JPRB)    ARR(MXDIM)
       REAL(KIND=JPRB)    XL

       INTEGER(KIND=JPIM) J
       REAL(KIND=JPRB)    PROD
!-----------------------------------------------------------------------
!
       PROD = 1.0
       DO 10 J=KS,KF
         IF (J .NE. K) THEN
           PROD = PROD*(ARG - ARR(J))
         ENDIF
   10  CONTINUE
!
       XL = PROD
!
!-----------------------------------------------------------------------
       RETURN
       END FUNCTION


! --------------------------------------------------------------------------
! The following subroutines have been taken from DISORT and CRTM
! --------------------------------------------------------------------------

! ASYMTX: from DISORT, style updated, the DOUBLE PREC arguments to have been
!         removed so that all real variables use JPRB.
! ASYMTX_TL/AD and matinv: from CRTM, only superficial changes to match RTTOV
!         coding style and use of RTTOV KINDs.

! ASYMTX called DISORT's D1MACH subroutine: this has been replaced by
! the d1mach4 module variable above.

!> Calculate eigenvalues and eigenvectors for a real asymmetric matrix for
!! which it is known a priori that the eigenvalues are real. This is taken
!! from DISORT with minor updates, mostly for RTTOV coding style, to use
!! RTTOV KINDs, and to normalise the eigenvectors.
!! @param[out]      err     status on exit (non-zero implies failure to converge)
!! @param[in,out]   aa      input asymmetric matrix, overwritten on exit
!! @param[out]      evec    normalised eigenvectors (column j corresponds to eval(j))
!! @param[out]      eval    unordered eigenvalues
!! @param[in]       m       order of aa
!! @param[in]       ia      first dimension of aa
!! @param[in]       ievec   first dimension of evec
SUBROUTINE ASYMTX(ERR, AA, EVEC, EVAL, M, IA, IEVEC)

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======
!
!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.
!
!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.
!
!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the Parlett-Reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.
!
!       There is a statement 'J  = WKD( I )' that converts a double
!       precision variable to an integer variable, that seems dangerous
!       to us in principle, but seems to work fine in practice.
!
!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford
!
!
!   I N P U T    V A R I A B L E S:
!
!       AA    :  input asymmetric matrix, destroyed after solved
!
!        M    :  order of  AA
!
!       IA    :  first dimension of  AA
!
!    IEVEC    :  first dimension of  EVEC
!
!
!   O U T P U T    V A R I A B L E S:
!
!       EVEC  :  (unnormalized) eigenvectors of  AA
!                   ( column J corresponds to EVAL(J) )
!
!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )
!
!       ERR   :  if  /=  0, signals that EVAL(ERR) failed to converge;
!                   in that case eigenvalues ERR+1,ERR+2,...,M  are
!                   correct but eigenvalues 1,...,ERR are set to zero.
!
! +-------------------------------------------------------------------+

  USE parkind1, ONLY : jpim, jprb, jplm

  IMPLICIT NONE

  ! Scalar Arguments

  INTEGER(jpim), INTENT(IN)  :: IA, IEVEC, M
  INTEGER(jpim), INTENT(OUT) :: ERR

  ! Array Arguments

  REAL(jprb), INTENT(INOUT)  :: AA(IA,M)
  REAL(jprb), INTENT(OUT)    :: EVAL(M), EVEC(IEVEC,M)

  ! Local Scalars

  LOGICAL(jplm) :: NOCONV, NOTLAS
  INTEGER(jpim) :: I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
  REAL(jprb)    :: C1, C2, C3, C4, C5, C6, COL, DISCRI, F, G, H, &
                   ONE, P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T, &
                   TOL, UU, VV, W, X, Y, Z, ZERO
  REAL(jprb)    :: C5_R, C6_R, F_R, H_R
  REAL(jprb)    :: WKD(2*M)

  DATA      C1 / 0.4375_jprb / , C2 / 0.5_jprb / , C3 / 0.75_jprb / , &
            C4 / 0.95_jprb / , C5 / 16._jprb / , C6 / 256._jprb / , &
            ZERO / 0._jprb / , ONE / 1._jprb /
  TRY

  ERR  = 0
  TOL  = D1MACH4

  C5_R = 1._jprb / C5
  C6_R = 1._jprb / C6

  IF( M < 1 .OR. IA < M .OR. IEVEC < M ) THEN
    ERR = ERRORSTATUS_FATAL
    THROWM(ERR.NE.0,"ASYMTX--bad input variable(s)")
  ENDIF

  ! Handle 1x1 and 2x2 special cases
  IF( M == 1 ) THEN
     EVAL( 1 )   = AA( 1,1 )
     EVEC( 1,1 ) = 1.0
     RETURN
  ELSE IF( M == 2 ) THEN
     DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 + 4.*AA( 1,2 )*AA( 2,1 )

     IF( DISCRI  <  0.0 ) THEN
       ERR = ERRORSTATUS_FATAL
       THROWM(ERR.NE.0,"ASYMTX--complex evals in 2x2 case")
     ENDIF

     SGN  = ONE
     IF( AA( 1,1 ) < AA( 2,2 ) ) SGN  = - ONE

     EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
     EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
     EVEC( 1,1 ) = 1.0
     EVEC( 2,2 ) = 1.0

     IF( AA( 1,1 )  ==  AA( 2,2 ) .AND. &
         ( AA( 2,1 ) == 0.0 .OR. AA( 1,2 ) == 0.0 ) ) THEN

        RNORM = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) + &
                ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
        W     = TOL * RNORM
        EVEC( 2,1 ) =   AA( 2,1 ) / W
        EVEC( 1,2 ) = - AA( 1,2 ) / W
     ELSE
        EVEC( 2,1 ) = AA( 2,1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
        EVEC( 1,2 ) = AA( 1,2 ) / ( EVAL( 2 ) - AA( 1,1 ) )
     ENDIF

     DO j = 1, M
        EVEC(:,J) = EVEC(:,J) / SQRT(DOT_PRODUCT(EVEC(:,J), EVEC(:,J)))
     ENDDO

     RETURN
  ENDIF

  ! Initialize output variables
  ERR  = 0
  EVAL(:) = ZERO
  EVEC(:,:) = ZERO
  DO I = 1, M
     EVEC( I, I ) = ONE
  ENDDO

  ! Balance the input matrix and reduce its norm by
  ! diagonal similarity transformation stored in WK;
  ! then search for rows isolating an eigenvalue
  ! and push them down
  RNORM = ZERO
  L = 1
  K = M

 50 CONTINUE
  KKK = K

  DO J = KKK, 1, -1
     ROW = ZERO
     DO I = 1, K
        IF( I /= J ) ROW  = ROW + ABS( AA( J,I ) )
     ENDDO
     IF( ROW == ZERO ) THEN
        WKD( K ) = J
        IF( J /= K ) THEN
           DO I = 1, K
              REPL       = AA( I, J )
              AA( I, J ) = AA( I, K )
              AA( I, K ) = REPL
           ENDDO
           DO I = L, M
              REPL       = AA( J, I )
              AA( J, I ) = AA( K, I )
              AA( K, I ) = REPL
           ENDDO
        ENDIF
        K  = K - 1
        GO TO  50
     ENDIF
  ENDDO

  ! Search for columns isolating an eigenvalue and push them left
100 CONTINUE
  LLL = L

  DO J = LLL, K
     COL  = ZERO
     DO I = L, K
        IF( I /= J ) COL  = COL + ABS( AA( I,J ) )
     ENDDO
     IF( COL == ZERO ) THEN
        WKD( L ) = J
        IF( J /= L ) THEN
           DO I = 1, K
              REPL       = AA( I, J )
              AA( I, J ) = AA( I, L )
              AA( I, L ) = REPL
           ENDDO
           DO I = L, M
              REPL       = AA( J, I )
              AA( J, I ) = AA( L, I )
              AA( L, I ) = REPL
           ENDDO
        ENDIF
        L  = L + 1
        GO TO  100
     ENDIF
  ENDDO

  ! Balance the submatrix in rows L through K
  WKD(L:K) = ONE

  NOCONV = .TRUE.
  DO WHILE (NOCONV)
    NOCONV = .FALSE.

    DO I = L, K
       COL  = ZERO
       ROW  = ZERO
       DO J = L, K
          IF( J /= I ) THEN
             COL  = COL + ABS( AA( J,I ) )
             ROW  = ROW + ABS( AA( I,J ) )
          ENDIF
       ENDDO

       F  = ONE
       G  = ROW * C5_R
       H  = COL + ROW
       DO WHILE (COL < G )
          F    = F * C5
          COL  = COL * C6
       ENDDO

       G  = ROW * C5
       DO WHILE (COL >= G)
          F    = F * C5_R
          COL  = COL * C6_R
       ENDDO

       ! Now balance
       F_R = 1._jprb / F
!        IF( ( COL + ROW ) * F_R < C4*H ) THEN
       IF( F_R < C4 ) THEN  ! H == COL + ROW (see above)
          WKD( I ) = WKD( I )*F
          NOCONV = .TRUE.
          AA(I,L:M) = AA(I,L:M) * F_R
          AA(1:K,I) = AA(1:K,I) * F
       ENDIF
    ENDDO
  ENDDO

  ! Is A already in Hessenberg form?
  IF( K-1 >= L+1 ) THEN

    ! Transfer A to a Hessenberg form
    DO N = L + 1, K - 1
       H  = ZERO
       WKD( N + M ) = ZERO

       ! Scale column
       SCALE = SUM(ABS(AA(N:K,N-1)))

       IF( SCALE /= ZERO ) THEN
          DO I = K, N, -1
             WKD( I + M ) = AA( I, N - 1 ) / SCALE
             H  = H + WKD( I + M )**2
          ENDDO

          G    = - SIGN( SQRT( H ), WKD( N + M ) )
          H    = H - WKD( N + M )*G
          H_R  = 1._jprb / H
          WKD( N + M ) = WKD( N + M ) - G

          ! Form (I-(U*UT)/H)*A
          DO J = N, M
             F  = ZERO
!              DO I = K, N, -1
!                 F  = F + WKD( I + M )*AA( I, J )
!              ENDDO
! No need for reversed loop:
             DO I = N, K
                F  = F + WKD( I + M )*AA( I, J )
             ENDDO

             DO I = N, K
                AA( I, J ) = AA( I, J ) - WKD( I + M )*F * H_R
             ENDDO
          ENDDO

          ! Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
          DO I = 1, K
             F  = ZERO
!              DO J = K, N, -1
!                 F  = F + WKD( J + M )*AA( I, J )
!              ENDDO
! No need for reversed loop:
             DO J = N, K
                F  = F + WKD( J + M )*AA( I, J )
             ENDDO

             DO J = N, K
                AA( I, J ) = AA( I, J ) - WKD( J + M )*F * H_R
             ENDDO
          ENDDO
          WKD( N + M ) = SCALE*WKD( N + M )
          AA( N, N - 1 ) = SCALE*G
       ENDIF
    ENDDO

    DO N = K - 2, L, -1
       N1   = N + 1
       N2   = N + 2
       F  = AA( N + 1, N )
       IF( F /= ZERO ) THEN
          F  = F*WKD( N + 1 + M )
          DO I = N + 2, K
             WKD( I + M ) = AA( I, N )
          ENDDO
          IF( N + 1 <= K ) THEN
             F_R = 1._jprb / F
             DO J = 1, M
!                 G  = ZERO
!                 DO I = N + 1, K
!                    G  = G + WKD( I + M )*EVEC( I, J )
!                 ENDDO
!                 G  = G * F_R
!                 DO I = N + 1, K
!                    EVEC( I, J ) = EVEC( I, J ) + G*WKD( I + M )
!                 ENDDO
! Equivalent:
                G = F_R * SUM(WKD(N+1+M:K+M) * EVEC(N+1:K,J))
                EVEC(N+1:K,J) = EVEC(N+1:K,J) + G * WKD(N+1+M:K+M)
             ENDDO
          ENDIF
       ENDIF
    ENDDO

  ENDIF

  N  = 1
  DO I = 1, M
     RNORM  = RNORM + SUM(ABS(AA(I,N:M)))
     N  = I
     IF( I < L .OR. I > K ) EVAL( I ) = AA( I, I )
  ENDDO

  N  = K
  T  = ZERO

  ! Search for next eigenvalues
400 CONTINUE
  IF( N < L ) GO TO  550

  IN  = 0
  N1  = N - 1
  N2  = N - 2

  ! Look for single small sub-diagonal element
410 CONTINUE

  DO I = L, N
     LB  = N + L - I
     IF( LB == L ) EXIT
     S  = ABS( AA( LB - 1,LB - 1 ) ) + ABS( AA( LB,LB ) )
     IF( S == ZERO ) S  = RNORM
     IF( ABS( AA( LB, LB-1 ) ) <= TOL*S ) EXIT
  ENDDO

  X  = AA( N, N )

  IF( LB == N ) THEN
     ! One eigenvalue found
     AA( N, N ) = X + T
     EVAL( N ) = AA( N, N )
     N  = N1
     GO TO  400
  ENDIF

  Y  = AA( N1, N1 )
  W  = AA( N, N1 )*AA( N1, N )

  IF( LB == N1 ) THEN
     ! Two eigenvalues found
     P  = ( Y - X )*C2
     Q  = P**2 + W
     Z  = SQRT( ABS( Q ) )
     AA( N, N ) = X + T
     X  = AA( N, N )
     AA( N1, N1 ) = Y + T

     ! Real pair
     Z  = P + SIGN( Z, P )
     EVAL( N1 ) = X + Z
     EVAL( N ) = EVAL( N1 )

     IF( Z /= ZERO ) EVAL( N ) = X - W / Z

     X  = AA( N, N1 )

     ! Employ scale factor in case X and Z are very small
     R  = SQRT( X*X + Z*Z )
     P  = X / R
     Q  = Z / R

     ! Row modification
     DO J = N1, M
        Z  = AA( N1, J )
        AA( N1, J ) = Q*Z + P*AA( N, J )
        AA( N, J ) = Q*AA( N, J ) - P*Z
     ENDDO

     ! Column modification
     DO I = 1, N
        Z  = AA( I, N1 )
        AA( I, N1 ) = Q*Z + P*AA( I, N )
        AA( I, N ) = Q*AA( I, N ) - P*Z
     ENDDO

     ! Accumulate transformations
     DO I = L, K
        Z  = EVEC( I, N1 )
        EVEC( I, N1 ) = Q*Z + P*EVEC( I, N )
        EVEC( I, N ) = Q*EVEC( I, N ) - P*Z
     ENDDO

     N  = N2
     GO TO  400
  ENDIF

  IF( IN == 30 ) THEN
     ! No convergence after 30 iterations; set error
     ! indicator to the index of the current eigenvalue
     ERR  = N
     RETURN
  ENDIF

  ! Form shift
  IF( IN == 10 .OR. IN == 20 ) THEN
     T  = T + X
     DO I = L, N
        AA( I, I ) = AA( I, I ) - X
     ENDDO
     S  = ABS( AA( N,N1 ) ) + ABS( AA( N1,N2 ) )
     X  = C3*S
     Y  = X
     W  = -C1*S**2
  ENDIF

  IN  = IN + 1

  ! Look for two consecutive small sub-diagonal elements
  DO J = LB, N2
     I  = N2 + LB - J
     Z  = AA( I, I )
     R  = X - Z
     S  = Y - Z
     P  = ( R*S - W ) / AA( I + 1, I ) + AA( I, I + 1 )
     Q  = AA( I + 1, I + 1 ) - Z - R - S
     R  = AA( I + 2, I + 1 )
     S  = ABS( P ) + ABS( Q ) + ABS( R )
     P  = P / S
     Q  = Q / S
     R  = R / S

     IF( I == LB ) EXIT

     UU   = ABS( AA( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
     VV   = ABS( P ) * ( ABS( AA( I-1, I-1 ) ) + ABS( Z ) + &
                         ABS( AA( I+1, I+1 ) ) )
     IF( UU <= TOL*VV ) EXIT
  ENDDO

  AA( I+2, I ) = ZERO
  DO J = I + 3, N
     AA( J, J - 2 ) = ZERO
     AA( J, J - 3 ) = ZERO
  ENDDO

  ! Double QR step involving rows K to N and columns M to N
  DO KA = I, N1
     NOTLAS = KA /= N1
     IF( KA == I ) THEN
        S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
        IF( LB /= I ) AA( KA, KA - 1 ) = -AA( KA, KA - 1 )
     ELSE
        P  = AA( KA, KA - 1 )
        Q  = AA( KA + 1, KA - 1 )
        R  = ZERO
        IF( NOTLAS ) R  = AA( KA + 2, KA - 1 )

        X  = ABS( P ) + ABS( Q ) + ABS( R )
        IF( X == ZERO ) EXIT
        P  = P / X
        Q  = Q / X
        R  = R / X
        S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
        AA( KA, KA - 1 ) = -S*X
     ENDIF

     P  = P + S
     X  = P / S
     Y  = Q / S
     Z  = R / S
     Q  = Q / P
     R  = R / P

     ! Row modification
     DO J = KA, M
        P  = AA( KA, J ) + Q*AA( KA + 1, J )
        IF( NOTLAS ) THEN
           P  = P + R*AA( KA + 2, J )
           AA( KA + 2, J ) = AA( KA + 2, J ) - P*Z
        ENDIF
        AA( KA + 1, J ) = AA( KA + 1, J ) - P*Y
        AA( KA, J ) = AA( KA, J ) - P*X
     ENDDO

     ! Column modification
     DO II = 1, MIN( N, KA + 3 )
        P  = X*AA( II, KA ) + Y*AA( II, KA + 1 )
        IF( NOTLAS ) THEN
           P  = P + Z*AA( II, KA + 2 )
           AA( II, KA + 2 ) = AA( II, KA + 2 ) - P*R
        ENDIF
        AA( II, KA + 1 ) = AA( II, KA + 1 ) - P*Q
        AA( II, KA ) = AA( II, KA ) - P
     ENDDO

     ! Accumulate transformations
     DO II = L, K
        P  = X*EVEC( II, KA ) + Y*EVEC( II, KA + 1 )
        IF( NOTLAS ) THEN
           P  = P + Z*EVEC( II, KA + 2 )
           EVEC( II, KA + 2 ) = EVEC( II, KA + 2 ) - P*R
        ENDIF
        EVEC( II, KA + 1 ) = EVEC( II, KA + 1 ) - P*Q
        EVEC( II, KA ) = EVEC( II, KA ) - P
     ENDDO
  ENDDO

  GO TO  410

    ! All evals found, now backsubstitute real vector
550 CONTINUE

  IF( RNORM /= ZERO ) THEN
     DO N = M, 1, -1
        N2   = N
        AA( N, N ) = ONE
        DO I = N - 1, 1, -1
           W  = AA( I, I ) - EVAL( N )
           IF( W == ZERO ) W  = TOL*RNORM
           R  = AA( I, N )
           R  = R + SUM(AA(I,N2:N-1)*AA(N2:N-1,N))
           AA( I, N ) = -R / W
           N2   = I
        ENDDO
     ENDDO

     ! End backsubstitution vectors of isolated evals
     DO I = 1, M
        IF( I < L .OR. I > K ) THEN
           EVEC(I,I:M) = AA(I,I:M)
        ENDIF
     ENDDO

     ! Multiply by transformation matrix
     IF( K /= 0 ) THEN
        DO J = M, L, -1
           DO I = L, K
              Z  = ZERO
              DO N = L, MIN( J, K )
                 Z  = Z + EVEC( I, N )*AA( N, J )
              ENDDO
              EVEC( I, J ) = Z
           ENDDO
        ENDDO
     ENDIF
  ENDIF

!   DO I = L, K
!      EVEC(I,1:M) = EVEC(I,1:M) * WKD(I)
!   ENDDO
! Equivalent loop with alternative indexing:
  DO I = 1, M
     EVEC(L:K,I) = EVEC(L:K,I) * WKD(L:K)
  ENDDO

  ! Interchange rows if permutations occurred
  DO I = L-1, 1, -1
     J  = WKD( I )
     IF( I /= J ) THEN
        DO N = 1, M
           REPL   = EVEC( I, N )
           EVEC( I, N ) = EVEC( J, N )
           EVEC( J, N ) = REPL
        ENDDO
     ENDIF
  ENDDO

  DO I = K + 1, M
     J  = WKD( I )
     IF( I /= J ) THEN
        DO N = 1, M
           REPL   = EVEC( I, N )
           EVEC( I, N ) = EVEC( J, N )
           EVEC( J, N ) = REPL
        ENDDO
     ENDIF
  ENDDO

!  normalise eigenvector
  DO J = 1, M
    EVEC(:,J) = EVEC(:,J) / SQRT(DOT_PRODUCT(EVEC(:,J), EVEC(:,J)))
  ENDDO
  CATCH
END SUBROUTINE ASYMTX

!> TL of ASYMTX - this code was taken from CRTM with only superficial modifications.
!! @param[out]      err     status on exit
!! @param[in]       n       size of matrix
!! @param[in]       v       normalised eigenvectors (column j corresponds to val(j))
!! @param[in]       val     eigenvalues
!! @param[in]       a_tl    input matrix perturbations
!! @param[in,out]   v_tl    eigenvector perturbations
!! @param[in,out]   val_tl  eigenvalue perturbations
SUBROUTINE ASYMTX_TL(err, n, v, val, a_tl, v_tl, val_tl)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: v(n,n), val(n), a_tl(n,n)
  REAL(jprb),    INTENT(INOUT) :: v_tl(n,n), val_tl(n)

  REAL(jprb)    :: v_int(n,n), a1_tl(n,n)
  REAL(jprb)    :: b_tl
  INTEGER(jpim) :: i, j

  TRY

  v_int = matinv(v(1:n,1:n), err)
  THROWM(err.NE.0,"Error in matrix inversion")

  ! part 1
  !!  a1_tl = 0._jprb
  a1_tl = MATMUL(v_int, MATMUL(a_tl, v))

  ! part 2
  !!  a1_tl = 0._jprb
  !!  v_tl = 0._jprb
  DO i = 1, n
    val_tl(i) = a1_tl(i,i)
    DO j = 1, n
     IF (i /= j) THEN
       IF (ABS(val(j) - val(i)) > eigen_threshold) THEN
         v_tl(i,j) = a1_tl(i,j) / (val(j) - val(i))
       ELSE
         v_tl(i,j) = 0._jprb
       ENDIF
     ELSE
       v_tl(i,j) = 0._jprb
     ENDIF
    ENDDO
  ENDDO

  ! part 3
  !!  a1_tl = 0._jprb
  a1_tl = MATMUL(v, v_tl)
  !!  v_tl = 0._jprb

  ! part 4
  DO i = 1, n
    b_tl = 0._jprb
    DO j = 1, n
      b_tl = b_tl - v(j,i) * a1_tl(j,i)
    ENDDO

    DO j = 1, n
      v_tl(j,i) = a1_tl(j,i) + v(j,i) * b_tl
    ENDDO
  ENDDO

  CATCH
END SUBROUTINE ASYMTX_TL

!> AD of ASYMTX - this is the AD written directly from ASYMTX_TL.
!! @param[out]      err     status on exit
!! @param[in]       n       size of matrix
!! @param[in]       v       normalised eigenvectors (column j corresponds to val(j))
!! @param[in]       val     eigenvalues
!! @param[in]       v_ad    eigenvector increments
!! @param[in]       val_ad  eigenvalue increments
!! @param[in,out]   a_ad    matrix increments
SUBROUTINE ASYMTX_AD(err, n, v, val, v_ad, val_ad, a_ad)

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: n
  REAL(jprb),    INTENT(IN)    :: v(:,:), val(:), v_ad(:,:), val_ad(:)
  REAL(jprb),    INTENT(INOUT) :: a_ad(:,:)

  INTEGER(jpim) :: i, j
  REAL(jprb)    :: v_int(n,n), a1_ad(n,n), vtmp_ad(n,n), b_ad

  TRY
  a1_ad = 0._jprb
!part4
  b_ad = 0._jprb
  DO i = n, 1, -1
    DO j = n, 1, -1
!       v_tl(j,i) = a1_tl(j,i) + v(j,i) * b_tl
      a1_ad(j,i) = a1_ad(j,i) + v_ad(j,i)
      b_ad = b_ad + v(j,i) * v_ad(j,i)
    ENDDO

    DO j = n, 1, -1
!       b_tl = b_tl - v(j,i) * a1_tl(j,i)
      a1_ad(j,i) = a1_ad(j,i) - v(j,i) * b_ad
    ENDDO

    b_ad = 0._jprb
  ENDDO

!part3
! a1_tl = MATMUL(v, v_tl)
  vtmp_ad = MATMUL(TRANSPOSE(v), a1_ad)
  a1_ad = 0._jprb
!part2
  DO i = n, 1, -1
    DO j = n, 1, -1
     IF (i /= j) THEN
       IF (ABS(val(j) - val(i)) > eigen_threshold) THEN
!          v_tl(i,j) = a1_tl(i,j) / (val(j) - val(i))
         a1_ad(i,j) = a1_ad(i,j) + vtmp_ad(i,j) / (val(j) - val(i))
       ELSE
         vtmp_ad(i,j) = 0._jprb
       ENDIF
     ELSE
       vtmp_ad(i,j) = 0._jprb
     ENDIF
    ENDDO
!     val_tl(i) = a1_tl(i,i)
    a1_ad(i,i) = a1_ad(i,i) + val_ad(i)
  ENDDO

  v_int = matinv(v(1:n,1:n), err)
  THROWM(err.NE.0,"Error in matrix inversion")

  a_ad = MATMUL(TRANSPOSE(v_int), MATMUL(a1_ad, TRANSPOSE(v)))
  CATCH
END SUBROUTINE ASYMTX_AD

! !> AD of ASYMTX - this code was taken from CRTM with only superficial modifications.
! !! @param[out]      err     status on exit
! !! @param[in]       n       size of matrix
! !! @param[in]       v       normalised eigenvectors (column j corresponds to val(j))
! !! @param[in]       val     eigenvalues
! !! @param[in]       v_ad    eigenvector increments
! !! @param[in]       val_ad  eigenvalue increments
! !! @param[in,out]   a_ad    matrix increments
! SUBROUTINE ASYMTX_AD(err, n, v, val, v_ad, val_ad, a_ad)
!
!   USE parkind1, ONLY : jpim, jprb
!   IMPLICIT NONE
!   INTEGER(jpim), INTENT(OUT)   :: err
!   INTEGER(jpim), INTENT(IN)    :: n
!   REAL(jprb),    INTENT(IN)    :: v(:,:), val(:), v_ad(:,:), val_ad(:)
!   REAL(jprb),    INTENT(INOUT) :: a_ad(:,:)
!
!   INTEGER(jpim) :: i, j, k
!   REAL(jprb)    :: v_int(n,n), a1_ad(n,n)
!
!   TRY
!   a1_ad = 0._jprb
!   DO i = 1, n
!     DO j = 1, n
!       DO k = 1, n
!         IF (k /= j) THEN
!           IF (ABS(val(j) - val(k)) > eigen_threshold) THEN
!             a1_ad(k,j) = a1_ad(k,j) + v(i,k) * v_ad(i,j) / (val(j) - val(k))
!           ENDIF
!         ENDIF
!       ENDDO
!     ENDDO
!     a1_ad(i,i) = a1_ad(i,i) + val_ad(i)
!   ENDDO
!   v_int = matinv(v(1:n,1:n), err)
!   THROWM(err.NE.0,"Error in matrix inversion")
!
!   a_ad = MATMUL(TRANSPOSE(v_int), MATMUL(a1_ad, TRANSPOSE(v)))
!   CATCH
! END SUBROUTINE ASYMTX_AD

!> Invert a square matrix by the Gauss method. Subroutine taken from CRTM.
!! @param[in]     a     square matrix
!! @param[out]    err   status on exit (i.e. flags singular matrix)
FUNCTION matinv(a, err)
  ! --------------------------------------------------------------------
  ! Compute the inversion of the matrix a
  ! Invert matrix by Gauss method
  ! --------------------------------------------------------------------
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  REAL(jprb),    INTENT(IN)  :: a(:,:)
  INTEGER(jpim), INTENT(OUT) :: err

  INTEGER(jpim) :: n
  REAL(jprb)    :: b(SIZE(a,1),SIZE(a,2))
  REAL(jprb)    :: matinv(SIZE(a,1),SIZE(a,2))
  REAL(jprb)    :: temp(SIZE(a,1))

  REAL(jprb)    :: c, d
  INTEGER(jpim) :: i, j, k, m, imax(1), ipvt(SIZE(a,1))
  ! --------------------------------------------------------------------
  TRY

  b = a
  n = SIZE(a,1)
  matinv = a
  ipvt = (/ (i, i = 1, n) /)
  ! Go into loop- b, the inverse, is initialized equal to a above
  DO k = 1, n
  ! Find the largest value and position of that value
    imax = MAXLOC(ABS(b(k:n,k)))
    m = k - 1 + imax(1)
  ! singular matrix check
    IF (ABS(b(m,k)) <= 1.E-40_jprb) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Singular matrix")
    ENDIF
  ! get the row beneath the current row if the current row will
  ! not compute
    IF (m /= k) THEN
      ipvt((/m,k/)) = ipvt((/k,m/))
      b((/m,k/),:)  = b((/k,m/),:)
    ENDIF
  ! d is a coefficient - brings the pivot value to one and then is applied
  ! to the rest of the row
    d = 1 / b(k,k)
    temp = b(:,k)
    DO j = 1, n
       c = b(k,j) * d
       b(:,j) = b(:,j) - temp * c
       b(k,j) = c
    ENDDO
    b(:,k) = temp * (-d)
    b(k,k) = d
  ENDDO
  matinv(:,ipvt) = b
  CATCH
END FUNCTION matinv

END MODULE rttov_scattering_mod
