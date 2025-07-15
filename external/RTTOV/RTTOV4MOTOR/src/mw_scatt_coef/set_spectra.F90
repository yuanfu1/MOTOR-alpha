subroutine set_spectra (i_type, lwc, temp, dia_froz, d_dia_froz, mass, d_equiv, &
  & nd, psd, a, b, ll_dsd, regime, &
  & max_renorm, mod_gamma_data, diag)

! Description:
!
!   Computes the particle size distribution appropriate to a given water content. Many 
!   options are available:
!
!      Field et al. (2005) - only for totalice (option 3)
!      Field et al. (2007) - only for totalice (option 4) or Liu DDA "snow"
!      Marshall-Palmer     - snow, rain, graupel, aggregate, totalice (option 2)
!      Modified Gamma      - CLW, CIW, totalice (option 1)
!
!   In most cases, the particle density relation (which may be constant or a function of 
!   particle size) m=D^b is a basic input, and is configurable for totalice.
!
!   IN: i_type     - hydrometeor type (see mod_scattering.F90)   
!       lwc        - liquid water content [kg m^-3]
!       temp       - temperature          [K]
!       dia_froz   - particle diameter or maximum dimension [m]
!       d_dia_froz - width of integration element [m]
!       mass       - particle mass [kg]
!       d_equiv    - particle equivalent diameter [m]
!       psd        - PSD parametrization
!       regime     - Field et al. (2007) PSD regime
!       ll_dsd     - Panegrossi et al. n0 vs T
!       i_scat     - scattering computation type (see mod_scattering.F90)
!       max_renorm - maximum renormalisation of PSD (e.g. where truncated)
!         
!   OUT: nd        - size distribution [m^-4]
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           12/03/2010  Put all size distributions here (Alan Geer)
!           09/02/2011  Simplify normalisation. Documentation. (Alan Geer)
!           15/03/2013  Fully-flexible PSD, density and shape options (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!           06/03/2020  Refactor / SI / more flexible / new PSDs (Alan Geer)

use parkind1, only: jprb, jpim, jplm
use mod_scattering,  only: n_dia, modified_gamma_data, diag_data
!INTF_OFF
use mod_scattering,  only: i_aggregate, &
               & psd_gamma, psd_marshall_palmer, psd_field_2005, psd_field_2007, &
               & modified_gamma_lambda, modified_gamma_n0, is_frozen, n0_frozen, n0_liquid, &
               & i_clw, i_ciw, modified_gamma_fixed, psd_modified_gamma, &
               & psd_heymsfield_2013, psd_mh97, t_0c
!INTF_ON

implicit none

!* common variables
real (kind=jprb),    intent ( in)   :: lwc, temp, max_renorm
real (kind=jprb),    intent ( in)   :: dia_froz(n_dia), d_dia_froz(n_dia)
real    (kind=jprb), intent ( in)   :: mass(n_dia), d_equiv(n_dia)
integer (kind=jpim), intent ( in)   :: psd, i_type
real    (kind=jprb), intent ( in)   :: a, b
logical (kind=jplm), intent ( in)   :: ll_dsd
character,           intent ( in)   :: regime
type (modified_gamma_data), intent(in) :: mod_gamma_data
real (kind=jprb),    intent (out)   :: nd(n_dia) 
type (diag_data), optional, intent(inout) :: diag

!INTF_END

!* local variables
real    (kind=jprb):: lwc_int, lwc_fac, n0, lambda
real    (kind=jprb):: mu, min_lambda, renorm_order
real    (kind=jprb):: integral(n_dia)
integer (kind=jpim):: i_row

#include "predict_psd.interface"
#include "predict_psd_F07.interface"
#include "modified_gamma.interface"
#include "n0_t.interface"
#include "make_psd_mh97.interface"

if (psd == psd_field_2005 ) then

  call predict_psd(lwc, temp, dia_froz, a, b, nd)

else if (psd == psd_field_2007) then

  call predict_psd_F07(lwc, temp, dia_froz, a, b, nd, regime)

elseif (psd == psd_marshall_palmer) then

  if (i_type /= i_aggregate .and. ll_dsd) then
    ! Panegrossi et al. parametrisation of N0 as function of T
    call n0_t (i_type, lwc, temp, n0)
  else
    if(is_frozen(i_type)) then
      n0 = n0_frozen
    else
      n0 = n0_liquid
    endif
  endif

  call modified_gamma(n_dia, dia_froz, lwc, modified_gamma_lambda, n0, 0.0_jprb, 0.0_jprb, 1.0_jprb, a, b, nd)

elseif (psd == psd_gamma) then

  if (i_type == i_clw) then
    lambda = 2.13e5_jprb
  elseif (i_type == i_ciw) then
    lambda = 2.05e5_jprb
  else
    write(0,*) 'No lambda parameter provided'
    stop
  endif

  ! For pre-v13 back-compatability, the approximate n0 setting used before
  ! (this is not quite what we'd get from PH'11 equation 31)
  n0 = 1.3655e27_jprb * lwc

  call modified_gamma(n_dia, dia_froz, lwc, modified_gamma_fixed, n0, 2.0_jprb, lambda, 1.0_jprb, a, b, nd)

elseif (psd == psd_modified_gamma) then

  i_row = mod_gamma_data%i_row(mod_gamma_data%i_mtype)
  if((i_row < 1) .or. (size(mod_gamma_data%i_free) < 1)) then
    write(0,*) 'Incorrect modified gamma table row, or table size zero'
    stop
  endif

  call modified_gamma(n_dia, dia_froz, lwc, mod_gamma_data%i_free(i_row), &
    & mod_gamma_data%params(1,i_row), mod_gamma_data%params(2,i_row), &
    & mod_gamma_data%params(3,i_row), mod_gamma_data%params(4,i_row), &
    & a, b, nd)

elseif (psd == psd_heymsfield_2013) then

  n0 = 0.0_jprb ! This is the free parameter

  if( regime == 'S' ) then

    ! Stratiform
    lambda     = 1530_jprb * exp (-0.053_jprb *(temp - t_0c))
    min_lambda = max(lambda,-0.0056_jprb) 
    mu         = 1.5_jprb * (min_lambda / 100_jprb) ** 0.223_jprb - 3.0_jprb

  elseif (regime == 'C') then

    ! Convective
    lambda     = 340_jprb * exp (-0.083_jprb *(temp - t_0c))
    min_lambda = max(lambda,-0.0146_jprb) 
    mu         = 0.6_jprb * (min_lambda / 100_jprb) ** 0.382_jprb - 3.0_jprb

  elseif (regime == 'A') then
   
    ! Lambda "composite" / mu "all"
    if (temp < (t_0c - 60.0_jprb)) then
      lambda   = 988_jprb * exp (-0.060_jprb *(temp - t_0c))
    else
      lambda   = 75_jprb * exp (-0.1057_jprb *(temp - t_0c))
    endif
    min_lambda = max(lambda,-0.0094_jprb) 
    mu         = 0.9_jprb * (min_lambda / 100_jprb) ** 0.308_jprb - 3.0_jprb

  else
    write(*,*) 'Unrecognised regime for the Heymsfield (2013) PSD: ', regime
    stop
  endif

  call modified_gamma(n_dia, dia_froz, lwc, modified_gamma_n0, &
    & n0, mu, lambda, 1.0_jprb, a, b, nd)

elseif (psd == psd_mh97) then

  call make_psd_mh97(lwc, temp, d_equiv, dia_froz, nd)

else
  write(*,*) 'Unrecognised PSD: ', psd
  stop
endif

! Normalise in case the range of dia_frozen omits part of the size distribution
! but double check the numerical integration is going to be good enough and/or
! the integration limits have been badly specified for this size distribution
integral = nd * mass * d_dia_froz
lwc_int = sum(integral)
lwc_fac = lwc / lwc_int

renorm_order = abs(dlog10(lwc_fac))
if ( renorm_order > max_renorm ) then
  write(*,'(A,E12.5,A,E12.5)') 'Renormalisation of ',renorm_order, &
     & ' orders of magnitude is greater than permitted ',max_renorm
  write(*,*) ' which suggests a fundamental size distribution problem: '
  write(*,*) minval(dia_froz), maxval(dia_froz), lwc, lwc_int, i_type
  stop 
endif 

nd(:) = lwc_fac * nd(:)

if(present(diag)) then
  diag%renorm  = lwc_fac - 1.0_jprb
endif

return
end subroutine set_spectra
