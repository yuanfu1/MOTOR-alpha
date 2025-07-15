subroutine make_psd_mh97 (lwc, temp, d_equiv, d_max, nd)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2020, EUMETSAT, All Rights Reserved.


! Implements the mass-based PSD of:

! McFarquhar, Greg M., and Andrew J. Heymsfield. Journal of the atmospheric sciences 54, no. 17 (1997): 2187-2200
! Parameterization of tropical cirrus ice crystal size distributions and implications for 
! radiative transfer: Results from CEPEX

! PSD is converted into a max-dimension basis before returning. Although this could be done analytically using
! the mass-dimension relations for most sizes, in the small-particle limit the particle is solid and these do
! not apply. Hence it is more flexible and robust to convert using numerical differentiation by calculating
! delta(d_e) / delta(d_g). Note that in creating most other PSDs the solid particle limit is ignored and an
! analytical calculation of LWC is done anyway.

! All in/out units are SI. 

! Note density of ice is set from mod_scattering, currently 917 kg/m3, not 910 as in the MH97 paper. This
! is used in both parts of the PSD (below) and in defining the mass equivalent diameter but the difference 
! between 910 and 917 does not appear to have much effect on the PSD in either basis.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           27/03/2020  First version (Alan Geer)

use parkind1, only: jprb
use mod_scattering, only: n_dia
!INTF_OFF
!use parkind1, only: jpim
use mod_scattering, only: pi, density_ice, t_0c
!INTF_ON

implicit none

real    (kind=jprb), intent ( in) :: lwc            ! water content [kg m^-3]
real    (kind=jprb), intent ( in) :: temp           ! temperature [K]
real    (kind=jprb), intent ( in) :: d_equiv(n_dia) ! Mass equivalent diameter [m]
real    (kind=jprb), intent ( in) :: d_max(n_dia)   ! Max diameter enclosing sphere [m]
real    (kind=jprb), intent (out) :: nd(n_dia)      ! D_g (max diameter) size distribution [m^-4]

!INTF_END
 
real (kind=jprb), parameter :: iwc_0 = 1e-3_jprb ! [kg m^-3]
real (kind=jprb), parameter :: d_0   = 1e-6_jprb ! [m]

! Lookup tables every 10 degrees from -70 to -20 and one bin -20 to 0 Celcius
! Actually these are superseded by Eq. 9-12, but having typed them in, they can stay for now.
!integer (kind=jpim), parameter :: ntbins = 6
!real (kind=jprb), parameter, dimension(ntbins) :: a_mu  = (/ 5.156, 5.131, 5.148, 5.148, 5.146, 5.203 /)
!real (kind=jprb), parameter, dimension(ntbins) :: b_mu  = (/ 0.091, 0.093, 0.089, 0.071, 0.027, 0.042 /)
!real (kind=jprb), parameter, dimension(ntbins) :: a_sig = (/ 0.370, 0.348, 0.396, 0.397, 0.440, 0.447 /)
!real (kind=jprb), parameter, dimension(ntbins) :: b_sig = (/ 0.030, 0.022, 0.044, 0.040, 0.038, 0.008 /)

real (kind=jprb) :: ne_small(n_dia), ne_large(n_dia) ! D_equiv size distribution, small and large [m^-4]
real (kind=jprb) :: iwc_small, iwc_large             ! Below and above 100 micron wc [kg m^-3]
real (kind=jprb) :: alpha_small, mu_large, sigma_large, iwc_ratio
real (kind=jprb) :: a_mu, b_mu, a_sig, b_sig, t_c
real (kind=jprb) :: part1(n_dia), part2(n_dia)
real (kind=jprb) :: de_dg(n_dia)

t_c = temp - t_0c

! Index in temperature bins 
!i_temp = min( ntbins, max( 1_jpim , 1_jpim + floor((temp - t_0c + 70.0_jprb) / 10_jpim)))

! Eq. 5
iwc_small = min(lwc, 0.252e-3_jprb * (lwc / iwc_0) ** 0.837_jprb)
iwc_large = lwc - iwc_small

! -- Small PSD (gamma)

! Eq. 6
! But limited to prevent it going negative for high water contents (-> -ve psd) - limited to min observed on Fig. 6
alpha_small = max( 5e4_jprb, -4.99e3_jprb -49.4e3_jprb * dlog10(iwc_small/iwc_0))

! Eq.  3 (note 6 / gamma(5) is 1/4 since gamma(5) = 24)
ne_small = 0.25_jprb / pi * iwc_small / density_ice * alpha_small ** 5.0_jprb & 
 & * d_equiv * exp ( -1.0_jprb * alpha_small * d_equiv)

! -- Large PSD (log normal)

iwc_ratio = dlog10(iwc_large/iwc_0)

! Eq. 9-12, 7 & 8
a_mu  = 5.200_jprb + 1.3e-3_jprb * t_c
b_mu  = 0.026_jprb - 1.2e-3_jprb * t_c 
a_sig = 0.470_jprb + 2.1e-3_jprb * t_c
b_sig = 0.018_jprb - 2.1e-4_jprb * t_c
mu_large = a_mu + b_mu * iwc_ratio
sigma_large = a_sig + b_sig * iwc_ratio

! Eq. 4 
part1 = pi**(1.5_jprb) * density_ice * sqrt(2.0_jprb) * &
  &     exp(3.0_jprb * mu_large + 4.5_jprb * sigma_large**2.0_jprb) * &
  &     d_equiv * sigma_large * d_0 ** 3.0_jprb

part2 = exp( -0.5_jprb * ( (dlog(d_equiv/d_0) - mu_large) / sigma_large) **2.0_jprb)

ne_large = 6.0_jprb * iwc_large / part1 * part2

! Combine both parts of PSD
nd = (ne_small + ne_large) 

! Numerical differentation of d_equiv vs d_max (d_e vs d_g in PH11 terms)
! Where possible the centred (average) differential is calculated, but at the edges uncentred.
de_dg(n_dia)     = 0.0_jprb
de_dg(1:n_dia-1) = (d_equiv(2:n_dia) - d_equiv(1:n_dia-1))/(d_max(2:n_dia) - d_max(1:n_dia-1))
de_dg(2:n_dia)   = de_dg(2:n_dia) + de_dg(1:n_dia-1)
de_dg(2:n_dia-1) = de_dg(2:n_dia-1)/2.0_jprb

! Convert to d_g (max diameter) basis using d(d_e)/d(d_g)
nd = nd * de_dg

return
end subroutine make_psd_mh97
