subroutine modified_gamma (n_dia, d, lwc, free_parameter, n0, mu, lambda, gamma_in, a, b, nd)

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

! This program calculates all size distributions fitting within the modified gamma
! family (e.g. this includes Marshall-Palmer) as described by Petty and Huang (2011, JAS, PH11).

! The size descriptor for the PSD is the geometric diameter of PH11, equivalently the
! maximum dimension of the particle (which is its diameter, if a sphere)

! Given a moment of the PSD (in this case, LWC) one free parameter of the PSD can be estimated.

! All units are in SI

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           06/03/2020   Initial version (Alan Geer)

use parkind1, only: jprb, jpim
!INTF_OFF
use mod_scattering, only: modified_gamma_n0, modified_gamma_lambda, modified_gamma_fixed
!use mod_scattering, only: not_available, modified_gamma_mu, modified_gamma_gamma
!INTF_ON

implicit none

!* common variables
integer (kind=jpim), intent ( in) :: n_dia
real    (kind=jprb), intent ( in) :: d(n_dia)                 ! Maximum dimension of particle [m]
real    (kind=jprb), intent ( in) :: lwc                      ! Water content [kg/m^3]
integer (kind=jpim), intent ( in) :: free_parameter           ! Free parameter of the PSD (see mod_scattering.F90)
real    (kind=jprb), intent ( in) :: n0, mu, lambda, gamma_in ! Parameters of the PSD [SI]
real    (kind=jprb), intent ( in) :: a, b                     ! Mass-size relation m=a.d^b
real    (kind=jprb), intent (out) :: nd(n_dia)                ! PSD [m^-4]
!INTF_END

!* local variables	
real    (kind=jprb) :: gamma_term, lambda_exponent, use_n0, use_lambda
integer (kind=jpim) :: i

use_n0     = n0
use_lambda = lambda

lambda_exponent = (mu + b + 1.0_jprb) / gamma_in
gamma_term      = gamma(lambda_exponent)

! Here we use the LWC (which is the bth moment of the PSD) to find a free parameter
! via PH11 equation 31. Note that only n0 and lambda are easy to get and hence mu and
! gamma are not yet allowed to be free parameters here.
if (free_parameter == modified_gamma_n0) then

  use_n0 = lwc * gamma_in * (lambda ** lambda_exponent) / (a * gamma_term)

elseif (free_parameter == modified_gamma_lambda) then

  use_lambda = ((a * n0 * gamma_term) / (gamma_in * lwc)) ** (1.0_jprb/lambda_exponent)

elseif (free_parameter == modified_gamma_fixed) then

  ! Do nothing - use all input settings

else
  write(*,*) 'Unsupported free parameter in modified gamma PSD'
  stop
endif

do i = 1, n_dia
  nd(i) = use_n0 * (d(i)**mu) * exp(-1.0_jprb * use_lambda * (d(i)**gamma_in))
enddo

end subroutine modified_gamma
