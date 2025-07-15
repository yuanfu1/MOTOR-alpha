subroutine perm_ice (f_ghz, t_k, m06, perm_re, perm_im)

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

! Calculates the permittivity of ice. Original (pre v13) formula is described as:
!   "Hufford (1991), Mishima et al. (1983) modified by Maetzler and Wegmueller (1987)"

! This was near-identical to the "Maetzler 2006" (M06) parametrisation used by e.g. ARTS & Liu database etc.
!   Mätzler, C.: Microwave dielectric properties of ice., 2006. Thermal microwave radiation - Applications for remote sensing, 
!   Electromagn., Waves Ser., vol. 52, edited by: C. Mätzler, Inst. Eng. Technol., Stevenage, UK, Sect. 5.3, 455-462

! From v13 an updated "dbeta" term and perm_re is available to make it precisely as Maetzler '06

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           27/03/2020   Minor update to Maetzler (2006) version (Alan Geer)

use parkind1, only: jprb, jplm

implicit none

!* common variables
real    (kind=jprb), intent (in)  :: f_ghz, t_k
logical (kind=jplm), intent (in)  :: m06
real    (kind=jprb), intent (out) :: perm_re, perm_im

!INTF_END

!*local variables
real (kind=jprb) :: tk, tc, theta, alpha, fac, betam, dbeta, beta

tk = min (t_k, 273.16_jprb)
tc = tk - 273.16_jprb

! M06 Eq. 5Q1
if(m06) then
  perm_re = 3.1884_jprb + 9.1e-04_jprb * (tk - 273.0_jprb)
else
  ! Note implied change to 273.16 not the 273 in M06 (above)
  perm_re = 3.1884_jprb + 9.1e-04_jprb * tc
endif

! M06 Eq. 5Q3
theta = 300.0_jprb / tk
alpha = 1.0e-04_jprb * (50.4_jprb + 62.0_jprb * (theta - 1.0_jprb)) * exp (-22.1_jprb * (theta - 1.0_jprb))

! M06 Eq. 5Q4
fac   = exp (335.0_jprb / tk)
betam = 0.0207_jprb * fac / (tk * (fac - 1.0_jprb) * (fac - 1.0_jprb)) + 1.16e-11_jprb * f_ghz  * f_ghz

if(m06) then
  ! M06 Eq. 5Q6
  dbeta = exp (-9.963_jprb + 0.0372_jprb * tc)
else
  ! Maetzler (1998) version as referred to in M06 (this used by pre-v13 RTTOV-SCATT)
  dbeta = exp (-10.02_jprb + 0.0364_jprb * tc)
endif

! M06 Eq. 5Q5
beta  = betam + dbeta

! M06 Eq. 5Q2
perm_im = alpha / f_ghz + beta * f_ghz

return
end subroutine perm_ice
