subroutine density_all (i_type, i_scat, particle_shape, dens, dia_froz, ll_ignore_density_floor, &
  & density, mass, d_equiv, a, b)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.


! Centralises all density-related computations. Outputs are (1) density
! itself and (2) the a and b coefficients of the mass-size relation m = ad^b.
!
! All units are SI
!
! The particle size/mass coordinate is assumed to be the maximum dimension (or diameter, if a sphere)
!
! All density and mass properties are calculated from a and b so density relations
! like Brown/Francis'95 differ numerically from their canonical forms (e.g. not exactly rho=0.035d^-1.1)
!
! One density relation is included that is not a power law. This cannot be combined with PSDs that
! are evaluated assuming m=aD^b power laws, but should work with with Mie Spheres.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           03/2013     First version (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!           06/03/2020  Refactor / SI / more flexible / new PSDs (Alan Geer)

use parkind1, only: jprb, jpim, jplm
use mod_scattering, only: n_dia
!INTF_OFF
use mod_scattering, only:  density_sphere, pi, dens_wilson_ballard_1999, dens_jones_1995, &
 & dens_brown_francis_1995, dens_fixed, i_clw, i_rain, not_available, i_liu, i_arts, &
 & density_ice, density_wat, is_frozen
use mod_arts, only: alpha_arts, beta_arts
!INTF_ON

implicit none

integer (kind=jpim), intent ( in) :: i_type            ! Hydrometeor type
integer (kind=jpim), intent ( in) :: i_scat            ! Scattering computation type
integer (kind=jpim), intent ( in) :: particle_shape    ! particle shape (e.g. Liu DDA habit number)
integer (kind=jpim), intent ( in) :: dens              ! Choice of density parametrisation
logical (kind=jplm), intent ( in) :: ll_ignore_density_floor ! Option for pre-v13 replication decreases 
                                                       ! ext for v low wcs by allowing super dense particles.
real    (kind=jprb), intent ( in) :: dia_froz(n_dia)   ! Diameter or maximum dimension [m]
real    (kind=jprb), intent (out) :: density(n_dia)    ! Density wrt to circumscribed sphere [kg/m3]
real    (kind=jprb), intent (out) :: mass(n_dia)       ! Particle mass [kg]
real    (kind=jprb), intent (out) :: d_equiv(n_dia)    ! Mass-equivalent diameter [m]
real    (kind=jprb), intent (out) :: a, b              ! mass = aD^b [SI]
!INTF_END

integer (kind=jpim) :: i_dia
real    (kind=jprb) :: x, y
real    (kind=jprb) :: material_density  ! Density of material composing particle [kg/m3]

real    (kind=jprb) :: density_ice_limit

#include "liu_density.interface"

x = not_available
y = not_available

! The formulae are all in SI units
if ( (i_scat == i_liu) .or. (i_scat == i_arts) ) then

  if (i_scat == i_liu) then
    ! Liu shapes (a,b are in SI.). See p. 27-139
    call liu_density(particle_shape, x, y)
  else
    x = alpha_arts(particle_shape)
    y = beta_arts(particle_shape)
  endif

elseif (dens == dens_wilson_ballard_1999) then

  ! set x and y to be values used in UM (density = 0.132*D-1)
  x=0.069
  y=2.0

elseif (dens == dens_jones_1995) then

  ! Power law mass-size (x, y) are not valid here!

elseif (dens == dens_brown_francis_1995) then

  ! set x and y to be the values equivalent to the Brown and Francis rho = 0.035*D-1.1
  x = 0.0185_jprb
  y = 1.9

elseif (dens == dens_fixed) then

  ! Fixed density sphere 
  x = density_sphere (i_type) * pi / 6.0_jprb
  y = 3.0_jprb

else

  write(*,*) 'Invalid density / mass-size relation / PSD'
  stop

endif

a=x
b=y

if(is_frozen(i_type)) then
  material_density=density_ice
else
  material_density=density_wat
endif

if (dens == dens_jones_1995) then

  ! density = 8.74E-4*exp(-0.625D2) + 4.5E-5
  density(:) = (0.000874_jprb*exp(-0.625_jprb*(dia_froz(:)*dia_froz(:))) + 0.000045_jprb)

else

  ! All power-law m-D relations
  do i_dia = 1, n_dia
    density(i_dia) = 6.0_JPRB/pi*x*(dia_froz(i_dia))**(y-3.0_JPRB)
  enddo

endif

! an upper limit so water-air fraction of the mie sphere cannot be > 1 (see e.g. Brown and Francis 1995 or PH'11 section 4b)
if(i_type /= i_clw .and. i_type /= i_rain) then
  ! density_ice_limit = density_ice ! AJGDB what it should be ! AJGDB do water too...
  density_ice_limit = 920_jprb ! AJGDB backward consistency
  where (density > density_ice_limit) density = density_ice_limit
endif

do i_dia = 1, n_dia

  ! This could be done with a and b more easily, but only if density did not have a floor
  ! To retain numerical equivalence with pre-v13 temporarily, this should be applied for F07 particles:
  if (ll_ignore_density_floor) then
    mass(i_dia) = x * dia_froz(i_dia) ** y
  else

    ! Particle mass
    mass(i_dia) = pi/6.0_jprb*density(i_dia)*dia_froz(i_dia)**3.0_jprb

  endif
enddo

! Mass equivalent sphere diameter
d_equiv = (6.0_jprb * mass / material_density / pi)**(1.0_jprb/3.0_jprb)

return
end subroutine density_all
