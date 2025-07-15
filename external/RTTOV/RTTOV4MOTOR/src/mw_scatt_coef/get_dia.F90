subroutine get_dia ( dia_froz, d_dia_froz, d_min, i_type, i_scat, particle_shape, psd, ll_extend, ll_new_integration)

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
!    Sets up the diameter range for calculations. Not appropriate in
!    the totalice case, where predict_psd_f07 provides it instead.
!
!    OUT: dia_froz   - particle diameter or maximum dimension [m]
!         d_dia_froz - integration weights == d_dia only for old rectangle rule integration [m]
!         d_min      - minimum diameter/size in database [m]
!    IN:  i_type     - hydrometeor type (see mod_scattering.F90)   
!         i_scat     - number indicates use of Mie, Liu (2008) DDA or other database scattering params 
!         particle_habit  - habit for above, as defined in Liu (2008) or other database
!         psd
!         ll_extend  - extend integration range below minimum of database particle
!         ll_new_integration - set up integration weights for trapezium rule on logarithmic size basis

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           02/03/2010  New function (Alan Geer)
!           15/03/2013  Fully-flexible Liu shapes (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok). d_dia bugfix.
!           26/03/2020  Allow extension below d_min for database particles, SI (Alan Geer)

use parkind1, only: jprb, jpim, jplm
use mod_scattering,  only: n_dia
!INTF_OFF
use mod_scattering,  only: d_min_sphere, d_max_sphere, d_liu_max, d_liu_min, &
 & i_liu, i_arts, d_field_min, psd_field_2007
use mod_arts, only: d_arts_max, d_arts_min
!INTF_ON
implicit none

! Interface
real (kind=jprb),    intent (  out) :: dia_froz(n_dia), d_dia_froz(n_dia), d_min
integer (kind=jpim), intent (in   ) :: i_type, i_scat, particle_shape
integer (kind=jpim), intent (in   ) :: psd
logical (kind=jplm), intent (in   ) :: ll_extend, ll_new_integration
!INTF_END

! Local variables
integer (kind=jpim) :: i_dia
real (kind=jprb)    :: d_dia_one, d_max_use, d_min_use

!* size ranges (computed later in the totalice case)
if ((i_scat == i_liu) .or. (i_scat == i_arts)) then
  if (i_scat == i_liu) then
    ! Apply size limits appropriate to Liu shapes 
    d_max_use = d_liu_max(particle_shape) 
    d_min     = d_liu_min(particle_shape)
  elseif (i_scat == i_arts) then
    ! ARTS database
    d_max_use = d_arts_max(particle_shape)
    d_min     = d_arts_min(particle_shape)
  endif
  if (ll_extend) then
    ! Extend below d_min using solid spheres
    d_min_use = d_min_sphere (i_type)
  else
    if (psd == psd_field_2007) then
      ! LWC -> moment conversion not valid below 100 micron
      d_min_use = max(d_field_min, d_min)
    else
      ! truncate PSD at bottom of allowed size range
      d_min_use = d_min
    endif
  endif
else
  ! Size ranges for Mie spheres
  d_max_use = d_max_sphere (i_type)
  d_min_use = d_min_sphere (i_type)
endif

if (.not. ll_new_integration) then

  ! Old regular-spaced rectangle-rule integration
  d_dia_one = (d_max_use - d_min_use) / n_dia 
  ! d_dia_froz = (d_max_use - d_min_use) / REAL(n_dia - 1) ! AJGDB Technically correct, but does not reproduce earlier versions

  do i_dia = 1, n_dia
    dia_froz (i_dia) = d_min_use + (i_dia - 1) * d_dia_one
    d_dia_froz(i_dia) = d_dia_one 
  enddo

else

  ! log spacing 
  d_dia_one = (dlog10(d_max_use) - dlog10(d_min_use)) / (n_dia - 1)
  dia_froz (1)     = d_min_use
  dia_froz (n_dia) = d_max_use
  do i_dia = 2, n_dia-1
    dia_froz (i_dia) = 10**(dlog10(d_min_use) + (i_dia - 1) * d_dia_one)
  enddo

  ! weights for trapezium rule integration
  d_dia_froz = 0.0_jprb
  do i_dia = 1, n_dia-1
    d_dia_one = dia_froz (i_dia+1) - dia_froz (i_dia)
    d_dia_froz(i_dia  ) = d_dia_froz(i_dia  ) + 0.5_jprb * d_dia_one
    d_dia_froz(i_dia+1) = d_dia_froz(i_dia+1) + 0.5_jprb * d_dia_one
  enddo

endif

return
end subroutine get_dia
