subroutine predict_psd_F07(iwcsi,tk_in,D,x,y,dN_dD,regime)

! Description:
!   Subroutine copied from Paul fields pvwave code compute_psd.pro
!   fortran code obtained here 
!   http://www.mmm.ucar.edu/people/field/files/web/Refereed%20publications.html
!   This routine uses the second moment and temperature to predict the
!   ice PSD.
!   The moment relations use the 2nd moment as a reference moment.If
!   IWC is computed using a different moment then the second moment 
!   needs to be estimated first.
!   Subroutine predict_mom07.f90 below carries out this step.
!   See 'Snow size distribution parameterization for midlatitude and 
!   tropical ice clouds' by Field et al JAS2007
!
! Copyright:
!   This software was developed within the context of
!   the EUMETSAT Satellite Application Facility on
!   Numerical Weather Prediction (NWP SAF), under the
!   Cooperation Agreement dated 7 December 2016, between
!   EUMETSAT and the Met Office, UK, by one or more partners
!   within the NWP SAF. The partners in the NWP SAF are
!   the Met Office, ECMWF, DWD and MeteoFrance.
!
!   Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008  transferred to fcm (Amy Doherty)
!           15/01/2009  updated to Field et al 2007 version (Amy Doherty)
!           02/03/2011  normalise (Alan Geer)
!           02/03/2011  use with Liu particle habits (Alan Geer)
!           12/05/2011  Corrected density parameters. Andrew Smith
!           15/03/2013  Fully-flexible PSD, density and shape (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!           20/03/2020  SI and more streamlined code (Alan Geer)
!
! Code Description: 
!   Language:           Fortran 90. 
!   Software Standards: "European Standards for Writing and 
!     Documenting Exchangeable Fortran 90 Code". 

use parkind1, only: jprb
use mod_scattering, only:  n_dia
!INTF_OFF
use parkind1, only: jpim
use mod_scattering, only:  not_available, lfield_ignore_temp, f07_fixed_temp
!INTF_ON

implicit none

real     (kind=jprb), intent(in)  :: iwcsi   ! IWC (kg/m^3)
real     (kind=jprb), intent(in)  :: tk_in   ! Temperature (K)
character,            intent(in)  :: regime  ! ('T' Tropical, 'M' mid-latitude)
real     (kind=jprb), intent(in)  :: D(n_dia) ! Particle diameter [m]
real     (kind=jprb), intent(in)  :: x,y      ! Refactor and exponent x, y in mass-size relation m=xD^y [SI]
real     (kind=jprb), intent(out) :: dN_dD(n_dia)   ! Size distribution

!INTF_END

#include "predict_mom07.interface"

! local variables
real     (kind=jprb) :: xx (n_dia)
real     (kind=jprb) :: phi (n_dia)
real     (kind=jprb) :: tk, tc, m3
real     (kind=jprb) :: My, M2, n
integer  (kind=jpim) :: i

!-------------------------------------------------------------------------------

if (x == not_available) then
  write(*,*) 'Invalid density treatment selected for use with Field PSD'
  stop
endif

! define constants for moment relations based on second moment as reference
!----------------------------------------------------------
!a1=5.065339
!a2=-0.062659
!a3=-3.032362
!a4=0.029469
!a5=-0.000285
!a6=0.31255
!a7=0.000204
!a8=0.003199
!a9=0.0
!a10=-0.015952
!
!b1=0.476221
!b2=-0.015896
!b3=0.165977
!b4=0.007468
!b5=-0.000141
!b6=0.060366
!b7=0.000079
!b8=0.000594
!b9=0.0
!b10=-0.003577
!
! check in cloud temperature
!----------------------------------
if(lfield_ignore_temp) then
  tk = f07_fixed_temp
else
  tk = tk_in
endif

! first change from Kelvin to degress C
tc = tk - 273.0_jprb

if (tc.gt.0) then
  write(*,*) 'in-cloud temperature was too warm (>0C)', tc, tk
  stop
endif
if (tc.lt.-70.0_jprb) then
  write(*,*) 'Warning: in-cloud temperature colder than included in original data'
endif

!------------------------------------

! Calculate second moment
My = iwcsi/x           !yth moment

! Use moment relations to get the other moment required for predicting PSD
M2 = My         ! Default for UM y=2 (do not need to do inversion in 
                !following if-block to avoid slight error introduced by plane fitting)

if (y.ne.2.0_jprb) then
  m2 = -9999.0_jprb
  call predict_mom07(m2,tc,y,My)
endif
  
! check M2 range (1E-5 to 3E-2 m-1)

if (M2.lt.0.0_jprb) then
  write(*,*) 'Negative M2'
  stop
endif

!if (m2.lt.1E-5) then
!  write(*,*)'Warning: M2 outside of range encountered in original data - too small', M2
!endif
!if (m2.gt.3e-2) then
!  write(*,*)'Warning: M2 outside of range encountered in original data - too large', M2
!endif

! Use 2nd and third moment to predict PSD
n=3.0_jprb
call predict_mom07(m2,tc,n,m3)

!define dimensionless sizes
do i = 1,n_dia
  xx(i) = (m2/m3) * D(i)
enddo

!define universal psd
do i = 1,n_dia

  if (xx(i) < 20.0_jprb) then

    if (regime.eq.'T') then
      phi(i)=152.0_jprb*exp(-12.4_jprb*xx(i))+3.28_jprb*xx(i)**(-0.78_jprb)*exp(-1.94_jprb*xx(i))
    else if (regime.eq.'M') then
      phi(i)=141.0_jprb*exp(-16.8_jprb*xx(i))+102.0_jprb*xx(i)**(2.07_jprb)*exp(-4.82_jprb*xx(i))
    else
      write(*,*) 'Unknown regime: ', regime
      stop
    endif

  else
    phi(i) = 0.0_jprb
  endif

  !;;scale universal psd;;;;;
  dN_dD(i)=phi(i)*m2**4/m3**3 !psd m^-4

enddo

end subroutine predict_psd_F07
