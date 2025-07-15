subroutine scattering_one_temp ( temp, wavelength, f_ghz, i_type, i_scat, dens, &
  & psd, regime, permwat, permice, ll_melt, ll_dsd, particle_habit, is_loaded, &
  & ll_extend_by_mie, ll_new_integration, max_renorm, &
  & ssp_arts, f_arts, T_arts, D_arts, &
  & mod_gamma_data, &
  & ext_tab, ssa_tab, asm_tab, zef_tab, &
  & diag)

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
! Computes scattering parameters for a particular frequency, 
! temperature, and hydrometeor type, as a function of lwp
!
!    IN: temp       - temperature       [K]
!        wavelength -                   [m]
!        f_ghz      - frequency         [GHz]
!        i_type     - hydrometeor type (see mod_scattering.F90)   
!        i_scat     - scattering computation type (see mod_scattering.F90)
!        ll_melt    - do melting layer of Bauer (2001) 
!        dens       - density parametrization 
!        psd        - PSD parametrization 
!        regime     - Field et al. (2007) PSD regime
!        permwat, permice - permittivity choices - see permittivity.F90
!        ll_melt
!        ll_dsd     - Panegrossi et al. n0 vs T
!        particle_shape  - shape ID for Liu (2008) DDA shapes or Baran ensemble etc.
!        is_loaded  - used by Liu DDA tables
!                     development work. Too slow for operational use.
!        ll_extend_by_mie - For database shapes, below D_min, extend using Mie
!        ll_new_integration - set up integration weights for trapezium rule on logarithmic size basis
!        max_renorm - Max PSD renormalisation factor
!        ssp_arts   - combined SSP from ARTS database for this i_type
!
!    OUT: ext_tab   - Bulk extinction               [m^-1]
!         ssa_tab   - Bulk single scattering albedo [ ]
!         asm_tab   - Bulk asymmetry paramter       [ ]
!         zef_tab   - Bulk radar reflectivity       [m^2] 
!

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/03/2010  Basic routine for parallelisation (Alan Geer)
!           13/03/2013  Fully-flexible PSD, density and shape (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!           10/01/2018  add option to select liquid water permittivity model (Katrin Lonitz) 
!           23/03/2020  SI, new PSDs, diagnostics mode (Alan Geer)

use parkind1, only: jprb, jpim, jplm
use mod_scattering, only: n_lwc, modified_gamma_data, diag_data
!INTF_OFF
use mod_scattering, only: n_dia, get_lwc, ctype, psd_field_2007
!INTF_ON
use mod_arts, only: nf_max_arts, nT_max_arts, nD_max_arts

implicit none

real (kind=jprb),    intent(in   ) :: temp, wavelength, f_ghz
integer (kind=jpim), intent(in   ) :: i_type, i_scat, particle_habit, dens, psd
integer (kind=jpim), intent(in   ) :: permwat, permice
logical,             intent(in   ) :: ll_melt, ll_dsd
character,           intent(in   ) :: regime
integer (kind=jpim), intent(inout) :: is_loaded
logical (kind=jplm), intent(in   ) :: ll_extend_by_mie, ll_new_integration
real (kind=jprb),    intent(in   ) :: max_renorm
real (kind=jprb),    intent(in   ) :: ssp_arts(nf_max_arts, nT_max_arts, nD_max_arts,4)
real (kind=jprb),    intent(in   ) :: f_arts(0:nf_max_arts), T_arts(0:nT_max_arts), D_arts(0:nD_max_arts) 
type (modified_gamma_data), intent(in) :: mod_gamma_data
real (kind=jprb),    intent(  out) :: ext_tab(n_lwc), ssa_tab(n_lwc), asm_tab(n_lwc), zef_tab(n_lwc)
type (diag_data), optional, intent(inout) :: diag
!INTF_END

integer (kind=jpim) :: i_lwc, iostat, i_dia, i_diag(1)
real (kind=jprb), dimension (n_dia) :: q_ext, q_sct, q_asm, q_bsct
real (kind=jprb), dimension (n_dia) :: nd, dia_froz, d_dia_froz, by_dia_density, mass, d_equiv
real (kind=jprb), dimension (n_lwc) :: lwc
real (kind=jprb) :: ext, sct, asm, bsct
real (kind=jprb) :: a, b, d_min

!* Interface blocks
#include "get_dia.interface"
#include "scattering.interface"
#include "set_spectra.interface"
#include "scattering_one_wc.interface"
#include "density_all.interface"

ext_tab = 0.0
ssa_tab = 0.0
asm_tab = 0.0
zef_tab = 0.0

! Diameter range
call get_dia( dia_froz, d_dia_froz, d_min, i_type, i_scat, particle_habit, psd, &
  & ll_extend_by_mie, ll_new_integration)

! Density by size, m=aDb relation
call density_all(i_type, i_scat, particle_habit, dens, dia_froz, &
  & psd == psd_field_2007, & ! AJGDB wrong but for pre-v13 snow compatibility
  & by_dia_density, mass, d_equiv, a, b)

! Compute scattering parameters as a function of diameter.
call scattering( dia_froz, d_min, by_dia_density, temp, wavelength, f_ghz, i_type, &
  & permwat, permice, &
  & i_scat, particle_habit, is_loaded, ll_extend_by_mie, &
  & ssp_arts, f_arts, T_arts, D_arts, &
  & q_ext, q_sct, q_asm, q_bsct)

do i_lwc = 1, n_lwc
  lwc(i_lwc) = get_lwc(i_lwc)
enddo
if (present(diag)) i_diag = minloc(abs(lwc - diag%iwc))

do i_lwc = 1, n_lwc

  ! Get size distribution for this LWC 
  call set_spectra (i_type, lwc(i_lwc), temp, dia_froz, d_dia_froz, &
    & mass, d_equiv, nd, psd, a, b, ll_dsd, regime, &
    & max_renorm, mod_gamma_data, diag=diag)

  ! Integrate over sizes to get bulk layer properties
  call scattering_one_wc( wavelength, f_ghz, i_type, permwat, permice, &
    & (ll_melt .and. nint(temp) == 273), dia_froz, d_dia_froz, nd, &
    & q_ext, q_sct, q_asm, q_bsct, ext, sct, asm, bsct, diag=diag)

  ssa_tab(i_lwc) = ssa_tab(i_lwc) + sct / ext 
  ext_tab(i_lwc) = ext_tab(i_lwc) + ext 
  asm_tab(i_lwc) = asm_tab(i_lwc) + asm
  zef_tab(i_lwc) = zef_tab(i_lwc) + bsct

  ! Diagnostic output
  if (present(diag)) then
    if ( i_lwc == i_diag(1) .and. i_type == diag%itype .and. &
       & abs(temp - diag%temp) < 1e-3_jprb .and. abs(f_ghz - diag%freq) < 1e-3_jprb) then

      open (20, file = diag%filename, iostat=iostat)
      if(iostat /= 0) then
        write(*,*)' Problem opening '//trim(diag%filename)
        exit
      endif

      write(20,*) 'Frequency [GHz]:', f_ghz
      write(20,*) 'Type: ',ctype(i_type)
      write(20,*) 'LWC: [kg m^-3]',lwc(i_lwc)
      write(20,*) 'T [K]', temp
      write(20,*) 'Renormalisation - 1: ',diag%renorm
      write(20,*) 'Bulk extinction, scattering, backscatter cross sections [m^2] / g [0-1]'
      write(20,'(20(e13.5e3))') ext, sct, bsct, asm
      write(20,*) ' D_g [m]       D_m [m]       Mass [kg]     Nd [m^-4]     Extinction ',&
                 &'   Scatter[m^2]  Backsc.[m^2]  Asymmetry     Extegrand'
      do i_dia = 1, n_dia
        write(20,'(20(e14.5e3))')  dia_froz(i_dia), d_equiv(i_dia), mass(i_dia), nd(i_dia), &
         & q_ext(i_dia), q_sct(i_dia), q_bsct(i_dia), q_asm(i_dia), diag%integrand(i_dia)
      enddo
      close(20)

    endif
  endif

end do

end subroutine scattering_one_temp
