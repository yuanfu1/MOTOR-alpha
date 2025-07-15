program rttov_scatt_make_coef

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

! Description:
! Main program for generating tables of scattering coefficients for use with RTTOV.
! Mie theory is used (spherical particles are assumed) as a default
! Coefficients (extinction coefficient, single scattering albedo and asymmetry parameter, backscatter)
! are generated for instruments efined in channels.dat
! at 70 temperatures, 401 water contents, and 100 particle diameters. 
! These tables are then read into rttovscatt
! (available with RTTOV version 9 onwards) in the subroutine rttov_readscattcoeffs.F90.
! 
! for more details on the instruments see channels.dat_all

! Method:
! 

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!   1.2     23/04/2009   changed field PSD treatement to Field '07 instead of Field '05
!                        giving choice between tropical and mid latitude regime (Amy Doherty)
!   1.3     12/03/2010   Optimisation and restructuring (Alan Geer)
!   1.4     02/03/2011   More of the above + extra choices for snow (Alan Geer)
!   1.5     15/03/2013   Fully-flexible PSD, shape and density options (Alan Geer)
!   ?.?     31/10/2017   Adapted for ARTS-SSDB use (Jana Mendrok)
!           10/01/2018   add option to select liquid water absorption model (Katrin Lonitz)
!   2.0     27/03/2020   Complete overhaul especially for PSDs, cloud ice, SI, integration (Alan Geer)

!**************************************************************************

use parkind1, only: jprb, jpim, jplm

use mod_scattering, only: n_temp, n_lwc, pi, i_liu, i_arts, &
                   temp_offset, lmax, i_rain, modified_gamma_data, diag_data, &
                   nmaxconfig, permwat_rosenkranz_15
use mod_arts, only: nf_max_arts, nT_max_arts, nD_max_arts

implicit none

real (kind=jprb)    :: f_ghz, fc_ghz, df_ghz
real (kind=jprb)    :: temp
real (kind=jprb)    :: wavelength, znorm
real (kind=jprb)    :: tq_ext, tq_sct, tq_asm, tq_bsct ! Not used except to help loading Liu(2008) parameters
integer :: nargs
integer (kind=jpim) :: i_freq, j_freq, i_type, i_mtype, i_temp, n_band, i_band, i_row, n_rows
integer (kind=jpim) :: i_foundrain, ipermwat
integer (kind=jpim) :: n_freq, n_sensor, nconfig
integer (kind=jpim) :: is_loaded
complex (kind=jprb) :: k, perm_wat
logical (kind=jplm) :: ll_melt, ll_dsd

integer (kind=jpim) :: mpol ! Polarisation options following RTTOV MPOL, plus -1 signifying polarisation ignored

character(len=lmax) :: line
character(len=lmax) :: output_dir
character(len=lmax) :: input_file

integer (kind=jpim) :: m_type 
integer (kind=jpim), allocatable, dimension (:) :: id_type, i_scat
character, allocatable, dimension (:) :: regime
integer (kind=jpim), allocatable, dimension (:) :: dens, psd, particle_shape
integer (kind=jpim), allocatable, dimension (:) :: iperm, wperm
logical (kind=jplm), allocatable, dimension (:) :: ll_extend_by_mie, ll_new_integration
real    (kind=jprb), allocatable, dimension (:) :: max_renorm

character (len=lmax) :: cconfig(nmaxconfig)
real (kind=jprb), allocatable, dimension (:) :: pol_alpha ! Scaling of v and h extinctions to crudely simulate polarised scattering

real (kind=jprb), dimension (:,:,:,:), allocatable :: ext_tab, ssa_tab, asm_tab, zef_tab
real (kind=jprb), dimension (:), allocatable, save :: ext_wc, ssa_wc, asm_wc, zef_wc

real (kind=jprb), dimension (:,:,:,:,:), allocatable :: ssp_arts
real (kind=jprb), dimension (:,:), allocatable :: f_arts, T_arts, D_arts

type(modified_gamma_data) :: mod_gamma_data

type(diag_data), pointer :: diag=>NULL(), use_diag=>NULL()

#ifndef RTTOV_NAG53
!$OMP THREADPRIVATE(ext_wc, ssa_wc, asm_wc, zef_wc)
#endif

!* Interface blocks
#include "convert_hydrotable.interface"
#include "scattering_one_temp.interface"
#include "permittivity.interface"
#include "liu_dda.interface"
#include "load_arts_ssp.interface"

!-------------------------------------------------------------------------------

!------------------------
! Read configuration file
!------------------------

! Use default file name 'channels.dat' unless one is specified on the command
! line.

nargs = command_argument_count()
if ( nargs > 0 ) then
  call get_command_argument(1,input_file)
else
  write(*,*) 'Please specify a channels.dat file'
  stop 
end if

open(10, file = input_file, status = 'old')

print '(a)', 'Reading configuration file ' // trim(input_file) // '...'

! Look for the header

do
  read (10,*) line

  if ( line == 'DIAG' ) then
    allocate(diag)
    read (10,*) diag%filename
    read (10,*) diag%itype, diag%iwc, diag%freq, diag%temp
    write(*,*) 'Diagnostic ouput to ', trim(diag%filename)
  endif

  if ( line == 'HEADER' ) exit
enddo

read (10,*)
read (10,'(a500)') output_dir
if( trim(output_dir) == '/enter/your/path/here/' ) stop 'Please set a valid directory name in channels.dat'
read (10,*) n_freq
read (10,*) n_sensor
read (10,*) m_type
allocate(id_type(m_type))
allocate(i_scat(m_type))
allocate(particle_shape(m_type))
allocate(dens(m_type))
allocate(psd(m_type))
allocate(regime(m_type))
allocate(wperm(m_type))
allocate(iperm(m_type))
allocate(ll_extend_by_mie(m_type))
allocate(ll_new_integration(m_type))
allocate(max_renorm(m_type))
allocate(mod_gamma_data%i_row(m_type))
allocate(pol_alpha(m_type))
read (10,*) id_type 
read (10,*)

nconfig = 0
do ! Find frequencies and read config lines
  read (10,'(a500)') line
  if (trim(line) == 'FREQUENCIES') then
    exit
  else
    nconfig = nconfig+1
    if (nconfig <= nmaxconfig) then
      cconfig(nconfig) = line
    else
      write(*,*) 'Too many config lines'
      stop
    endif
  endif
enddo

read (cconfig(2),*) ll_melt 
read (cconfig(3),*) ll_dsd
read (cconfig(6),*) i_scat
read (cconfig(8),*) particle_shape
read (cconfig(12),*) dens
read (cconfig(16),*) psd
read (cconfig(18),*) mod_gamma_data%i_row
read (cconfig(21),*) n_rows
if(any(mod_gamma_data%i_row > n_rows)) then
  write(*,*) 'Requested row from modified gamma table does not exist'
  stop
endif
allocate(mod_gamma_data%params(4,n_rows))
allocate(mod_gamma_data%i_free(n_rows))
do i_row=1,n_rows
  read (cconfig(21+i_row),*) mod_gamma_data%params(:,i_row), mod_gamma_data%i_free(i_row)
enddo
read (cconfig(24+n_rows),*) regime
read (cconfig(29+n_rows),*) wperm
read (cconfig(33+n_rows),*) iperm 
read (cconfig(35+n_rows),*) ll_extend_by_mie
read (cconfig(37+n_rows),*) ll_new_integration
read (cconfig(40+n_rows),*) max_renorm
read (cconfig(42+n_rows),*) pol_alpha

print '(a)', 'Config:'
print '(a)', ' Output dir = ' // trim(output_dir)
print '(a,i0)', ' Number of frequencies = ', n_freq
print '(a,i0)', ' Number of sensors = ', n_sensor
print '(a,i0)', ' Number of hydrometeor types = ', m_type
print '(a,l1)', ' Melting layer = ', ll_melt
print '(a,l1)', ' Parameterized n0 = ', ll_dsd
print '(a,12(i0,1x))', ' Hydrometeor IDs = ', id_type
print '(a,12(i0,1x))', ' Scattering computation = ', i_scat
print '(a,12(i0,1x))', ' Particle shapes = ', particle_shape
print '(a,12(i0,1x))', ' Density scheme = ', dens
print '(a,12(i0,1x))', ' PSD = ', psd
print '(a,12(i0,1x))', ' Modified gamma row = ', mod_gamma_data%i_row
print '(a)', 'Row  N0           Mu           Lambda       Gamma       Free parameter'
do i_row=1,n_rows
  print '(i0,3x,4(e12.5,1x),i0)', i_row, mod_gamma_data%params(:,i_row), mod_gamma_data%i_free(i_row)
enddo
print '(a,12(a1,1x))', ' PSD Regime = ', regime
print '(a,12(i0,1x))', ' Liquid water permittivity model = ', wperm
print '(a,12(i0,1x))', ' Ice water permittivity model = ', iperm
print '(a,12(l1,1x))', ' Extend small size regime using solid Mie Spheres = ', ll_extend_by_mie
print '(a,12(l1,1x))', ' Use improved integration basis = ', ll_new_integration
print '(a,12(f9.3))',  ' Max allowed PSD renormalisation (order of magnitude) = ', max_renorm
print '(a,12(f5.2,1x))', ' v/h extinction scaling alpha = ', pol_alpha

! Start of frequencies

read (10,*)
read (10,*)
read (10,*)
read (10,*)
read (10,*)

!---------------------
! Calculate hydrotables
!---------------------

print '(a)', 'Calculating hydrotables...'

! Initialize hydrotable arrays

allocate(ext_tab(n_freq,m_type,n_temp,n_lwc)) 
allocate(ssa_tab(n_freq,m_type,n_temp,n_lwc)) 
allocate(asm_tab(n_freq,m_type,n_temp,n_lwc))
allocate(zef_tab(n_freq,m_type,n_temp,n_lwc))
ssa_tab = 0.0
ext_tab = 0.0
asm_tab = 0.0
zef_tab = 0.0

! Initialize arrays over water content
#ifndef RTTOV_NAG53
!$OMP PARALLEL DEFAULT(SHARED)
#endif
allocate(ext_wc(n_lwc)) 
allocate(ssa_wc(n_lwc)) 
allocate(asm_wc(n_lwc))
allocate(zef_wc(n_lwc))
#ifndef RTTOV_NAG53
!$OMP END PARALLEL
#endif

is_loaded = 0_jpim

! Force Liu DDA parameters to be read into memory (not threadsafe, so should not be done later)
if (any(i_scat == i_liu)) then
  call liu_dda(100.0_jprb, 250.0_jprb, 0_jpim, 1e-3_jprb, tq_ext, tq_sct, tq_asm, tq_bsct, is_loaded)
endif

! Identify the rain hydrometeor (needed for radar reflectivity factor later)
do i_foundrain=1,m_type
  if(id_type(i_foundrain) == i_rain) exit
enddo
if( i_foundrain <= m_type ) then
  ipermwat = iperm(i_foundrain)
else
  ! No rain hydrometoer type included
  ipermwat = permwat_rosenkranz_15
endif

! Load SSP from ARTS SSP database habit files
! For the moment we read all data of all habits used at once and store them in a
! "big" array.
! Up to ECMWF to make that smarter in future.
allocate(f_arts(0:nf_max_arts,m_type))
allocate(T_arts(0:nT_max_arts,m_type))
allocate(D_arts(0:nD_max_arts,m_type))
if (any(i_scat == i_arts)) then
  allocate(ssp_arts(nf_max_arts,nT_max_arts,nD_max_arts,4,m_type))
  ssp_arts = -1.
  f_arts = -1.
  T_arts = -1.
  D_arts = -1.
  call load_arts_ssp( m_type, i_scat, particle_shape, &
                      ssp_arts, f_arts, T_arts, D_arts )
else
  allocate(ssp_arts(nf_max_arts,nT_max_arts,nD_max_arts,0,m_type))
endif

! Loop (1) over frequencies

do i_freq = 1, n_freq

  ! read next line of frequency data from configuration file
  read (10,*) j_freq, fc_ghz, n_band, df_ghz, mpol

  write(*,'(a,i0,a,i0,a,f8.3,a,i2)') &
    'Freq ', i_freq, ' of ', n_freq, ': ', fc_ghz, ' GHz, polarisation option: ', mpol

  wavelength = 0.29997_jprb / fc_ghz  

  !* radar reflectivity norm. factor
  perm_wat = permittivity (i_rain, fc_ghz, 273.0_JPRB, ipermwat=ipermwat)
  k        = (perm_wat - 1.0_jprb) / (perm_wat + 2.0_jprb)
  znorm    = 1.0e+18 * wavelength ** 4.0_jprb / (pi ** 5.0_jprb * abs (k) ** 2.0_jprb)

  ! Loop (2) over bands

  do i_band = 1, n_band

    ! If the instrument channel has side bands, calculate freqs of these based
    ! on info read in from channels.dat then calculate the scattering coefs for
    ! each band add them together and divide by number of bands to get coefs
    ! for channel.

    if (n_band > 1) then
      f_ghz = fc_ghz + (2 * i_band - 3) * df_ghz
      write(*,'(2(a,i0),a,f8.3)') ' Band ', i_band, ' of ', n_band, ': ', f_ghz
    else
      f_ghz = fc_ghz
      write(*,'(a)') ' Band 1 of 1'
    end if

    ! Loop (3) over hydrometeor types

    do i_mtype = 1, m_type

      i_type = id_type(i_mtype)
      mod_gamma_data%i_mtype = i_mtype

      write(*,'(3(a,i0))') '  Hydrometeor ', i_mtype, ' of ', m_type, ', type: ', i_type

      ! Loop (4) over temperatures. Allow this loop to use multiple threads to
      ! speed computation.

#ifndef RTTOV_NAG53
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_temp, temp)
#endif
      do i_temp = 1, n_temp

        temp = temp_offset(i_type) + i_temp

        ! Diagnostics are not otherwise threadsafe, so just activate for the required temperature
        if(associated(diag)) then
          if (abs(temp - diag%temp) < 1e-3_jprb) then
            use_diag=>diag
          else
            use_diag=>NULL()
          endif
        endif

        ! Compute scattering parameters for this temperature, hydrometeor type and frequency
        call scattering_one_temp (temp, wavelength, f_ghz, i_type, i_scat(i_mtype), &
          & dens(i_mtype), psd(i_mtype), regime(i_mtype), wperm(m_type), iperm(m_type), &
          & ll_melt, ll_dsd, particle_shape(i_mtype), is_loaded, &
          & ll_extend_by_mie(i_mtype), ll_new_integration(i_mtype), max_renorm(i_mtype), &
          & ssp_arts(:,:,:,:,i_mtype), f_arts(:,i_mtype), T_arts(:,i_mtype), D_arts(:,i_mtype), &
          & mod_gamma_data, &
          & ext_wc, ssa_wc, asm_wc, zef_wc, &
          & diag=use_diag)

        ! Simple way to get polarisation
        if (mpol /= -1 .and. pol_alpha(i_mtype) /= 0.0_jprb) then
          if (mpol == 4) then ! h
            ext_wc = ext_wc * (1.0_jprb + pol_alpha(i_mtype))
          elseif (mpol == 3) then ! v
            ext_wc = ext_wc * (1.0_jprb - pol_alpha(i_mtype))
          endif
        endif

        ssa_tab(i_freq,i_mtype,i_temp,:) = ssa_tab(i_freq,i_mtype,i_temp,:) + ssa_wc
        ext_tab(i_freq,i_mtype,i_temp,:) = ext_tab(i_freq,i_mtype,i_temp,:) + ext_wc
        asm_tab(i_freq,i_mtype,i_temp,:) = asm_tab(i_freq,i_mtype,i_temp,:) + asm_wc
        zef_tab(i_freq,i_mtype,i_temp,:) = zef_tab(i_freq,i_mtype,i_temp,:) + zef_wc

      end do !(4)
#ifndef RTTOV_NAG53
!$OMP END PARALLEL DO
#endif

    end do !(3)

  end do !(2)

  ! Divide arrays by number of sidebands and convert to output (non-SI) units

  ssa_tab(i_freq,:,:,:) =  ssa_tab(i_freq,:,:,:) / n_band
  ext_tab(i_freq,:,:,:) =  ext_tab(i_freq,:,:,:) / n_band * 1e3_jprb ! to km-1
  asm_tab(i_freq,:,:,:) =  asm_tab(i_freq,:,:,:) / n_band
  where( zef_tab(i_freq,:,:,:) > 0.0 )
    zef_tab(i_freq,:,:,:) =  znorm * zef_tab(i_freq,:,:,:) / n_band
  end where

end do !(1)

close(10)

if(associated(diag)) write (*,*) 'Written: ', trim(diag%filename)

!---------------------------
! Write hydrotable data files
!---------------------------

call convert_hydrotable(input_file, cconfig, nconfig, ext_tab, ssa_tab, asm_tab, zef_tab)

deallocate(asm_tab, ssa_tab, ext_tab, zef_tab, id_type, i_scat, particle_shape, psd, dens, regime)
deallocate(iperm,wperm,max_renorm,ll_extend_by_mie,ll_new_integration)
deallocate(mod_gamma_data%i_row, mod_gamma_data%params, mod_gamma_data%i_free)
if(associated(diag)) deallocate(diag)

#ifndef RTTOV_NAG53
!$OMP PARALLEL DEFAULT(SHARED)
#endif
deallocate(ext_wc, ssa_wc, asm_wc, zef_wc)
#ifndef RTTOV_NAG53
!$OMP END PARALLEL
#endif

end program rttov_scatt_make_coef
