program example_rttovscatt
! Copyright:
!
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
! Example call of RTTOV_SCATT forward operator
!
! This should be run in the same directory as:
!
!   (a) the input profiles, e.g. example_rttovscatt.asc
!   (b) a soft link to the relevant rttov coefficient file, e.g.
!            rtcoef_noaa_15_amsua.dat
!   (c) a soft link to the relevant scattering coefficients, e.g.
!            hydrotable_noaa_amsua.dat
!
! If the rttovscatt_test.sh script has been run succesfully,
! these files will exist in rttov_test/test_rttovscatt.1, so
! you can run the executable from there, e.g.:
!   
!   cd <rttov root>/rttov_test/test_rttovscatt.1
!   ../../bin/example_rttovscatt.exe
!
!
! 28/04/2010   Alan Geer    First version 
! 29/09/2010   Alan Geer    Added a Jacobian calculation
!
!-------------------------------------------------------
   
Use rttov_types, only :   &
  & rttov_coefs          ,&
  & rttov_scatt_coef     ,&
  & rttov_chanprof       ,&
  & rttov_options        ,&
  & rttov_options_scatt  ,&
  & rttov_profile        ,&
  & rttov_profile_cloud  ,&
  & rttov_radiance       ,&
  & rttov_emissivity


Use parkind1, only: jpim, jprb, jplm

    IMPLICIT NONE
  
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_scatt.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_init_rad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"      

logical (kind=jplm)       , allocatable :: calcemiss   (:)
type (rttov_emissivity)   , allocatable :: emissivity  (:), emissivity_k (:)
integer (kind=jpim)       , allocatable :: frequencies (:)
type (rttov_chanprof)     , allocatable :: chanprof    (:) ! Channel and profile indices
type (rttov_profile)      , allocatable :: profiles    (:), profiles_k    (:)
type (rttov_profile_cloud), allocatable :: cld_profiles(:), cld_profiles_k(:)

integer (kind=jpim)        :: errorstatus
type (rttov_radiance)      :: radiance, radiance_k
type (rttov_options)       :: opts     ! Defaults to everything optional switched off
type (rttov_options_scatt) :: opts_scatt
type (rttov_coefs     )    :: coef_rttov
type (rttov_scatt_coef)    :: coef_scatt

integer (kind=jpim) :: instrument (3)
integer (kind=jpim) :: ilev, iprof, ichan, nprof, nchan, nlev, nchannels, err
real    (kind=jprb) :: zenangle
integer (kind=jpim), parameter :: fin = 10
integer (kind=jpim), parameter :: nhydro_frac = 1
character (len=256) :: outstring

! -------------------------------------------------------------------------------------
! Forward calculation 
! -------------------------------------------------------------------------------------

! Open input file
open (fin, file = 'example_rttovscatt.asc', action='READ')
 
! Read instrument and platform specification
read (fin,*) instrument (1)   ! Satellite series (see rttov_const.F90)
read (fin,*) instrument (2)   ! Satellite ID     (see rttov_const.F90)
read (fin,*) instrument (3)   ! Sensor ID        (see rttov_const.F90)
read (fin,*) zenangle         ! Zenith angle [degrees]
 
! Read number of profiles and levels
read (fin,*) nprof
read (fin,*) nlev

! Read / initialise RTTOV coefficients for this instrument
write(*,*) 'Reading coefficients' 
call rttov_read_coefs      (err, coef_rttov, opts, instrument=instrument)
call rttov_read_scattcoeffs (err, opts_scatt, coef_rttov, coef_scatt)
  
nchan     = coef_rttov%coef%fmv_chn   ! number of channels on instrument
nchannels = nprof * nchan             ! total channels to simulate

allocate (chanprof     (nchannels))
allocate (frequencies  (nchannels))
allocate (emissivity   (nchannels))
allocate (calcemiss    (nchannels))
allocate (profiles     (nprof))
allocate (cld_profiles (nprof))

! Request RTTOV / FASTEM to calculate surface emissivity
calcemiss  = .true.
emissivity % emis_in = 0.0_JPRB

! Setup indices
call rttov_scatt_setupindex ( &
 & errorstatus,     & ! out
 & nprof,           & ! in  
 & nchan,           & ! in 
 & coef_rttov,      & ! in
 & coef_scatt,      & ! in
 & nchannels,       & ! in
 & chanprof,        & ! out
 & frequencies)       ! out

! Allocate profiles (input) and radiance (output) structures
call rttov_alloc_prof      ( err, nprof,     profiles, nlev, opts, 1_jpim, init = .true._jplm)
call rttov_alloc_scatt_prof( err, nprof, cld_profiles, nlev, coef_scatt%nhydro, nhydro_frac, 1_jpim, init = .true._jplm)
call rttov_alloc_rad       ( err, nchannels, radiance, nlev, 1_jpim)

write(*,*) 'Reading profiles' 

profiles(:) % gas_units = 2 ! ppmv over moist air

! Read input profiles
do iprof = 1, nprof

  read (fin,*) ! Skip line

  ! Surface
  read (fin,*) profiles (iprof) % skin % t        ! surface skin temperature (K)
  read (fin,*) profiles (iprof) % s2m % t         ! 2-meter temperature (K)
  read (fin,*) profiles (iprof) % s2m % q         ! 2-meter specific humidity (ppmv)
  read (fin,*) profiles (iprof) % s2m % u         ! 10-meter wind u (m/s)
  read (fin,*) profiles (iprof) % s2m % v         ! 10-meter wind v (m/s) 
  read (fin,*) cld_profiles (iprof) % ph (nlev+1) ! lowest half-level pressure (hPa)

  ! Surface and other initialisations
  profiles (iprof) % skin % surftype = 1  ! Ocean
  profiles (iprof) % skin % watertype = 1 ! Ocean water
  profiles (iprof) % ctp = 500.0_JPRB     ! Not used but still required by RTTOV

  ! Surface pressure is lowest half-level pressure
  profiles (iprof) % s2m % p  = cld_profiles (iprof) % ph (nlev+1) 

  ! Zenith angle
  profiles (iprof) % zenangle = zenangle

  read (fin,*) ! Skip line
  read (fin,*) ! Skip line
  read (fin,*) ! Skip line

  ! Levels
  do ilev = 1,nlev

    read (fin,'(9e10.3)') &
      & profiles     (iprof) % p    (ilev), &   ! full level pressure (hPa)
      & cld_profiles (iprof) % ph   (ilev), &   ! half level pressure (hPa)
      & profiles     (iprof) % t    (ilev), &   ! temperature (K)
      & profiles     (iprof) % q    (ilev), &   ! specific humidity (ppmv)
      & cld_profiles (iprof) % hydro_frac (ilev,1), & ! cloud cover    (0-1)
      & cld_profiles (iprof) % hydro (ilev,4), &   ! liquid water   (kg/kg)
      & cld_profiles (iprof) % hydro (ilev,5), &   ! ice water      (kg/kg)
      & cld_profiles (iprof) % hydro (ilev,1), &   ! rain           (kg/kg)
      & cld_profiles (iprof) % hydro (ilev,2)      ! frozen precip. (kg/kg)
  
    cld_profiles (iprof) % hydro (ilev,3) = 0.0_JPRB  ! Graupel (kg/kg) not yet provided in data file

  enddo 
enddo  
close (fin)  

write(*,*) 'Calling forward model' 

! Forward model run
call rttov_scatt ( &
  & errorstatus,         &! out
  & opts_scatt,          &! in
  & nlev,                &! in
  & chanprof,            &! in
  & frequencies,         &! in
  & profiles,            &! in  
  & cld_profiles,        &! in
  & coef_rttov,          &! in
  & coef_scatt,          &! in
  & calcemiss,           &! in
  & emissivity,          &! in
  & radiance)             ! out

! Write output
write(*,*) ' '
write(*,*) 'Channel Profile          Tb '

do ichan = 1, nchannels
  write (*,'(2i8,3x,f10.3)') chanprof(ichan), radiance % bt (ichan)
enddo

! -------------------------------------------------------------------------------------
! Jacobian calculation - achieved by calling the adjoint code in "K" mode
! -------------------------------------------------------------------------------------

write(*,*) ' '
write(*,*) 'Calling Jacobian model' 

! Allocate & initialiase K profiles (input/output) and K radiance (input) structures
allocate (profiles_k     (nchannels))
allocate (cld_profiles_k (nchannels))
call rttov_alloc_prof      ( err, nchannels,     profiles_k, nlev, opts, 1_jpim, init = .true._jplm)
call rttov_alloc_scatt_prof( err, nchannels, cld_profiles_k, nlev, coef_scatt%nhydro, nhydro_frac, 1_jpim, init = .true._jplm)
call rttov_alloc_rad       ( err, nchannels,     radiance_k, nlev, 1_jpim)
call rttov_init_rad        ( radiance_k)

radiance_k % bt(:)       = 1.0_JPRB ! (ACTIVATES K behaviour in adjoint code)
  
allocate (emissivity_k(nchannels))
emissivity_k(1:nchannels) % emis_in  = 0.0_JPRB
emissivity_k(1:nchannels) % emis_out = 0.0_JPRB

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & nlev,               &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles,           &! inout  
  & cld_profiles,       &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity,         &! inout
  & profiles_k,         &! in (ACTIVATES K behaviour in adjoint code: note dimension)
  & cld_profiles_k,     &! in
  & emissivity_k,       &! inout
  & radiance,           &! inout
  & radiance_k)          ! inout 

! Write output
write(*,*) ' '
write(*,*) 'Jacobian wrt T [1E-3 K / K] for profile 1'

write(*,*) ' '
write(*,*) 'Level                               Channel'
write (*,'(6x,30i4)') chanprof (1:nchan) % chan 
write(*,*) ' '
do ilev = 1, nlev
  write (outstring,'(i3,3x)') ilev
  do ichan = 1, nchan
    write (outstring, '(a,i4)') trim(outstring), int(1000.*profiles_k(ichan) % t (ilev))
  end do
  write (*,*) trim(outstring)
enddo

! Deallocate all storage
call rttov_dealloc_coefs      ( err,          coef_rttov )
call rttov_dealloc_scattcoeffs(               coef_scatt )
call rttov_alloc_prof         ( err, nprof,     profiles, nlev, opts,    0_jpim)
call rttov_alloc_scatt_prof   ( err, nprof, cld_profiles, nlev, coef_scatt%nhydro, nhydro_frac, 0_jpim)
call rttov_alloc_rad          ( err, nchannels, radiance, nlev,          0_jpim)
call rttov_alloc_prof         ( err, nchannels, profiles_k, nlev, opts,    0_jpim)
call rttov_alloc_scatt_prof   ( err, nchannels, cld_profiles_k, nlev, coef_scatt%nhydro, nhydro_frac, 0_jpim)
call rttov_alloc_rad          ( err, nchannels, radiance_k, nlev,          0_jpim)
deallocate (chanprof, frequencies, emissivity, emissivity_k, calcemiss)
deallocate (profiles, profiles_k, cld_profiles, cld_profiles_k)

end program example_rttovscatt
