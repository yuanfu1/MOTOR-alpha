subroutine rttovscatt_test_one ( nchannels, opts_scatt, coef_rttov, coef_scatt, &
                        & chanprof, &
                        & frequencies,  &
                        & emissivity, lradar)
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

!
! Main part of the RTTOV-SCATT test routines
!
!       2005   Peter Bauer        First version
! 13/12/2007   Alan Geer          RTTOV9 version; streamlined
! 02/01/2019   Alan Geer          RTTOV13 version with 5 hydrometeors and radar
!
!-------------------------------------------------------

Use rttov_types, only :   &
  & rttov_coefs          ,&
  & rttov_scatt_coef     ,&
  & rttov_chanprof       ,&
  & rttov_options_scatt

Use parkind1, only: jpim, jprb, jplm
!INTF_OFF
Use mod_rttovscatt_test, only: kflevg, kproma, fastem_land_coeff, ioin, zenangle, iooutstreams

Use rttov_types, only :   &
  & rttov_options        ,&
  & rttov_profile        ,&
  & rttov_profile_cloud  ,&
  & rttov_radiance       ,&
  & rttov_reflectivity   ,&
  & rttov_emissivity

Use rttov_const, only :   &
  & errorstatus_fatal,    &
  & errorstatus_success,  &
  & default_err_unit,     &
  & fastem_sp

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!INTF_ON
    IMPLICIT NONE

integer (kind=jpim),  intent (in) :: nchannels
type(rttov_options_scatt), intent(in) :: opts_scatt
real    (kind=jprb),  intent (in) , dimension (nchannels) :: emissivity    
integer (kind=jpim),  intent (in) , dimension (nchannels) :: frequencies
Type(rttov_chanprof), Intent (in) , dimension (nchannels) :: chanprof 

type (rttov_coefs     ), intent (inout) :: coef_rttov        
type (rttov_scatt_coef), intent (inout) :: coef_scatt  

logical(kind=jplm), intent(in) :: lradar

!INTF_END

#include "rttov_scatt.interface"
#include "rttov_scatt_tl.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_scatt_prof.interface"

!* FORWARD   
type (rttov_profile)        :: profiles_d1     (kproma)
type (rttov_profile_cloud)  :: cld_profiles_d1 (kproma)
type (rttov_radiance)       :: radiance_d1
type (rttov_reflectivity), pointer :: reflectivity_d1 => null()
type (rttov_emissivity), dimension (nchannels) :: emissivity_d1

!* TL   
type (rttov_profile)       :: profiles_d2     (kproma)
type (rttov_profile)       :: profiles_tl     (kproma)
type (rttov_profile)       :: profiles_tl2     (kproma)
type (rttov_profile_cloud) :: cld_profiles_d2 (kproma)
type (rttov_profile_cloud) :: cld_profiles_tl (kproma)
type (rttov_profile_cloud) :: cld_profiles_tl2 (kproma)
type (rttov_radiance)      :: radiance_d2
type (rttov_radiance)      :: radiance_d3
type (rttov_radiance)      :: radiance_tl
type (rttov_radiance)      :: radiance_tl2
type (rttov_reflectivity), pointer :: reflectivity_d2 => null()
type (rttov_reflectivity), pointer :: reflectivity_d3 => null()
type (rttov_reflectivity), pointer :: reflectivity_tl => null()
type (rttov_reflectivity), pointer :: reflectivity_tl2 => null()


type (rttov_emissivity), dimension (nchannels) :: emissivity_d2
type (rttov_emissivity), dimension (nchannels) :: emissivity_tl
type (rttov_emissivity), dimension (nchannels) :: emissivity_tl2
           
!* AD   
type (rttov_profile)       :: profiles_ad     (kproma)
type (rttov_profile_cloud) :: cld_profiles_ad (kproma)
type (rttov_radiance)      :: radiance_ad
type (rttov_reflectivity), pointer :: reflectivity_ad => null()

type (rttov_profile)       :: profiles_ad2     (kproma)
type (rttov_profile_cloud) :: cld_profiles_ad2 (kproma)
type (rttov_radiance)      :: radiance_ad2
type (rttov_reflectivity), pointer :: reflectivity_ad2 => null()

type (rttov_emissivity), dimension (nchannels) :: emissivity_ad
    
type (rttov_emissivity), dimension (nchannels) :: emissivity_ad2
       
!* K   
type (rttov_profile)       :: profiles_k     (nchannels)
type (rttov_profile_cloud) :: cld_profiles_k (nchannels)
type (rttov_radiance)      :: radiance_k
type (rttov_reflectivity), pointer :: reflectivity_k => null()

type (rttov_emissivity), dimension (nchannels) :: emissivity_k

!* OTHER                      
logical (kind=jplm) :: calcemiss (nchannels)

integer (kind=jpim) :: errorstatus, erroralloc
integer (kind=jpim) :: i_lev, i_proma, i_chan, i_btout, i_lambda, i_fast, i_output, ioout, i_type

real    (kind=jprb) :: lambda, zeps, zdelta1, zdelta2, threshold, z
Real    (Kind=jprb)       :: ratio(2)
Real(Kind=jprb),    Allocatable :: radiance_total_ref (:)
Real(Kind=jprb),    Allocatable :: reflectivity_total_ref (:,:)
REAL(KIND=jprb) :: ZHOOK_HANDLE
real    (kind=jprb), parameter :: imposed_cfrac    = 0.74_jprb  
real    (kind=jprb), dimension (kproma) :: cfrac

integer (kind=jpim), parameter :: nhydro = 5

integer (kind=jpim), parameter :: iallocate = 1, ideallocate = 0

character(len=50), dimension(2) :: format_out = (/'(i4,3x,30f24.16)','(i4,3x,30f14.6) '/)
character(len=50), dimension(2) :: format_radar = (/'(i4,3x,i4,3x,2e24.16)','(i4,3x,i4,3x,2e18.10)'/)

!- End of header ------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TEST',0_jpim,ZHOOK_HANDLE)

if (lradar) then
  ! AJGDB for whatever reason, possibly the surfeit of logarithms, the radar simulator doesn't achieve quite the same TL/AD accuracy
  threshold = 1.0E-6_jprb
else
  threshold = 1.0E-8_jprb
endif


!* FORWARD-MODEL TEST ***********************************************************************************
!* Set-up
errorstatus = errorstatus_success
emissivity_d1 (1:nchannels) % emis_in = emissivity    (1:nchannels)
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) % emis_in < 0.01_jprb

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_d1, cld_profiles_d1, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d1, kflevg, iallocate)
Allocate ( radiance_total_ref  ( nchannels ) )
if (lradar) then
  allocate ( reflectivity_total_ref  ( kflevg, nchannels ) )
  allocate(reflectivity_d1, reflectivity_d2, reflectivity_d3, reflectivity_tl, reflectivity_tl2)
  allocate(reflectivity_ad, reflectivity_ad2, reflectivity_k)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d1, kflevg, iallocate)
endif


!* Read profiles
open (ioin,  file = 'example_rttovscatt.asc', status = 'old')

! Skip things that are provided from elsewhere
do i_proma = 1, 6
  read (ioin,*) ! Skip line
enddo

do i_proma = 1, kproma

  read (ioin,*) ! Skip line

  ! Surface
  read (ioin,*) profiles_d1 (i_proma) % skin % t ! surface skin temperature (K)
  read (ioin,*) profiles_d1 (i_proma) % s2m % t  ! 2-meter temperature (K)
  read (ioin,*) profiles_d1 (i_proma) % s2m % q  ! 2-meter specific humidity (ppmv)
  read (ioin,*) profiles_d1 (i_proma) % s2m % u  ! 10-meter wind u (m/s)
  read (ioin,*) profiles_d1 (i_proma) % s2m % v  ! 10-meter wind v (m/s) 
  read (ioin,*) profiles_d1 (i_proma) % s2m % p  ! lowest half-level pressure (hPa)
     
  profiles_d1 (i_proma) % s2m % o = 0._jprb
  profiles_d1 (i_proma) % s2m % wfetc = 0._jprb ! Unused (solar computations only)     
  profiles_d1 (i_proma) % skin % surftype  = 1  ! Ocean  
  profiles_d1 (i_proma) % skin % watertype = 0  ! Ocean water   
  profiles_d1 (i_proma) % skin % fastem(:) = fastem_land_coeff (:)

  profiles_d1 (i_proma) % zenangle   = zenangle
  profiles_d1 (i_proma) % azangle    = 0._jprb   ! default value ! FASTEM-3 will require an actual value here
  profiles_d1 (i_proma) % ctp        = 500.0_jprb ! default value
  profiles_d1 (i_proma) % cfraction  = 0._jprb   ! default value
  profiles_d1 (i_proma) % elevation  = 0._jprb   ! default value
  profiles_d1 (i_proma) % sunzenangle  = 0._jprb   ! default value
  profiles_d1 (i_proma) % sunazangle   = 0._jprb   ! default value
  profiles_d1 (i_proma) % latitude     = 0._jprb   ! default value
  profiles_d1 (i_proma) % longitude    = 0._jprb   ! default value
  profiles_d1 (i_proma) % Be           = 0._jprb   ! default value
  profiles_d1 (i_proma) % cosbk        = 0._jprb   ! default value
 
  if(opts_scatt%lusercfrac) cld_profiles_d1 (i_proma) % cfrac = imposed_cfrac

  ! Surface pressure is also lowest half-level pressure
  cld_profiles_d1 (i_proma) % ph (kflevg+1) = profiles_d1 (i_proma) % s2m % p

  read (ioin,*) ! Skip line
  read (ioin,*) ! Skip line
  read (ioin,*) ! Skip line

  ! Levels
  do i_lev = 1, kflevg

    read (ioin,'(9e10.3)') &
      & profiles_d1     (i_proma) % p    (i_lev), &   ! full level pressure (hPa)
      & cld_profiles_d1 (i_proma) % ph   (i_lev), &   ! half level pressure (hPa)
      & profiles_d1     (i_proma) % t    (i_lev), &   ! temperature (K)
      & profiles_d1     (i_proma) % q    (i_lev), &   ! specific humidity (ppmv)
      & cld_profiles_d1 (i_proma) % hydro_frac (i_lev,1), &   ! cloud cover    (0-1)
      & cld_profiles_d1 (i_proma) % hydro (i_lev,4), &  ! liquid water   (kg/kg)
      & cld_profiles_d1 (i_proma) % hydro (i_lev,5), &  ! ice water      (kg/kg)
      & cld_profiles_d1 (i_proma) % hydro (i_lev,1), &  ! rain           (kg/kg)
      & cld_profiles_d1 (i_proma) % hydro (i_lev,2)     ! frozen precip. (kg/kg)

    ! Graupel (kg/kg) not yet provided in data file
    cld_profiles_d1 (i_proma) % hydro (i_lev,3) = 0.0_JPRB  

    ! Simulate use of a hydrometeor fraction per hydrometeor
    do i_type = 2,nhydro
      cld_profiles_d1 (i_proma) % hydro_frac (:,i_type) = cld_profiles_d1 (i_proma) % hydro_frac (:,1)
    enddo

  enddo 
enddo  

close (ioin)  

!* Reference forward model run
call rttov_scatt ( &
  & errorstatus,         &! out
  & opts_scatt,          &! in
  & kflevg,              &! in
  & chanprof,            &! in
  & frequencies,         &! in
  & profiles_d1,         &! inout  
  & cld_profiles_d1,     &! in
  & coef_rttov,          &! in
  & coef_scatt,          &! in
  & calcemiss,           &! in
  & emissivity_d1,       &! inout
  & radiance_d1,         &! inout
  & cfrac,               &! out, diagnostic only
  & reflectivity=reflectivity_d1)! inout

! main output:
! radiance_d1%tb        = cloud-affected Tbs
! radiance_d1%tb_clear  = clear-sky Tbs
  
! Full and reduced output
do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)
  
  write(ioout,'(A10,i4)') 'nchan ', nchannels

  write(ioout,*) '--------------------------------------------------------------------------'
  write(ioout,*) '--------------------------------------------------------------------------'
  write(ioout,*) 'This dataset is made of ',kproma,' ECMWF model profiles'
  write(ioout,*)
  write(ioout,*)
  write(ioout,*) 'Call to RTTOV_SCATT'
  write(ioout,*) '-------------------'
  write(ioout,*)
  write(ioout,*) 'Channel  cloudy Tb '

  do i_chan = 1, nchannels
    write (ioout,format_out(i_output)) i_chan, radiance_d1 % bt (i_chan)
  enddo

  write(ioout,*)
  write(ioout,*) 'Channel  clear Tb '

  do i_chan = 1, nchannels
    write (ioout,format_out(i_output)) i_chan, radiance_d1 % bt_clear (i_chan)
  enddo

  if (lradar) then
    write(ioout,*)
    write(ioout,*) 'Channel Level reflectivity (corrected)   reflectivity (attenuated) '

    do i_chan = 1, nchannels
      do i_lev = 1, kflevg
        write (ioout,format_radar(i_output)) i_chan, i_lev, reflectivity_d1 % zef (i_lev,i_chan), &
                                                            reflectivity_d1 % azef (i_lev,i_chan)
      enddo
    enddo
  endif

enddo

!* TANGENT-LINEAR TEST ***********************************************************************************

zeps = 0.01_jprb

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_tl, cld_profiles_tl, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d3, kflevg, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl, kflevg, iallocate)
if (lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d3, kflevg, iallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_tl, kflevg, iallocate)
endif

call set_perturbation ( kproma, profiles_tl, cld_profiles_tl, profiles_d1, cld_profiles_d1, zeps)

emissivity_tl (1:nchannels) % emis_in = emissivity_d1 (1:nchannels) % emis_in * zeps
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) % emis_in < 0.01_jprb

call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! in  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl,        &! in
  & cld_profiles_tl,    &! in
  & emissivity_tl,      &! inout
  & radiance_d3,        &! inout
  & radiance_tl,        &! inout
  & reflectivity = reflectivity_d3, & ! inout
  & reflectivity_tl =  reflectivity_tl) ! inout

If ( Any( abs(radiance_d3%bt(:) - radiance_d1%bt(:)) > threshold * radiance_d1%bt(:)  ))  Then 
  write(default_err_unit,*) 'wrong forward model in TL'
  write(default_err_unit,*) radiance_d3%bt(:)
  write(default_err_unit,*) abs(radiance_d3%bt(:)-radiance_d1%bt(:)) / (threshold * radiance_d1%bt(:))
  Stop 1
Endif

if(lradar) then
  do i_lev = 1, kflevg
    If ( Any( abs(reflectivity_d3%zef(i_lev,:) - reflectivity_d1%zef(i_lev,:)) > &
              abs(threshold * reflectivity_d1%zef(i_lev,:))))  Then
      write(default_err_unit,*) 'wrong forward model in TL- reflectivity'
      write(default_err_unit,*) reflectivity_d1%zef(i_lev,:)
      write(default_err_unit,*) reflectivity_d3%zef(i_lev,:)
      write(default_err_unit,*) abs(reflectivity_d3%zef(i_lev,:)-reflectivity_d1%zef(i_lev,:))
      Stop 1
    Endif
  enddo
endif

! Save radiance as a reference for the trajectory
! TL is used instead of rttov_scatt because
! calcemis = F and reflectivities have not been saved
radiance_total_ref(:) = radiance_d1%bt(:)
if (lradar) then
  reflectivity_total_ref(:,:) = reflectivity_d1%azef(:,:)
endif

do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)
 
  write(ioout,*)
  write(ioout,*) 'Call to RTTOV_SCATT_TL'
  write(ioout,*) '----------------------'
  write(ioout,*)
  write(ioout,*) 'Channel  cloudy Tb TL '

  do i_chan = 1, nchannels
    write (ioout,format_out(i_output)) i_chan, radiance_tl % bt (i_chan)
  enddo

  write(ioout,*)
  write(ioout,*) 'Channel  clear Tb TL '

  do i_chan = 1, nchannels
    write (ioout,format_out(i_output)) i_chan, radiance_tl % bt_clear (i_chan)
  enddo

  if (lradar) then
    write(ioout,*)
    write(ioout,*) 'Channel Level reflectivity TL (corrected)   reflectivity TL (attenuated) '

    do i_chan = 1, nchannels
      do i_lev = 1, kflevg
        write (ioout,format_radar(i_output)) i_chan, i_lev, reflectivity_tl % zef (i_lev,i_chan), &
                                                            reflectivity_tl % azef (i_lev,i_chan)
      enddo
    enddo
  endif

enddo

do i_output=1,1
  ioout = iooutstreams(i_output)
  write(ioout,*)
  write(ioout,*) 'Test TL'
  write(ioout,*) '-------'
  write(ioout,*)
enddo
!---------------------------
! second run of TL
!---------------------------
lambda = 0.5_jprb
  
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_tl2, cld_profiles_tl2, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl2, kflevg, iallocate)
if (lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_tl2, kflevg, iallocate)
endif

call set_perturbation ( kproma, profiles_tl2, cld_profiles_tl2, profiles_tl, cld_profiles_tl, lambda)
  
emissivity_tl2 (1:nchannels) % emis_in = emissivity_tl (1:nchannels) % emis_in *  lambda
calcemiss     (1:nchannels) = emissivity_tl (1:nchannels) % emis_in < 0.01_jprb

call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl2,       &! in
  & cld_profiles_tl2,   &! in
  & emissivity_tl2,     &! inout
  & radiance_d1,        &! inout
  & radiance_tl2,       &! inout
  & reflectivity = reflectivity_d1,     &! inout
  & reflectivity_tl = reflectivity_tl2)  ! inout

!---------------------------

do i_chan = 1, nchannels
  if( abs(lambda * radiance_tl%bt_clear(i_chan) - radiance_tl2%bt_clear(i_chan)) > threshold ) &
    call error_stop( 'TL test fails for radiance_tl%bt_clear for channel ', i_chan )
  if( abs(lambda * radiance_tl%bt(i_chan) - radiance_tl2%bt(i_chan)) > threshold ) &
    call error_stop( 'TL test fails for radiance_tl%bt for channel ', i_chan )
end do

if (lradar) then
  do i_chan = 1, nchannels
    do i_lev = 1, kflevg
      if( abs(lambda * reflectivity_tl % zef (i_lev,i_chan) - reflectivity_tl2 % zef (i_lev,i_chan)) > threshold ) &
        call error_stop( 'TL test fails for reflectivity_tl % zef for channel, level ', i_chan, inum2=i_lev )
      if( abs(lambda * reflectivity_tl % azef (i_lev,i_chan) - reflectivity_tl2 % azef (i_lev,i_chan)) > threshold ) &
        call error_stop( 'TL test fails for reflectivity_tl % azef for channel, level ', i_chan, inum2=i_lev )
    enddo
  enddo
endif

! Now run the Taylor test
!-------------------------

! AJGDB no radar yet, but allocate reflectivity array for re-use later.

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_d2, cld_profiles_d2, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d2, kflevg, iallocate)
if(lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d2, kflevg, iallocate)
endif

do i_output=1,1
  ioout = iooutstreams(i_output)
  write(ioout,*) '  Unless the LGRADP option is set to true to activate the RTTOV '
  write(ioout,*) ' TL/adjoint pressure sensitivity, the following test will show very non-linear '
  write(ioout,*) ' behaviour from the operator, i.e (h*-h)/H very different from 1.'
enddo

do i_lambda = 10, 1, -1
  lambda = 10.0_jprb ** (-1.0_jprb * i_lambda) 

  errorstatus = errorstatus_success

  emissivity_d2 (1:nchannels) % emis_in = emissivity_d1 (1:nchannels) % emis_in + &
                                        & emissivity_tl (1:nchannels) % emis_in * lambda
  calcemiss     (1:nchannels) = emissivity_d2 (1:nchannels) % emis_in < 0.01_jprb

  do i_proma = 1, kproma    
  
    !* Add perturbations
    profiles_d2 (i_proma) % p (:) = profiles_d1 (i_proma) % p (:) + profiles_tl (i_proma) % p (:) * lambda
    profiles_d2 (i_proma) % t (:) = profiles_d1 (i_proma) % t (:) + profiles_tl (i_proma) % t (:) * lambda 
    profiles_d2 (i_proma) % q (:) = profiles_d1 (i_proma) % q (:) + profiles_tl (i_proma) % q (:) * lambda 

    cld_profiles_d2 (i_proma) % cfrac   = cld_profiles_d1 (i_proma) % cfrac   + cld_profiles_tl (i_proma) % cfrac   * lambda
    cld_profiles_d2 (i_proma) % ph  (:) = cld_profiles_d1 (i_proma) % ph  (:) + cld_profiles_tl (i_proma) % ph  (:) * lambda
    cld_profiles_d2 (i_proma) % hydro (:,:) = cld_profiles_d1 (i_proma) % hydro (:,:) + &
       & cld_profiles_tl (i_proma) % hydro (:,:) * lambda
    cld_profiles_d2 (i_proma) % hydro_frac (:,:) = cld_profiles_d1 (i_proma) % hydro_frac (:,:) + &
       & cld_profiles_tl (i_proma) % hydro_frac (:,:) * lambda

!* Fill in RTTOV/RTTOVSCATT arrays once per profile
    profiles_d2 (i_proma) % s2m % p = profiles_d1 (i_proma) % s2m % p + profiles_tl (i_proma) % s2m % p * lambda
    profiles_d2 (i_proma) % s2m % q = profiles_d1 (i_proma) % s2m % q + profiles_tl (i_proma) % s2m % q * lambda
    profiles_d2 (i_proma) % s2m % o = profiles_d1 (i_proma) % s2m % o + profiles_tl (i_proma) % s2m % o * lambda
    profiles_d2 (i_proma) % s2m % t = profiles_d1 (i_proma) % s2m % t + profiles_tl (i_proma) % s2m % t * lambda
    profiles_d2 (i_proma) % s2m % u = profiles_d1 (i_proma) % s2m % u + profiles_tl (i_proma) % s2m % u * lambda
    profiles_d2 (i_proma) % s2m % v = profiles_d1 (i_proma) % s2m % v + profiles_tl (i_proma) % s2m % v * lambda
    profiles_d2 (i_proma) % s2m % wfetc = 0._jprb
    profiles_d2 (i_proma) % skin % surftype   = profiles_d1 (i_proma) % skin % surftype
    profiles_d2 (i_proma) % skin % t = profiles_d1 (i_proma) % skin % t + profiles_tl (i_proma) % skin % t * lambda 
    profiles_d2 (i_proma) % skin % fastem (:) = profiles_d1 (i_proma) % skin % fastem (:) + &
                                              & profiles_tl (i_proma) % skin % fastem (:) * lambda 
    profiles_d2 (i_proma) % skin % watertype = 0     ! Ocean water
    
    profiles_d2 (i_proma) % zenangle   = zenangle
    profiles_d2 (i_proma) % azangle    = 0._jprb     ! default value
    profiles_d2 (i_proma) % ctp        = 500._jprb   ! default value
    profiles_d2 (i_proma) % cfraction  = 0._jprb     ! default value
    profiles_d2 (i_proma) % elevation  = 0._jprb     ! default value
    profiles_d2 (i_proma) % sunzenangle  = 0._jprb   ! default value
    profiles_d2 (i_proma) % sunazangle   = 0._jprb   ! default value
    profiles_d2 (i_proma) % latitude     = 0._jprb   ! default value
    profiles_d2 (i_proma) % longitude    = 0._jprb   ! default value
    profiles_d2 (i_proma) % Be           = 0._jprb   ! default value
    profiles_d2 (i_proma) % cosbk        = 0._jprb   ! default value

  end do    

  !* Reference forward model run
  call rttov_scatt ( &
    & errorstatus,         &! out
    & opts_scatt,          &! in
    & kflevg,              &! in
    & chanprof,            &! in
    & frequencies,         &! in
    & profiles_d2,         &! inout
    & cld_profiles_d2,     &! in
    & coef_rttov,          &! in
    & coef_scatt,          &! in
    & calcemiss,           &! in
    & emissivity_d2,       &! inout
    & radiance_d2,         &! inout
    & cfrac)                ! out, diagnostic only

  do i_output=1,1
    ioout = iooutstreams(i_output)

    write(ioout,*)
    write(ioout,*) 'Chan      Lambda         Cloudy  [h(x*)-h(x)]/H(x*-x)  Clear '

    do i_chan = 1, nchannels
      if((radiance_tl % bt(i_chan) /= 0) .and. (radiance_tl % bt_clear(i_chan) /= 0)) then
        ratio(1) = (radiance_d2 % bt(i_chan) - radiance_d1 % bt(i_chan)) / (lambda * radiance_tl % bt(i_chan))
        ratio(2) = (radiance_d2 % bt_clear(i_chan) - radiance_d1 % bt_clear(i_chan)) / (lambda * radiance_tl % bt_clear(i_chan))
        write (ioout,'(i4,3x,1e9.2,2f25.16)') i_chan,  lambda, ratio(1), ratio(2)
      endif
    enddo
  enddo
enddo
 
!* ADJOINT TEST ***********************************************************************************

do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)

  write(ioout,*)
  write(ioout,*) 'Test AD'
  write(ioout,*) '-------'
  write(ioout,*)

  write(ioout,*) '1 - Test Linearity'
  write(ioout,*)
enddo

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_ad, cld_profiles_ad, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad, kflevg, iallocate)
if(lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_ad, kflevg, iallocate)
endif

call rttov_init_prof(profiles_ad)
call rttov_init_scatt_prof(cld_profiles_ad)

emissivity_ad (1:nchannels) % emis_in  = 0._jprb
emissivity_ad (1:nchannels) % emis_out = 0._jprb
   
! Set perturbations
!
call rttov_init_rad(radiance_ad)
radiance_ad % bt_clear(:)   = 0.05_jprb * radiance_d1 % bt_clear(:)
radiance_ad % bt(:)         = 0.05_jprb * radiance_d1 % bt(:)
if (lradar) then
  reflectivity_ad % zef(:,:)  = 0.05_jprb * reflectivity_d1 % zef(:,:)
  reflectivity_ad % azef(:,:) = 0.05_jprb * reflectivity_d1 % azef(:,:)
endif

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad,        &! in
  & cld_profiles_ad,    &! in
  & emissivity_ad,      &! inout
  & radiance_d2,        &! inout
  & radiance_ad,        &! inout
  & reflectivity=reflectivity_d2,    &! inout
  & reflectivity_ad=reflectivity_ad)  ! inout

If ( errorstatus == errorstatus_fatal ) Then
  call error_stop2('rttov_scatt_ad error 1')
End If

If ( Any( abs(radiance_total_ref(:) - radiance_d2%bt(:)) > threshold * radiance_total_ref(:)  ))  Then 
  write(default_err_unit,*) 'wrong forward model in AD'
  write(default_err_unit,*) radiance_total_ref(:)
  write(default_err_unit,*) abs(radiance_total_ref(:)-radiance_d2%bt(:)) / (threshold * radiance_total_ref(:))
  Stop 1
Endif

if(lradar) then
  do i_lev = 1, kflevg
    If ( Any( abs(reflectivity_total_ref(i_lev,:) - reflectivity_d2%azef(i_lev,:)) > &
              abs(threshold * reflectivity_d2%azef(i_lev,:))))  Then
      write(default_err_unit,*) 'wrong forward model in AD- reflectivity: i_lev ', i_lev
      write(default_err_unit,*) reflectivity_total_ref(i_lev,:)
      write(default_err_unit,*) reflectivity_d2%azef(i_lev,:)
      write(default_err_unit,*) abs(reflectivity_total_ref(i_lev,:)-reflectivity_d2%azef(i_lev,:))
      Stop 1
    Endif
  enddo
endif

!---------------------------
! Second run of AD

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_ad2, cld_profiles_ad2, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad2, kflevg, iallocate)
if(lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_ad2, kflevg, iallocate)
endif

call rttov_init_prof(profiles_ad2)
call rttov_init_scatt_prof(cld_profiles_ad2)

emissivity_ad2 (1:nchannels) % emis_in  = 0._jprb
emissivity_ad2 (1:nchannels) % emis_out = 0._jprb
     
! Set perturbations
!
call rttov_init_rad(radiance_ad2)
radiance_ad2 % bt_clear(:)   = 0.05_jprb * radiance_d1 % bt_clear(:) * lambda
radiance_ad2 % bt(:)         = 0.05_jprb * radiance_d1 % bt(:) * lambda
if (lradar) then
  reflectivity_ad2 % zef(:,:)  = 0.05_jprb * reflectivity_d1 % zef(:,:) * lambda
  reflectivity_ad2 % azef(:,:) = 0.05_jprb * reflectivity_d1 % azef(:,:) * lambda
endif

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad2,       &! in
  & cld_profiles_ad2,   &! in
  & emissivity_ad2,     &! inout
  & radiance_d2,        &! inout
  & radiance_ad2,       &! inout
  & reflectivity=reflectivity_d2,    &! inout
  & reflectivity_ad=reflectivity_ad2)  ! inout

If ( errorstatus == errorstatus_fatal ) Then
  call error_stop2('rttov_scatt_ad error 2')
End If

do i_proma = 1, kproma
  do i_lev = 1, profiles_ad (i_proma) % nlevels
    if ( abs(lambda * profiles_ad (i_proma) % t (i_lev) - profiles_ad2 (i_proma) % t (i_lev)) > threshold ) &
      call error_stop( 'test AD 1 fails', i_lev ) 
    if ( abs(lambda * profiles_ad (i_proma) % q (i_lev) - profiles_ad2 (i_proma) % q (i_lev)) > threshold ) &
      call error_stop( 'test AD 2 fails', i_lev )
    if ( abs(lambda * profiles_ad (i_proma) % p (i_lev) - profiles_ad2 (i_proma) % p (i_lev)) > threshold ) &
       call error_stop( 'test AD 3 fails', i_lev )
  enddo
enddo


do i_proma = 1, kproma
  if ( opts_scatt%lusercfrac .and. &
     & ((abs(lambda * cld_profiles_ad (i_proma) % cfrac - cld_profiles_ad2 (i_proma) % cfrac)) > threshold)) &
    call error_stop( 'test AD 3b fails', i_lev )
  do i_lev = 1, cld_profiles_ad (i_proma) % nlevels
    if ( abs(lambda * cld_profiles_ad (i_proma) % ph (i_lev) - cld_profiles_ad2 (i_proma) % ph (i_lev)) > threshold ) &
      call error_stop( 'test AD 4 fails', i_lev )
    do i_type = 1, nhydro
      if ( abs(lambda * cld_profiles_ad (i_proma) % hydro (i_lev,i_type) - &
             & cld_profiles_ad2 (i_proma) % hydro (i_lev,i_type)) > threshold ) then
        write(default_err_unit,*) lambda, cld_profiles_ad (i_proma) % hydro (i_lev,i_type), &
                                  cld_profiles_ad2 (i_proma) % hydro (i_lev,i_type), threshold
        call error_stop( 'test AD 5 fails', i_lev, inum2=i_type )
      endif
      if ( abs(lambda * cld_profiles_ad (i_proma) % hydro_frac (i_lev,i_type) - &
             & cld_profiles_ad2 (i_proma) % hydro_frac (i_lev,i_type)) > threshold ) &
        call error_stop( 'test AD 6 fails', i_lev, inum2=i_type )
    enddo
  enddo
enddo

do i_proma = 1, kproma
  if ( abs(lambda * profiles_ad (i_proma) % s2m % t - profiles_ad2 (i_proma) % s2m % t) > threshold ) &
    call error_stop( 'test AD 12 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % q - profiles_ad2 (i_proma) % s2m % q) > threshold ) &
    call error_stop( 'test AD 13 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % p - profiles_ad2 (i_proma) % s2m % p) > threshold ) &
    call error_stop( 'test AD 14 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % u - profiles_ad2 (i_proma) % s2m % u) > threshold ) &
    call error_stop( 'test AD 15 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % v - profiles_ad2 (i_proma) % s2m % v) > threshold ) &
    call error_stop( 'test AD 16 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % skin % t - profiles_ad2 (i_proma) % skin % t) > threshold ) &
    call error_stop( 'test AD 17 fails', i_proma ) 
enddo

do i_chan = 1, nchannels
  if ( abs(lambda * emissivity_ad (i_chan) % emis_in - emissivity_ad2 (i_chan) % emis_in ) > threshold ) &
    call error_stop( 'test AD 18 fails', i_chan )
enddo

do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)
  write(ioout,*) '2 - Test Equality of Norms'
  write(ioout,*)
enddo

call set_perturbation ( kproma, profiles_tl, cld_profiles_tl, profiles_d1, cld_profiles_d1, zeps)

call rttov_init_prof(profiles_ad)
call rttov_init_scatt_prof(cld_profiles_ad)

emissivity_d1 (1:nchannels) % emis_in = 0._jprb
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) % emis_in < 0.01_jprb

emissivity_tl (1:nchannels) % emis_in  = emissivity_d1 (1:nchannels) % emis_in * zeps
emissivity_ad (1:nchannels) % emis_in  = 0._jprb
emissivity_ad (1:nchannels) % emis_out = 0._jprb

radiance_tl % bt_clear(:) = 0._jprb
radiance_tl % bt(:)       = 0._jprb
if(lradar) then
  reflectivity_tl % zef(:,:) = 0.0_JPRB
  reflectivity_tl % azef(:,:) = 0.0_JPRB
endif
call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl,        &! in
  & cld_profiles_tl,    &! in
  & emissivity_tl,      &! inout
  & radiance_d3,        &! inout
  & radiance_tl,        &! inout
  & reflectivity=reflectivity_d3,    &! inout
  & reflectivity_tl=reflectivity_tl)  ! inout

If ( errorstatus == errorstatus_fatal ) Then
  call error_stop2('rttov_scatt_tl error')
End If

!* compute <subtl(delta_x),delta_z>

zdelta1 = 0._jprb

do i_chan = 1, nchannels
  zdelta1 = zdelta1 + (radiance_tl % bt(i_chan)) ** 2.0_jprb 
  zdelta1 = zdelta1 + (radiance_tl % bt_clear(i_chan)) ** 2.0_jprb  
  if(lradar) then
    zdelta1 = zdelta1 + sum(reflectivity_tl % zef(:,i_chan) ** 2.0_jprb)
    zdelta1 = zdelta1 + sum(reflectivity_tl % azef(:,i_chan) ** 2.0_jprb)
  endif
enddo
   
!* Initialize     
call rttov_init_rad(radiance_ad)
call rttov_init_rad(radiance_d1)

radiance_ad % bt_clear(:)   = radiance_tl % bt_clear (:) 
radiance_ad % bt(:)         = radiance_tl % bt (:) 
if(lradar) then
  reflectivity_ad % zef(:,:)  = reflectivity_tl % zef(:,:)
  reflectivity_ad % azef(:,:) = reflectivity_tl % azef(:,:)
endif

!---------------------------
! Now run AD code with TL radiances in input
! move TL results to AD radiance increments
  
call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad,        &! in
  & cld_profiles_ad,    &! in
  & emissivity_ad,      &! inout
  & radiance_d2,        &! inout
  & radiance_ad,        &! inout
  & reflectivity=reflectivity_d2,    &! inout
  & reflectivity_ad=reflectivity_ad)  ! inout

If ( errorstatus == errorstatus_fatal ) Then
  call error_stop2('rttov_scatt_ad error 3')
End If

!* compute <delta_x,subad(delta_z)>

zdelta2 = 0._jprb

do i_proma = 1, kproma

  if(opts_scatt%lusercfrac) zdelta2 = zdelta2 &
    & + cld_profiles_tl (i_proma) % cfrac * cld_profiles_ad (i_proma) % cfrac 

  zdelta2 = zdelta2 + sum( cld_profiles_tl (i_proma) % hydro (:,:) * cld_profiles_ad (i_proma) % hydro(:,:) )
  zdelta2 = zdelta2 + sum( cld_profiles_tl (i_proma) % hydro_frac (:,:) * cld_profiles_ad (i_proma) % hydro_frac (:,:) )

  do i_lev = 1, kflevg + 1
    zdelta2 = zdelta2 &
      & + cld_profiles_tl (i_proma) % ph (i_lev) * cld_profiles_ad (i_proma) % ph (i_lev)
  enddo

  do i_lev = 1, kflevg  
    zdelta2 = zdelta2  &
      & + profiles_tl (i_proma) % p   (i_lev) * profiles_ad (i_proma) % p   (i_lev) &
      & + profiles_tl (i_proma) % t   (i_lev) * profiles_ad (i_proma) % t   (i_lev) &
      & + profiles_tl (i_proma) % q   (i_lev) * profiles_ad (i_proma) % q   (i_lev)  
  enddo

  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % p * profiles_ad (i_proma) % s2m % p
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % q * profiles_ad (i_proma) % s2m % q
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % o * profiles_ad (i_proma) % s2m % o
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % t * profiles_ad (i_proma) % s2m % t
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % u * profiles_ad (i_proma) % s2m % u
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % v * profiles_ad (i_proma) % s2m % v

  zdelta2 = zdelta2 + profiles_tl (i_proma) % skin % t * profiles_ad (i_proma) % skin % t  

  do i_fast = 1, fastem_sp     
    zdelta2 = zdelta2 + profiles_tl (i_proma) % skin % fastem(i_fast)*profiles_ad(i_proma)%skin%fastem(i_fast)
  enddo         

  do i_chan = 1, nchannels 
    zdelta2 = zdelta2 + emissivity_tl (i_chan) % emis_in * emissivity_ad (i_chan) % emis_in
  enddo
enddo 

if (.not. abs(zdelta2) > 0._jprb) then
  z = 1._jprb
else
  z = zdelta2
endif

do i_output=1,1
  ioout = iooutstreams(i_output)
  write (ioout,'(A10,f23.16)') 'delta1 = ', zdelta1
  write (ioout,'(A10,f23.16)') 'delta2 = ', zdelta2    
  write (ioout,fmt= &
    & '('' The difference is '',f22.1, '' times the zero of the machine '')') &
    & abs((zdelta2-zdelta1)/epsilon(z)/z) 
enddo

if( abs((zdelta2-zdelta1)/epsilon(z)/z) < 50) then
  do i_output=1,size(iooutstreams)
    ioout = iooutstreams(i_output)
    write (ioout,*) 'AD is OK'
  enddo
else
  call error_stop( 'Adjoint test fails', 0_jpim)
endif

!* K-TEST ***********************************************************************************

do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)

  write(ioout,*)
  write(ioout,*) 'Test K'
  write(ioout,*) '------'
  write(ioout,*)
enddo

call allocate_profs( kflevg, nchannels, nhydro, nhydro, profiles_k, cld_profiles_k, iallocate)
call rttov_init_prof(profiles_k)
call rttov_init_scatt_prof(cld_profiles_k)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_k, kflevg, iallocate)
if(lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_k, kflevg, iallocate)
endif

emissivity_d1 (1:nchannels) % emis_in = 0._jprb 
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) % emis_in < 0.01_jprb      
emissivity_k  (1:nchannels) % emis_in  = 0._jprb
emissivity_k  (1:nchannels) % emis_out = 0._jprb

call rttov_init_rad(radiance_k)
radiance_k % bt(:) = 1.0_jprb ! (ACTIVATES K behaviour in adjoint code)
if (lradar) then
  radiance_k % bt(:) = 0.0_jprb
  reflectivity_k % azef(:,:)  = 0.0_jprb
  reflectivity_k % azef(50,:) = 1.0_jprb ! AJGDB still need to add the nlevels dimension in the K, so test a single level for now
  reflectivity_k % zef(:,:)  = 0.0_jprb
endif

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & opts_scatt,         &! in
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_k,         &! in (ACTIVATES K behaviour in adjoint code: note dimension)
  & cld_profiles_k,     &! in
  & emissivity_k,       &! inout
  & radiance_d1,        &! inout
  & radiance_k,         &! inout
  & reflectivity=reflectivity_d1, &! inout
  & reflectivity_ad=reflectivity_k) ! Inout

If ( Any( abs(radiance_total_ref(:) - radiance_d1 % bt(:)) > threshold * radiance_total_ref(:) ))  Then
  do i_output=1,size(iooutstreams)
    write(default_err_unit,*) 'wrong forward model in K'
    write(default_err_unit,*) radiance_total_ref(:)
    write(default_err_unit,*) radiance_d1 % bt(:)
    write(default_err_unit,*) abs(radiance_total_ref(:)-radiance_d1 %bt(:)) / ( threshold * radiance_total_ref(:))
  enddo
  Stop 1
Endif

if(lradar) then
  do i_lev = 1, kflevg
    If ( Any( abs(reflectivity_total_ref(i_lev,:) - reflectivity_d1%azef(i_lev,:)) > &
              abs(threshold * reflectivity_d1%azef(i_lev,:))))  Then
      write(default_err_unit,*) 'Wrong forward model in K- reflectivity'
      write(default_err_unit,*) reflectivity_total_ref(i_lev,:)
      write(default_err_unit,*) reflectivity_d1%azef(i_lev,:)
      write(default_err_unit,*) abs(reflectivity_total_ref(i_lev,:)-reflectivity_d1%azef(i_lev,:))
      Stop 1
    Endif
  enddo
endif

!---------------------------
! Compares K to AD

do i_btout = 1, nchannels   

  call rttov_init_rad(radiance_ad)
  radiance_ad % bt (i_btout) = 1.0_jprb
  if (lradar) then
    radiance_ad % bt(:) = 0.0_jprb
    reflectivity_ad % azef(:,:)        = 0.0_jprb
    reflectivity_ad % azef(50,i_btout) = 1.0_jprb ! AJGDB still need to add the nlevels dimension 
                                                  !    in the K, so test a single level for now
    reflectivity_ad % zef(:,:)         = 0.0_jprb
  endif

  call rttov_init_prof(profiles_ad)
  call rttov_init_scatt_prof(cld_profiles_ad)

  emissivity_ad (1:nchannels) % emis_in  = 0._jprb
  emissivity_ad (1:nchannels) % emis_out = 0._jprb

  call rttov_scatt_ad ( &
    & errorstatus,        &! out
    & opts_scatt,         &! in
    & kflevg,             &! in
    & chanprof,           &! in
    & frequencies,        &! in
    & profiles_d1,        &! inout
    & cld_profiles_d1,    &! in
    & coef_rttov,         &! in
    & coef_scatt,         &! in
    & calcemiss,          &! in
    & emissivity_d1,      &! inout
    & profiles_ad,        &! in
    & cld_profiles_ad,    &! in
    & emissivity_ad,      &! inout
    & radiance_d2,        &! inout
    & radiance_ad,        &! inout
    & reflectivity=reflectivity_d2, & ! inout
    & reflectivity_ad=reflectivity_ad) ! Inout

  i_proma = chanprof(i_btout)%prof

  do i_lev = 1, profiles_ad(i_proma) % nlevels
    if ( abs (profiles_ad (i_proma) % p   (i_lev) - profiles_k (i_btout) % p   (i_lev)) > threshold ) &
      call error_stop( 'test K 1 fails',i_lev)
    if ( abs (profiles_ad (i_proma) % t   (i_lev) - profiles_k (i_btout) % t   (i_lev)) > threshold ) &
      call error_stop( 'test K 2 fails',i_lev)
    if ( abs (profiles_ad (i_proma) % q   (i_lev) - profiles_k (i_btout) % q  (i_lev)) > threshold ) &
      call error_stop( 'test K 3 fails',i_lev)    
  End Do

  if( opts_scatt%lusercfrac .and. &
    & (abs (cld_profiles_ad (i_proma) % cfrac - cld_profiles_k (i_btout) % cfrac) > threshold))&
     call error_stop( 'test K 4 fails',i_lev)  

  do i_lev = 1, cld_profiles_ad(i_proma) % nlevels
    if ( abs (cld_profiles_ad (i_proma) % ph   (i_lev) - cld_profiles_k (i_btout) % ph   (i_lev)) > threshold)&
      call error_stop( 'test K 6 fails',i_lev)
    do i_type = 1, nhydro
      if ( abs (cld_profiles_ad (i_proma) % hydro (i_lev,i_type) - &
        & cld_profiles_k (i_btout) % hydro (i_lev,i_type))>threshold)&
            call error_stop( 'test K 8 fails', i_lev, inum2=i_type)
      if ( abs (cld_profiles_ad (i_proma) % hydro_frac (i_lev,i_type) - &
        & cld_profiles_k (i_btout) % hydro_frac (i_lev,i_type)) > threshold)&
            call error_stop( 'test K 9 fails',i_lev, inum2=i_type)
    enddo
  enddo

  if ( abs (profiles_ad (i_proma)  % s2m % p - profiles_k (i_btout) %  s2m % p) > threshold ) &
    call error_stop( 'test K 13 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % q - profiles_k (i_btout) %  s2m % q) > threshold ) &
    call error_stop( 'test K 14 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % o - profiles_k (i_btout) %  s2m % o) > threshold ) &
    call error_stop( 'test K 15 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % t - profiles_k (i_btout) %  s2m % t) > threshold ) &
    call error_stop( 'test K 16 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % u - profiles_k (i_btout) %  s2m % u) > threshold ) &
    call error_stop( 'test K 17 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % v - profiles_k (i_btout) %  s2m % v) > threshold ) &
    call error_stop( 'test K 18 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % skin % t - profiles_k (i_btout) %  skin % t) > threshold ) &
    call error_stop( 'test K 19 fails',i_lev )
  if ( abs (emissivity_ad (i_btout) % emis_in - emissivity_k (i_btout) % emis_in ) > threshold ) &
    call error_stop('test K 23 fails',i_lev )
          
enddo

do i_output=1,size(iooutstreams)
  ioout = iooutstreams(i_output)
  write(ioout,*) 'K is ok'
  write(ioout,*)
  write(ioout,*) 'End of RTTOVSCATT tests'
  write(ioout,*)
enddo

call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_d1, cld_profiles_d1, ideallocate)
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_tl, cld_profiles_tl, ideallocate)
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_tl2, cld_profiles_tl2, ideallocate)
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_ad, cld_profiles_ad, ideallocate)
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_d2, cld_profiles_d2, ideallocate)
call allocate_profs( kflevg, kproma, nhydro, nhydro, profiles_ad2, cld_profiles_ad2, ideallocate)
call allocate_profs( kflevg, nchannels, nhydro, nhydro, profiles_k, cld_profiles_k, ideallocate)

call rttov_alloc_rad ( erroralloc, nchannels, radiance_d1, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d3, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl2, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d2, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad2, kflevg, ideallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_k, kflevg, ideallocate)

if (lradar) then
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d1, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d3, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_tl, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_tl2, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_d2, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_ad, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_ad2, kflevg, ideallocate)
  call rttov_alloc_reflectivity ( erroralloc, nchannels, reflectivity_k, kflevg, ideallocate)

  deallocate(reflectivity_d1, reflectivity_d2, reflectivity_d3, reflectivity_tl, reflectivity_tl2)
  deallocate(reflectivity_ad, reflectivity_ad2, reflectivity_k)
endif

IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TEST',1_jpim,ZHOOK_HANDLE)

!*******

contains

!*******

!* Allocate and prepare / deallocate profile structures
subroutine allocate_profs( kflevg, kproma, khydro, khydro_frac, profiles, cld_profiles, asw)

Integer(Kind=jpim),  Intent(in) :: kflevg       ! number of levels
Integer(Kind=jpim),  Intent(in) :: kproma       ! number of profiles
Integer(Kind=jpim),  Intent(in) :: khydro       ! number of hydrometeors
Integer(Kind=jpim),  Intent(in) :: khydro_frac  ! number of hydrometeor fractions
Integer(Kind=jpim),  Intent(in) :: asw          ! 1=allocate, 0=deallocate
Type(rttov_profile), Intent(inout)       :: profiles (kproma)
Type(rttov_profile_cloud), Intent(inout) :: cld_profiles (kproma)

Type(rttov_options) :: opts ! Defaults to everything switched off
integer(kind=jpim)  :: err
   
call rttov_alloc_prof( err, kproma, profiles, kflevg, opts, asw, init = .true._jplm)
call rttov_alloc_scatt_prof ( err, kproma, cld_profiles, kflevg, khydro, khydro_frac, asw, &
  init = .true._jplm)

! Test input data gas units are ppmv over moist air
if (asw == 1) profiles(:) % gas_units = 2

end subroutine

!*******

subroutine set_perturbation( kproma, profiles_tl, cld_profiles_tl, &
  & profiles_d, cld_profiles_d, zeps)

Integer(Kind=jpim), Intent(in) :: kproma       ! number of profiles (or channels in case of K/AD)
Type(rttov_profile),       Intent (inout) :: profiles_tl  (kproma)
Type(rttov_profile_cloud), Intent (inout) :: cld_profiles_tl (kproma)
Type(rttov_profile),       Intent (in)    :: profiles_d   (kproma)
Type(rttov_profile_cloud), Intent (in)    :: cld_profiles_d (kproma)
Real(Kind=jprb),           Intent (in)    :: zeps

do i_proma = 1, kproma

  !* Set perturbation
  
  profiles_tl (i_proma) % p (:)        = profiles_d (i_proma) % p (:)        * zeps
  profiles_tl (i_proma) % t (1:kflevg) = profiles_d (i_proma) % t (1:kflevg) * zeps
  profiles_tl (i_proma) % q (1:kflevg) = profiles_d (i_proma) % q (1:kflevg) * zeps

  cld_profiles_tl (i_proma) % cfrac           = cld_profiles_d (i_proma) % cfrac  * zeps

  cld_profiles_tl (i_proma) % ph (:)          = cld_profiles_d (i_proma) % ph (:) * zeps
  cld_profiles_tl (i_proma) % hydro(1:kflevg,:) = cld_profiles_d (i_proma) % hydro (1:kflevg,:) * zeps
  cld_profiles_tl (i_proma) % hydro_frac (1:kflevg,:) = cld_profiles_d (i_proma) % hydro_frac (1:kflevg,:) * zeps

  profiles_tl (i_proma) % s2m % p = profiles_d (i_proma) % s2m % p * zeps
  profiles_tl (i_proma) % s2m % q = profiles_d (i_proma) % s2m % q * zeps
  profiles_tl (i_proma) % s2m % o = profiles_d (i_proma) % s2m % o * zeps
  profiles_tl (i_proma) % s2m % t = profiles_d (i_proma) % s2m % t * zeps
  profiles_tl (i_proma) % s2m % u = profiles_d (i_proma) % s2m % u * zeps
  profiles_tl (i_proma) % s2m % v = profiles_d (i_proma) % s2m % v * zeps


  profiles_tl (i_proma) % s2m % wfetc = 0._jprb ! Unused (solar computations only)

  profiles_tl (i_proma) % skin % surftype = -1 
  profiles_tl (i_proma) % skin % watertype = -1   

  profiles_tl (i_proma) % skin % t          = profiles_d (i_proma) % skin % t          * zeps
  profiles_tl (i_proma) % skin % fastem (:) = profiles_d (i_proma) % skin % fastem (:) * zeps

  profiles_tl (i_proma) % zenangle   = -1._jprb
  profiles_tl (i_proma) % azangle    = -1._jprb
  profiles_tl (i_proma) % ctp    = 0._jprb
  profiles_tl (i_proma) % cfraction    = 0._jprb
  profiles_tl (i_proma) % sunzenangle   = 0._jprb
  profiles_tl (i_proma) % sunazangle    = 0._jprb
  profiles_tl (i_proma) % latitude      = 0._jprb
  profiles_tl (i_proma) % longitude     = 0._jprb
  profiles_tl (i_proma) % elevation     = 0._jprb
  profiles_tl (i_proma) % Be            = 0._jprb
  profiles_tl (i_proma) % cosbk         = 0._jprb

enddo 

end subroutine

!*******

subroutine error_stop( text, inum, inum2)

character (len=*), intent(in)   :: text
integer( kind=jpim), intent(in) :: inum
integer( kind=jpim), optional, intent(in) :: inum2

write(default_err_unit,*) ' -------------------------------------------------- '
write(default_err_unit,*) ' ------------------- TEST HAS FAILED -------------- '
write(default_err_unit,*)
if (present(inum2)) then
  write(default_err_unit,*) text, inum, inum2
else
  write(default_err_unit,*) text, inum
endif
write(default_err_unit,*)
write(default_err_unit,*) ' -------------------------------------------------- '
write(default_err_unit,*) ' -------------------------------------------------- '

stop 1

end subroutine

subroutine error_stop2( text)
  character (len=*), intent(in)   :: text
  do i_output=1,size(iooutstreams)
    ioout = iooutstreams(i_output)
    write (ioout, * ) text
  enddo
  stop 1
end subroutine

End subroutine rttovscatt_test_one
