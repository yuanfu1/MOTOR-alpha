program rttovscatt_test
!
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
! Top level control of the RTTOV-SCATT test routines
!
!       2005   Peter Bauer        First version
! 13/12/2007   Alan Geer          RTTOV9 version
!
!-------------------------------------------------------
use parkind1, only: jpim, jprb, jplm
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

Use rttov_types, only : rttov_coefs, rttov_scatt_coef, rttov_options, rttov_chanprof, rttov_options_scatt

use mod_rttovscatt_test, only: kproma, zenangle, iooutstreams

implicit none

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttovscatt_test_one.interface"

integer (kind=jpim) :: nchannels,n_chan 
real    (kind=jprb), allocatable :: emissivity    (:), surfem (:)
integer (kind=jpim), allocatable :: frequencies   (:)
Type(rttov_chanprof), allocatable :: chanprof(:) ! Channel and profile indices

type (rttov_coefs     ) :: coef_rttov
type (rttov_scatt_coef) :: coef_scatt

Type(rttov_options)       :: opts ! defaults to everything off
Type(rttov_options_scatt) :: opts_scatt

character(len=256)  :: coef_filename
integer (kind=jpim) :: errorstatus
logical(kind=jplm) :: lradar
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

!* Read coefficient filename
IF (LHOOK) CALL DR_HOOK('RTTOVSCATT_TEST',0_jpim,ZHOOK_HANDLE)

read (*,'(a)') coef_filename
read (*,*) zenangle
read (*,*) lradar


! Test is configured for hydrometeor TL/AD sensitivity only on cloudy levels
opts_scatt%zero_hydro_tlad = .false.


!* Read coefficients
call rttov_read_coefs (errorstatus,coef_rttov, opts, file_coef = coef_filename)

call rttov_read_scattcoeffs (errorstatus, opts_scatt, coef_rttov, coef_scatt)
  
n_chan=coef_rttov%coef%fmv_chn   ! channels on  instrument
nchannels = kproma * n_chan ! total channels to simulate

allocate (chanprof     (nchannels))
allocate (frequencies  (nchannels))
allocate (emissivity   (nchannels))
allocate (surfem       (nchannels))

surfem (:) = 0.0_JPRB
emissivity = 0.0_JPRB

call rttov_scatt_setupindex ( &
        & errorstatus,  & ! out
        & kproma,       & ! in  
        & n_chan,       & ! in 
        & coef_rttov,   & ! in
        & coef_scatt,   & ! in
        & nchannels,    & ! in
        & chanprof,     & ! out
        & frequencies)    ! out

open (iooutstreams(1), file='outputscatt_full.ascii',form='formatted')
open (iooutstreams(2), file='outputscatt.ascii',form='formatted')

! Test previous basic config (no ice polarisation)
opts_scatt%lusercfrac = .FALSE._jplm
opts_scatt%ice_polarisation = -1.0_jprb
call rttovscatt_test_one ( nchannels,         &
          & opts_scatt,                       &
          & coef_rttov, coef_scatt         ,  &
          & chanprof                       ,  &
          & frequencies                    ,  &
          & emissivity                     ,  &
          & lradar)

! Test with imposed cloud fraction
opts_scatt%lusercfrac = .TRUE._jplm
call rttovscatt_test_one ( nchannels,         &
          & opts_scatt,                       &
          & coef_rttov, coef_scatt         ,  &
          & chanprof                       ,  &
          & frequencies                    ,  &
          & emissivity                     ,  &
          & lradar)

! Test with polarisation active (final v13 configuration)
opts_scatt%lusercfrac = .FALSE._jplm
opts_scatt%ice_polarisation = 1.3_jprb
call rttovscatt_test_one ( nchannels,         &
          & opts_scatt,                       &
          & coef_rttov, coef_scatt         ,  &
          & chanprof                       ,  &
          & frequencies                    ,  &
          & emissivity                     ,  &
          & lradar)

close (iooutstreams(1))
close (iooutstreams(2))

deallocate (chanprof, frequencies, emissivity, surfem)

Call rttov_dealloc_coefs( errorstatus, coef_rttov )
Call rttov_dealloc_scattcoeffs( coef_scatt )

IF (LHOOK) CALL DR_HOOK('RTTOVSCATT_TEST',1_jpim,ZHOOK_HANDLE)
 
end program rttovscatt_test
