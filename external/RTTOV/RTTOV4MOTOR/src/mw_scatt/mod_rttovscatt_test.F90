  module mod_rttovscatt_test  

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
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------

  Use parkind1, only: jpim  ,jprb
  
  Use rttov_types, only : rttov_coef, rttov_scatt_coef 

  IMPLICIT NONE

  integer (kind=jpim), parameter :: nulout = 5
  integer (kind=jpim), parameter :: nulerr = 6
  integer (kind=jpim), parameter :: kproma = 3
  integer (kind=jpim), parameter :: kflevg = 61
  
  integer (kind=jpim), parameter :: ioin   = 10
  
  ! Output streams: first is verbose output, second essential
  integer (kind=jpim), dimension(2), parameter :: iooutstreams  = (/20,21/)
  
  integer (kind=jpim), parameter :: inproc  = 1
  integer (kind=jpim), parameter :: imyproc = 1
  integer (kind=jpim), parameter :: iioproc = 1
  
  
  !* from module /satrad/module/onedvar_const.F90
  real    (kind=jprb), parameter :: fastem_land_coeff (5) = (/ 3.0_JPRB, 5.0_JPRB, 15.0_JPRB, 0.1_JPRB, 0.3_JPRB /)
  real    (kind=jprb), parameter :: fastem_ocean          = 0.0_JPRB
  
  real    (kind=jprb) :: zenangle
  
  end module mod_rttovscatt_test
