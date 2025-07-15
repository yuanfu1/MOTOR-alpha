! Description:
!> @file
!!   Copy an internal raytracing structure.
!
!> @brief
!!   Copy an internal raytracing structure.
!!
!! @param[in]     addsolar       flag indicating whether solar simulations are being performed
!! @param[in,out] raytracing1    copy of raytracing structure
!! @param[in]     raytracing2    input raytracing structure
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_copy_raytracing(addsolar, raytracing1, raytracing2)

  USE rttov_types, ONLY : rttov_raytracing
  USE parkind1, ONLY : jplm
  IMPLICIT NONE
  LOGICAL(jplm),          INTENT(IN)    :: addsolar
  TYPE(rttov_raytracing), INTENT(INOUT) :: raytracing1
  TYPE(rttov_raytracing), INTENT(IN)    :: raytracing2
!INTF_END

  ! Only copy those members which are used outside of locpat
  raytracing1%zasat    = raytracing2%zasat
  raytracing1%pathsat  = raytracing2%pathsat
  raytracing1%ltick    = raytracing2%ltick
  raytracing1%hgpl     = raytracing2%hgpl
  raytracing1%dmair    = raytracing2%dmair
  IF (addsolar) THEN
    raytracing1%zasun    = raytracing2%zasun
    raytracing1%pathsun  = raytracing2%pathsun
    raytracing1%patheff  = raytracing2%patheff
  ENDIF
END SUBROUTINE 
