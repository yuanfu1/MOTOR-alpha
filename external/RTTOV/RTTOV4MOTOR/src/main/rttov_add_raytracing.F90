! Description:
!> @file
!!   Adds two internal raytracing structures.
!
!> @brief
!!   Adds two internal raytracing structures.
!!
!! @param[in]     addsolar       flag indicating whether solar simulations are being performed
!! @param[in,out] raytracing     output summed raytracing structure
!! @param[in]     raytracing1    first input raytracing structure
!! @param[in]     raytracing2    second input raytracing structure
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
SUBROUTINE rttov_add_raytracing(addsolar, raytracing, raytracing1, raytracing2)

  USE rttov_types, ONLY : rttov_raytracing
  USE parkind1, ONLY : jplm
  IMPLICIT NONE
  LOGICAL(jplm),          INTENT(IN)    :: addsolar
  TYPE(rttov_raytracing), INTENT(INOUT) :: raytracing
  TYPE(rttov_raytracing), INTENT(IN)    :: raytracing1
  TYPE(rttov_raytracing), INTENT(IN)    :: raytracing2
!INTF_END

  ! Only add those members which are used outside of locpat
  raytracing%zasat    = raytracing1%zasat + raytracing2%zasat
  raytracing%pathsat  = raytracing1%pathsat + raytracing2%pathsat
  raytracing%ltick    = raytracing1%ltick + raytracing2%ltick
  raytracing%hgpl     = raytracing1%hgpl + raytracing2%hgpl
  raytracing%dmair    = raytracing1%dmair + raytracing2%dmair
  IF (addsolar) THEN
    raytracing%zasun    = raytracing1%zasun + raytracing2%zasun
    raytracing%pathsun  = raytracing1%pathsun + raytracing2%pathsun
    raytracing%patheff  = raytracing1%patheff + raytracing2%patheff
  ENDIF
END SUBROUTINE 
