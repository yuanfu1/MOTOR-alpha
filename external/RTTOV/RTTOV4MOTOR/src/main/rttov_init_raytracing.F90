! Description:
!> @file
!!   Initialise internal raytracing structure.
!
!> @brief
!!   Initialise internal raytracing structure.
!!
!! @details
!!   This subroutine is primarily intended for use by the AD/K models
!!   and as such only the members which are required to be initialised
!!   by these models are set to zero.
!!
!! @param[in]     addsolar       flag indicating whether solar simulations are being performed
!! @param[in,out] raytracings    raytracing structure
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
SUBROUTINE rttov_init_raytracing(addsolar, raytracings)

  USE rttov_types, ONLY : rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(jplm),          INTENT(IN)    :: addsolar
  TYPE(rttov_raytracing), INTENT(INOUT) :: raytracings
!INTF_END

  raytracings%ltick        = 0._jprb
  raytracings%hgpl         = 0._jprb
  raytracings%dmair        = 0._jprb
  raytracings%zasat        = 0._jprb
  raytracings%pathsat      = 0._jprb

  IF (addsolar) THEN
    raytracings%zasun        = 0._jprb
    raytracings%pathsun      = 0._jprb
    raytracings%patheff      = 0._jprb
  ENDIF
END SUBROUTINE 
