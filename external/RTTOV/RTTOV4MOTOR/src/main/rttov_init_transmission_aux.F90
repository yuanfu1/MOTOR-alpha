! Description:
!> @file
!!   Initialise internal transmission_aux structure.
!
!> @brief
!!   Initialise internal transmission_aux structure.
!!
!! @param[in]     opts              options to configure the simulations
!! @param[in,out] transmission_aux  transmission_aux structure to initialise
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_init_transmission_aux(opts, transmission_aux)

  USE rttov_types, ONLY : rttov_options, rttov_transmission_aux
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),          INTENT(IN)    :: opts
  TYPE(rttov_transmission_aux), INTENT(INOUT) :: transmission_aux
!INTF_END

!- End of header --------------------------------------------------------

  ! Don't need to initialise these:
!   IF (ASSOCIATED (transmission_aux%fac1))            transmission_aux%fac1            = 0._jprb

  ! Thermal path1
  IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf))       transmission_aux%thermal_path1%tau_surf       = 0._jprb
  IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level))      transmission_aux%thermal_path1%tau_level      = 0._jprb
  IF (ASSOCIATED(transmission_aux%thermal_path1%od_singlelayer)) transmission_aux%thermal_path1%od_singlelayer = 0._jprb
  IF (ASSOCIATED(transmission_aux%thermal_path1%od_sfrac))       transmission_aux%thermal_path1%od_sfrac       = 0._jprb
  IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_ac))    transmission_aux%thermal_path1%tau_surf_ac    = 0._jprb

  IF (ASSOCIATED(transmission_aux%solar_path2)) THEN
    IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf))       transmission_aux%solar_path2%tau_surf       = 0._jprb
    IF (ASSOCIATED(transmission_aux%solar_path2%tau_level))      transmission_aux%solar_path2%tau_level      = 0._jprb
    IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf_ac))    transmission_aux%solar_path2%tau_surf_ac    = 0._jprb

    IF (ASSOCIATED(transmission_aux%solar_path1%tau_level))      transmission_aux%solar_path1%tau_level      = 0._jprb
    IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf))       transmission_aux%solar_path1%tau_surf       = 0._jprb
    IF (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) THEN
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_level_r))    transmission_aux%solar_path1%tau_level_r    = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf_ac))    transmission_aux%solar_path1%tau_surf_ac    = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path1%od_singlelayer)) transmission_aux%solar_path1%od_singlelayer = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path2%od_singlelayer)) transmission_aux%solar_path2%od_singlelayer = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path1%od_sfrac))       transmission_aux%solar_path1%od_sfrac       = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path2%od_sfrac))       transmission_aux%solar_path2%od_sfrac       = 0._jprb
      IF (ASSOCIATED(transmission_aux%solar_path2%od_frac_ac))     transmission_aux%solar_path2%od_frac_ac     = 0._jprb
    ENDIF
  ENDIF
END SUBROUTINE 
