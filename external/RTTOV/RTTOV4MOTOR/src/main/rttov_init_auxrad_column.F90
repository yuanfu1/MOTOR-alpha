! Description:
!> @file
!!   Initialise internal dyanmically-sized auxiliary radiance
!!   structure.
!
!> @brief
!!   Initialise internal dyanmically-sized auxiliary radiance
!!   structure.
!!
!! @param[in,out] auxrad_column          auxiliary radiance structure to initialise
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
SUBROUTINE rttov_init_auxrad_column(auxrad_column)

  USE rttov_types, ONLY : rttov_radiance_aux
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_radiance_aux), INTENT(INOUT) :: auxrad_column
!INTF_END

  IF (ASSOCIATED(auxrad_column%up))                 auxrad_column%up                 = 0._jprb
  IF (ASSOCIATED(auxrad_column%down))               auxrad_column%down               = 0._jprb
  IF (ASSOCIATED(auxrad_column%down_p))             auxrad_column%down_p             = 0._jprb
  IF (ASSOCIATED(auxrad_column%up_solar))           auxrad_column%up_solar           = 0._jprb
  IF (ASSOCIATED(auxrad_column%down_solar))         auxrad_column%down_solar         = 0._jprb
  IF (ASSOCIATED(auxrad_column%down_ref))           auxrad_column%down_ref           = 0._jprb
  IF (ASSOCIATED(auxrad_column%down_p_ref))         auxrad_column%down_p_ref         = 0._jprb
  IF (ASSOCIATED(auxrad_column%down_ref_solar))     auxrad_column%down_ref_solar     = 0._jprb
  IF (ASSOCIATED(auxrad_column%meanrad_up))         auxrad_column%meanrad_up         = 0._jprb
  IF (ASSOCIATED(auxrad_column%meanrad_down))       auxrad_column%meanrad_down       = 0._jprb
  IF (ASSOCIATED(auxrad_column%meanrad_down_p))     auxrad_column%meanrad_down_p     = 0._jprb
  IF (ASSOCIATED(auxrad_column%meanrad_up_solar))   auxrad_column%meanrad_up_solar   = 0._jprb
  IF (ASSOCIATED(auxrad_column%meanrad_down_solar)) auxrad_column%meanrad_down_solar = 0._jprb
  IF (ASSOCIATED(auxrad_column%cloudy))             auxrad_column%cloudy             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac1_2))             auxrad_column%fac1_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac3_2))             auxrad_column%fac3_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac4_2))             auxrad_column%fac4_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac5_2))             auxrad_column%fac5_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac6_2))             auxrad_column%fac6_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac7_2))             auxrad_column%fac7_2             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac4_3))             auxrad_column%fac4_3             = 0._jprb
  IF (ASSOCIATED(auxrad_column%fac5_3))             auxrad_column%fac5_3             = 0._jprb

END SUBROUTINE rttov_init_auxrad_column
