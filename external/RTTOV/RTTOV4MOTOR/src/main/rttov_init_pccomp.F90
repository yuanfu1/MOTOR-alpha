! Description:
!> @file
!!   Initialise a pccomp structure.
!
!> @brief
!!   Initialise a pccomp structure.
!!
!! @param[in,out]  pccomp   PC-RTTOV pccomp structure to initialise
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
SUBROUTINE rttov_init_pccomp(pccomp)

  USE rttov_types, ONLY : rttov_pccomp
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp
!INTF_END
  IF (ASSOCIATED(pccomp%clear_pcscores))    pccomp%clear_pcscores    = 0._jprb
  IF (ASSOCIATED(pccomp%total_pcscores))    pccomp%total_pcscores    = 0._jprb
  IF (ASSOCIATED(pccomp%overcast_pcscores)) pccomp%overcast_pcscores = 0._jprb
  IF (ASSOCIATED(pccomp%cloudy_pcscores))   pccomp%cloudy_pcscores   = 0._jprb

  IF (ASSOCIATED(pccomp%clear_pccomp))      pccomp%clear_pccomp      = 0._jprb
  IF (ASSOCIATED(pccomp%total_pccomp))      pccomp%total_pccomp      = 0._jprb
  IF (ASSOCIATED(pccomp%overcast_pccomp))   pccomp%overcast_pccomp   = 0._jprb
  IF (ASSOCIATED(pccomp%cloudy_pccomp))     pccomp%cloudy_pccomp     = 0._jprb

  IF (ASSOCIATED(pccomp%bt_pccomp))         pccomp%bt_pccomp         = 0._jprb
  IF (ASSOCIATED(pccomp%bt_clear_pccomp))   pccomp%bt_clear_pccomp   = 0._jprb
END SUBROUTINE 
