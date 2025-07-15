! Description:
!> @file
!!   Copy a PC components structure.
!
!> @brief
!!   Copy a PC components structure.
!!
!! @param[in,out] pccomp1    copy of PC components structure
!! @param[in]     pccomp2    input PC components structure
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
SUBROUTINE rttov_copy_pccomp(pccomp1, pccomp2)

  USE rttov_types, ONLY : rttov_pccomp
  IMPLICIT NONE
  TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp1
  TYPE(rttov_pccomp), INTENT(IN)    :: pccomp2
!INTF_END
    IF (ASSOCIATED(pccomp1%total_pcscores) .AND. ASSOCIATED(pccomp2%total_pcscores)) THEN
      pccomp1%total_pcscores = pccomp2%total_pcscores
    ENDIF
    IF (ASSOCIATED(pccomp1%total_pccomp) .AND. ASSOCIATED(pccomp2%total_pccomp)) THEN
      pccomp1%total_pccomp = pccomp2%total_pccomp
    ENDIF
    IF (ASSOCIATED(pccomp1%bt_pccomp) .AND. ASSOCIATED(pccomp2%bt_pccomp)) THEN
      pccomp1%bt_pccomp = pccomp2%bt_pccomp
    ENDIF
    IF (ASSOCIATED(pccomp1%clear_pcscores) .AND. ASSOCIATED(pccomp2%clear_pcscores)) THEN
      pccomp1%clear_pcscores = pccomp2%clear_pcscores
    ENDIF
    IF (ASSOCIATED(pccomp1%clear_pccomp) .AND. ASSOCIATED(pccomp2%clear_pccomp)) THEN
      pccomp1%clear_pccomp = pccomp2%clear_pccomp
    ENDIF
    IF (ASSOCIATED(pccomp1%bt_clear_pccomp) .AND. ASSOCIATED(pccomp2%bt_clear_pccomp)) THEN
      pccomp1%bt_clear_pccomp = pccomp2%bt_clear_pccomp
    ENDIF
    IF (ASSOCIATED(pccomp1%overcast_pcscores) .AND. ASSOCIATED(pccomp2%overcast_pcscores)) THEN
      pccomp1%overcast_pcscores = pccomp2%overcast_pcscores
    ENDIF
    IF (ASSOCIATED(pccomp1%cloudy_pcscores) .AND. ASSOCIATED(pccomp2%cloudy_pcscores)) THEN
      pccomp1%cloudy_pcscores = pccomp2%cloudy_pcscores
    ENDIF
    IF (ASSOCIATED(pccomp1%overcast_pccomp) .AND. ASSOCIATED(pccomp2%overcast_pccomp)) THEN
      pccomp1%overcast_pccomp = pccomp2%overcast_pccomp
    ENDIF
    IF (ASSOCIATED(pccomp1%cloudy_pccomp) .AND. ASSOCIATED(pccomp2%cloudy_pccomp)) THEN
      pccomp1%cloudy_pccomp = pccomp2%cloudy_pccomp
    ENDIF
END SUBROUTINE 
