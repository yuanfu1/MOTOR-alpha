! Description:
!> @file
!!   Compute PC score or reconstructed radiance/BT increments suitable for use
!!   in AD
!
!> @brief
!!   Compute PC score or reconstructed radiance/BT increments suitable for use
!!   in AD
!!
!! @details
!!   The increment is calculated in PC scores if addradrec is false,
!!   otherwise it is calculated in radiances if swtchrad is false or BTs
!!   if switchrad is true.
!!
!! @param[in,out] pccomp_inc       computed PC score/radiance/BT increments
!! @param[in]     opts             options to configure the simulations
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_make_pccomp_inc(pccomp_inc, opts)

  USE rttov_types, ONLY : rttov_pccomp, rttov_options
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_pccomp),  INTENT(INOUT) :: pccomp_inc
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END
  INTEGER(KIND=jpim) :: i

  IF (opts%rt_ir%pc%addradrec) THEN
    IF (opts%rt_all%switchrad) THEN
      DO i = 1, SIZE(pccomp_inc%bt_pccomp)
        pccomp_inc%bt_pccomp(i) = 0.01_jprb * (MODULO(i, 113_jpim) + 1)
      ENDDO
    ELSE
      DO i = 1, SIZE(pccomp_inc%total_pccomp)
        pccomp_inc%total_pccomp(i) = 0.01_jprb * (MODULO(i, 113_jpim) + 1)
      ENDDO
    ENDIF
  ELSE
    DO i = 1, SIZE(pccomp_inc%total_pcscores)
      pccomp_inc%total_pcscores(i) = 0.01_jprb * (MODULO(i, 113_jpim) + 1)
    ENDDO
  ENDIF
END SUBROUTINE rttov_make_pccomp_inc
