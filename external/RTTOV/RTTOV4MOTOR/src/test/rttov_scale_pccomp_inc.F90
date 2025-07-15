! Description:
!> @file
!!   Scales PC score or reconstructed radiance or brightness temperature
!!   increment
!
!> @brief
!!   Scales PC score or reconstructed radiance or brightness temperature
!!   increment
!!
!! @details
!!   The scale factor is applied to PC scores if addradrec is false,
!!   otherwise it is applied to radiances if swtchrad is false or BTs
!!   if switchrad is true.
!!
!! @param[in,out] pccomp_inc       PC score/radiance structure to scale
!! @param[in]     factor           scale factor
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
SUBROUTINE rttov_scale_pccomp_inc(pccomp_inc, factor, opts)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_pccomp, rttov_options
  IMPLICIT NONE
  TYPE(rttov_pccomp),  INTENT(INOUT) :: pccomp_inc
  REAL(jprb),          INTENT(IN)    :: factor
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END

  IF (opts%rt_ir%pc%addradrec) THEN
    IF (opts%rt_all%switchrad) THEN
      pccomp_inc%bt_pccomp = pccomp_inc%bt_pccomp  * factor
    ELSE
      pccomp_inc%total_pccomp = pccomp_inc%total_pccomp * factor
    ENDIF
  ELSE
    pccomp_inc%total_pcscores = pccomp_inc%total_pcscores * factor
  ENDIF

END SUBROUTINE rttov_scale_pccomp_inc
