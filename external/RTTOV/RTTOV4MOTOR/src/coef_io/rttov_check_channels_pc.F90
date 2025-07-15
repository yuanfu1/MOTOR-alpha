! Description:
!> @file
!!   Check specified channel list matches the PC predictor channel set
!!   specified in the options.
!
!> @brief
!!   Check specified channel list matches the PC predictor channel set
!!   specified in the options.
!!
!! @details
!!   This is called by rttov_read_coefs to determine if the specified
!!   channel subset to be read in matches the specified regression set in
!!   the options.
!!
!!   First the channel list is compared to the predictor sets in the PC
!!   coefficient structure to find a match. The matched predictor set is
!!   then compared to the specified predictor set in options.
!!
!!   An error is returned if no matching predictor set is found or if the
!!   selection in options does not match the supplied channel list.
!!
!! @param[out]     err              status on exit
!! @param[in]      opts             RTTOV options structure
!! @param[in]      coefs            RTTOV coefficients structure
!! @param[in]      channels         list of channels
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
SUBROUTINE rttov_check_channels_pc(err, opts, coefs, channels)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_options, rttov_coefs

  IMPLICIT NONE

  INTEGER(jpim),       INTENT(INOUT) :: err
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(rttov_coefs),   INTENT(IN)    :: coefs
  INTEGER(jpim),       INTENT(IN)    :: channels(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim) :: ipcbnd, ipcreg, m, n
  !----------------------------------------------------------------------------
TRY

  ipcbnd = 0_jpim
  err = errorstatus_fatal
  DO m = 1, coefs%coef_pccomp%fmv_pc_bands
    ipcbnd = ipcbnd + 1
    ipcreg = 0_jpim
    DO n = 1, coefs%coef_pccomp%fmv_pc_sets(m)
      ipcreg = ipcreg + 1
      IF (SIZE(channels) == coefs%coef_pccomp%pcreg(m,n)%fmv_pc_npred) THEN
        IF (ALL((coefs%coef%ff_ori_chn(:) - coefs%coef_pccomp%pcreg(m,n)%predictindex(:)) == 0_jpim)) THEN
          err = 0_jpim
          EXIT
        ENDIF
      ENDIF
    ENDDO
    IF (err == 0) EXIT
  ENDDO
  THROWM(err .NE. 0, "invalid regression channel set")

  IF (opts%rt_ir%pc%ipcreg > 0 .AND. opts%rt_ir%pc%ipcreg /= ipcreg) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "mismatch with regression channel set selected in opts%rt_ir%pc%ipcreg")
  ENDIF

  IF (opts%rt_ir%pc%ipcbnd > 0 .AND. opts%rt_ir%pc%ipcbnd /= ipcbnd) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, "mismatch with spectral band selected in opts%rt_ir%pc%ipcbnd")
  ENDIF

CATCH
END SUBROUTINE rttov_check_channels_pc
