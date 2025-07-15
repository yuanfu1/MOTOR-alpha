! Description:
!> @file
!!   Check input channel list against PC regression channel set
!
!> @brief
!!   Check input channel list against PC regression channel set
!!
!! @param[in]     nprofiles            number of input profiles
!! @param[in]     nchanprof            number of channels being simulated
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     chanprof             specifies channels and profiles to simulate
!! @param[in]     coefs                coefficients structure for instrument to simulate
!! @param[out]    err                  status on exit
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
SUBROUTINE rttov_checkpcchan( &
         nprofiles,         &
         nchanprof,         &
         opts,              &
         chanprof,          &
         coefs,             &
         err)
!INTF_OFF
#include "throw.h"
!INTF_ON
USE parkind1, ONLY : jpim
USE rttov_types, ONLY: rttov_options, rttov_chanprof, rttov_coefs
!INTF_OFF
USE parkind1, ONLY : jpim
!INTF_ON

  IMPLICIT NONE
  INTEGER(KIND=jpim),   INTENT(IN)  :: nprofiles
  INTEGER(KIND=jpim),   INTENT(IN)  :: nchanprof
  TYPE(rttov_options),  INTENT(IN)  :: opts
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof(:)
  TYPE(rttov_coefs),    INTENT(IN)  :: coefs
  INTEGER(KIND=jpim),   INTENT(OUT) :: err
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: i
  INTEGER(KIND=jpim) :: dc(nchanprof/nprofiles)
!-----End of header-------------------------------------------------------------

  TRY

  ! Test if the number of channels per profile is compatible with the regression set
  IF (nchanprof/nprofiles .NE. coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%fmv_pc_npred) THEN
    err = 1_jpim
  ENDIF
  THROWM(err.NE.0, "PC calc; invalid number of regression channels")

  IF (chanprof(1)%chan == 1 .AND. chanprof(nchanprof)%chan == nchanprof/nprofiles) THEN

    ! Test if channels are [1..N]
    DO prof = 1, nprofiles
      i = 1 + (prof-1) * (nchanprof/nprofiles)
      dc = coefs%coef%ff_ori_chn(chanprof(i:i+(nchanprof/nprofiles)-1)%chan) - &
           coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%predictindex(:)
      IF (ANY(dc .NE. 0_jpim)) err = 1_jpim
    ENDDO

  ELSE

    ! Test if channels are exactly the same as the predictors
    DO prof = 1, nprofiles
      i = 1 + (prof-1) * (nchanprof/nprofiles)
      dc = chanprof(i:i+(nchanprof/nprofiles)-1)%chan - &
           coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%predictindex(:)
      IF (ANY(dc .NE. 0_jpim)) err = 1_jpim
    ENDDO

  ENDIF
  THROWM(err.NE.0, "PC calc; invalid regression channels indices")

  CATCH
END SUBROUTINE
