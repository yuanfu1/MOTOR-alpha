! Description:
!> @file
!!   Compute PC scores for PC-RTTOV simulations
!
!> @brief
!!   Compute PC scores for PC-RTTOV simulations
!!
!! @details
!!   PC scores are computed by multiplying the predictor channel radiances by
!!   the PC-RTTOV regression coefficients.
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanprof_pc     internal chanprof structure for PC predictor radiances
!! @param[in,out] pccomp          output PC scores
!! @param[in]     coef_pccomp     PC-RTTOV coefficient structure
!! @param[in]     radiance        output radiance structure (for predictor channels)
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
SUBROUTINE rttov_pcscores( &
             opts,         &
             chanprof,     &
             chanprof_pc,  &
             pccomp,       &
             coef_pccomp,  &
             radiance)

  USE rttov_types, ONLY :  &
        rttov_options,     &
        rttov_chanprof,    &
        rttov_pccomp,      &
        rttov_coef_pccomp, &
        rttov_radiance
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp),      INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_radiance),    INTENT(IN)    :: radiance
!INTF_END
  INTEGER(KIND=jpim) :: i, j, chan, chan_pc, prof
  INTEGER(KIND=jpim) :: nchannels, nchannels_p
  INTEGER(KIND=jpim) :: npcscores, npcscores_p
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: noise_r(SIZE(coef_pccomp%noise_r))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  npcscores = SIZE(chanprof_pc)
  nprofiles = chanprof(nchannels)%prof

  nchannels_p = nchannels / nprofiles
  npcscores_p = npcscores / nprofiles

! DAR: pre-calculate required channel noise
  DO i = 1, nchannels_p
    chan = chanprof(i)%chan
    noise_r(i) = coef_pccomp%noise_r(chan)
  ENDDO

! DAR: RTTOV 11.2 - Improve pcscores performance by reducing memory access and
! by improving use of vectorisation
  DO prof = 1, nprofiles
    DO j = (prof-1) * npcscores_p + 1, prof * npcscores_p
      chan_pc = chanprof_pc(j)%chan
      pccomp%total_pcscores(j) = DOT_PRODUCT( &
        coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg)%coefficients(:, chan_pc), &
        radiance%total((prof-1) * nchannels_p + 1 : prof * nchannels_p) * noise_r(1:nchannels_p))
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores
