! Description:
!> @file
!!   TL of reconstructed radiance calculation for PC-RTTOV simulations
!
!> @brief
!!   TL of reconstructed radiance calculation for PC-RTTOV simulations
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof_in     internal chanprof structure for reconstructed radiances
!! @param[in]     chanprof_pc     internal chanprof structure for PC predictor radiances
!! @param[in,out] pccomp_tl       PC scores and reconstructed radiance perturbations
!! @param[in]     coef_pccomp     PC-RTTOV coefficient structure
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
SUBROUTINE rttov_reconstruct_tl( &
             opts,        &
             chanprof_in, &
             chanprof_pc, &
             pccomp_tl,   &
             coef_pccomp)

  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp, rttov_options
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp),      INTENT(INOUT) :: pccomp_tl
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
  INTEGER(KIND=jpim) :: i, j, k, prof, chan
  INTEGER(KIND=jpim) :: nchannels_rec, nchannels_rec_p
  INTEGER(KIND=jpim) :: npcscores, npcscores_p
  INTEGER(KIND=jpim) :: nprofiles
  REAL(KIND=jprb)    :: rad_tl
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_TL', 0_jpim, ZHOOK_HANDLE)

  nchannels_rec = SIZE(chanprof_in)
  npcscores     = SIZE(chanprof_pc)
  nprofiles = chanprof_in(nchannels_rec)%prof

  nchannels_rec_p = nchannels_rec / nprofiles
  npcscores_p = npcscores / nprofiles

! DAR: RTTOV 11.2 - use pre-calculated transpose to improve memory locality
  k = 1
  DO prof = 1, nprofiles
    DO i = 1, nchannels_rec_p
      chan = chanprof_in(i)%chan
      rad_tl = 0._jprb
      DO j = 1, npcscores_p
        rad_tl = rad_tl + &
          coef_pccomp%eigen(opts%rt_ir%pc%ipcbnd)%eigenvectors_t(j, chan) * &
          pccomp_tl%total_pcscores((prof-1) * npcscores_p + j)
      ENDDO

      pccomp_tl%total_pccomp(k) = rad_tl * coef_pccomp%noise_in(chan)
      k = k + 1
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_reconstruct_tl
