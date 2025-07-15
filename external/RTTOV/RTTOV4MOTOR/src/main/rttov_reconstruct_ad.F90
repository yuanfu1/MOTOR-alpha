! Description:
!> @file
!!   AD of reconstructed radiance calculation for PC-RTTOV simulations
!
!> @brief
!!   AD of reconstructed radiance calculation for PC-RTTOV simulations
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof_in     internal chanprof structure for reconstructed radiances
!! @param[in]     chanprof_pc     internal chanprof structure for PC predictor radiances
!! @param[in,out] pccomp_ad       PC scores and reconstructed radiance increments
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
SUBROUTINE rttov_reconstruct_ad( &
             opts,        &
             chanprof_in, &
             chanprof_pc, &
             pccomp_ad,   &
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
  TYPE(rttov_pccomp),      INTENT(INOUT) :: pccomp_ad
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
  INTEGER(KIND=jpim) :: nchannels_rec, nchannels_rec_p
  INTEGER(KIND=jpim) :: npcscores, npcscores_p
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: i, j, prof, chan
  REAL(KIND=jprb)    :: rad_ad_array(SIZE(chanprof_in))
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels_rec = SIZE(chanprof_in)
  npcscores     = SIZE(chanprof_pc)
  nprofiles = chanprof_in(nchannels_rec)%prof

  nchannels_rec_p = nchannels_rec / nprofiles
  npcscores_p = npcscores / nprofiles

  DO i=1, nchannels_rec
    chan = chanprof_in(i)%chan
    rad_ad_array(i) = pccomp_ad%total_pccomp(i) * coef_pccomp%noise_in(chan)
  ENDDO

! DAR: RTTOV 11.2 - swap loops to reduce memory accesses and pre-calculate rad_ad
  DO prof = 1, nprofiles
    DO j = 1, npcscores_p
      DO i = 1, nchannels_rec_p
        chan = chanprof_in(i)%chan
        pccomp_ad%total_pcscores((prof-1) * npcscores_p + j) = &
            pccomp_ad%total_pcscores((prof-1) * npcscores_p + j) + &
            coef_pccomp%eigen(opts%rt_ir%pc%ipcbnd)%eigenvectors(chan, j) * &
            rad_ad_array((prof-1) * nchannels_rec_p + i)
      ENDDO
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_reconstruct_ad
