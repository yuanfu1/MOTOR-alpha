! Description:
!> @file
!!   K of reconstructed radiance calculation for PC-RTTOV simulations
!
!> @brief
!!   K of reconstructed radiance calculation for PC-RTTOV simulations
!!
!! @details
!!   This differs to the AD because of the way the PC-RTTOV Jacobians are
!!   calculated (see rttov_mult_profiles_k). This is used to compute the
!!   Jacobians of reconstructed radiances wrt PC scores.
!!
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof_in     internal chanprof structure for reconstructed radiances
!! @param[in]     chanprof_pc     internal chanprof structure for PC predictor radiances
!! @param[in]     pccomp_k        PC scores and reconstructed radiance increments
!! @param[in,out] pcscores_k      output Jacobian of reconstructed radiances wrt PC scores
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

SUBROUTINE rttov_reconstruct_k( &
             opts,        &
             chanprof_in, &
             chanprof_pc, &
             pccomp_k,    &
             pcscores_k,  &
             coef_pccomp)

  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp, rttov_options
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_pc(:)
  REAL(KIND=jprb),         INTENT(INOUT) :: pcscores_k(:,:,:)
  TYPE(rttov_pccomp),      INTENT(INOUT) :: pccomp_k
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels_rec, nchannels_rec_p
  INTEGER(KIND=jpim) :: npcscores, npcscores_p
  INTEGER(KIND=jpim) :: i, prof, chan
  REAL(KIND=jprb)    :: rad_k
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_K', 0_jpim, ZHOOK_HANDLE)
  nchannels_rec = SIZE(chanprof_in)
  npcscores     = SIZE(chanprof_pc)
  nprofiles     = SIZE(pcscores_k, 3)

  nchannels_rec_p = nchannels_rec / nprofiles
  npcscores_p = npcscores / nprofiles

! DAR: RTTOV 11.2 - use transposed pcscores_k array to improve memory locality
  DO prof = 1, nprofiles
    DO i = 1, nchannels_rec_p
      chan = chanprof_in(i)%chan
      rad_k = pccomp_k%total_pccomp((prof - 1) * nchannels_rec_p + i) * &
              coef_pccomp%noise_in(chan)
      pcscores_k(1:npcscores_p, i, prof) = rad_k * &
        coef_pccomp%eigen(opts%rt_ir%pc%ipcbnd)% &
        eigenvectors_t(1:npcscores_p, chan)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_RECONSTRUCT_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_reconstruct_k
