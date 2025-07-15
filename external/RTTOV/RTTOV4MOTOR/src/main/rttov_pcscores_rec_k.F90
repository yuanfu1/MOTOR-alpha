! Description:
!> @file
!!   K of calculation of PC scores for PC-RTTOV simulations
!
!> @brief
!!   K of calculation of PC scores for PC-RTTOV simulations
!!
!! @details
!!   This is the PC scores Jacobian calculation for the case where
!!   reconstructed radiances are calculated (cf rttov_pcscores_k).
!!
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     pcscores_k      Jacobian of reconstructed radiances wrt PC scores
!! @param[in]     coef_pccomp     PC-RTTOV coefficient structure
!! @param[in,out] total_k_pc      Jacobians of reconstructed radiances wrt predictor channel radiances
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
SUBROUTINE rttov_pcscores_rec_k( &
             opts,        &
             chanprof,    &
             pcscores_k,  &
             coef_pccomp, &
             total_k_pc)

  USE rttov_types, ONLY :  &
        rttov_options,     &
        rttov_chanprof,    &
        rttov_coef_pccomp
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  REAL(KIND=jprb),         INTENT(IN)    :: pcscores_k(:,:,:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  REAL(KIND=jprb),         INTENT(INOUT) :: total_k_pc(:,:,:)
!INTF_END
  INTEGER(KIND=jpim) :: i, k, chan, prof
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels_p
  INTEGER(KIND=jpim) :: nchannels_rec_p
  INTEGER(KIND=jpim) :: npcscores_p
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_REC_K', 0_jpim, ZHOOK_HANDLE)
  nchannels_p     = SIZE(total_k_pc, 1)
  nchannels_rec_p = SIZE(total_k_pc, 2)
  nprofiles       = SIZE(total_k_pc, 3)
  npcscores_p     = SIZE(pcscores_k, 1)

! DAR: Update for RTTOV 11.2
! DAR: RTTOV 11.2 - Improve pcscores performance by using transposed total_k_pc
! array facilitating improved use of vectorisation. 
  DO prof = 1, nprofiles
    DO i = 1, nchannels_rec_p
      DO k = 1, nchannels_p
        chan = chanprof(k)%chan

        total_k_pc(k, i, prof) = &
          DOT_PRODUCT( &
            pcscores_k(1:npcscores_p, i, prof), &
            coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg)% &
            coefficients_t(1:npcscores_p, k)) * &
          coef_pccomp%noise_r(chan)
      ENDDO
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_REC_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores_rec_k
