! Description:
!> @file
!!   K of calculation of PC scores for PC-RTTOV simulations
!
!> @brief
!!   K of calculation of PC scores for PC-RTTOV simulations
!!
!! @details
!!   This differs to the AD because of the way the PC-RTTOV Jacobians are
!!   calculated (see rttov_mult_profiles_k). This subroutines calculates
!!   PC scores Jacobians in the case were reconstructed radiances are not
!!   computed (cf rttov_pcscores_rec_k).
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanprof_pc     internal chanprof structure for PC predictor radiances
!! @param[in]     pccomp_k        PC scores increments
!! @param[in]     coef_pccomp     PC-RTTOV coefficient structure
!! @param[in,out] total_k_pc      Jacobians of PC scores wrt predictor channel radiances
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
SUBROUTINE rttov_pcscores_k( &
             opts,        &
             chanprof,    &
             chanprof_pc, &
             pccomp_k,    &
             coef_pccomp, &
             total_k_pc)

  USE rttov_types, ONLY :  &
        rttov_options,     &
        rttov_chanprof,    &
        rttov_pccomp,      &
        rttov_coef_pccomp
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp),      INTENT(IN)    :: pccomp_k
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  REAL(KIND=jprb),         INTENT(INOUT) :: total_k_pc(:,:,:)
!INTF_END
  INTEGER(KIND=jpim) :: i, j, chan, chan_pc, prof
  INTEGER(KIND=jpim) :: nprofiles, nchannels
  INTEGER(KIND=jpim) :: npcscores_p, nchannels_p
  REAL   (KIND=jprb) :: noise_r(SIZE(coef_pccomp%noise_r))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 0_jpim, ZHOOK_HANDLE)

  npcscores_p = SIZE(total_k_pc, 2)
  nprofiles   = SIZE(total_k_pc, 3)
  nchannels   = SIZE(chanprof)

  nchannels_p = nchannels / nprofiles 

  DO i = 1, nchannels_p
    chan = chanprof(i)%chan
    noise_r(i) = coef_pccomp%noise_r(chan)
  ENDDO

! DAR: Update for RTTOV 11.2
! DAR: RTTOV 11.2 - Improve pcscores performance by using transposed total_k_pc
! array facilitating improved use of vectorisation
  DO prof = 1, nprofiles
    DO j = 1, npcscores_p
      chan_pc = chanprof_pc(j)%chan

      total_k_pc(1:nchannels_p, j, prof) = &
        pccomp_k%total_pcscores((prof - 1) * npcscores_p + j) * &
        coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg)% &
        coefficients(1:nchannels_p, chan_pc) * &
        noise_r(1:nchannels_p)
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores_k
