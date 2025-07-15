! Description:
!> @file
!!   AD of PC-RTTOV brightness temperature calculation
!
!> @brief
!!   AD of PC-RTTOV brightness temperature calculation
!!
!! @param[in]     chanprof_in      specifies channels and profiles for reconstructed radiances
!! @param[in]     coef_pccomp      PC-RTTOV coefficient structure
!! @param[in]     pccomp           PC-RTTOV radiance structure
!! @param[in,out] pccomp_ad        PC-RTTOV radiance increments
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
SUBROUTINE rttov_calcbt_pc_ad(chanprof_in, coef_pccomp, pccomp, pccomp_ad)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef_pccomp, rttov_pccomp
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_pccomp     ), INTENT(IN)    :: pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_ad
!INTF_END
  REAL   (KIND=jprb) :: tstar
  REAL   (KIND=jprb) :: tstar_ad
  REAL   (KIND=jprb) :: radtotal
  INTEGER(KIND=jpim) :: chan, i
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels_rec = SIZE(chanprof_in)
  DO i = 1, nchannels_rec
    chan = chanprof_in(i)%chan
! clear radiance
! T star for direct model
    tstar = coef_pccomp%ff_bco_in(chan) + coef_pccomp%ff_bcs_in(chan) * pccomp%bt_pccomp(i)
    radtotal = pccomp%total_pccomp(i)
! AD
    tstar_ad = pccomp_ad%bt_pccomp(i) / coef_pccomp%ff_bcs_in(chan)
    pccomp_ad%total_pccomp(i) = coef_pccomp%planck1_in(chan) * tstar ** 2 / &
        (coef_pccomp%planck2_in(chan) * radtotal * (radtotal + coef_pccomp%planck1_in(chan))) * tstar_ad
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcbt_pc_ad
