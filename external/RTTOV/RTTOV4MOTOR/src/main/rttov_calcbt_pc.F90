! Description:
!> @file
!!   Calculate brightness temperatures from radiances for PC-RTTOV
!
!> @brief
!!   Calculate brightness temperatures from radiances for PC-RTTOV
!!
!! @details
!!   The Planck function is calculated using the channel central wavenumbers
!!   as specified in the coefficient file. To adjust for the finite spectral
!!   bandwidth the temperatures are modified using band correction coefficients
!!   also stored in the coefficient file.
!!
!!   Radiances units are mW/cm-1/ster/m2 and temperature units are Kelvin.
!!
!!   Reference:
!!   Sharp, J.C, 1983: A comparison of approximate methods for converting
!!   between radiance and equivalent black body temperature for a radiometer
!!   channel, Met Office technical report (see docs/ directory)
!!
!! @param[in]     chanprof_in      specifies channels and profiles for reconstructed radiances
!! @param[in]     coef_pccomp      PC-RTTOV coefficient structure
!! @param[in,out] pccomp           PC-RTTOV radiance structure
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
SUBROUTINE rttov_calcbt_pc(chanprof_in, coef_pccomp, pccomp)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef_pccomp, rttov_pccomp
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
!INTF_END
  REAL   (KIND=jprb) :: tstore2
  INTEGER(KIND=jpim) :: chan, i
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header ------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC', 0_jpim, ZHOOK_HANDLE)
  nchannels_rec = SIZE(chanprof_in)
  DO i = 1, nchannels_rec
    chan = chanprof_in(i)%chan
    tstore2 = coef_pccomp%planck2_in(chan) / LOG(1 + coef_pccomp%planck1_in(chan) / pccomp%total_pccomp(i))
    pccomp%bt_pccomp(i) = (tstore2 - coef_pccomp%ff_bco_in(chan)) / coef_pccomp%ff_bcs_in(chan)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCBT_PC', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcbt_pc
