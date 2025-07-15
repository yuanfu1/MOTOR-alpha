! Description:
!> @file
!!   Print out human-readable interpretation of a radiance%quality output.
!
!> @brief
!!   Print out human-readable interpretation of a radiance%quality output.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   quality   integer quality output from radiance structure
!! @param[in]   lu        logical unit for output, optional
!! @param[in]   text      additional text to print, optional
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_print_radiance_quality(quality, lu, text)

  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
  USE rttov_const, ONLY : &
    qflag_reg_limits,             &
    qflag_pc_aer_reg_limits,      &
    qflag_mfasis_zenangle,        &
    qflag_mfasis_sumzenangle,     &
    qflag_mfasis_geometry_bounds, &
    qflag_mfasis_opdpedia_bounds
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)           :: quality
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lu
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: text
!INTF_END

  INTEGER(KIND=jpim)  :: iu       ! logical unit for print

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (PRESENT(text)) WRITE(iu,'(a)') TRIM(text)

  IF (quality == 0) THEN
    WRITE(iu,'(a)') 'Quality OK - no bits set'
  ELSE
    IF (BTEST(quality, qflag_reg_limits)) THEN
      WRITE(iu,'(a)') 'Gas optical depth regression limits exceeded for some variable'
    ENDIF

    IF (BTEST(quality, qflag_pc_aer_reg_limits)) THEN
      WRITE(iu,'(a)') 'PC-RTTOV aerosol regression limits exceeded'
    ENDIF

    IF (BTEST(quality, qflag_mfasis_zenangle)) THEN
      WRITE(iu,'(a)') 'MFASIS warning: zenith angle greater than max valid value'
    ENDIF
    IF (BTEST(quality, qflag_mfasis_sumzenangle)) THEN
      WRITE(iu,'(a)') 'MFASIS warning: sum of zenith angles greater than max valid value'
    ENDIF
    IF (BTEST(quality, qflag_mfasis_geometry_bounds)) THEN
      WRITE(iu,'(a)') 'MFASIS warning: scattering angle out of LUT bounds'
    ENDIF
    IF (BTEST(quality, qflag_mfasis_opdpedia_bounds)) THEN
      WRITE(iu,'(a)') 'MFASIS warning: some optical depth/effective diamater value out of LUT bounds'
    ENDIF
  ENDIF

END SUBROUTINE rttov_print_radiance_quality
