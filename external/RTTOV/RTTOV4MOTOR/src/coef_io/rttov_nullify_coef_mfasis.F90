! Description:
!> @file
!!   Nullify/zero an MFASIS LUT structure.
!
!> @brief
!!   Nullify/zero an MFASIS LUT structure.
!!
!!
!! @param[in,out]  coef_mfasis     the MFASIS LUT structure to nullify/zero
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
SUBROUTINE rttov_nullify_coef_mfasis(coef_mfasis)
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_mfasis
  IMPLICIT NONE
  TYPE(rttov_coef_mfasis), INTENT(INOUT) :: coef_mfasis
!INTF_END

  coef_mfasis%file_type      = 0_jpim
  coef_mfasis%version        = 0_jpim
  coef_mfasis%readme_lut     = 'xxxx'
  coef_mfasis%nparticles     = 0_jpim
  coef_mfasis%ndims          = 0_jpim
  coef_mfasis%nchannels      = 0_jpim
  coef_mfasis%nchannels_coef = 0_jpim
  coef_mfasis%clw_scheme     = 0_jpim
  coef_mfasis%ice_scheme     = 0_jpim
  coef_mfasis%maxzenangle    = 0._jprb

  NULLIFY(coef_mfasis%aer_types)
  NULLIFY(coef_mfasis%channel_list)
  NULLIFY(coef_mfasis%channel_lut_index)
  NULLIFY(coef_mfasis%lut_axes)
  NULLIFY(coef_mfasis%lut)
END SUBROUTINE rttov_nullify_coef_mfasis
