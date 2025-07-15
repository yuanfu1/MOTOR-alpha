! Description:
!> @file
!!   Nullify scattering coefficients.
!
!> @brief
!!   Nullify scattering coefficients.
!!
!! @details
!!   Nuliifies the pointers in the hydrotable structure passed as an argument.
!!
!! @param[in,out]  coef_scatt   hydrotable structure to nullify
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
SUBROUTINE rttov_nullify_scattcoeffs(coef_scatt)

  USE rttov_types, ONLY : rttov_scatt_coef

  IMPLICIT NONE

  TYPE(rttov_scatt_coef), INTENT(INOUT) :: coef_scatt
!INTF_END

  NULLIFY(coef_scatt%is_frozen,   &
          coef_scatt%freq,        &
          coef_scatt%mpol,        &
          coef_scatt%offset_temp, &
          coef_scatt%ext,         &
          coef_scatt%ssa,         &
          coef_scatt%asp,         &
          coef_scatt%zef)

END SUBROUTINE rttov_nullify_scattcoeffs
