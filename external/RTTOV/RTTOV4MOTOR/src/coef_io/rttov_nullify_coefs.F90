! Description:
!> @file
!!   Nullify/zero the RTTOV coefficients structure.
!
!> @brief
!!   Nullify/zero the RTTOV coefficients structure.
!!
!! @details
!!   This includes the optical depth coefficients, PC
!!   coefficients and cloud/aerosol scattering coefficients.
!!
!! @param[in,out]  coefs     the coefficient structure to nullify/zero
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
SUBROUTINE rttov_nullify_coefs(coefs)

  USE rttov_types, ONLY : rttov_coefs
  IMPLICIT NONE
  TYPE(rttov_coefs), INTENT(INOUT) :: coefs
!INTF_END

#include "rttov_nullify_coef.interface"
#include "rttov_nullify_coef_pccomp.interface"
#include "rttov_nullify_coef_mfasis.interface"
#include "rttov_nullify_coef_scatt.interface"
#include "rttov_nullify_coef_htfrtc.interface"

  CALL rttov_nullify_coef(coefs%coef)
  CALL rttov_nullify_coef_pccomp(coefs%coef_pccomp)
  CALL rttov_nullify_coef_mfasis(coefs%coef_mfasis_cld)
  CALL rttov_nullify_coef_mfasis(coefs%coef_mfasis_aer)
  CALL rttov_nullify_coef_scatt(coefs%coef_scatt)
  CALL rttov_nullify_coef_htfrtc(coefs%coef_htfrtc)

END SUBROUTINE 
