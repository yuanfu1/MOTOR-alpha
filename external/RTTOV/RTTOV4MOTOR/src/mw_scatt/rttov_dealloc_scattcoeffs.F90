! Description:
!> @file
!!   Deallocate scattering coefficients.
!
!> @brief
!!   Deallocate scattering coefficients.
!!
!! @details
!!   Deallocates the hydrotable structure passed as an argument.
!!
!! @param[in,out]  coef_scatt   hydrotable structure to deallocate
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_dealloc_scattcoeffs(coef_scatt)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_scatt_coef
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_scatt_coef), INTENT(INOUT) :: coef_scatt
!INTF_END

#include "rttov_nullify_scattcoeffs.interface"
#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: err
  REAL(KIND=jprb) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_DEALLOC_SCATTCOEFFS',0_jpim,ZHOOK_HANDLE)

  IF (ASSOCIATED(coef_scatt%offset_temp)) DEALLOCATE(coef_scatt%offset_temp, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%freq)) DEALLOCATE(coef_scatt%freq, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%mpol)) DEALLOCATE(coef_scatt%mpol, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%is_frozen)) DEALLOCATE(coef_scatt%is_frozen, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%ext)) DEALLOCATE(coef_scatt%ext, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%ssa)) DEALLOCATE(coef_scatt%ssa, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%asp)) DEALLOCATE(coef_scatt%asp, STAT = err)
  THROW(err .NE. 0)
  IF (ASSOCIATED(coef_scatt%zef)) DEALLOCATE(coef_scatt%zef, STAT = err)
  THROW(err .NE. 0)

  CALL rttov_nullify_scattcoeffs(coef_scatt)
  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_DEALLOC_SCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
END SUBROUTINE rttov_dealloc_scattcoeffs
