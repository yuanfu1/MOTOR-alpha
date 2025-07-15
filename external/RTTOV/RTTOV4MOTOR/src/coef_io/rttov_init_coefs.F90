! Description:
!> @file
!!   Initialise the coefficients structure.
!
!> @brief
!!   Initialise the coefficients structure.
!!
!! @details
!!   This is automatically called after the coefficients are read by
!!   rttov_read_coefs. It allocates some additional arrays and carries
!!   out some useful pre-calculations.
!!
!! @param[out]     err      status on exit
!! @param[in]      opts     RTTOV options structure
!! @param[in,out]  coefs    the coefficients structure to initialise
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
SUBROUTINE rttov_init_coefs(err, opts, coefs)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),  INTENT(OUT)   :: err
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(rttov_coefs),   INTENT(INOUT) :: coefs
!INTF_END
  REAL(KIND=jprb) :: ZHOOK_HANDLE
#include "rttov_errorreport.interface"
#include "rttov_init_coef.interface"
#include "rttov_init_coef_scatt.interface"
#include "rttov_init_coef_pccomp.interface"
!- End of header --------------------------------------------------------
  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 0_jpim, ZHOOK_HANDLE)

  IF (.NOT. coefs%initialised) THEN

    CALL rttov_init_coef(err, coefs%coef)
    THROW(err.NE.0)

    IF (opts%rt_ir%addclouds) THEN
      CALL rttov_init_coef_scatt(err, coefs%coef, coefs%coef_scatt)
      THROW(err.NE.0)
    ENDIF

    IF (opts%rt_ir%pc%addpc) THEN
      CALL rttov_init_coef_pccomp(err, coefs%coef, coefs%coef_pccomp)
      THROW(err.NE.0)
    ENDIF

    coefs%initialised = .TRUE.
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_INIT_COEFS', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE
