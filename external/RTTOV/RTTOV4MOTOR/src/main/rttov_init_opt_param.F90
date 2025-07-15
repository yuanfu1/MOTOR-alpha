! Description:
!> @file
!!   Initialise auxiliary arrays in the user optical parameter structure for
!!   aerosols or clouds.
!
!> @brief
!!   Initialise auxiliary arrays in the user optical parameter structure for
!!   aerosols or clouds.
!!
!! @details
!!   This subroutine serves two purposes depending on the zero_only argument.
!!
!!   If zero_only is FALSE or is not present:
!!
!!     This subroutine precalculates some data used for interpolating the phase
!!     functions. The phangle array MUST FIRST have been populated with the phase
!!     angle grid. This subroutine MUST be called before calling RTTOV when running
!!     solar scattering simulations with explicit input optical properties.
!!     It must be called for the input aerosol and/or cloud optical parameter
!!     structures separately, but only for the *direct* model inputs.
!!
!!     If you are making multiple calls to RTTOV this subroutine only needs
!!     to be called again if the angle grid on which the phase functions are
!!     specified changes.
!!
!!     The opts argument is mandatory in this case.
!!
!!   If the zero_only argument is TRUE:
!!
!!     This subroutine zeros the contents of the opt_param structure. This is
!!     useful for AD/K simulations where the AD/K opt_param structures should
!!     be zero before calling rttov_ad/k.
!!
!!     The opts argument may be omitted in this case.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           RTTOV options structure, optional (required if zero_only omitted or false)
!! @param[in,out] opt_param      optical property structure
!! @param[in]     zero_only      flag to request zeroing of arrays only, optional
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
SUBROUTINE rttov_init_opt_param(err, opts, opt_param, zero_only)

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_options, rttov_opt_param
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : errorstatus_fatal
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),    INTENT(INOUT)          :: err
  TYPE(rttov_options),   INTENT(IN),   OPTIONAL :: opts
  TYPE(rttov_opt_param), INTENT(INOUT)          :: opt_param
  LOGICAL(KIND=jplm),    INTENT(IN),   OPTIONAL :: zero_only
!INTF_END

  LOGICAL(KIND=jplm) :: lzero_only

#include "rttov_errorreport.interface"
#include "rttov_alloc_phfn_int.interface"
  TRY

  lzero_only = .FALSE.
  IF (PRESENT(zero_only)) lzero_only = zero_only

  IF (lzero_only) THEN

    opt_param%nmom = 0
    IF (ASSOCIATED(opt_param%abs)) opt_param%abs = 0._jprb
    IF (ASSOCIATED(opt_param%sca)) opt_param%sca = 0._jprb
    IF (ASSOCIATED(opt_param%bpr)) opt_param%bpr = 0._jprb
    IF (ASSOCIATED(opt_param%legcoef)) opt_param%legcoef = 0._jprb
    IF (ASSOCIATED(opt_param%phangle)) opt_param%phangle = 0._jprb
    IF (ASSOCIATED(opt_param%pha)) opt_param%pha = 0._jprb

  ELSE

    IF (.NOT. PRESENT(opts)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "opts argument is mandatory")
    ENDIF

    IF (opts%rt_ir%addsolar) THEN

      IF (.NOT. ASSOCIATED(opt_param%phangle)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "opt_param%phangle must be allocated: call rttov_alloc_opt_param first")
      ENDIF

      IF (ASSOCIATED(opt_param%phasefn_int%cosphangle)) DEALLOCATE(opt_param%phasefn_int%cosphangle, STAT = err)
      THROWM(err.NE.0, "deallocation of aux phase array cosphangle")
      IF (ASSOCIATED(opt_param%phasefn_int%iphangle)) DEALLOCATE(opt_param%phasefn_int%iphangle, STAT = err)
      THROWM(err.NE.0, "deallocation of aux phase array iphangle")

      CALL rttov_alloc_phfn_int(err, opt_param%phangle, opt_param%phasefn_int, 1_jpim)
      THROW(err.NE.0)

    ENDIF

  ENDIF

  CATCH
END SUBROUTINE rttov_init_opt_param
