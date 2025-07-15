! Description:
!> @file
!!   Allocate/deallocate the user optical parameter structure for
!!   aerosols or clouds.
!
!> @brief
!!   Allocate/deallocate the user optical parameter structure for
!!   aerosols or clouds.
!!
!! @details
!!   This structure is used to pass in profiles of aerosol or cloud
!!   optical properties for each channel being simulated. For thermal
!!   channels these consist of absorption and scattering coefficients
!!   and the back-scattering parameter. Phase functions are only
!!   required for solar single-scattering simulations or for calculating
!!   Legendre coefficients or back-scattering parameters for thermal
!!   scattering calculations.
!!
!!   The init argument can be used to zero the newly allocated structure.
!!
!!   Note that for solar simulations, you must call the rttov_init_opt_param
!!   subroutine for the direct model optical property structure without the
!!   zero_only argument to pre-compute some data related to phase functions.
!!   This must be done after you store the phase function in the structure.
!!
!! @param[out]    err            status on exit
!! @param[in,out] opt_param      optical property structure
!! @param[in]     nchanprof      total number of channels being simulated
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     nmom           number of Legendre coefficients for each phase fn
!!                                  NB indexed 1:nmom+1, first coefficient is always 1.
!!                                  (for non-DOM simulations this can be 0)
!! @param[in]     nphangle       number of angles over which phase functions are
!!                                 specified (for non-solar simulations this can be 0)
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to zero newly allocated structure, optional
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
SUBROUTINE rttov_alloc_opt_param( &
              err,       &
              opt_param, &
              nchanprof, &
              nlevels,   &
              nmom,      &
              nphangle,  &
              asw,       &
              init)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_opt_param

  IMPLICIT NONE
  INTEGER(KIND=jpim),    INTENT(INOUT)        :: err
  TYPE(rttov_opt_param), INTENT(INOUT)        :: opt_param
  INTEGER(KIND=jpim),    INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),    INTENT(IN)           :: nlevels
  INTEGER(KIND=jpim),    INTENT(IN)           :: nmom
  INTEGER(KIND=jpim),    INTENT(IN)           :: nphangle
  INTEGER(KIND=jpim),    INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),    INTENT(IN), OPTIONAL :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_alloc_phfn_int.interface"
#include "rttov_init_opt_param.interface"

  LOGICAL(KIND=jplm) :: init1

  TRY
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1_jpim) THEN
    NULLIFY (opt_param%abs,                    &
             opt_param%sca,                    &
             opt_param%bpr,                    &
             opt_param%legcoef,                &
             opt_param%pha,                    &
             opt_param%phangle,                &
             opt_param%phasefn_int%cosphangle, &
             opt_param%phasefn_int%iphangle)

    ALLOCATE (opt_param%abs(nlevels-1,nchanprof),              &
              opt_param%sca(nlevels-1,nchanprof),              &
              opt_param%bpr(nlevels-1,nchanprof),              &
              opt_param%legcoef(1:nmom+1,nlevels-1,nchanprof), &
              opt_param%pha(nphangle,nlevels-1,nchanprof),     &
              opt_param%phangle(nphangle),                   &
              STAT = err) ! arrays in phasefn_int are allocated by calling rttov_init_opt_param
    THROWM(err.NE.0,"Allocation of opt_param failed")

    IF (init1) CALL rttov_init_opt_param(err, opt_param=opt_param, zero_only=.TRUE._jplm)

    opt_param%nmom = nmom
  ENDIF

  IF (asw == 0_jpim) THEN
    CALL rttov_alloc_phfn_int(err, opt_param%phangle, opt_param%phasefn_int, asw)
    THROW(err.NE.0)

    IF (ASSOCIATED(opt_param%abs))     DEALLOCATE(opt_param%abs)
    IF (ASSOCIATED(opt_param%sca))     DEALLOCATE(opt_param%sca)
    IF (ASSOCIATED(opt_param%bpr))     DEALLOCATE(opt_param%bpr)
    IF (ASSOCIATED(opt_param%legcoef)) DEALLOCATE(opt_param%legcoef)
    IF (ASSOCIATED(opt_param%pha))     DEALLOCATE(opt_param%pha)
    IF (ASSOCIATED(opt_param%phangle)) DEALLOCATE(opt_param%phangle)
    NULLIFY (opt_param%abs,     &
             opt_param%sca,     &
             opt_param%bpr,     &
             opt_param%legcoef, &
             opt_param%pha,     &
             opt_param%phangle)
  ENDIF

  CATCH
END SUBROUTINE rttov_alloc_opt_param
