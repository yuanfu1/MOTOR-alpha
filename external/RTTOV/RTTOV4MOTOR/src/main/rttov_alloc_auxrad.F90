! Description:
!> @file
!!   Allocate/deallocate internal auxiliary radiance structure containing
!!   Planck radiance calculations.
!
!> @brief
!!   Allocate/deallocate internal auxiliary radiance structure containing
!!   Planck radiance calculations.
!!
!! @details
!!   This structure holds calculated Planck radiances.
!!
!!   Some arrays are only required by the direct model instance
!!   of this structure: these are allocated if the optional "direct"
!!   argument is TRUE. If omitted, it is assumed to be FALSE.
!!
!! @param[out]    err                    status on exit
!! @param[in,out] auxrad                 auxiliary radiance structure to (de)allocate
!! @param[in]     nlevels                number of input profile levels
!! @param[in]     nchanprof              total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw                    1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     direct                 if FALSE then direct-model-only arrays are not allocated, optional
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
SUBROUTINE rttov_alloc_auxrad( &
              err,       &
              auxrad,    &
              nlevels,   &
              nchanprof, &
              asw,       &
              direct)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_radiance_aux
  IMPLICIT NONE
  INTEGER(KIND=jpim),       INTENT(OUT)          :: err
  TYPE(rttov_radiance_aux), INTENT(INOUT)        :: auxrad
  INTEGER(KIND=jpim),       INTENT(IN)           :: nlevels
  INTEGER(KIND=jpim),       INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),       INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),       INTENT(IN), OPTIONAL :: direct
!INTF_END

  LOGICAL(KIND=jplm) :: direct1

#include "rttov_errorreport.interface"

  TRY
  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (asw == 1) THEN
    CALL nullify_struct()

    IF (direct1) THEN
      ALLOCATE( &
        auxrad%air_t_eff(nlevels, nchanprof), &
        auxrad%surf_t_eff(nchanprof), &
        auxrad%skin_t_eff(nchanprof), &
        auxrad%cosmic_t_eff(nchanprof), &
        auxrad%cosmic(nchanprof), STAT = err)
      THROWM(err.NE.0,"Allocation of auxrad failed")
    ENDIF

    ALLOCATE( &
      auxrad%air(nlevels, nchanprof), &
      auxrad%surfair(nchanprof), &
      auxrad%skin(nchanprof), STAT = err)
    THROWM(err.NE.0,"Allocation of auxrad failed")
  ENDIF

  IF (asw == 0) THEN
    IF (ASSOCIATED(auxrad%air_t_eff))    DEALLOCATE(auxrad%air_t_eff)
    IF (ASSOCIATED(auxrad%surf_t_eff))   DEALLOCATE(auxrad%surf_t_eff)
    IF (ASSOCIATED(auxrad%skin_t_eff))   DEALLOCATE(auxrad%skin_t_eff)
    IF (ASSOCIATED(auxrad%cosmic_t_eff)) DEALLOCATE(auxrad%cosmic_t_eff)
    IF (ASSOCIATED(auxrad%cosmic))       DEALLOCATE(auxrad%cosmic)
    IF (ASSOCIATED(auxrad%air))          DEALLOCATE(auxrad%air)
    IF (ASSOCIATED(auxrad%surfair))      DEALLOCATE(auxrad%surfair)
    IF (ASSOCIATED(auxrad%skin))         DEALLOCATE(auxrad%skin)
    CALL nullify_struct()
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY(auxrad%air_t_eff,    &
            auxrad%surf_t_eff,   &
            auxrad%skin_t_eff,   &
            auxrad%cosmic_t_eff, &
            auxrad%air,          &
            auxrad%surfair,      &
            auxrad%skin,         &
            auxrad%cosmic)
  END SUBROUTINE nullify_struct
END SUBROUTINE 
