! Description:
!> @file
!!   Allocate an internal optical depth path structure.
!
!> @brief
!!   Allocate an internal optical depth path structure.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           RTTOV options structure
!! @param[in,out] opdp_path      opdep_path structure to allocate
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     nchanprof      total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_opdp_path( &
              err,       &
              opts,      &
              opdp_path, &
              nlevels,   &
              nchanprof, &
              asw,       &
              init)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_options, rttov_opdp_path
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim),    INTENT(OUT)            :: err
  TYPE(rttov_options),   INTENT(IN)             :: opts
  INTEGER(KIND=jpim),    INTENT(IN)             :: nlevels
  INTEGER(KIND=jpim),    INTENT(IN)             :: nchanprof
  TYPE(rttov_opdp_path), INTENT(INOUT)          :: opdp_path
  INTEGER(KIND=jpim),    INTENT(IN)             :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm),    INTENT(IN),   OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_opdp_path.interface"

  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1) THEN
    NULLIFY(opdp_path%atm_level)
    NULLIFY(opdp_path%sun_level_path2)

    ALLOCATE (opdp_path%atm_level(nlevels, nchanprof), STAT = err)
    THROWM(err.NE.0, "allocation of opdp_path%atm_level")
    IF (opts%rt_ir%addsolar) THEN
      ALLOCATE (opdp_path%sun_level_path2(nlevels, nchanprof), STAT = err)
      THROWM(err.NE.0, "allocation of opdp_path%sun_level_path2")
    ENDIF

    IF (init1) CALL rttov_init_opdp_path(opts, opdp_path)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE (opdp_path%atm_level, STAT = err)
    THROWM(err.NE.0, "deallocation of opdp_path%atm_level")
    IF (opts%rt_ir%addsolar) THEN
      DEALLOCATE (opdp_path%sun_level_path2, STAT = err)
      THROWM(err.NE.0, "deallocation of opdp_path%sun_level_path2")
    ENDIF

    NULLIFY(opdp_path%atm_level)
    NULLIFY(opdp_path%sun_level_path2)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_opdp_path
