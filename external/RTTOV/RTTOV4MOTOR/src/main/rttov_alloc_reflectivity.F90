! Description:
!> @file
!!   Allocate/deallocate reflectivity structure.
!
!> @brief
!!   Allocate/deallocate reflectivity structure.
!!
!! @details
!!
!! @param[out]    err            status on exit
!! @param[in]     nchanprof      total number of channels being simulated (SIZE(chanprof))
!! @param[in,out] reflectivity   reflectivities
!! @param[in]     nlevels        number of nlevels in input profiles
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structure, optional
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
!    Copyright 2020, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_alloc_reflectivity( &
              err,          &
              nchanprof,    &
              reflectivity, &
              nlevels,      &
              asw,          &
              init)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_reflectivity
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),           INTENT(OUT)   :: err
  INTEGER(KIND=jpim),           INTENT(IN)    :: nchanprof
  TYPE(rttov_reflectivity),     INTENT(INOUT) :: reflectivity
  INTEGER(KIND=jpim),           INTENT(IN)    :: nlevels
  INTEGER(KIND=jpim),           INTENT(IN)    :: asw
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(IN)    :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_reflectivity.interface"

  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY

  IF (asw == 1) THEN
    init1   = .FALSE.
    IF (PRESENT(init)) init1 = init

    ALLOCATE (reflectivity%zef(nlevels, nchanprof), &
              reflectivity%azef(nlevels, nchanprof), STAT = err)
    THROWM(err.NE.0, "allocation of radar reflectivities in reflectivity type")
    
    IF (init1) CALL rttov_init_reflectivity(reflectivity)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE (reflectivity%zef, reflectivity%azef, STAT = err)
    THROWM(err.NE.0, "deallocation of radar reflectivities in reflectivity type")

    NULLIFY (reflectivity%zef, reflectivity%azef)
  ENDIF

  CATCH
END SUBROUTINE rttov_alloc_reflectivity
