! Description:
!> @file
!!   Allocate/deallocate internal raytracing structure.
!
!> @brief
!!   Allocate/deallocate internal raytracing structure.
!!
!! @details
!!   The raytracing structure contains information related to the
!!   calculation of the local angles of the radiation path through
!!   the atmosphere.
!!
!! @param[out]    err            status on exit
!! @param[in]     nprofiles      number of profiles
!! @param[in]     addsolar       flag indicating whether solar simulations are being performed
!! @param[in,out] raytracings    raytracing structure
!! @param[in]     nlevels        number of nlevels in input profiles
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_alloc_raytracing( &
              err,          &
              nprofiles,    &
              addsolar,     &
              raytracings,  &
              nlevels,      &
              asw,          &
              init)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_raytracing
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim),     INTENT(OUT)          :: err
  INTEGER(KIND=jpim),     INTENT(IN)           :: nprofiles
  LOGICAL(KIND=jplm),     INTENT(IN)           :: addsolar
  INTEGER(KIND=jpim),     INTENT(IN)           :: nlevels
  TYPE(rttov_raytracing), INTENT(INOUT)        :: raytracings
  INTEGER(KIND=jpim),     INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),     INTENT(IN), OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_raytracing.interface"

  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1) THEN
    ALLOCATE ( &
        raytracings%hgpl(nlevels, nprofiles), raytracings%dmair(nlevels, nprofiles),     &
        raytracings%r(nlevels, nprofiles), raytracings%refractivity(nlevels, nprofiles), &
        raytracings%r_r(nlevels, nprofiles), raytracings%z_r(nlayers, nprofiles),        &
        raytracings%ratoesat(nlayers, nprofiles), raytracings%zasat(nlayers, nprofiles), &
        raytracings%int(nlevels, nprofiles), raytracings%ztemp(nlevels, nprofiles),      &
        raytracings%ppw(nlevels, nprofiles), raytracings%dispco2(nlevels, nprofiles),    &
        raytracings%pathsat(nlayers, nprofiles), raytracings%ltick(nlayers, nprofiles), STAT = err)
    THROWM(err .NE. 0,"Allocation of raytracing failed")

    IF (addsolar) THEN
      ALLOCATE (                                    &
          raytracings%ratoesun(nlayers, nprofiles), &
          raytracings%zasun(nlayers, nprofiles),    &
          raytracings%pathsun(nlayers, nprofiles),  &
          raytracings%patheff(nlayers, nprofiles),  &
          STAT = err)
      THROWM(err .NE. 0,"Allocation of solar raytracing failed")
    ENDIF

    IF (init1) CALL rttov_init_raytracing(addsolar, raytracings)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE (                                        &
        raytracings%hgpl, raytracings%dmair,            &
        raytracings%r, raytracings%refractivity,        &
        raytracings%r_r, raytracings%z_r,               &
        raytracings%ratoesat, raytracings%zasat,        &
        raytracings%int, raytracings%ztemp,             &
        raytracings%ppw, raytracings%dispco2,           &
        raytracings%pathsat, raytracings%ltick, STAT = err)
    THROWM(err .NE. 0,"Deallocation of raytracing failed")

    IF (addsolar) THEN
      DEALLOCATE (              &
          raytracings%ratoesun, &
          raytracings%zasun,    &
          raytracings%pathsun,  &
          raytracings%patheff, STAT = err)
          THROWM(err .NE. 0,"Deallocation of solar raytracing failed")
    ENDIF

    NULLIFY (raytracings%hgpl)
    NULLIFY (raytracings%dmair)
    NULLIFY (raytracings%r)
    NULLIFY (raytracings%refractivity)
    NULLIFY (raytracings%r_r)
    NULLIFY (raytracings%z_r)
    NULLIFY (raytracings%ratoesun)
    NULLIFY (raytracings%ratoesat)
    NULLIFY (raytracings%zasun)
    NULLIFY (raytracings%zasat)
    NULLIFY (raytracings%int)
    NULLIFY (raytracings%ztemp)
    NULLIFY (raytracings%ppw)
    NULLIFY (raytracings%dispco2)
    NULLIFY (raytracings%pathsat)
    NULLIFY (raytracings%pathsun)
    NULLIFY (raytracings%patheff)
    NULLIFY (raytracings%ltick)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_raytracing
