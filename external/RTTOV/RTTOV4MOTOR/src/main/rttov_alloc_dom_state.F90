! Description:
!> @file
!!   Allocate/deallocate structure for DOM direct model internal state.
!
!> @brief
!!   Allocate/deallocate structure for DOM direct model internal state.
!!
!! @details
!!   The DOM state structure is only used by the TL/AD/K to store
!!   some of the intermediate calculation results from the direct model
!!   run of the Discrete Ordinates Method algorithm. This enables the
!!   TL/AD/K to avoid recomputing these results which is more computationally
!!   efficient, but requires a certain amount of memory.
!!
!!   The DOM state is stored in an array of structures. The array is of
!!   size nchanprof and each structure is only allocated if the DOM algorithm
!!   is being used for that channel. There are separate structures for "solar"
!!   and "thermal" channels as the DOM algorithm is called separately to treat
!!   the thermal emission and solar source terms. The dosolar argument indicates
!!   whether the allocation is for the solar (true) or thermal (false) structures.
!!
!!   If the direct model is called on its own this internal state is not
!!   stored and this structure is not allocated.
!!
!! @param[out]    err              status on exit
!! @param         dom_state        pointer to RTTOV array of rttov_dom_state structures
!! @param[in]     nchanprof        total number of channels being simulated (SIZE(chanprof))
!! @param[in]     nlayers          number of input profile layers (nlevels-1)
!! @param[in]     ncolumns         number of cloud columns as calculated by rttov_cloud_overlap
!! @param[in]     dom_nstr         number of DOM streams (discrete ordinates)
!! @param[in]     dosolar          flag to indicate if this allocation is for solar-enabled channels
!! @param[in]     chanflag         flags to indicate which channels should be allocated (either channels with
!!                                 significant thermal or solar contributions)
!! @param[in]     asw              1_jpim => allocate; 0_jpim => deallocate
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
SUBROUTINE rttov_alloc_dom_state( &
              err,          &
              dom_state,    &
              nchanprof,    &
              nlayers,      &
              ncolumns,     &
              dom_nstr,     &
              dosolar,      &
              chanflag,     &
              asw)

#include "throw.h"
  USE rttov_types, ONLY : rttov_dom_state
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)          :: err
  TYPE(rttov_dom_state),   POINTER              :: dom_state(:)
  INTEGER(KIND=jpim),      INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),      INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim),      INTENT(IN)           :: ncolumns
  INTEGER(KIND=jpim),      INTENT(IN)           :: dom_nstr
  LOGICAL(KIND=jplm),      INTENT(IN)           :: dosolar
  LOGICAL(KIND=jplm),      INTENT(IN)           :: chanflag(nchanprof)
  INTEGER(KIND=jpim),      INTENT(IN)           :: asw             ! 1=allocate, 0=deallocate
!INTF_END
#include "rttov_errorreport.interface"

  INTEGER(jpim) :: i, maxnlayers, naz

  TRY

  IF (asw == 1_jpim) THEN

    maxnlayers = nlayers

    ALLOCATE(dom_state(nchanprof), STAT=err)
    THROWM(err.NE.0,"allocation of dom_state")

    CALL nullify_struct()

    DO i = 1, nchanprof
      IF (chanflag(i)) THEN
        IF (dosolar) THEN
          naz = dom_nstr - 1
          ALLOCATE(dom_state(i)%nazloops(0:ncolumns),            &
                   dom_state(i)%z(dom_nstr,0:nlayers,0:1,0:naz), &
                   STAT = err)
        ELSE ! dothermal
          naz = 0
          ALLOCATE(dom_state(i)%y0(dom_nstr,nlayers,0:ncolumns), &
                   dom_state(i)%y1(dom_nstr,nlayers,0:ncolumns), &
                   STAT = err)
        ENDIF
        THROWM(err.NE.0,"allocation of thermal/solar dom_state")

        ALLOCATE(dom_state(i)%thisrad(0:ncolumns),       &
                 dom_state(i)%radsurfup(0:ncolumns),     &  ! also 0:naz if surface not Lambertian
                 dom_state(i)%radsurfup_sum(0:ncolumns), &  ! also 0:naz if surface not Lambertian
                 dom_state(i)%eval(dom_nstr/2,0:nlayers,0:1,0:naz),         &
                 dom_state(i)%xp(dom_nstr/2,dom_nstr/2,nlayers,0:1,0:naz),  &
                 dom_state(i)%x(dom_nstr*maxnlayers,0:naz,0:ncolumns),      &
                 STAT = err)
        THROWM(err.NE.0,"allocation of dom_state")
      ENDIF
    ENDDO

  ENDIF

  IF (asw == 0_jpim .AND. ASSOCIATED(dom_state)) THEN

    DO i = 1, nchanprof
      IF (ASSOCIATED(dom_state(i)%nazloops))      DEALLOCATE(dom_state(i)%nazloops)
      IF (ASSOCIATED(dom_state(i)%z))             DEALLOCATE(dom_state(i)%z)
      IF (ASSOCIATED(dom_state(i)%y0))            DEALLOCATE(dom_state(i)%y0)
      IF (ASSOCIATED(dom_state(i)%y1))            DEALLOCATE(dom_state(i)%y1)
      IF (ASSOCIATED(dom_state(i)%thisrad))       DEALLOCATE(dom_state(i)%thisrad)
      IF (ASSOCIATED(dom_state(i)%radsurfup))     DEALLOCATE(dom_state(i)%radsurfup)
      IF (ASSOCIATED(dom_state(i)%radsurfup_sum)) DEALLOCATE(dom_state(i)%radsurfup_sum)
      IF (ASSOCIATED(dom_state(i)%eval))          DEALLOCATE(dom_state(i)%eval)
      IF (ASSOCIATED(dom_state(i)%xp))            DEALLOCATE(dom_state(i)%xp)
      IF (ASSOCIATED(dom_state(i)%x))             DEALLOCATE(dom_state(i)%x)
    ENDDO
    CALL nullify_struct()

    DEALLOCATE(dom_state, STAT=err)
    THROWM(err.NE.0,"deallocation of dom_state")

    NULLIFY(dom_state)

  ENDIF

  CATCH
CONTAINS

  SUBROUTINE nullify_struct()
    DO i = 1, nchanprof
      NULLIFY(dom_state(i)%nazloops,      &
              dom_state(i)%z,             &
              dom_state(i)%y0,            &
              dom_state(i)%y1,            &
              dom_state(i)%thisrad,       &
              dom_state(i)%radsurfup,     &
              dom_state(i)%radsurfup_sum, &
              dom_state(i)%x,             &
              dom_state(i)%xp,            &
              dom_state(i)%eval)
    ENDDO
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_dom_state
