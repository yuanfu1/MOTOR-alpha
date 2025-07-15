! Description:
!> @file
!!   Allocate/deallocate internal ircld structure.
!
!> @brief
!!   Allocate/deallocate internal ircld structure.
!!
!! @details
!!   The ircld structure contains information about the calculated
!!   cloud columns for cloudy IR scattering calculations.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           options to configure the simulations
!! @param[in]     nprofiles      number of profiles being simulated
!! @param[in,out] ircld          ircld structure to (de)allocate
!! @param[in]     nlayers        number of layers in input profiles
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     direct         set .TRUE. if allocating direct model ircld structure, optional
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
SUBROUTINE rttov_alloc_ircld( &
              err,       &
              opts,      &
              nprofiles, &
              ircld,     &
              nlayers,   &
              asw,       &
              init,      &
              direct)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_ircld, rttov_options
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)          :: err
  TYPE(rttov_options), INTENT(IN)           :: opts
  INTEGER(KIND=jpim),  INTENT(IN)           :: nprofiles
  INTEGER(KIND=jpim),  INTENT(IN)           :: nlayers
  TYPE(rttov_ircld),   INTENT(INOUT)        :: ircld
  INTEGER(KIND=jpim),  INTENT(IN)           :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_ircld.interface"

  LOGICAL(KIND=jplm) :: init1, direct1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (asw == 1) THEN
    CALL nullify_struct()

    ! This is used as a lookup: for each channel/layer/column it tells us
    ! whether the layer is cloudy or non-cloudy
    IF (direct1) THEN
      IF (opts%rt_ir%addclouds) THEN
        ALLOCATE (ircld%icldarr(0:2*nlayers, nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%icldarr")
      ELSE
        ALLOCATE (ircld%icldarr(0:0, nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%icldarr")
      ENDIF
    ENDIF

    IF (opts%rt_ir%addclouds) THEN
      ALLOCATE (ircld%xcol(2*nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%xcol")
      ALLOCATE (ircld%cldcfr(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%cldcfr")
      ALLOCATE (ircld%maxcov(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%maxcov")
      ALLOCATE (ircld%xcolmax(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%xcolmax")
      ALLOCATE (ircld%xcolmin(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%xcolmin")
      ALLOCATE (ircld%a(2*nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%a")

      IF (direct1) THEN
        ALLOCATE (ircld%xcolref1(2*nlayers, 2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%xcolref1")
        ALLOCATE (ircld%xcolref2(2*nlayers, 2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%xcolref2")
        ALLOCATE (ircld%xcolref(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%xcolref")
        ALLOCATE (ircld%xcolminref(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%xcolminref")
        ALLOCATE (ircld%ntotref(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%ntotref")
        ALLOCATE (ircld%indexcol(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%indexcol")
        ALLOCATE (ircld%icount1ref(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%icount1ref")
        ALLOCATE (ircld%iloopin(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%iloopin")
        ALLOCATE (ircld%flag(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%flag")
        ALLOCATE (ircld%iflag(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld%iflag")
      ENDIF

      IF (direct1) THEN
        ALLOCATE (ircld%ncolumnref(nprofiles), &
                  ircld%iloop(nprofiles),      &
                  ircld%icount(nprofiles),     &
                  ircld%icouncol(nprofiles),   &
                  ircld%icount1(nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of ircld")
      ENDIF
    ENDIF

    ALLOCATE (ircld%xcolclr(nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of ircld%xcolclr")
    ircld%xcolclr = 1._jprb

    IF (direct1) THEN
      ALLOCATE (ircld%ncolumn(nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of ircld%ncolumn")
      ircld%ncolumn = 0_jpim
    ENDIF

    IF (init1) THEN
      IF (direct1 .AND. opts%rt_ir%addclouds) THEN
        ircld%ncolumnref = 0_jpim
        ircld%iloop      = 0_jpim
        ircld%icount     = 0_jpim
        ircld%icouncol   = 0_jpim
        ircld%icount1    = 0_jpim
      ENDIF
      CALL rttov_init_ircld(ircld)
    ENDIF
  ENDIF

  IF (asw == 0) THEN
    IF (ASSOCIATED(ircld%icldarr)) THEN
      DEALLOCATE (ircld%icldarr, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%icldarr")
    ENDIF
    IF (ASSOCIATED(ircld%xcolref1)) THEN
      DEALLOCATE (ircld%xcolref1, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolref1")
    ENDIF
    IF (ASSOCIATED(ircld%xcolref2)) THEN
      DEALLOCATE (ircld%xcolref2, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolref2")
    ENDIF
    IF (ASSOCIATED(ircld%xcol)) THEN
      DEALLOCATE (ircld%xcol, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcol")
    ENDIF
    IF (ASSOCIATED(ircld%xcolminref)) THEN
      DEALLOCATE (ircld%xcolminref, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolminref")
    ENDIF
    IF (ASSOCIATED(ircld%xcolref)) THEN
      DEALLOCATE (ircld%xcolref, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolref")
    ENDIF
    IF (ASSOCIATED(ircld%cldcfr)) THEN
      DEALLOCATE (ircld%cldcfr, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%cldcfr")
    ENDIF
    IF (ASSOCIATED(ircld%maxcov)) THEN
      DEALLOCATE (ircld%maxcov, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%maxcov")
    ENDIF
    IF (ASSOCIATED(ircld%xcolmax)) THEN
      DEALLOCATE (ircld%xcolmax, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolmax")
    ENDIF
    IF (ASSOCIATED(ircld%xcolmin)) THEN
      DEALLOCATE (ircld%xcolmin, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolmin")
    ENDIF
    IF (ASSOCIATED(ircld%a)) THEN
      DEALLOCATE (ircld%a, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%a")
    ENDIF
    IF (ASSOCIATED(ircld%ntotref)) THEN
      DEALLOCATE (ircld%ntotref, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%ntotref")
    ENDIF
    IF (ASSOCIATED(ircld%indexcol)) THEN
      DEALLOCATE (ircld%indexcol, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%indexcol")
    ENDIF
    IF (ASSOCIATED(ircld%icount1ref)) THEN
      DEALLOCATE (ircld%icount1ref, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%icount1ref")
    ENDIF
    IF (ASSOCIATED(ircld%iloopin)) THEN
      DEALLOCATE (ircld%iloopin, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%iloopin")
    ENDIF
    IF (ASSOCIATED(ircld%flag)) THEN
      DEALLOCATE (ircld%flag, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%flag")
    ENDIF
    IF (ASSOCIATED(ircld%iflag)) THEN
      DEALLOCATE (ircld%iflag, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%iflag")
    ENDIF

    IF (ASSOCIATED(ircld%ncolumn)) THEN
      DEALLOCATE (ircld%ncolumn, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%ncolumn")
    ENDIF
    IF (ASSOCIATED(ircld%ncolumnref)) THEN
      DEALLOCATE (ircld%ncolumnref, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%ncolumnref")
    ENDIF
    IF (ASSOCIATED(ircld%iloop)) THEN
      DEALLOCATE (ircld%iloop, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%iloop")
    ENDIF
    IF (ASSOCIATED(ircld%icount)) THEN
      DEALLOCATE (ircld%icount, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%icount")
    ENDIF
    IF (ASSOCIATED(ircld%icouncol)) THEN
      DEALLOCATE (ircld%icouncol, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%icouncol")
    ENDIF
    IF (ASSOCIATED(ircld%icount1)) THEN
      DEALLOCATE (ircld%icount1, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%icount1")
    ENDIF
    IF (ASSOCIATED(ircld%xcolclr)) THEN
      DEALLOCATE (ircld%xcolclr, STAT = err)
      THROWM(err .NE. 0, "deallocation of ircld%xcolclr")
    ENDIF

    CALL nullify_struct()

  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    NULLIFY (ircld%icldarr)
    NULLIFY (ircld%xcolref1)
    NULLIFY (ircld%xcolref2)
    NULLIFY (ircld%xcol)
    NULLIFY (ircld%xcolminref)
    NULLIFY (ircld%xcolref)
    NULLIFY (ircld%cldcfr)
    NULLIFY (ircld%maxcov)
    NULLIFY (ircld%xcolmax)
    NULLIFY (ircld%xcolmin)
    NULLIFY (ircld%a)
    NULLIFY (ircld%ntotref)
    NULLIFY (ircld%indexcol)
    NULLIFY (ircld%icount1ref)
    NULLIFY (ircld%iloopin)
    NULLIFY (ircld%flag)
    NULLIFY (ircld%iflag)
    NULLIFY (ircld%ncolumn)
    NULLIFY (ircld%ncolumnref)
    NULLIFY (ircld%iloop)
    NULLIFY (ircld%icount)
    NULLIFY (ircld%icouncol)
    NULLIFY (ircld%icount1)
    NULLIFY (ircld%xcolclr)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_ircld
