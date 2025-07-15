! Description:
!> @file
!!   Initialise internal ircld structure.
!
!> @brief
!!   Initialise internal ircld structure.
!!
!! @param[in,out] ircld          ircld structure to initialise
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
SUBROUTINE rttov_init_ircld(ircld)

  USE rttov_types, ONLY : rttov_ircld
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_ircld), INTENT(INOUT) :: ircld
!INTF_END
  ircld%xcolclr = 0._jprb
  IF (ASSOCIATED(ircld%icldarr)   ) ircld%icldarr    = 0_jpim
  IF (ASSOCIATED(ircld%xcolref1)  ) ircld%xcolref1   = 0._jprb
  IF (ASSOCIATED(ircld%xcolref2)  ) ircld%xcolref2   = 0._jprb
  IF (ASSOCIATED(ircld%indexcol)  ) ircld%indexcol   = 0_jpim
  IF (ASSOCIATED(ircld%icount1ref)) ircld%icount1ref = 0_jpim
  IF (ASSOCIATED(ircld%iloopin)   ) ircld%iloopin    = 0_jpim
  IF (ASSOCIATED(ircld%iflag)     ) ircld%iflag      = 0_jpim
  IF (ASSOCIATED(ircld%xcol)      ) ircld%xcol       = 0._jprb
  IF (ASSOCIATED(ircld%xcolminref)) ircld%xcolminref = 0._jprb
  IF (ASSOCIATED(ircld%xcolref)   ) ircld%xcolref    = 0._jprb
  IF (ASSOCIATED(ircld%cldcfr)    ) ircld%cldcfr     = 0._jprb
  IF (ASSOCIATED(ircld%maxcov)    ) ircld%maxcov     = 0._jprb
  IF (ASSOCIATED(ircld%xcolmax)   ) ircld%xcolmax    = 0._jprb
  IF (ASSOCIATED(ircld%xcolmin)   ) ircld%xcolmin    = 0._jprb
  IF (ASSOCIATED(ircld%a)         ) ircld%a          = 0._jprb
  IF (ASSOCIATED(ircld%ntotref)   ) ircld%ntotref    = 0._jprb
  IF (ASSOCIATED(ircld%flag)      ) ircld%flag       = .FALSE.
END SUBROUTINE rttov_init_ircld
