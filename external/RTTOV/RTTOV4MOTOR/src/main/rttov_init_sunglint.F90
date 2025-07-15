! Description:
!> @file
!!   Initialise internal sunglint structure.
!
!> @brief
!!   Initialise internal sunglint structure.
!!
!! @details
!!   This subroutine is primarily intended for use by the AD/K models
!!   and as such only the members which are required to be initialised
!!   by these models are set to zero. The arrays for wave spectrum
!!   parameterisations are only allocated in the direct model-only and
!!   so are not initialised here for efficiency.
!!
!! @param[in,out] sunglint    sunglint structure
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
SUBROUTINE rttov_init_sunglint(sunglint)

  USE rttov_types, ONLY : rttov_sunglint
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_sunglint), INTENT(INOUT) :: sunglint
!INTF_END
  sunglint%s%windsp      = 0._jprb
  sunglint%s%wangl       = 0._jprb
  sunglint%s%dazng       = 0._jprb
  sunglint%s%zensat      = 0._jprb
  sunglint%s%zensun      = 0._jprb
  sunglint%s%gamma_sq    = 0._jprb
  sunglint%s%gamma_o     = 0._jprb
  sunglint%s%gamma_p     = 0._jprb
  sunglint%s%csi         = 0._jprb
  sunglint%s%alfa        = 0._jprb
  sunglint%s%omega       = 0._jprb
  sunglint%s%gammax      = 0._jprb
  sunglint%s%q_shad      = 0._jprb
  sunglint%s%a_shad      = 0._jprb
  sunglint%s%b_shad      = 0._jprb
  sunglint%s%lambda_a    = 0._jprb
  sunglint%s%lambda_b    = 0._jprb
  sunglint%s%c_shad      = 0._jprb
  sunglint%s%p_prime     = 0._jprb
  sunglint%s%g_shad      = 0._jprb
  sunglint%s%fac1        = 0._jprb
  sunglint%s%pxy_gammaxy = 0._jprb
  sunglint%s%glint       = 0._jprb
  sunglint%s%x_u         = 0._jprb
  sunglint%s%alfa1       = 0._jprb
  sunglint%s%omega_m     = 0._jprb
  sunglint%s%omega_c     = 0._jprb
END SUBROUTINE
