! Description:
!> @file
!!   Initialise a radiance structure.
!
!> @brief
!!   Initialise a radiance structure and optionally also a
!!   secondary radiance structure.
!!
!! @param[in,out]  rad   Radiance structure to initialise
!! @param[in,out]  rad2  Secondary radiance structure to initialise, optional
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
SUBROUTINE rttov_init_rad(rad, rad2)

  USE rttov_types, ONLY : rttov_radiance, rttov_radiance2
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_radiance),  INTENT(INOUT)           :: rad
  TYPE(rttov_radiance2), INTENT(INOUT), OPTIONAL :: rad2
!INTF_END

  rad%clear            = 0._jprb
  rad%total            = 0._jprb
  rad%bt_clear         = 0._jprb
  rad%bt               = 0._jprb
  rad%refl_clear       = 0._jprb
  rad%refl             = 0._jprb
  rad%cloudy           = 0._jprb
  rad%overcast         = 0._jprb
  rad%quality          = 0_jpim
  rad%geometric_height = 0._jprb

  IF (PRESENT(rad2)) THEN
    rad2%upclear     = 0._jprb
    rad2%dnclear     = 0._jprb
    rad2%refldnclear = 0._jprb
    rad2%up          = 0._jprb
    rad2%down        = 0._jprb
    rad2%surf        = 0._jprb
  ENDIF
END SUBROUTINE
