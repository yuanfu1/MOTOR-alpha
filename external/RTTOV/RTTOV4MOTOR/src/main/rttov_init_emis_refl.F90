! Description:
!> @file
!!   Initialise an emissivity and/or a reflectance structure.
!
!> @brief
!!   Initialise an emissivity and/or a reflectance structure.
!!
!! @param[in,out]  emis   Surface emissivity structure to initialise, optional
!! @param[in,out]  refl   Surface reflectance structure to initialise, optional
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
SUBROUTINE rttov_init_emis_refl(emis, refl)

  USE rttov_types, ONLY : rttov_emissivity, rttov_reflectance
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_emissivity),  INTENT(INOUT), OPTIONAL :: emis(:)
  TYPE(rttov_reflectance), INTENT(INOUT), OPTIONAL :: refl(:)
!INTF_END

  IF (PRESENT(emis)) THEN
    emis%emis_in          = 0._jprb
    emis%emis_out         = 0._jprb
    emis%specularity      = 0._jprb
  ENDIF

  IF (PRESENT(refl)) THEN
    refl%refl_in          = 0._jprb
    refl%refl_out         = 0._jprb
    refl%diffuse_refl_in  = 0._jprb
    refl%diffuse_refl_out = 0._jprb
    refl%refl_cloud_top   = 0._jprb
  ENDIF
END SUBROUTINE rttov_init_emis_refl
