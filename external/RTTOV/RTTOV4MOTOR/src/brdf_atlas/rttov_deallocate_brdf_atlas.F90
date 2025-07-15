! Description:
!> @file
!!   Deallocate memory for BRDF atlas.
!
!> @brief
!!   Deallocate memory for BRDF atlas.
!!
!! @details
!!   This subroutine is used to deallocate an BRDF atlas
!!   data structure. This subroutine should be called for every
!!   atlas data structure that was initialised.
!!
!! @param[in,out]  atlas  Atlas data structure to deallocate
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
SUBROUTINE rttov_deallocate_brdf_atlas(atlas)

  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data
!INTF_OFF
  USE mod_rttov_brdf_atlas, ONLY : brdf_atlas_id

  USE mod_brdf_atlas, ONLY : rttov_visnirbrdf_close_atlas
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_brdf_atlas_data), INTENT(INOUT) :: atlas
!INTF_END

  IF (atlas%init) THEN
    IF (atlas%atlas_id == brdf_atlas_id) CALL rttov_visnirbrdf_close_atlas(atlas%brdf_atlas)
    atlas%init = .FALSE.
  ENDIF


END SUBROUTINE rttov_deallocate_brdf_atlas
