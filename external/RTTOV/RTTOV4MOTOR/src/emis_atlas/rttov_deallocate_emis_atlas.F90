! Description:
!> @file
!!   Deallocate memory for emissivity atlas.
!
!> @brief
!!   Deallocate memory for emissivity atlas.
!!
!! @details
!!   This subroutine is used to deallocate an emissivity atlas
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
SUBROUTINE rttov_deallocate_emis_atlas(atlas)

  USE mod_rttov_emis_atlas, ONLY : rttov_emis_atlas_data
!INTF_OFF
  USE mod_rttov_emis_atlas, ONLY : &
    uwiremis_atlas_id, camel_atlas_id, camel_clim_atlas_id, &
    telsem2_atlas_id, cnrm_mw_atlas_id

  USE mod_uwiremis_atlas,   ONLY : rttov_uwiremis_close_atlas
  USE mod_camel_atlas,      ONLY : rttov_camel_close_atlas
  USE mod_camel_clim_atlas, ONLY : rttov_camel_clim_close_atlas
  USE mod_mwatlas_m2,       ONLY : rttov_close_telsem2_atlas => rttov_closemw_atlas
  USE mod_cnrm_mw_atlas,    ONLY : rttov_cnrm_close_atlas
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_emis_atlas_data), INTENT(INOUT) :: atlas
!INTF_END

  IF (atlas%init) THEN
    IF (atlas%is_mw) THEN
      ! MW atlas
      IF (atlas%atlas_id == telsem2_atlas_id) CALL rttov_close_telsem2_atlas(atlas%telsem2_atlas)
      IF (atlas%atlas_id == cnrm_mw_atlas_id) CALL rttov_cnrm_close_atlas(atlas%cnrm_mw_atlas)
    ELSE
      ! IR atlas
      IF (atlas%atlas_id == uwiremis_atlas_id) CALL rttov_uwiremis_close_atlas(atlas%uwiremis_atlas)
      IF (atlas%atlas_id == camel_atlas_id) CALL rttov_camel_close_atlas(atlas%camel_atlas)
      IF (atlas%atlas_id == camel_clim_atlas_id) CALL rttov_camel_clim_close_atlas(atlas%camel_clim_atlas)
    ENDIF
    atlas%init = .FALSE.
  ENDIF

END SUBROUTINE rttov_deallocate_emis_atlas
