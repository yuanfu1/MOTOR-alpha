! Description:
!> @file
!!   Global module for RTTOV emissivity atlases
!
!> @brief
!!   Global module for RTTOV emissivity atlases
!!
!! @details
!!   This module defines the rttov_emis_atlas_data type
!!   which contains data for a loaded emissivity atlas.
!!
!!   A variable of this type should be defined for each
!!   set of atlas data the user requires. The variable
!!   is passed into subroutines to initialise, read and
!!   deallocate the atlases.
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
MODULE mod_rttov_emis_atlas

  USE parkind1, Only : jpim, jplm

  USE mod_uwiremis_atlas,   ONLY : uwiremis_atlas_data
  USE mod_camel_atlas,      ONLY : camel_atlas_data
  USE mod_camel_clim_atlas, ONLY : camel_clim_atlas_data
  USE mod_mwatlas_m2,       ONLY : telsem2_atlas_data
  USE mod_cnrm_mw_atlas,    ONLY : cnrm_mw_atlas_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atlas_type_mw,        &
            atlas_type_ir,        &
            uwiremis_atlas_id,    &
            camel_atlas_id,       &
            camel_clim_atlas_id,  &
            telsem2_atlas_id,     &
            cnrm_mw_atlas_id,     &
            rttov_emis_atlas_data

  INTEGER(jpim), PARAMETER :: atlas_type_mw = 1       !< ID to specify MW atlases
  INTEGER(jpim), PARAMETER :: atlas_type_ir = 2       !< ID to specify IR atlases

  INTEGER(jpim), PARAMETER :: uwiremis_atlas_id   = 1 !< ID for UWIRemis IR atlas
  INTEGER(jpim), PARAMETER :: camel_atlas_id      = 2 !< ID for CAMEL 2007 IR atlas
  INTEGER(jpim), PARAMETER :: camel_clim_atlas_id = 3 !< ID for CAMEL climatology IR atlas
  INTEGER(jpim), PARAMETER :: telsem2_atlas_id    = 1 !< ID for TELSEM2 MW atlas
  INTEGER(jpim), PARAMETER :: cnrm_mw_atlas_id    = 2 !< ID for CNRM MW atlas

  !> Structure to hold data from one emissivity atlas for one month (and, 
  !! in some cases, for one particular instrument)
  TYPE rttov_emis_atlas_data
    LOGICAL(jplm) :: init = .FALSE.
    LOGICAL(jplm) :: is_mw
    INTEGER(jpim) :: atlas_id = 1
    TYPE(uwiremis_atlas_data)   :: uwiremis_atlas
    TYPE(camel_atlas_data)      :: camel_atlas
    TYPE(camel_clim_atlas_data) :: camel_clim_atlas
    TYPE(telsem2_atlas_data)    :: telsem2_atlas
    TYPE(cnrm_mw_atlas_data)    :: cnrm_mw_atlas
  ENDTYPE rttov_emis_atlas_data

END MODULE  mod_rttov_emis_atlas
