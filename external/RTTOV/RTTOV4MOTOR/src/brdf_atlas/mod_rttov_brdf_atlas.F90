! Description:
!> @file
!!   Global module for RTTOV BRDF atlas
!
!> @brief
!!   Global module for RTTOV BRDF atlas
!!
!! @details
!!   This module defines the rttov_brdf_atlas_data type
!!   which contains data for a loaded BRDF atlas.
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
MODULE mod_rttov_brdf_atlas

  USE parkind1, ONLY : jpim, jplm

  USE mod_brdf_atlas, ONLY : brdf_atlas_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: brdf_atlas_id, rttov_brdf_atlas_data

  INTEGER(jpim), PARAMETER :: brdf_atlas_id = 1  !< ID of BRDF atlas (currently there is only one)

  !> Structure to hold BRDF atlas data for one month (and optionally for one
  !! particular instrument)
  TYPE rttov_brdf_atlas_data
    LOGICAL(jplm) :: init = .FALSE.
    INTEGER(jpim) :: atlas_id = 1
    TYPE(brdf_atlas_data) :: brdf_atlas
  ENDTYPE rttov_brdf_atlas_data

END MODULE  mod_rttov_brdf_atlas
