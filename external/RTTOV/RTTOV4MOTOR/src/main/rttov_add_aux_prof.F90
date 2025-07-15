! Description:
!> @file
!!   Add two arrays of auxiliary profile nearest-level structures.
!
!> @brief
!!   Add two arrays of auxiliary profile nearest-level structures.
!!
!! @details
!!   The nearest-level data are calculated on user levels and coef levels,
!!   but if the interpolator is not used then the coef-level values are
!!   taken directly from those calculated from the user levels. This routine
!!   is called by the AD/K as the counterpart to the rttov_copy_aux_prof
!!   subroutine called by the direct and TL models.
!!
!! @param[in,out] aux_prof_s          output summed auxiliary profile structure
!! @param[in]     aux_prof_s1         first input auxiliary profile structure
!! @param[in]     aux_prof_s2         second input auxiliary profile structure
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
SUBROUTINE rttov_add_aux_prof(aux_prof_s, aux_prof_s1, aux_prof_s2)

  USE rttov_types, ONLY : rttov_profile_aux_s
  IMPLICIT NONE
  TYPE(rttov_profile_aux_s), INTENT(INOUT) :: aux_prof_s(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: aux_prof_s1(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: aux_prof_s2(:)
!INTF_END
  aux_prof_s%pfraction_surf = aux_prof_s1%pfraction_surf + aux_prof_s2%pfraction_surf
  aux_prof_s%pfraction_ctp  = aux_prof_s1%pfraction_ctp  + aux_prof_s2%pfraction_ctp
  aux_prof_s%cfraction      = aux_prof_s1%cfraction      + aux_prof_s2%cfraction
END SUBROUTINE
