! Description:
!> @file
!!   Copy an array of auxiliary profile nearest-level structures.
!
!> @brief
!!   Copy an array of auxiliary profile nearest-level structures.
!!
!! @details
!!   The nearest-level data are calculated on user levels and coef levels,
!!   but if the interpolator is not used then the coef-level values are
!!   taken directly from those calculated from the user levels.
!!
!! @param[in,out] aux_prof_s1         copy of auxiliary profile structure
!! @param[in]     aux_prof_s2         input auxiliary profile structure
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
SUBROUTINE rttov_copy_aux_prof(aux_prof_s1, aux_prof_s2)

  USE rttov_types, ONLY : rttov_profile_aux_s
  IMPLICIT NONE
  TYPE(rttov_profile_aux_s), INTENT(INOUT) :: aux_prof_s1(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: aux_prof_s2(:)
!INTF_END
  aux_prof_s1 = aux_prof_s2
END SUBROUTINE
