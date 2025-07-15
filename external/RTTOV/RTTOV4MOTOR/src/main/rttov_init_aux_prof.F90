! Description:
!> @file
!!   Initialise internal auxiliary profile structure.
!
!> @brief
!!   Initialise internal auxiliary profile structure.
!!
!! @param[in,out] aux_prof          auxiliary profile structure to initialise
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
SUBROUTINE rttov_init_aux_prof(aux_prof)

  USE rttov_types, ONLY : rttov_profile_aux
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_prof
!INTF_END
  aux_prof%s%nearestlev_surf = 0_jpim
  aux_prof%s%pfraction_surf  = 0._jprb
  aux_prof%s%nearestlev_ctp  = 0_jpim
  aux_prof%s%pfraction_ctp   = 0._jprb
  aux_prof%s%cfraction       = 0._jprb
  aux_prof%ros_z2            = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  IF (ASSOCIATED(aux_prof%clw))               aux_prof%clw               = 0._jprb
  IF (ASSOCIATED(aux_prof%debye_prof))        aux_prof%debye_prof        = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_eps_s))         aux_prof%ros_eps_s         = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_dr))            aux_prof%ros_dr            = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_gammar))        aux_prof%ros_gammar        = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_db))            aux_prof%ros_db            = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_nu1))           aux_prof%ros_nu1           = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_z1))            aux_prof%ros_z1            = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  IF (ASSOCIATED(aux_prof%ros_log_abs_z1_sq)) aux_prof%ros_log_abs_z1_sq = 0._jprb
  IF (ASSOCIATED(aux_prof%ros_log_z1star_z2)) aux_prof%ros_log_z1star_z2 = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  IF (ASSOCIATED(aux_prof%ros_div1))          aux_prof%ros_div1          = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  IF (ASSOCIATED(aux_prof%ros_div2))          aux_prof%ros_div2          = CMPLX(0._jprb, 0._jprb, KIND=jprb)
  IF (ASSOCIATED(aux_prof%tkc_eps_s))         aux_prof%tkc_eps_s         = 0._jprb
  IF (ASSOCIATED(aux_prof%tkc_delta_1))       aux_prof%tkc_delta_1       = 0._jprb
  IF (ASSOCIATED(aux_prof%tkc_delta_2))       aux_prof%tkc_delta_2       = 0._jprb
  IF (ASSOCIATED(aux_prof%tkc_tau_1))         aux_prof%tkc_tau_1         = 0._jprb
  IF (ASSOCIATED(aux_prof%tkc_tau_2))         aux_prof%tkc_tau_2         = 0._jprb
  IF (ASSOCIATED(aux_prof%clw_dg))            aux_prof%clw_dg            = 0._jprb
  IF (ASSOCIATED(aux_prof%clw_dg_ref))        aux_prof%clw_dg_ref        = 0._jprb
  IF (ASSOCIATED(aux_prof%ice_dg))            aux_prof%ice_dg            = 0._jprb
  IF (ASSOCIATED(aux_prof%ice_dg_ref))        aux_prof%ice_dg_ref        = 0._jprb
  IF (ASSOCIATED(aux_prof%fac1_ice_dg))       aux_prof%fac1_ice_dg       = 0._jprb
  IF (ASSOCIATED(aux_prof%relhum))            aux_prof%relhum            = 0._jprb
  IF (ASSOCIATED(aux_prof%tave))              aux_prof%tave              = 0._jprb
  IF (ASSOCIATED(aux_prof%wmixave))           aux_prof%wmixave           = 0._jprb
  IF (ASSOCIATED(aux_prof%xpresave))          aux_prof%xpresave          = 0._jprb
  IF (ASSOCIATED(aux_prof%ppv))               aux_prof%ppv               = 0._jprb
  IF (ASSOCIATED(aux_prof%esw))               aux_prof%esw               = 0._jprb
  IF (ASSOCIATED(aux_prof%esi))               aux_prof%esi               = 0._jprb
END SUBROUTINE 
