! Description:
!> @file
!!   Initialise internal coef-level auxiliary profile structure.
!
!> @brief
!!   Initialise internal coef-level auxiliary profile structure.
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
SUBROUTINE rttov_init_aux_prof_coef(aux_prof)

  USE rttov_types, ONLY : rttov_profile_aux_coef
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_prof
!INTF_END
  aux_prof%s%nearestlev_surf = 0_jpim
  aux_prof%s%pfraction_surf  = 0._jprb
  aux_prof%s%nearestlev_ctp  = 0_jpim
  aux_prof%s%pfraction_ctp   = 0._jprb
  aux_prof%s%cfraction       = 0._jprb

  aux_prof%pathsat_sqrt = 0._jprb
  aux_prof%pathsat_4rt  = 0._jprb

  ! Initialise only those predictor variables which are used in rttov_setpredictors
  IF (ASSOCIATED(aux_prof%t_layer)) THEN
    ! v789/13
    aux_prof%t_layer = 0._jprb
    aux_prof%dt      = 0._jprb
    aux_prof%dtabsdt = 0._jprb
    aux_prof%tr      = 0._jprb
    aux_prof%tr2     = 0._jprb
    aux_prof%tr_r    = 0._jprb
    aux_prof%tr_sqrt = 0._jprb
    aux_prof%wr_sqrt = 0._jprb
    aux_prof%wrw_r   = 0._jprb
    aux_prof%wr_4rt  = 0._jprb
    aux_prof%wr      = 0._jprb
    aux_prof%ww      = 0._jprb
    aux_prof%co2_cm  = 0._jprb

    ! v78
    IF (ASSOCIATED(aux_prof%tw)) THEN
      aux_prof%tw      = 0._jprb
      aux_prof%tw_4rt  = 0._jprb
    ENDIF

    ! v789/13
    IF (ASSOCIATED(aux_prof%or)) THEN
      aux_prof%dto     = 0._jprb
      aux_prof%or      = 0._jprb
      aux_prof%ow      = 0._jprb
      aux_prof%ow_r    = 0._jprb
      aux_prof%or_sqrt = 0._jprb
      aux_prof%ow_sqrt = 0._jprb
    ENDIF

    ! v89/13
    IF (ASSOCIATED(aux_prof%wwr)) THEN
      aux_prof%twr     = 0._jprb
      aux_prof%wrwr_r  = 0._jprb

      IF (ASSOCIATED(aux_prof%co2_layer)) THEN
        aux_prof%co2r      = 0._jprb
        aux_prof%co2r_sqrt = 0._jprb
        aux_prof%co2w      = 0._jprb
      ENDIF

      ! v9/13
      IF (ASSOCIATED(aux_prof%ww_4rt)) THEN

        IF (ASSOCIATED(aux_prof%tuw)) aux_prof%tuw = 0._jprb !v9 only

        aux_prof%ww_sqrt = 0._jprb
        aux_prof%ww_4rt  = 0._jprb
        aux_prof%tro     = 0._jprb

        IF (ASSOCIATED(aux_prof%n2o_layer)) THEN
          aux_prof%n2or      = 0._jprb
          aux_prof%n2or_sqrt = 0._jprb
          aux_prof%n2or_4rt  = 0._jprb
          aux_prof%n2ow      = 0._jprb
          aux_prof%n2owr     = 0._jprb
          aux_prof%n2ow_r    = 0._jprb
        ENDIF

        IF (ASSOCIATED(aux_prof%co_layer)) THEN
          aux_prof%cor_sqrt   = 0._jprb
          aux_prof%cor_4rt    = 0._jprb
          aux_prof%corw_rsqrt = 0._jprb
          aux_prof%corw_r     = 0._jprb
          aux_prof%corw_r4rt  = 0._jprb
          aux_prof%cow        = 0._jprb
          aux_prof%cowr       = 0._jprb
          aux_prof%cowr_4rt   = 0._jprb
        ENDIF

        IF (ASSOCIATED(aux_prof%ch4_layer)) THEN
          aux_prof%ch4r      = 0._jprb
          aux_prof%ch4r_sqrt = 0._jprb
          aux_prof%ch4r_4rt  = 0._jprb
          aux_prof%ch4w      = 0._jprb
          aux_prof%ch4w_4rt  = 0._jprb
          aux_prof%ch4wr     = 0._jprb
          aux_prof%ch4rw_r   = 0._jprb
        ENDIF

        IF (ASSOCIATED(aux_prof%so2_layer)) THEN
          aux_prof%so2r      = 0._jprb
          aux_prof%so2r_sqrt = 0._jprb
          aux_prof%so2r_4rt  = 0._jprb
          aux_prof%so2w_sqrt = 0._jprb
          aux_prof%so2w_4rt  = 0._jprb
          aux_prof%so2rw_r   = 0._jprb
          aux_prof%so2rwr_r  = 0._jprb
        ENDIF
      ENDIF
    ENDIF
  ENDIF

END SUBROUTINE rttov_init_aux_prof_coef
