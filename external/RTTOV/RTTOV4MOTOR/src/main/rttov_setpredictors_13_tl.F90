! Description:
!> @file
!!   TL of v13 predictor calculation.
!!
!> @brief
!!   TL of v13 predictor calculation.
!!
!! @details
!!   This subroutine operates on coefficient layers/levels.
!!
!! @param[in]     opts            RTTOV options
!! @param[in]     prof            profiles on coefficient levels
!! @param[in]     coef            rttov_coef structure
!! @param[in]     aux             coef level auxiliary profile data structure
!! @param[in]     aux_tl          coef level auxiliary profile data perturbations
!! @param[in]     predictors      calculated predictors
!! @param[in,out] predictors_tl   calculated predictor perturbations
!! @param[in]     raytracing      RTTOV raytracing structure
!! @param[in]     raytracing_tl   raytracing structure perturbations
!! @param[in]     raypath         integer indicating local path: 1 => pathsat, 2 => patheff
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
!    Copyright 2019, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_setpredictors_13_tl( &
             opts,          &
             prof,          &
             coef,          &
             aux,           &
             aux_tl,        &
             predictors,    &
             predictors_tl, &
             raytracing,    &
             raytracing_tl, &
             raypath)

  USE rttov_types, ONLY :  &
        rttov_coef,             &
        rttov_options,          &
        rttov_profile,          &
        rttov_profile_aux_coef, &
        rttov_path_pred,        &
        rttov_raytracing
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  USE parkind1, ONLY : jpim
  IMPLICIT NONE

  TYPE(rttov_options),          INTENT(IN)    :: opts
  TYPE(rttov_profile),          INTENT(IN)    :: prof(:)
  TYPE(rttov_coef),             INTENT(IN)    :: coef
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux_tl
  TYPE(rttov_path_pred),        INTENT(IN)    :: predictors(:)
  TYPE(rttov_path_pred),        INTENT(INOUT) :: predictors_tl(:)
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing_tl
  INTEGER(jpim),                INTENT(IN)    :: raypath
!INTF_END
 
  INTEGER(jpim) :: lay, i, nlayers, nprofiles

  REAL(jprb), POINTER :: path(:,:), path_sqrt(:,:), path_4rt(:,:)
  REAL(jprb), POINTER :: path_tl(:,:), path_sqrt_tl(:,:), path_4rt_tl(:,:)

  REAL(jprb) :: sec_tr(prof(1)%nlayers, SIZE(prof)), sec_twr(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r(prof(1)%nlayers), sec_wrwrtr_r(prof(1)%nlayers)
  REAL(jprb) :: sec_gr(prof(1)%nlayers), sec_gr_sqrt(prof(1)%nlayers)
  REAL(jprb) :: sec_gw(prof(1)%nlayers), sec_gw_sqrt(prof(1)%nlayers), sec_gw_4rt(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr(prof(1)%nlayers)

  REAL(jprb) :: sec_tr_tl(prof(1)%nlayers, SIZE(prof)), sec_twr_tl(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r_tl(prof(1)%nlayers), sec_wrwrtr_r_tl(prof(1)%nlayers)
  REAL(jprb) :: sec_gr_tl(prof(1)%nlayers), sec_gr_sqrt_tl(prof(1)%nlayers), sec_gr_4rt_tl(prof(1)%nlayers)
  REAL(jprb) :: sec_gw_tl(prof(1)%nlayers), sec_gw_sqrt_tl(prof(1)%nlayers), sec_gw_4rt_tl(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr_tl(prof(1)%nlayers)

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13_TL', 0_jpim, ZHOOK_HANDLE)

  IF (raypath == 1) THEN
    path => raytracing%pathsat
    path_sqrt => aux%pathsat_sqrt
    path_4rt => aux%pathsat_4rt
    path_tl => raytracing_tl%pathsat
    path_sqrt_tl => aux_tl%pathsat_sqrt
    path_4rt_tl => aux_tl%pathsat_4rt
  ELSEIF (raypath == 2) THEN
    path => raytracing%patheff
    path_sqrt => aux%patheff_sqrt
    path_4rt => aux%patheff_4rt
    path_tl => raytracing_tl%patheff
    path_sqrt_tl => aux_tl%patheff_sqrt
    path_4rt_tl => aux_tl%patheff_4rt
  ENDIF

  nlayers = prof(1)%nlayers
  nprofiles = SIZE(prof)

  ! ---------------------------------------------------------------------------
  ! Mixed gases
  ! ---------------------------------------------------------------------------

  sec_tr = path * aux%tr
  sec_tr_tl = path_tl * aux%tr + path * aux_tl%tr
  
  DO i = 1, nprofiles
    DO lay = 1, nlayers
      predictors_tl(i)%mixedgas(1,lay)  = path_tl(lay,i)
      predictors_tl(i)%mixedgas(2,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
      predictors_tl(i)%mixedgas(3,lay)  = sec_tr_tl(lay,i)
      predictors_tl(i)%mixedgas(4,lay)  = sec_tr_tl(lay,i) * aux%tr(lay,i) + sec_tr(lay,i) * aux_tl%tr(lay,i)
      predictors_tl(i)%mixedgas(5,lay)  = aux_tl%tr(lay,i)
      predictors_tl(i)%mixedgas(6,lay)  = aux_tl%tr2(lay,i)
      predictors_tl(i)%mixedgas(7,lay)  = path_tl(lay,i) * aux%twr(lay,i) + path(lay,i) * aux_tl%twr(lay,i)
      predictors_tl(i)%mixedgas(8,lay)  = sec_tr_tl(lay,i) * aux%tr2(lay,i) + sec_tr(lay,i) * aux_tl%tr2(lay,i)
      predictors_tl(i)%mixedgas(9,lay)  = path_sqrt(lay,i) * (1.5_jprb * path_tl(lay,i) * aux%tr_sqrt(lay,i) + &
                                                              path(lay,i) * aux_tl%tr_sqrt(lay,i))
      predictors_tl(i)%mixedgas(10,lay) = 0._jprb
    ENDDO
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Water vapour
  ! ---------------------------------------------------------------------------

  DO i = 1, nprofiles
    sec_gr_sqrt    = path_sqrt(:,i) * aux%wr_sqrt(:,i)
    sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%wr_sqrt(:,i) + &
                     path_sqrt(:,i) * aux_tl%wr_sqrt(:,i)

    sec_gr_4rt_tl  = path_4rt_tl(:,i) * aux%wr_4rt(:,i) + &
                     path_4rt(:,i) * aux_tl%wr_4rt(:,i)

    sec_gr         = sec_gr_sqrt**2
    sec_gr_tl      = 2._jprb * sec_gr_sqrt * sec_gr_sqrt_tl

    sec_gw_sqrt    = path_sqrt(:,i) * aux%ww_sqrt(:,i)
    sec_gw_sqrt_tl = path_sqrt_tl(:,i) * aux%ww_sqrt(:,i) + &
                     path_sqrt(:,i) * aux_tl%ww_sqrt(:,i)

    sec_gw_4rt     = path_4rt(:,i) * aux%ww_4rt(:,i)
    sec_gw_4rt_tl  = path_4rt_tl(:,i) * aux%ww_4rt(:,i) + &
                     path_4rt(:,i) * aux_tl%ww_4rt(:,i)

    sec_gw         = path(:,i) * aux%ww(:,i)
    sec_gw_tl      = path_tl(:,i) * aux%ww(:,i) + &
                     path(:,i) * aux_tl%ww(:,i)

    DO lay = 1, nlayers
      predictors_tl(i)%watervapour(1,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
      predictors_tl(i)%watervapour(2,lay)  = sec_gw_tl(lay)
      predictors_tl(i)%watervapour(3,lay)  = 2._jprb * sec_gw(lay) * sec_gw_tl(lay)
      predictors_tl(i)%watervapour(4,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
      predictors_tl(i)%watervapour(5,lay)  = sec_gr_sqrt_tl(lay)
      predictors_tl(i)%watervapour(6,lay)  = sec_gr_4rt_tl(lay)
      predictors_tl(i)%watervapour(7,lay)  = sec_gr_tl(lay)
      predictors_tl(i)%watervapour(8,lay)  = sec_gw(lay) * sec_gw_sqrt_tl(lay) + sec_gw_tl(lay) * sec_gw_sqrt(lay)
      predictors_tl(i)%watervapour(9,lay)  = sec_gr(lay) * sec_gr_sqrt_tl(lay) + sec_gr_tl(lay) * sec_gr_sqrt(lay)
      predictors_tl(i)%watervapour(10,lay) = sec_gr(lay) * (sec_gr_sqrt(lay) * aux_tl%dt(lay,i)  + &
                                                            sec_gr_sqrt_tl(lay) * aux%dt(lay,i)) + &
                                             sec_gr_tl(lay) * sec_gr_sqrt(lay) * aux%dt(lay,i)
      predictors_tl(i)%watervapour(11,lay) = sec_gr_sqrt_tl(lay) * aux%dt(lay,i) + sec_gr_sqrt(lay) * aux_tl%dt(lay,i)
      predictors_tl(i)%watervapour(12,lay) = sec_gw(lay) * sec_gw_4rt_tl(lay) + sec_gw_tl(lay) * sec_gw_4rt(lay)
      predictors_tl(i)%watervapour(13,lay) = sec_gr_tl(lay) * aux%wrw_r(lay,i) + sec_gr(lay) * aux_tl%wrw_r(lay,i)
      predictors_tl(i)%watervapour(14,lay) = sec_gr_sqrt_tl(lay) * aux%wrwr_r(lay,i) + &
                                             sec_gr_sqrt(lay) * aux_tl%wrwr_r(lay,i)
      predictors_tl(i)%watervapour(15,lay) = path_sqrt_tl(lay,i) * sec_gw_sqrt(lay) + path_sqrt(lay,i) * sec_gw_sqrt_tl(lay)
    ENDDO

    IF (coef%nwvcont > 0) THEN
      sec_wrtr_r    = sec_gr * aux%tr_r(:,i)
      sec_wrtr_r_tl = sec_gr_tl * aux%tr_r(:,i) + &
                      sec_gr * aux_tl%tr_r(:,i)

      sec_wrwrtr_r    = sec_wrtr_r * aux%wr(:,i)
      sec_wrwrtr_r_tl = sec_wrtr_r * aux_tl%wr(:,i) + &
                        sec_wrtr_r_tl * aux%wr(:,i)

      DO lay = 1, nlayers
        predictors_tl(i)%wvcont(1,lay) = sec_wrwrtr_r_tl(lay)
        predictors_tl(i)%wvcont(2,lay) = sec_wrtr_r_tl(lay)
        predictors_tl(i)%wvcont(3,lay) = aux%tr_r(lay,i)**2 * &
                                         (sec_wrwrtr_r(lay) * 3._jprb * aux_tl%tr_r(lay,i) + &
                                          sec_wrwrtr_r_tl(lay) * aux%tr_r(lay,i))
        predictors_tl(i)%wvcont(4,lay) = sec_wrtr_r_tl(lay) * aux%tr_r(lay,i) + &
                                         sec_wrtr_r(lay) * aux_tl%tr_r(lay,i)
      ENDDO
    ENDIF
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Ozone
  ! ---------------------------------------------------------------------------

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nprofiles
        sec_gr_sqrt    = path_sqrt(:,i) * aux%or_sqrt(:,i)
        sec_gr         = sec_gr_sqrt**2
        sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%or_sqrt(:,i) + &
                         path_sqrt(:,i) * aux_tl%or_sqrt(:,i)
        sec_gr_tl      = 2._jprb * sec_gr_sqrt * sec_gr_sqrt_tl

        sec_gw_sqrt    = path_sqrt(:,i) * aux%ow_sqrt(:,i)
        sec_gw         = sec_gw_sqrt**2
        sec_gw_sqrt_tl = path_sqrt_tl(:,i) * aux%ow_sqrt(:,i) + &
                         path_sqrt(:,i) * aux_tl%ow_sqrt(:,i)
        sec_gw_tl      = 2._jprb * sec_gw_sqrt * sec_gw_sqrt_tl
        sec_gw_4rt     = path_4rt(:,i) * aux%ow_4rt(:,i)

        DO lay = 1, nlayers
          predictors_tl(i)%ozone(1,lay)  = sec_gr_tl(lay)
          predictors_tl(i)%ozone(2,lay)  = sec_gr_sqrt_tl(lay)
          predictors_tl(i)%ozone(3,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
          predictors_tl(i)%ozone(4,lay)  = sec_gr_tl(lay) * aux%ow_r(lay,i) + sec_gr(lay) * aux_tl%ow_r(lay,i)
          predictors_tl(i)%ozone(5,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
          predictors_tl(i)%ozone(6,lay)  = sec_gr_tl(lay) * aux%or(lay,i) * aux%ow(lay,i) + &
                                           sec_gr(lay) * (aux_tl%or(lay,i) * aux%ow(lay,i) + &
                                                          aux%or(lay,i) * aux_tl%ow(lay,i))
          predictors_tl(i)%ozone(7,lay)  = sec_gr_sqrt_tl(lay) * aux%or(lay,i) * aux%ow_r(lay,i) + &
                                           sec_gr_sqrt(lay) * (aux_tl%or(lay,i) * aux%ow_r(lay,i) + &
                                                               aux%or(lay,i) * aux_tl%ow_r(lay,i))
          predictors_tl(i)%ozone(8,lay)  = sec_gr_tl(lay) * aux%ow(lay,i) + sec_gr(lay) * aux_tl%ow(lay,i)
          predictors_tl(i)%ozone(9,lay)  = 1.75_jprb * sec_gw_sqrt(lay) * sec_gw_4rt(lay) * sec_gw_tl(lay)
          predictors_tl(i)%ozone(10,lay) = sec_gr(lay) * sec_gw_sqrt_tl(lay) + sec_gr_tl(lay) * sec_gw_sqrt(lay)
          predictors_tl(i)%ozone(11,lay) = 2._jprb * sec_gw(lay) * sec_gw_tl(lay)
          predictors_tl(i)%ozone(12,lay) = aux%ow(lay,i) * (2._jprb * aux_tl%ow(lay,i) * &
                                           path_sqrt(lay,i) * aux%dt(lay,i) + &
                                           aux%ow(lay,i) * (path_sqrt(lay,i) * aux_tl%dt(lay,i) + &
                                                            path_sqrt_tl(lay,i) * aux%dt(lay,i)))
          predictors_tl(i)%ozone(13,lay) = sec_gw_tl(lay)
        ENDDO
      ENDDO
    ELSE ! no user ozone data
      DO i = 1, nprofiles   
        DO lay = 1, nlayers
          predictors_tl(i)%ozone(1,lay)  = path_tl(lay,i)
          predictors_tl(i)%ozone(2,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%ozone(3,lay)  = path_tl(lay,i) * aux%dt(lay,i) + path(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%ozone(4,lay)  = path_tl(lay,i)
          predictors_tl(i)%ozone(5,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ozone(6,lay)  = path_tl(lay,i)
          predictors_tl(i)%ozone(7,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%ozone(8,lay)  = path_tl(lay,i)
          predictors_tl(i)%ozone(9,lay)  = 1.75_jprb * path_sqrt(lay,i) * path_4rt(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ozone(10,lay) = 1.5_jprb * path_sqrt(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ozone(11,lay) = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ozone(12,lay) = path_sqrt_tl(lay,i) * aux%dt(lay,i) + path_sqrt(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%ozone(13,lay) = path_tl(lay,i)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO i = 1, nprofiles
    ! -----------------------------------------------------------------------
    ! CO2
    ! -----------------------------------------------------------------------

    IF (coef%nco2 > 0) THEN
      sec_twr    = path(:,i) * aux%twr(:,i)
      sec_twr_tl = path_tl(:,i) * aux%twr(:,i) + path(:,i) * aux_tl%twr(:,i)

      IF (opts%rt_all%co2_data) THEN
        DO lay = 1, nlayers
          predictors_tl(i)%co2(1,lay)  = path_tl(lay,i) * aux%co2r(lay,i) + path(lay,i) * aux_tl%co2r(lay,i)
          predictors_tl(i)%co2(2,lay)  = aux_tl%tr2(lay,i)
          predictors_tl(i)%co2(3,lay)  = sec_tr_tl(lay,i)
          predictors_tl(i)%co2(4,lay)  = sec_tr_tl(lay,i) * aux%tr(lay,i) + sec_tr(lay,i) * aux_tl%tr(lay,i)
          predictors_tl(i)%co2(5,lay)  = aux_tl%tr(lay,i)
          predictors_tl(i)%co2(6,lay)  = sec_twr_tl(lay)
          predictors_tl(i)%co2(7,lay)  = 2._jprb * path(lay,i) * aux%co2w(lay,i) * &
                                                   (path_tl(lay,i) * aux%co2w(lay,i) + &
                                                    path(lay,i) * aux_tl%co2w(lay,i))
          predictors_tl(i)%co2(8,lay)  = sec_twr_tl(lay) * aux%tr_sqrt(lay,i) + sec_twr(lay) * aux_tl%tr_sqrt(lay,i)
          predictors_tl(i)%co2(9,lay)  = path_sqrt_tl(lay,i) * aux%co2r_sqrt(lay,i) + path_sqrt(lay,i) * aux_tl%co2r_sqrt(lay,i)
          predictors_tl(i)%co2(10,lay) = 3._jprb * aux%tr2(lay,i) * aux_tl%tr(lay,i)
          predictors_tl(i)%co2(11,lay) = aux%tr(lay,i) * (sec_tr(lay,i) * 2._jprb * aux_tl%tr(lay,i) + & 
                                                          sec_tr_tl(lay,i) * aux%tr(lay,i))
          predictors_tl(i)%co2(12,lay) = aux%twr(lay,i)**2 * &
                                         (3._jprb * path_sqrt(lay,i) * aux%tr2(lay,i) * aux_tl%twr(lay,i) + &
                                          path_sqrt(lay,i) * aux%twr(lay,i) * aux_tl%tr2(lay,i) + &
                                          aux%tr2(lay,i) * aux%twr(lay,i) * path_sqrt_tl(lay,i))
          predictors_tl(i)%co2(13,lay) = aux%twr(lay,i) * (aux%tr2(lay,i) * 2._jprb * aux_tl%twr(lay,i) + &
                                                           aux_tl%tr2(lay,i) * aux%twr(lay,i))
          predictors_tl(i)%co2(14,lay) = path_tl(lay,i) * aux%co2w(lay,i) + path(lay,i) * aux_tl%co2w(lay,i)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors_tl(i)%co2(1,lay)  = path_tl(lay,i)
          predictors_tl(i)%co2(2,lay)  = aux_tl%tr2(lay,i)
          predictors_tl(i)%co2(3,lay)  = sec_tr_tl(lay,i)
          predictors_tl(i)%co2(4,lay)  = sec_tr_tl(lay,i) * aux%tr(lay,i) + sec_tr(lay,i) * aux_tl%tr(lay,i)
          predictors_tl(i)%co2(5,lay)  = aux_tl%tr(lay,i)
          predictors_tl(i)%co2(6,lay)  = sec_twr_tl(lay)
          predictors_tl(i)%co2(7,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%co2(8,lay)  = sec_twr_tl(lay) * aux%tr_sqrt(lay,i) + sec_twr(lay) * aux_tl%tr_sqrt(lay,i)
          predictors_tl(i)%co2(9,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%co2(10,lay) = 3._jprb * aux%tr2(lay,i) * aux_tl%tr(lay,i)
          predictors_tl(i)%co2(11,lay) = aux%tr(lay,i) * (sec_tr(lay,i) * 2._jprb * aux_tl%tr(lay,i) + & 
                                                          sec_tr_tl(lay,i) * aux%tr(lay,i))
          predictors_tl(i)%co2(12,lay) = aux%twr(lay,i)**2 * &
                                         (3._jprb * path_sqrt(lay,i) * aux%tr2(lay,i) * aux_tl%twr(lay,i) + &
                                          path_sqrt(lay,i) * aux%twr(lay,i) * aux_tl%tr2(lay,i) + &
                                          aux%tr2(lay,i) * aux%twr(lay,i) * path_sqrt_tl(lay,i))
          predictors_tl(i)%co2(13,lay) = aux%twr(lay,i) * (aux%tr2(lay,i) * 2._jprb * aux_tl%twr(lay,i) + &
                                                           aux_tl%tr2(lay,i) * aux%twr(lay,i))
          predictors_tl(i)%co2(14,lay) = path_tl(lay,i)
        ENDDO
      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------
    ! N2O
    ! ---------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN

      sec_gwr          = path(:,i) * aux%n2owr(:,i)
      sec_gwr_tl       = path_tl(:,i) * aux%n2owr(:,i) + path(:,i) * aux_tl%n2owr(:,i)

      IF (opts%rt_all%n2o_data) THEN
        sec_gr_4rt_tl  = path_4rt_tl(:,i) * aux%n2or_4rt(:,i) + path_4rt(:,i) * aux_tl%n2or_4rt(:,i)

        sec_gr_sqrt    = path_sqrt(:,i) * aux%n2or_sqrt(:,i)
        sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%n2or_sqrt(:,i) +  path_sqrt(:,i) * aux_tl%n2or_sqrt(:,i)

        sec_gr         = sec_gr_sqrt**2
        sec_gr_tl      = 2._jprb * sec_gr_sqrt * sec_gr_sqrt_tl

        DO lay = 1, nlayers
          predictors_tl(i)%n2o(1,lay)  = sec_gr_tl(lay)
          predictors_tl(i)%n2o(2,lay)  = sec_gr_sqrt_tl(lay)
          predictors_tl(i)%n2o(3,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
          predictors_tl(i)%n2o(4,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
          predictors_tl(i)%n2o(5,lay)  = aux_tl%n2or(lay,i) * aux%dt(lay,i) + aux%n2or(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%n2o(6,lay)  = sec_gr_4rt_tl(lay)
          predictors_tl(i)%n2o(7,lay)  = path_tl(lay,i) * aux%n2ow(lay,i) + path(lay,i) * aux_tl%n2ow(lay,i)
          predictors_tl(i)%n2o(8,lay)  = sec_gwr_tl(lay)
          predictors_tl(i)%n2o(9,lay)  = sec_gr_sqrt(lay) * (aux%n2or(lay,i) * aux_tl%n2ow_r(lay,i) + &
                                                             aux_tl%n2or(lay,i) * aux%n2ow_r(lay,i)) + &
                                         sec_gr_sqrt_tl(lay) * aux%n2or(lay,i) * aux%n2ow_r(lay,i)
          predictors_tl(i)%n2o(10,lay) = 2._jprb * sec_gwr(lay) * sec_gwr_tl(lay)
          predictors_tl(i)%n2o(11,lay) = 3._jprb * sec_gwr(lay)**2 * sec_gwr_tl(lay)
          predictors_tl(i)%n2o(12,lay) = path(lay,i) * (sec_gwr(lay) * aux_tl%dt(lay,i) + &
                                                        sec_gwr_tl(lay) * aux%dt(lay,i)) + &
                                         path_tl(lay,i) * sec_gwr(lay) * aux%dt(lay,i)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors_tl(i)%n2o(1,lay)  = path_tl(lay,i)
          predictors_tl(i)%n2o(2,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%n2o(3,lay)  = path_tl(lay,i) * aux%dt(lay,i) + path(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%n2o(4,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%n2o(5,lay)  = aux_tl%dt(lay,i)
          predictors_tl(i)%n2o(6,lay)  = path_4rt_tl(lay,i)
          predictors_tl(i)%n2o(7,lay)  = path_tl(lay,i)
          predictors_tl(i)%n2o(8,lay)  = sec_gwr_tl(lay)
          predictors_tl(i)%n2o(9,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%n2o(10,lay) = 2._jprb * sec_gwr(lay) * sec_gwr_tl(lay)
          predictors_tl(i)%n2o(11,lay) = 3._jprb * sec_gwr(lay)**2 * sec_gwr_tl(lay)
          predictors_tl(i)%n2o(12,lay) = path(lay,i) * (sec_gwr(lay) * aux_tl%dt(lay,i) + &
                                                        sec_gwr_tl(lay) * aux%dt(lay,i)) + &
                                         path_tl(lay,i) * sec_gwr(lay) * aux%dt(lay,i)
        ENDDO
      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------
    ! CO
    ! ---------------------------------------------------------------------

    IF (coef%nco > 0) THEN

      sec_gwr          = path(:,i) * aux%cowr(:,i)
      sec_gwr_tl       = path_tl(:,i) * aux%cowr(:,i) + path(:,i) * aux_tl%cowr(:,i)

      IF (opts%rt_all%co_data) THEN
        sec_gr_4rt_tl  = path_4rt_tl(:,i) * aux%cor_4rt(:,i) + path_4rt(:,i) * aux_tl%cor_4rt(:,i)

        sec_gr_sqrt    = path_sqrt(:,i) * aux%cor_sqrt(:,i)
        sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%cor_sqrt(:,i) +  path_sqrt(:,i) * aux_tl%cor_sqrt(:,i)

        sec_gr         = sec_gr_sqrt**2
        sec_gr_tl      = 2._jprb * sec_gr_sqrt * sec_gr_sqrt_tl
        sec_gw         = path(:,i) * aux%cow(:,i)
        sec_gw_tl      = path_tl(:,i) * aux%cow(:,i) + path(:,i) * aux_tl%cow(:,i)

        DO lay = 1, nlayers
          predictors_tl(i)%co(1,lay)  = sec_gr_tl(lay)
          predictors_tl(i)%co(2,lay)  = sec_gr_sqrt_tl(lay)
          predictors_tl(i)%co(3,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
          predictors_tl(i)%co(4,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
          predictors_tl(i)%co(5,lay)  = sec_gr_sqrt_tl(lay) * aux%dt(lay,i) + sec_gr_sqrt(lay) * aux_tl%dt(lay,i)
          predictors_tl(i)%co(6,lay)  = sec_gr_4rt_tl(lay)
          predictors_tl(i)%co(7,lay)  = sec_gr_tl(lay) * aux%dtabsdt(lay,i) + sec_gr(lay) * aux_tl%dtabsdt(lay,i)
          predictors_tl(i)%co(8,lay)  = sec_gr_tl(lay) * aux%corw_r(lay,i) + sec_gr(lay) * aux_tl%corw_r(lay,i)
          predictors_tl(i)%co(9,lay)  = sec_gr_sqrt_tl(lay) * aux%corw_r(lay,i) + sec_gr_sqrt(lay) * aux_tl%corw_r(lay,i)
          predictors_tl(i)%co(10,lay) = sec_gr_tl(lay) * aux%corw_rsqrt(lay,i) + sec_gr(lay) * aux_tl%corw_rsqrt(lay,i)
          predictors_tl(i)%co(12,lay) = path_4rt_tl(lay,i) * aux%cowr_4rt(lay,i) + path_4rt(lay,i) * aux_tl%cowr_4rt(lay,i)
          predictors_tl(i)%co(13,lay) = sec_gr_tl(lay) * sec_gw(lay) + sec_gr(lay) * sec_gw_tl(lay)
          predictors_tl(i)%co(14,lay) = sec_gw_tl(lay)
          predictors_tl(i)%co(15,lay) = sec_gwr_tl(lay)
          predictors_tl(i)%co(16,lay) = 2._jprb * sec_gw(lay) * sec_gw_tl(lay)
        ENDDO
        predictors_tl(i)%co(11,:) = 0.4_jprb * predictors(i)%co(11,:) * sec_gw_tl / sec_gw
      ELSE
        DO lay = 1, nlayers
          predictors_tl(i)%co(1,lay)  = path_tl(lay,i)
          predictors_tl(i)%co(2,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%co(3,lay)  = path_tl(lay,i) * aux%dt(lay,i) + path(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%co(4,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%co(5,lay)  = path_sqrt_tl(lay,i) * aux%dt(lay,i) + path_sqrt(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%co(6,lay)  = path_4rt_tl(lay,i)
          predictors_tl(i)%co(7,lay)  = path_tl(lay,i) * aux%dtabsdt(lay,i) + path(lay,i) * aux_tl%dtabsdt(lay,i)
          predictors_tl(i)%co(8,lay)  = path_tl(lay,i)
          predictors_tl(i)%co(9,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%co(10,lay) = path_tl(lay,i)
          predictors_tl(i)%co(12,lay) = path_4rt_tl(lay,i) * aux%cowr_4rt(lay,i) + path_4rt(lay,i) * aux_tl%cowr_4rt(lay,i)
          predictors_tl(i)%co(13,lay) = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%co(14,lay) = path_tl(lay,i)
          predictors_tl(i)%co(15,lay) = sec_gwr_tl(lay)
          predictors_tl(i)%co(16,lay) = 2._jprb * path(lay,i) * path_tl(lay,i)
        ENDDO
        predictors_tl(i)%co(11,:) = 0.4_jprb * predictors(i)%co(11,:) * path_tl(:,i) / path(:,i)
      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------
    ! CH4
    ! ---------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN

      IF (opts%rt_all%ch4_data) THEN
        sec_gr_4rt_tl  = path_4rt_tl(:,i) * aux%ch4r_4rt(:,i) + path_4rt(:,i) * aux_tl%ch4r_4rt(:,i)

        sec_gr_sqrt    = path_sqrt(:,i) * aux%ch4r_sqrt(:,i)
        sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%ch4r_sqrt(:,i) +  path_sqrt(:,i) * aux_tl%ch4r_sqrt(:,i)

        sec_gr         = path(:,i) * aux%ch4r(:,i)
        sec_gr_tl      = path_tl(:,i) * aux%ch4r(:,i) + path(:,i) * aux_tl%ch4r(:,i)

        sec_gw         = path(:,i) * aux%ch4w(:,i)
        sec_gw_tl      = path_tl(:,i) * aux%ch4w(:,i) + path(:,i) * aux_tl%ch4w(:,i)

        sec_gw_4rt     = path_4rt(:,i) * aux%ch4w_4rt(:,i)
        sec_gw_4rt_tl  = path_4rt_tl(:,i) * aux%ch4w_4rt(:,i) + path_4rt(:,i) * aux_tl%ch4w_4rt(:,i)
        
        DO lay = 1, nlayers
          predictors_tl(i)%ch4(1,lay)  = sec_gr_tl(lay)
          predictors_tl(i)%ch4(2,lay)  = sec_gr_sqrt_tl(lay)
          predictors_tl(i)%ch4(3,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
          predictors_tl(i)%ch4(4,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
          predictors_tl(i)%ch4(5,lay)  = aux_tl%ch4r(lay,i) * aux%dt(lay,i) + aux%ch4r(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%ch4(6,lay)  = sec_gr_4rt_tl(lay)
          predictors_tl(i)%ch4(7,lay)  = path_tl(lay,i) * aux%ch4wr(lay,i) + path(lay,i) * aux_tl%ch4wr(lay,i)
          predictors_tl(i)%ch4(8,lay)  = aux_tl%ch4wr(lay,i)
          predictors_tl(i)%ch4(9,lay)  = 2._jprb * sec_gw(lay) * sec_gw_tl(lay)
          predictors_tl(i)%ch4(10,lay) = sec_gw_tl(lay)
          predictors_tl(i)%ch4(11,lay) = sec_gr_sqrt_tl(lay) * aux%ch4rw_r(lay,i) + &
                                         sec_gr_sqrt(lay) * aux_tl%ch4rw_r(lay,i)
          predictors_tl(i)%ch4(12,lay) = sec_gw_tl(lay) * sec_gw_4rt(lay) + sec_gw(lay) * sec_gw_4rt_tl(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors_tl(i)%ch4(1,lay)  = path_tl(lay,i)
          predictors_tl(i)%ch4(2,lay)  = path_sqrt_tl(lay,i)
          predictors_tl(i)%ch4(3,lay)  = path_tl(lay,i) * aux%dt(lay,i) + path(lay,i) * aux_tl%dt(lay,i)
          predictors_tl(i)%ch4(4,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ch4(5,lay)  = aux_tl%dt(lay,i)
          predictors_tl(i)%ch4(6,lay)  = path_4rt_tl(lay,i)
          predictors_tl(i)%ch4(7,lay)  = path_tl(lay,i) * aux%ch4wr(lay,i) + path(lay,i) * aux_tl%ch4wr(lay,i)
          predictors_tl(i)%ch4(8,lay)  = aux_tl%ch4wr(lay,i)
          predictors_tl(i)%ch4(9,lay)  = 2._jprb * path(lay,i) * path_tl(lay,i)
          predictors_tl(i)%ch4(10,lay) = path_tl(lay,i)
          predictors_tl(i)%ch4(11,lay) = path_sqrt_tl(lay,i)
          predictors_tl(i)%ch4(12,lay) = path_tl(lay,i) * path_4rt(lay,i) + path(lay,i) * path_4rt_tl(lay,i)
        ENDDO
      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------
    ! SO2
    ! ---------------------------------------------------------------------

    IF (coef%nso2 > 0) THEN
      sec_gr_sqrt    = path_sqrt(:,i) * aux%so2r_sqrt(:,i)
      sec_gr_sqrt_tl = path_sqrt_tl(:,i) * aux%so2r_sqrt(:,i) + &
                       path_sqrt(:,i) * aux_tl%so2r_sqrt(:,i)

      sec_gr_4rt_tl  = path_4rt_tl(:,i) * aux%so2r_4rt(:,i) + &
                       path_4rt(:,i) * aux_tl%so2r_4rt(:,i)

      sec_gr         = sec_gr_sqrt**2
      sec_gr_tl      = 2._jprb * sec_gr_sqrt * sec_gr_sqrt_tl

      sec_gw_sqrt    = path_sqrt(:,i) * aux%so2w_sqrt(:,i)
      sec_gw_sqrt_tl = path_sqrt_tl(:,i) * aux%so2w_sqrt(:,i) + &
                       path_sqrt(:,i) * aux_tl%so2w_sqrt(:,i)

      sec_gw_4rt     = path_4rt(:,i) * aux%so2w_4rt(:,i)
      sec_gw_4rt_tl  = path_4rt_tl(:,i) * aux%so2w_4rt(:,i) + &
                       path_4rt(:,i) * aux_tl%so2w_4rt(:,i)

      sec_gw         = sec_gw_sqrt**2
      sec_gw_tl      = 2._jprb * sec_gw_sqrt * sec_gw_sqrt_tl

      DO lay = 1, nlayers
        predictors_tl(i)%so2(1,lay)  = 2._jprb * sec_gr(lay) * sec_gr_tl(lay)
        predictors_tl(i)%so2(2,lay)  = sec_gw_tl(lay)
        predictors_tl(i)%so2(3,lay)  = 2._jprb * sec_gw(lay) * sec_gw_tl(lay)
        predictors_tl(i)%so2(4,lay)  = sec_gr_tl(lay) * aux%dt(lay,i) + sec_gr(lay) * aux_tl%dt(lay,i)
        predictors_tl(i)%so2(5,lay)  = sec_gr_sqrt_tl(lay)
        predictors_tl(i)%so2(6,lay)  = sec_gr_4rt_tl(lay)
        predictors_tl(i)%so2(7,lay)  = sec_gr_tl(lay)
        predictors_tl(i)%so2(8,lay)  = sec_gr_sqrt_tl(lay) * aux%so2rwr_r(lay,i) + &
                                       sec_gr_sqrt(lay) * aux_tl%so2rwr_r(lay,i)
        predictors_tl(i)%so2(9,lay)  = sec_gw(lay) * sec_gw_sqrt_tl(lay) + sec_gw_tl(lay) * sec_gw_sqrt(lay)
        predictors_tl(i)%so2(10,lay) = sec_gr(lay) * sec_gr_sqrt_tl(lay) + sec_gr_tl(lay) * sec_gr_sqrt(lay)
        predictors_tl(i)%so2(11,lay) = sec_gr(lay) * (sec_gr_sqrt(lay) * aux_tl%dt(lay,i)  + &
                                                      sec_gr_sqrt_tl(lay) * aux%dt(lay,i)) + &
                                       sec_gr_tl(lay) * sec_gr_sqrt(lay) * aux%dt(lay,i)
        predictors_tl(i)%so2(12,lay) = sec_gr_sqrt_tl(lay) * aux%dt(lay,i) + sec_gr_sqrt(lay) * aux_tl%dt(lay,i)
        predictors_tl(i)%so2(13,lay) = sec_gw(lay) * sec_gw_4rt_tl(lay) + sec_gw_tl(lay) * sec_gw_4rt(lay)
        predictors_tl(i)%so2(14,lay) = sec_gr_tl(lay) * aux%so2rw_r(lay,i) + sec_gr(lay) * aux_tl%so2rw_r(lay,i)
        predictors_tl(i)%so2(15,lay) = path_sqrt_tl(lay,i) * sec_gw_sqrt(lay) + path_sqrt(lay,i) * sec_gw_sqrt_tl(lay)
      ENDDO
    ENDIF
  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_13_tl
