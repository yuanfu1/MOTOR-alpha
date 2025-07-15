! Description:
!> @file
!!   Calculates v13 predictors.
!!
!> @brief
!!   Calculates v13 predictors.
!!
!! @details
!!   Various quantities used in the calculations were precalculated
!!   in rttov_predictor_precalc_13.
!!
!!   This subroutine operates on coefficient layers/levels.
!!
!! @param[in]     opts            RTTOV options
!! @param[in]     prof            profiles on coefficient levels
!! @param[in]     coef            rttov_coef structure
!! @param[in]     aux             coef level auxiliary profile data structure
!! @param[in,out] predictors      calculated predictors
!! @param[in]     raytracing      RTTOV raytracing structure
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
SUBROUTINE rttov_setpredictors_13( &
             opts,       &
             prof,       &
             coef,       &
             aux,        &
             predictors, &
             raytracing, &
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
  TYPE(rttov_path_pred),        INTENT(INOUT) :: predictors(:)
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  INTEGER(jpim),                INTENT(IN)    :: raypath
!INTF_END
 
  INTEGER(jpim) :: lay, i, nlayers, nprofiles

  REAL(jprb), POINTER :: path(:,:), path_sqrt(:,:), path_4rt(:,:)

  REAL(jprb) :: sec_tr(prof(1)%nlayers, SIZE(prof)), sec_twr(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r(prof(1)%nlayers), sec_wrwrtr_r(prof(1)%nlayers)
  REAL(jprb) :: sec_gr(prof(1)%nlayers), sec_gr_sqrt(prof(1)%nlayers), sec_gr_4rt(prof(1)%nlayers)
  REAL(jprb) :: sec_gw(prof(1)%nlayers), sec_gw_sqrt(prof(1)%nlayers), sec_gw_4rt(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr(prof(1)%nlayers)

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13', 0_jpim, ZHOOK_HANDLE)

  IF (raypath == 1) THEN
    path => raytracing%pathsat
    path_sqrt => aux%pathsat_sqrt
    path_4rt => aux%pathsat_4rt
  ELSEIF (raypath == 2) THEN
    path => raytracing%patheff
    path_sqrt => aux%patheff_sqrt
    path_4rt => aux%patheff_4rt
  ENDIF

  nlayers = prof(1)%nlayers
  nprofiles = SIZE(prof)

  ! aux% variables are calculated in rttov_predictor_precalc_13

  ! ---------------------------------------------------------------------------
  ! Mixed gases
  ! ---------------------------------------------------------------------------

  sec_tr = path * aux%tr

  DO i = 1, nprofiles
    DO lay = 1, nlayers
      predictors(i)%mixedgas(1,lay)  = path(lay,i)
      predictors(i)%mixedgas(2,lay)  = path(lay,i)**2
      predictors(i)%mixedgas(3,lay)  = sec_tr(lay,i)
      predictors(i)%mixedgas(4,lay)  = sec_tr(lay,i) * aux%tr(lay,i)
      predictors(i)%mixedgas(5,lay)  = aux%tr(lay,i)
      predictors(i)%mixedgas(6,lay)  = aux%tr2(lay,i)
      predictors(i)%mixedgas(7,lay)  = path(lay,i) * aux%twr(lay,i)
      predictors(i)%mixedgas(8,lay)  = sec_tr(lay,i) * aux%tr2(lay,i)
      predictors(i)%mixedgas(9,lay)  = path(lay,i) * path_sqrt(lay,i) * aux%tr_sqrt(lay,i)
      predictors(i)%mixedgas(10,lay) = 1._jprb
    ENDDO
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Water vapour
  ! ---------------------------------------------------------------------------

  DO i = 1, nprofiles
    sec_gr_4rt  = path_4rt(:,i) * aux%wr_4rt(:,i)
    sec_gr_sqrt = path_sqrt(:,i) * aux%wr_sqrt(:,i)
    sec_gr      = sec_gr_sqrt**2

    sec_gw_4rt  = path_4rt(:,i) * aux%ww_4rt(:,i)
    sec_gw_sqrt = path_sqrt(:,i) * aux%ww_sqrt(:,i)
    sec_gw      = sec_gw_sqrt**2

    DO lay = 1, nlayers
      predictors(i)%watervapour(1,lay)  = sec_gr(lay)**2
      predictors(i)%watervapour(2,lay)  = sec_gw(lay)
      predictors(i)%watervapour(3,lay)  = sec_gw(lay)**2
      predictors(i)%watervapour(4,lay)  = sec_gr(lay) * aux%dt(lay,i)
      predictors(i)%watervapour(5,lay)  = sec_gr_sqrt(lay)
      predictors(i)%watervapour(6,lay)  = sec_gr_4rt(lay)
      predictors(i)%watervapour(7,lay)  = sec_gr(lay)
      predictors(i)%watervapour(8,lay)  = sec_gw(lay) * sec_gw_sqrt(lay)
      predictors(i)%watervapour(9,lay)  = sec_gr(lay) * sec_gr_sqrt(lay)
      predictors(i)%watervapour(10,lay) = sec_gr(lay) * sec_gr_sqrt(lay) * aux%dt(lay,i)
      predictors(i)%watervapour(11,lay) = sec_gr_sqrt(lay) * aux%dt(lay,i)
      predictors(i)%watervapour(12,lay) = sec_gw(lay) * sec_gw_4rt(lay) 
      predictors(i)%watervapour(13,lay) = sec_gr(lay) * aux%wrw_r(lay,i)
      predictors(i)%watervapour(14,lay) = sec_gr_sqrt(lay) * aux%wrwr_r(lay,i)
      predictors(i)%watervapour(15,lay) = path_sqrt(lay,i) * sec_gw_sqrt(lay)
    ENDDO

    IF (coef%nwvcont > 0) THEN
      sec_wrtr_r   = sec_gr * aux%tr_r(:,i)
      sec_wrwrtr_r = sec_wrtr_r * aux%wr(:,i) 

      DO lay = 1, nlayers
        predictors(i)%wvcont(1,lay) = sec_wrwrtr_r(lay)
        predictors(i)%wvcont(2,lay) = sec_wrtr_r(lay)
        predictors(i)%wvcont(3,lay) = sec_wrwrtr_r(lay) * aux%tr_r(lay,i)**3
        predictors(i)%wvcont(4,lay) = sec_wrtr_r(lay) * aux%tr_r(lay,i)
      ENDDO
    ENDIF
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Ozone
  ! ---------------------------------------------------------------------------

  ! If there us no input O3 profile the predictor calculations are simplified

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nprofiles
        sec_gr_sqrt = path_sqrt(:,i) * aux%or_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2
        sec_gw_sqrt = path_sqrt(:,i) * aux%ow_sqrt(:,i)
        sec_gw      = sec_gw_sqrt**2
        sec_gw_4rt  = path_4rt(:,i) * aux%ow_4rt(:,i)

        DO lay = 1, nlayers
          predictors(i)%ozone(1,lay)  = sec_gr(lay)
          predictors(i)%ozone(2,lay)  = sec_gr_sqrt(lay)
          predictors(i)%ozone(3,lay)  = sec_gr(lay) * aux%dt(lay,i)
          predictors(i)%ozone(4,lay)  = sec_gr(lay) * aux%ow_r(lay,i)
          predictors(i)%ozone(5,lay)  = sec_gr(lay)**2
          predictors(i)%ozone(6,lay)  = sec_gr(lay) * aux%or(lay,i) * aux%ow(lay,i)
          predictors(i)%ozone(7,lay)  = sec_gr_sqrt(lay) * aux%or(lay,i) * aux%ow_r(lay,i)
          predictors(i)%ozone(8,lay)  = sec_gr(lay) * aux%ow(lay,i)
          predictors(i)%ozone(9,lay)  = sec_gw(lay) * sec_gw_sqrt(lay) * sec_gw_4rt(lay) ! x**1.75
          predictors(i)%ozone(10,lay) = sec_gr(lay) * sec_gw_sqrt(lay)
          predictors(i)%ozone(11,lay) = sec_gw(lay)**2
          predictors(i)%ozone(12,lay) = path_sqrt(lay,i) * aux%ow(lay,i)**2 * aux%dt(lay,i)
          predictors(i)%ozone(13,lay) = sec_gw(lay)
        ENDDO
      ENDDO
    ELSE ! no user ozone data
      DO i = 1, nprofiles
        DO lay = 1, nlayers
          predictors(i)%ozone(1,lay)  = path(lay,i)
          predictors(i)%ozone(2,lay)  = path_sqrt(lay,i)
          predictors(i)%ozone(3,lay)  = path(lay,i) * aux%dt(lay,i)
          predictors(i)%ozone(4,lay)  = path(lay,i)
          predictors(i)%ozone(5,lay)  = path(lay,i)**2
          predictors(i)%ozone(6,lay)  = path(lay,i)
          predictors(i)%ozone(7,lay)  = path_sqrt(lay,i)
          predictors(i)%ozone(8,lay)  = path(lay,i)
          predictors(i)%ozone(9,lay)  = path(lay,i) * path_sqrt(lay,i) * path_4rt(lay,i)
          predictors(i)%ozone(10,lay) = path(lay,i) * path_sqrt(lay,i)
          predictors(i)%ozone(11,lay) = path(lay,i)**2
          predictors(i)%ozone(12,lay) = path_sqrt(lay,i) * aux%dt(lay,i)
          predictors(i)%ozone(13,lay) = path(lay,i)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! -----------------------------------------------------------------------
  ! CO2
  ! -----------------------------------------------------------------------

  IF (coef%nco2 > 0) THEN
    DO i = 1, nprofiles
      sec_twr = path(:,i) * aux%twr(:,i)
      IF (opts%rt_all%co2_data) THEN
        DO lay = 1, nlayers
          predictors(i)%co2(1,lay)  = path(lay,i) * aux%co2r(lay,i)
          predictors(i)%co2(2,lay)  = aux%tr2(lay,i)
          predictors(i)%co2(3,lay)  = sec_tr(lay,i)
          predictors(i)%co2(4,lay)  = sec_tr(lay,i) * aux%tr(lay,i)
          predictors(i)%co2(5,lay)  = aux%tr(lay,i)
          predictors(i)%co2(6,lay)  = sec_twr(lay)
          predictors(i)%co2(7,lay)  = (path(lay,i) * aux%co2w(lay,i))**2
          predictors(i)%co2(8,lay)  = sec_twr(lay) * aux%tr_sqrt(lay,i)
          predictors(i)%co2(9,lay)  = path_sqrt(lay,i) * aux%co2r_sqrt(lay,i)
          predictors(i)%co2(10,lay) = aux%tr2(lay,i) * aux%tr(lay,i)
          predictors(i)%co2(11,lay) = sec_tr(lay,i) * aux%tr2(lay,i)
          predictors(i)%co2(12,lay) = path_sqrt(lay,i) * aux%tr2(lay,i) * aux%twr(lay,i)**3
          predictors(i)%co2(13,lay) = aux%tr2(lay,i) * aux%twr(lay,i)**2
          predictors(i)%co2(14,lay) = path(lay,i) * aux%co2w(lay,i)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors(i)%co2(1,lay)  = path(lay,i)
          predictors(i)%co2(2,lay)  = aux%tr2(lay,i)
          predictors(i)%co2(3,lay)  = sec_tr(lay,i)
          predictors(i)%co2(4,lay)  = sec_tr(lay,i) * aux%tr(lay,i)
          predictors(i)%co2(5,lay)  = aux%tr(lay,i)
          predictors(i)%co2(6,lay)  = sec_twr(lay)
          predictors(i)%co2(7,lay)  = path(lay,i)**2
          predictors(i)%co2(8,lay)  = sec_twr(lay) * aux%tr_sqrt(lay,i)
          predictors(i)%co2(9,lay)  = path_sqrt(lay,i)
          predictors(i)%co2(10,lay) = aux%tr2(lay,i) * aux%tr(lay,i)
          predictors(i)%co2(11,lay) = sec_tr(lay,i) * aux%tr(lay,i)**2
          predictors(i)%co2(12,lay) = path_sqrt(lay,i) * aux%tr2(lay,i) * aux%twr(lay,i)**3
          predictors(i)%co2(13,lay) = aux%tr2(lay,i) * aux%twr(lay,i)**2
          predictors(i)%co2(14,lay) = path(lay,i)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! -----------------------------------------------------------------------
  ! N2O
  ! -----------------------------------------------------------------------

  IF (coef%nn2o > 0) THEN
    DO i = 1, nprofiles
      sec_gwr = path(:,i) * aux%n2owr(:,i)
      
      IF (opts%rt_all%n2o_data) THEN
        sec_gr_4rt  = path_4rt(:,i) * aux%n2or_4rt(:,i) 
        sec_gr_sqrt = path_sqrt(:,i) * aux%n2or_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2

        DO lay = 1, nlayers
          predictors(i)%n2o(1,lay)  = sec_gr(lay)
          predictors(i)%n2o(2,lay)  = sec_gr_sqrt(lay)
          predictors(i)%n2o(3,lay)  = sec_gr(lay) * aux%dt(lay,i)
          predictors(i)%n2o(4,lay)  = sec_gr(lay)**2
          predictors(i)%n2o(5,lay)  = aux%n2or(lay,i) * aux%dt(lay,i)
          predictors(i)%n2o(6,lay)  = sec_gr_4rt(lay)
          predictors(i)%n2o(7,lay)  = path(lay,i) * aux%n2ow(lay,i)
          predictors(i)%n2o(8,lay)  = sec_gwr(lay)
          predictors(i)%n2o(9,lay)  = sec_gr_sqrt(lay) * aux%n2or(lay,i) * aux%n2ow_r(lay,i)
          predictors(i)%n2o(10,lay) = sec_gwr(lay)**2
          predictors(i)%n2o(11,lay) = sec_gwr(lay)**3
          predictors(i)%n2o(12,lay) = path(lay,i) * sec_gwr(lay) * aux%dt(lay,i)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors(i)%n2o(1,lay)  = path(lay,i)
          predictors(i)%n2o(2,lay)  = path_sqrt(lay,i)
          predictors(i)%n2o(3,lay)  = path(lay,i) * aux%dt(lay,i)
          predictors(i)%n2o(4,lay)  = path(lay,i)**2
          predictors(i)%n2o(5,lay)  = aux%dt(lay,i)
          predictors(i)%n2o(6,lay)  = path_4rt(lay,i)
          predictors(i)%n2o(7,lay)  = path(lay,i)
          predictors(i)%n2o(8,lay)  = sec_gwr(lay)
          predictors(i)%n2o(9,lay)  = path_sqrt(lay,i)
          predictors(i)%n2o(10,lay) = sec_gwr(lay)**2
          predictors(i)%n2o(11,lay) = sec_gwr(lay)**3
          predictors(i)%n2o(12,lay) = path(lay,i) * sec_gwr(lay) * aux%dt(lay,i)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! -----------------------------------------------------------------------
  ! CO
  ! -----------------------------------------------------------------------

  IF (coef%nco > 0) THEN

    DO i = 1, nprofiles
      sec_gwr = path(:,i) * aux%cowr(:,i)

      IF (opts%rt_all%co_data) THEN
        sec_gr_4rt  = path_4rt(:,i) * aux%cor_4rt(:,i)
        sec_gr_sqrt = path_sqrt(:,i) * aux%cor_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2
        sec_gw      = path(:,i) * aux%cow(:,i)

        DO lay = 1, nlayers
          predictors(i)%co(1,lay)  = sec_gr(lay)
          predictors(i)%co(2,lay)  = sec_gr_sqrt(lay)
          predictors(i)%co(3,lay)  = sec_gr(lay) * aux%dt(lay,i)
          predictors(i)%co(4,lay)  = sec_gr(lay)**2
          predictors(i)%co(5,lay)  = sec_gr_sqrt(lay) * aux%dt(lay,i)
          predictors(i)%co(6,lay)  = sec_gr_4rt(lay)
          predictors(i)%co(7,lay)  = sec_gr(lay) * aux%dtabsdt(lay,i)
          predictors(i)%co(8,lay)  = sec_gr(lay) * aux%corw_r(lay,i)
          predictors(i)%co(9,lay)  = sec_gr_sqrt(lay) * aux%corw_r(lay,i)
          predictors(i)%co(10,lay) = sec_gr(lay) * aux%corw_rsqrt(lay,i)
          predictors(i)%co(12,lay) = path_4rt(lay,i) * aux%cowr_4rt(lay,i)
          predictors(i)%co(13,lay) = sec_gr(lay) * sec_gw(lay)
          predictors(i)%co(14,lay) = sec_gw(lay)
          predictors(i)%co(15,lay) = sec_gwr(lay)
          predictors(i)%co(16,lay) = sec_gw(lay)**2
        ENDDO
        predictors(i)%co(11,:) = (sec_gw)**0.4_jprb
      ELSE
        DO lay = 1, nlayers
          predictors(i)%co(1,lay)  = path(lay,i)
          predictors(i)%co(2,lay)  = path_sqrt(lay,i)
          predictors(i)%co(3,lay)  = path(lay,i) * aux%dt(lay,i)
          predictors(i)%co(4,lay)  = path(lay,i)**2
          predictors(i)%co(5,lay)  = path_sqrt(lay,i) * aux%dt(lay,i)
          predictors(i)%co(6,lay)  = path_4rt(lay,i)
          predictors(i)%co(7,lay)  = path(lay,i) * aux%dtabsdt(lay,i)
          predictors(i)%co(8,lay)  = path(lay,i)
          predictors(i)%co(9,lay)  = path_sqrt(lay,i)
          predictors(i)%co(10,lay) = path(lay,i)
          predictors(i)%co(12,lay) = path_4rt(lay,i) * aux%cowr_4rt(lay,i)
          predictors(i)%co(13,lay) = path(lay,i)**2
          predictors(i)%co(14,lay) = path(lay,i)
          predictors(i)%co(15,lay) = sec_gwr(lay)
          predictors(i)%co(16,lay) = path(lay,i)**2
        ENDDO
        predictors(i)%co(11,:) = (path(:,i))**0.4_jprb
      ENDIF
    ENDDO
  ENDIF

  ! -----------------------------------------------------------------------
  ! CH4
  ! -----------------------------------------------------------------------

  IF (coef%nch4 > 0) THEN

    DO i = 1, nprofiles
      IF (opts%rt_all%ch4_data) THEN
        sec_gr_4rt  = path_4rt(:,i) * aux%ch4r_4rt(:,i)
        sec_gr_sqrt = path_sqrt(:,i) * aux%ch4r_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2
        sec_gw      = path(:,i) * aux%ch4w(:,i)
        sec_gw_4rt  = path_4rt(:,i) * aux%ch4w_4rt(:,i)

        DO lay = 1, nlayers
          predictors(i)%ch4(1,lay)  = sec_gr(lay)
          predictors(i)%ch4(2,lay)  = sec_gr_sqrt(lay)
          predictors(i)%ch4(3,lay)  = sec_gr(lay) * aux%dt(lay,i)
          predictors(i)%ch4(4,lay)  = sec_gr(lay)**2
          predictors(i)%ch4(5,lay)  = aux%ch4r(lay,i) * aux%dt(lay,i)
          predictors(i)%ch4(6,lay)  = sec_gr_4rt(lay)
          predictors(i)%ch4(7,lay)  = path(lay,i) * aux%ch4wr(lay,i)
          predictors(i)%ch4(8,lay)  = aux%ch4wr(lay,i)
          predictors(i)%ch4(9,lay)  = sec_gw(lay)**2
          predictors(i)%ch4(10,lay) = sec_gw(lay)
          predictors(i)%ch4(11,lay) = sec_gr_sqrt(lay) * aux%ch4rw_r(lay,i)
          predictors(i)%ch4(12,lay) = sec_gw(lay) * sec_gw_4rt(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          predictors(i)%ch4(1,lay)  = path(lay,i)
          predictors(i)%ch4(2,lay)  = path_sqrt(lay,i)
          predictors(i)%ch4(3,lay)  = path(lay,i) * aux%dt(lay,i)
          predictors(i)%ch4(4,lay)  = path(lay,i)**2
          predictors(i)%ch4(5,lay)  = aux%dt(lay,i)
          predictors(i)%ch4(6,lay)  = path_4rt(lay,i)
          predictors(i)%ch4(7,lay)  = path(lay,i) * aux%ch4wr(lay,i)
          predictors(i)%ch4(8,lay)  = aux%ch4wr(lay,i)
          predictors(i)%ch4(9,lay)  = path(lay,i)**2
          predictors(i)%ch4(10,lay) = path(lay,i)
          predictors(i)%ch4(11,lay) = path_sqrt(lay,i)
          predictors(i)%ch4(12,lay) = path(lay,i) * path_4rt(lay,i)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! -----------------------------------------------------------------------
  ! SO2
  ! -----------------------------------------------------------------------

  IF (coef%nso2 > 0) THEN
    DO i = 1, nprofiles
      sec_gr_4rt  = path_4rt(:,i) * aux%so2r_4rt(:,i)
      sec_gr_sqrt = path_sqrt(:,i) * aux%so2r_sqrt(:,i)
      sec_gr      = sec_gr_sqrt**2
      sec_gw_4rt  = path_4rt(:,i) * aux%so2w_4rt(:,i)
      sec_gw_sqrt = path_sqrt(:,i) * aux%so2w_sqrt(:,i)
      sec_gw      = sec_gw_sqrt**2
    
      DO lay = 1, nlayers
        predictors(i)%so2(1,lay)  = sec_gr(lay)**2
        predictors(i)%so2(2,lay)  = sec_gw(lay)
        predictors(i)%so2(3,lay)  = sec_gw(lay)**2
        predictors(i)%so2(4,lay)  = sec_gr(lay) * aux%dt(lay,i)
        predictors(i)%so2(5,lay)  = sec_gr_sqrt(lay)
        predictors(i)%so2(6,lay)  = sec_gr_4rt(lay)
        predictors(i)%so2(7,lay)  = sec_gr(lay)
        predictors(i)%so2(8,lay)  = sec_gr_sqrt(lay) * aux%so2rwr_r(lay,i)
        predictors(i)%so2(9,lay)  = sec_gw(lay) * sec_gw_sqrt(lay)
        predictors(i)%so2(10,lay) = sec_gr(lay) * sec_gr_sqrt(lay)
        predictors(i)%so2(11,lay) = sec_gr(lay) * sec_gr_sqrt(lay) * aux%dt(lay,i)
        predictors(i)%so2(12,lay) = sec_gr_sqrt(lay) * aux%dt(lay,i)
        predictors(i)%so2(13,lay) = sec_gw(lay) * sec_gw_4rt(lay)
        predictors(i)%so2(14,lay) = sec_gr(lay) * aux%so2rw_r(lay,i)
        predictors(i)%so2(15,lay) = path(lay,i) * aux%so2w_sqrt(lay,i)
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_13
