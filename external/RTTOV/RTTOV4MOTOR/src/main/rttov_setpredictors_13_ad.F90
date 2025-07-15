! Description:
!> @file
!!   AD of v13 predictor calculation.
!!
!> @brief
!!   AD of v13 predictor calculation.
!!
!! @details
!!   This subroutine operates on coefficient layers/levels.
!!
!! @param[in]     opts            RTTOV options
!! @param[in]     prof            profiles on coefficient levels
!! @param[in]     coef            rttov_coef structure
!! @param[in]     aux             coef level auxiliary profile data structure
!! @param[in,out] aux_ad          coef level auxiliary profile data perturbations
!! @param[in]     predictors      calculated predictors
!! @param[in,out] predictors_ad   calculated predictor perturbations
!! @param[in]     raytracing      RTTOV raytracing structure
!! @param[in,out] raytracing_ad   raytracing structure perturbations
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
SUBROUTINE rttov_setpredictors_13_ad( &
             opts,          &
             prof,          &
             coef,          &
             aux,           &
             aux_ad,        &
             predictors,    &
             predictors_ad, &
             raytracing,    &
             raytracing_ad, &
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
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_ad
  TYPE(rttov_path_pred),        INTENT(IN)    :: predictors(:)
  TYPE(rttov_path_pred),        INTENT(INOUT) :: predictors_ad(:)
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_ad
  INTEGER(jpim),                INTENT(IN)    :: raypath
!INTF_END
 
  INTEGER(jpim) :: lay, i, nlayers, nprofiles

  REAL(jprb), POINTER :: path(:,:), path_sqrt(:,:), path_4rt(:,:)
  REAL(jprb), POINTER :: path_ad(:,:), path_sqrt_ad(:,:), path_4rt_ad(:,:)

  REAL(jprb) :: ztemp_1d(prof(1)%nlayers)

  REAL(jprb) :: sec_tr(prof(1)%nlayers, SIZE(prof)), sec_twr(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r(prof(1)%nlayers), sec_wrwrtr_r(prof(1)%nlayers)
  REAL(jprb) :: sec_gr(prof(1)%nlayers), sec_gr_sqrt(prof(1)%nlayers)
  REAL(jprb) :: sec_gw(prof(1)%nlayers), sec_gw_sqrt(prof(1)%nlayers), sec_gw_4rt(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr(prof(1)%nlayers)

  REAL(jprb) :: sec_tr_ad(prof(1)%nlayers, SIZE(prof)), sec_twr_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r_ad(prof(1)%nlayers), sec_wrwrtr_r_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gr_ad(prof(1)%nlayers), sec_gr_sqrt_ad(prof(1)%nlayers), sec_gr_4rt_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gw_ad(prof(1)%nlayers), sec_gw_sqrt_ad(prof(1)%nlayers), sec_gw_4rt_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr_ad(prof(1)%nlayers)

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13_AD', 0_jpim, ZHOOK_HANDLE)

  IF (raypath == 1) THEN
    path => raytracing%pathsat
    path_sqrt => aux%pathsat_sqrt
    path_4rt => aux%pathsat_4rt
    path_ad => raytracing_ad%pathsat
    path_sqrt_ad => aux_ad%pathsat_sqrt
    path_4rt_ad => aux_ad%pathsat_4rt
  ELSEIF (raypath == 2) THEN
    path => raytracing%patheff
    path_sqrt => aux%patheff_sqrt
    path_4rt => aux%patheff_4rt
    path_ad => raytracing_ad%patheff
    path_sqrt_ad => aux_ad%patheff_sqrt
    path_4rt_ad => aux_ad%patheff_4rt
  ENDIF

  nlayers = prof(1)%nlayers
  nprofiles = SIZE(prof)

  ! ---------------------------------------------------------------------------
  ! Mixed gases
  ! ---------------------------------------------------------------------------

  sec_tr = path * aux%tr

  DO i = 1, nprofiles
    DO lay = 1, nlayers
      path_ad(lay,i) = path_ad(lay,i) + &
                                                           predictors_ad(i)%mixedgas(1,lay) + &
        2._jprb * path(lay,i) *                            predictors_ad(i)%mixedgas(2,lay) + &
        aux%twr(lay,i) *                                   predictors_ad(i)%mixedgas(7,lay) + &
        1.5_jprb * aux%tr_sqrt(lay,i) * path_sqrt(lay,i) * predictors_ad(i)%mixedgas(9,lay)

      sec_tr_ad(lay,i) = &!sec_tr_ad(lay,i) + &
                                                           predictors_ad(i)%mixedgas(3,lay) + &
        aux%tr(lay,i) *                                    predictors_ad(i)%mixedgas(4,lay) + &
        aux%tr2(lay,i) *                                   predictors_ad(i)%mixedgas(8,lay)

      aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
        sec_tr(lay,i) *                                    predictors_ad(i)%mixedgas(4,lay) + &
                                                           predictors_ad(i)%mixedgas(5,lay)

      aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
                                                           predictors_ad(i)%mixedgas(6,lay) + &
        sec_tr(lay,i) *                                    predictors_ad(i)%mixedgas(8,lay)

      aux_ad%twr(lay,i) = aux_ad%twr(lay,i) + &
        path(lay,i) *                                      predictors_ad(i)%mixedgas(7,lay)

      aux_ad%tr_sqrt(lay,i) = aux_ad%tr_sqrt(lay,i) + &
        path_sqrt(lay,i) * path(lay,i) *                   predictors_ad(i)%mixedgas(9,lay)

      predictors_ad(i)%mixedgas(10,lay) = 0._jprb
    ENDDO
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Water vapour
  ! ---------------------------------------------------------------------------

  DO i = 1, nprofiles
    sec_gr_sqrt = path_sqrt(:,i) * aux%wr_sqrt(:,i)
    sec_gr      = sec_gr_sqrt**2

    sec_gw_4rt  = path_4rt(:,i) * aux%ww_4rt(:,i)
    sec_gw_sqrt = path_sqrt(:,i) * aux%ww_sqrt(:,i)
    sec_gw      = path(:,i) * aux%ww(:,i)

    IF (coef%nwvcont > 0) THEN
      sec_wrtr_r   = sec_gr * aux%tr_r(:,i)
      sec_wrwrtr_r = sec_wrtr_r * aux%wr(:,i)

      DO lay = 1, nlayers
        sec_wrwrtr_r_ad(lay) = &!sec_wrwrtr_r_ad(lay) + &
                                                             predictors_ad(i)%wvcont(1,lay) + &
          aux%tr_r(lay,i)**3 *                               predictors_ad(i)%wvcont(3,lay)

        aux_ad%tr_r(lay,i) = aux_ad%tr_r(lay,i) + &
          3._jprb * aux%tr_r(lay,i)**2 * sec_wrwrtr_r(lay) * predictors_ad(i)%wvcont(3,lay) + &
          sec_wrtr_r(lay) *                                  predictors_ad(i)%wvcont(4,lay)

        sec_wrtr_r_ad(lay) = &!sec_wrtr_r_ad(lay) + &
                                                             predictors_ad(i)%wvcont(2,lay) + &
          aux%tr_r(lay,i) *                                  predictors_ad(i)%wvcont(4,lay)
      ENDDO

      aux_ad%wr(:,i) = aux_ad%wr(:,i) + &
        sec_wrtr_r * sec_wrwrtr_r_ad

      sec_wrtr_r_ad = sec_wrtr_r_ad + &
        aux%wr(:,i) * sec_wrwrtr_r_ad

      sec_gr_ad = &!sec_gr_ad + &
         aux%tr_r(:,i) * sec_wrtr_r_ad

      aux_ad%tr_r(:,i) = aux_ad%tr_r(:,i) + &
         sec_gr * sec_wrtr_r_ad
    ELSE
      sec_gr_ad = 0._jprb ! need sec_gr_ad to be initialised
    ENDIF

    DO lay = 1, nlayers
      sec_gr_ad(lay) = sec_gr_ad(lay) + &
        2._jprb * sec_gr(lay) *            predictors_ad(i)%watervapour(1,lay) + &
        aux%dt(lay,i) *                    predictors_ad(i)%watervapour(4,lay) + &
                                           predictors_ad(i)%watervapour(7,lay) + &
        sec_gr_sqrt(lay) *                 predictors_ad(i)%watervapour(9,lay) + &
        sec_gr_sqrt(lay) * aux%dt(lay,i) * predictors_ad(i)%watervapour(10,lay) + &
        aux%wrw_r(lay,i) *                 predictors_ad(i)%watervapour(13,lay)

      sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
                                           predictors_ad(i)%watervapour(2,lay) + &
        2._jprb * sec_gw(lay) *            predictors_ad(i)%watervapour(3,lay) + &
        sec_gw_sqrt(lay) *                 predictors_ad(i)%watervapour(8,lay) + &
        sec_gw_4rt(lay) *                  predictors_ad(i)%watervapour(12,lay)

      sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                           predictors_ad(i)%watervapour(5,lay) + &
        sec_gr(lay) *                      predictors_ad(i)%watervapour(9,lay) + &
        sec_gr(lay) * aux%dt(lay,i) *      predictors_ad(i)%watervapour(10,lay) + &
        aux%dt(lay,i) *                    predictors_ad(i)%watervapour(11,lay) + &
        aux%wrwr_r(lay,i) *                predictors_ad(i)%watervapour(14,lay)

      sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
                                           predictors_ad(i)%watervapour(6,lay)
      
      aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
        sec_gr(lay) *                      predictors_ad(i)%watervapour(4,lay) + &
        sec_gr(lay) * sec_gr_sqrt(lay) *   predictors_ad(i)%watervapour(10,lay) + &
        sec_gr_sqrt(lay) *                 predictors_ad(i)%watervapour(11,lay)

      sec_gw_sqrt_ad(lay) = &!sec_gw_sqrt_ad(lay) + &
        sec_gw(lay) *                      predictors_ad(i)%watervapour(8,lay) + &
        path_sqrt(lay,i) *                 predictors_ad(i)%watervapour(15,lay)

      sec_gw_4rt_ad(lay) = &!sec_gw_4rt_ad(lay) + &
        sec_gw(lay) *                      predictors_ad(i)%watervapour(12,lay)

      aux_ad%wrw_r(lay,i) = aux_ad%wrw_r(lay,i) + &
        sec_gr(lay) *                      predictors_ad(i)%watervapour(13,lay)

      aux_ad%wrwr_r(lay,i) = aux_ad%wrwr_r(lay,i) + &
        sec_gr_sqrt(lay) *                 predictors_ad(i)%watervapour(14,lay)

      path_sqrt_ad(lay,i) = &!path_sqrt_ad(lay,i) + &
        sec_gw_sqrt(lay) *                 predictors_ad(i)%watervapour(15,lay)
    ENDDO

    sec_gr_sqrt_ad = sec_gr_sqrt_ad + &
      2._jprb * sec_gr_sqrt * sec_gr_ad

    aux_ad%wr_sqrt(:,i) = aux_ad%wr_sqrt(:,i) + &
      path_sqrt(:,i) * sec_gr_sqrt_ad

    aux_ad%wr_4rt(:,i) = aux_ad%wr_4rt(:,i) + &
      path_4rt(:,i) * sec_gr_4rt_ad

    aux_ad%ww(:,i) = aux_ad%ww(:,i) + &
      path(:,i) * sec_gw_ad

    aux_ad%ww_sqrt(:,i) = aux_ad%ww_sqrt(:,i) + &
      path_sqrt(:,i) * sec_gw_sqrt_ad

    aux_ad%ww_4rt(:,i) = aux_ad%ww_4rt(:,i) + &
      path_4rt(:,i) * sec_gw_4rt_ad

    path_4rt_ad(:,i) = &!path_4rt_ad(:,i) + &
      aux%ww_4rt(:,i) * sec_gw_4rt_ad + &
      aux%wr_4rt(:,i) * sec_gr_4rt_ad

    path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
      aux%wr_sqrt(:,i) * sec_gr_sqrt_ad + &
      aux%ww_sqrt(:,i) * sec_gw_sqrt_ad

    path_ad(:,i) = path_ad(:,i) + & 
      aux%ww(:,i) * sec_gw_ad
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Ozone
  ! ---------------------------------------------------------------------------

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nprofiles
        sec_gr_sqrt = path_sqrt(:,i) * aux%or_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2
        sec_gw_sqrt = path_sqrt(:,i) * aux%ow_sqrt(:,i)
        sec_gw      = sec_gw_sqrt**2
        sec_gw_4rt  = path_4rt(:,i) * aux%ow_4rt(:,i)

        DO lay = 1, nlayers
          sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
                                                             predictors_ad(i)%ozone(1,lay) + &
            aux%dt(lay,i) *                                  predictors_ad(i)%ozone(3,lay) + &
            aux%ow_r(lay,i) *                                predictors_ad(i)%ozone(4,lay) + &
            2._jprb * sec_gr(lay) *                          predictors_ad(i)%ozone(5,lay) + &
            aux%or(lay,i) * aux%ow(lay,i) *                  predictors_ad(i)%ozone(6,lay) + &
            aux%ow(lay,i) *                                  predictors_ad(i)%ozone(8,lay) + &
            sec_gw_sqrt(lay) *                               predictors_ad(i)%ozone(10,lay)

          sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                                             predictors_ad(i)%ozone(2,lay) + &
            aux%or(lay,i) * aux%ow_r(lay,i) *                predictors_ad(i)%ozone(7,lay)

          aux_ad%or(lay,i) = aux_ad%or(lay,i) + &
            sec_gr(lay) * aux%ow(lay,i) *                    predictors_ad(i)%ozone(6,lay) + &
            sec_gr_sqrt(lay) * aux%ow_r(lay,i) *             predictors_ad(i)%ozone(7,lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
            sec_gr(lay) *                                    predictors_ad(i)%ozone(3,lay) + &
            path_sqrt(lay,i) * aux%ow(lay,i)**2 *            predictors_ad(i)%ozone(12,lay)

          sec_gw_sqrt_ad(lay) = &! sec_gw_sqrt_ad(lay) + &
            sec_gr(lay) *                                    predictors_ad(i)%ozone(10,lay)

          aux_ad%ow(lay,i) = aux_ad%ow(lay,i) + &
            sec_gr(lay) * aux%or(lay,i) *                    predictors_ad(i)%ozone(6,lay) + &
            sec_gr(lay) *                                    predictors_ad(i)%ozone(8,lay) + &
            2._jprb * aux%ow(lay,i) * &
              path_sqrt(lay,i) * aux%dt(lay,i) *             predictors_ad(i)%ozone(12,lay)

          aux_ad%ow_r(lay,i) = aux_ad%ow_r(lay,i) + &
            sec_gr(lay) *                                    predictors_ad(i)%ozone(4,lay) + &
            sec_gr_sqrt(lay) * aux%or(lay,i) *               predictors_ad(i)%ozone(7,lay)

          sec_gw_ad(lay) = &! sec_gw_ad(lay) + &
            1.75_jprb * sec_gw_sqrt(lay) * sec_gw_4rt(lay) * predictors_ad(i)%ozone(9,lay) + &
            2._jprb * sec_gw(lay) *                          predictors_ad(i)%ozone(11,lay) + &
                                                             predictors_ad(i)%ozone(13,lay)

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
            aux%dt(lay,i) * aux%ow(lay,i)**2 *               predictors_ad(i)%ozone(12,lay)
        ENDDO

        sec_gr_sqrt_ad = sec_gr_sqrt_ad + &
          2._jprb * sec_gr_sqrt * sec_gr_ad
        sec_gw_sqrt_ad = sec_gw_sqrt_ad + &
          2._jprb * sec_gw_sqrt * sec_gw_ad

        aux_ad%or_sqrt(:,i) = aux_ad%or_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gr_sqrt_ad
        aux_ad%ow_sqrt(:,i) = aux_ad%ow_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gw_sqrt_ad

        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
          aux%ow_sqrt(:,i) * sec_gw_sqrt_ad + &
          aux%or_sqrt(:,i) * sec_gr_sqrt_ad
      ENDDO
    ELSE ! no user ozone data
      DO i = 1, nprofiles   
        DO lay = 1, nlayers
          path_ad(lay,i) = path_ad(lay,i) + &
                                                             predictors_ad(i)%ozone(1,lay) + &
            aux%dt(lay,i) *                                  predictors_ad(i)%ozone(3,lay) + &
                                                             predictors_ad(i)%ozone(4,lay) + &
                                                             predictors_ad(i)%ozone(6,lay) + &
                                                             predictors_ad(i)%ozone(8,lay) + &
            1.75_jprb * path_sqrt(lay,i) * path_4rt(lay,i) * predictors_ad(i)%ozone(9,lay) + &
            1.5_jprb * path_sqrt(lay,i) *                    predictors_ad(i)%ozone(10,lay) + &
            2._jprb * path(lay,i) *                         (predictors_ad(i)%ozone(5,lay) + &
                                                             predictors_ad(i)%ozone(11,lay)) + &
                                                             predictors_ad(i)%ozone(13,lay)

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                                                             predictors_ad(i)%ozone(2,lay) + &
                                                             predictors_ad(i)%ozone(7,lay) + &
            aux%dt(lay,i) *                                  predictors_ad(i)%ozone(12,lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
            path(lay,i) *                                    predictors_ad(i)%ozone(3,lay) + &
            path_sqrt(lay,i) *                               predictors_ad(i)%ozone(12,lay)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO i = 1, nprofiles

    ! -----------------------------------------------------------------------
    ! CO2
    ! -----------------------------------------------------------------------

    IF (coef%nco2 > 0) THEN
      sec_twr = path(:,i) * aux%twr(:,i)

      DO lay = 1, nlayers
        aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
                                                           predictors_ad(i)%co2(2,lay) + &
          aux%twr(lay,i)**2 * &
            (aux%twr(lay,i) * path_sqrt(lay,i) *           predictors_ad(i)%co2(12,lay) + &
                                                           predictors_ad(i)%co2(13,lay))

        sec_tr_ad(lay,i) = sec_tr_ad(lay,i) + &
                                                           predictors_ad(i)%co2(3,lay) + &
          aux%tr(lay,i) *                                  predictors_ad(i)%co2(4,lay) + &
          aux%tr(lay,i)**2 *                               predictors_ad(i)%co2(11,lay)

        aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
          sec_tr(lay,i) *                                  predictors_ad(i)%co2(4,lay) + &
                                                           predictors_ad(i)%co2(5,lay) + &
          3._jprb * aux%tr2(lay,i) *                       predictors_ad(i)%co2(10,lay) + &
          aux%tr(lay,i) * sec_tr(lay,i) * 2._jprb *        predictors_ad(i)%co2(11,lay)

        sec_twr_ad(lay) = &!sec_twr_ad(lay) + &
                                                           predictors_ad(i)%co2(6,lay) + &
          aux%tr_sqrt(lay,i) *                             predictors_ad(i)%co2(8,lay)

        aux_ad%twr(lay,i) = aux_ad%twr(lay,i) + &
          aux%twr(lay,i) * aux%tr2(lay,i) * &
            (3._jprb * path_sqrt(lay,i) * aux%twr(lay,i) * predictors_ad(i)%co2(12,lay) + &
                                                 2._jprb * predictors_ad(i)%co2(13,lay))

        aux_ad%tr_sqrt(lay,i) = aux_ad%tr_sqrt(lay,i) + &
          sec_twr(lay) *                                   predictors_ad(i)%co2(8,lay)
          
        path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
          aux%tr2(lay,i) * aux%twr(lay,i)**3 *             predictors_ad(i)%co2(12,lay)

        IF (opts%rt_all%co2_data) THEN
          path_ad(lay,i) = path_ad(lay,i) + &
            aux%co2r(lay,i) *                              predictors_ad(i)%co2(1,lay) + &
            2._jprb * path(lay,i) * aux%co2w(lay,i)**2 *   predictors_ad(i)%co2(7,lay) + &
            aux%co2w(lay,i) *                              predictors_ad(i)%co2(14,lay)

          aux_ad%co2r(lay,i) = aux_ad%co2r(lay,i) + &
            path(lay,i) *                                  predictors_ad(i)%co2(1,lay)

          aux_ad%co2w(lay,i) = aux_ad%co2w(lay,i) + &
            2._jprb * path(lay,i)**2 * aux%co2w(lay,i) *   predictors_ad(i)%co2(7,lay) + &
            path(lay,i) *                                  predictors_ad(i)%co2(14,lay)

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
            aux%co2r_sqrt(lay,i) *                         predictors_ad(i)%co2(9,lay)

          aux_ad%co2r_sqrt(lay,i) = aux_ad%co2r_sqrt(lay,i) + &
            path_sqrt(lay,i) *                             predictors_ad(i)%co2(9,lay)
        ELSE
          path_ad(lay,i) = path_ad(lay,i) + &
                                                           predictors_ad(i)%co2(1,lay) + &
            2._jprb * path(lay,i) *                        predictors_ad(i)%co2(7,lay) + &
                                                           predictors_ad(i)%co2(14,lay)

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                                                           predictors_ad(i)%co2(9,lay)
        ENDIF
      ENDDO

      aux_ad%twr(:,i) = aux_ad%twr(:,i) + &
        path(:,i) * sec_twr_ad

      path_ad(:,i) = path_ad(:,i) + &
        aux%twr(:,i) * sec_twr_ad
    ENDIF

    ! ---------------------------------------------------------------------
    ! N2O
    ! ---------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN
      sec_gwr = path(:,i) * aux%n2owr(:,i)

      IF (opts%rt_all%n2o_data) THEN
        sec_gr_sqrt = path_sqrt(:,i) * aux%n2or_sqrt(:,i)
        sec_gr      = path(:,i) * aux%n2or(:,i)

        DO lay = 1, nlayers
          aux_ad%n2or(lay,i) = aux_ad%n2or(lay,i) + &
            aux%dt(lay,i) *                        predictors_ad(i)%n2o(5,lay) + &
            sec_gr_sqrt(lay) * aux%n2ow_r(lay,i) * predictors_ad(i)%n2o(9,lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                              aux%n2or(lay,i) *    predictors_ad(i)%n2o(5,lay)

          aux_ad%n2ow(lay,i) = aux_ad%n2ow(lay,i) + &
            path(lay,i) *                          predictors_ad(i)%n2o(7,lay)

          path_ad(lay,i) = path_ad(lay,i) + &
            aux%n2ow(lay,i) *                      predictors_ad(i)%n2o(7,lay)

          aux_ad%n2ow_r(lay,i) = aux_ad%n2ow_r(lay,i) + &
            sec_gr_sqrt(lay) * aux%n2or(lay,i) *   predictors_ad(i)%n2o(9,lay)

          sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
            aux%n2or(lay,i) * aux%n2ow_r(lay,i) *  predictors_ad(i)%n2o(9,lay)
        ENDDO
      ELSE
        sec_gr_sqrt    = path_sqrt(:,i)
        sec_gr         = path(:,i)

        aux_ad%dt(:,i) = aux_ad%dt(:,i) +          predictors_ad(i)%n2o(5,:)

        path_ad(:,i) = path_ad(:,i) + &
                                                   predictors_ad(i)%n2o(7,:)

        sec_gr_sqrt_ad = &!sec_gr_sqrt_ad + &
                                                   predictors_ad(i)%n2o(9,:)
      ENDIF

      DO lay = 1, nlayers
        path_ad(lay,i) = path_ad(lay,i) + &
          sec_gwr(lay) * aux%dt(lay,i) *           predictors_ad(i)%n2o(12,lay)

        sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
          predictors_ad(i)%n2o(1,lay) + &
          aux%dt(lay,i) *                          predictors_ad(i)%n2o(3,lay) + &
          2._jprb * sec_gr(lay) *                  predictors_ad(i)%n2o(4,lay)

        sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
                                                   predictors_ad(i)%n2o(2,lay)

        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
                                                   predictors_ad(i)%n2o(6,lay)

        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
          sec_gr(lay) *                            predictors_ad(i)%n2o(3,lay) + &
          path(lay,i) * sec_gwr(lay)  *            predictors_ad(i)%n2o(12,lay)

        sec_gwr_ad(lay) = &!sec_gwr_ad(lay) + &
                                                   predictors_ad(i)%n2o(8,lay) + &
          2._jprb * sec_gwr(lay) *                 predictors_ad(i)%n2o(10,lay) + &
          3._jprb * sec_gwr(lay)**2 *              predictors_ad(i)%n2o(11,lay) + &
          path(lay,i) * aux%dt(lay,i) *            predictors_ad(i)%n2o(12,lay)
      ENDDO

      IF (opts%rt_all%n2o_data) THEN
        sec_gr_sqrt_ad = sec_gr_sqrt_ad + &
          2._jprb * sec_gr_sqrt * sec_gr_ad

        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
          aux%n2or_sqrt(:,i) * sec_gr_sqrt_ad

        aux_ad%n2or_sqrt(:,i) = aux_ad%n2or_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gr_sqrt_ad

        path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
          aux%n2or_4rt(:,i) * sec_gr_4rt_ad

        aux_ad%n2or_4rt(:,i) = aux_ad%n2or_4rt(:,i) + &
          path_4rt(:,i) * sec_gr_4rt_ad
      ELSE
        path_ad(:,i) = path_ad(:,i) + sec_gr_ad
        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
        path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
      ENDIF

      path_ad(:,i) = path_ad(:,i) + &
        aux%n2owr(:,i) * sec_gwr_ad

      aux_ad%n2owr(:,i) = aux_ad%n2owr(:,i) + &
        path(:,i) * sec_gwr_ad
    ENDIF

    ! ---------------------------------------------------------------------
    ! CO
    ! ---------------------------------------------------------------------

    IF (coef%nco > 0) THEN
      sec_gwr       = path(:,i) * aux%cowr(:,i)

      IF (opts%rt_all%co_data) THEN
        sec_gr_sqrt = path_sqrt(:,i) * aux%cor_sqrt(:,i)
        sec_gr      = sec_gr_sqrt**2
        sec_gw      = path(:,i) * aux%cow(:,i)
        ztemp_1d    = 0.4_jprb * predictors(i)%co(11,:) / sec_gw

        DO lay = 1, nlayers
          sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
            aux%corw_r(lay,i) *            predictors_ad(i)%co(8,lay) + &
            aux%corw_rsqrt(lay,i) *        predictors_ad(i)%co(10,lay) + &
            sec_gw(lay) *                  predictors_ad(i)%co(13,lay)

          sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
            aux%corw_r(lay,i) *            predictors_ad(i)%co(9,lay)

          sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
            ztemp_1d(lay) *                predictors_ad(i)%co(11,lay) + &
            sec_gr(lay) *                  predictors_ad(i)%co(13,lay) + &
                                           predictors_ad(i)%co(14,lay) + &
            2._jprb * sec_gw(lay) *        predictors_ad(i)%co(16,lay)

          aux_ad%corw_r(lay,i) = aux_ad%corw_r(lay,i) + &
             sec_gr(lay) *                 predictors_ad(i)%co(8,lay) + &
             sec_gr_sqrt(lay) *            predictors_ad(i)%co(9,lay)

          aux_ad%corw_rsqrt(lay,i) = aux_ad%corw_rsqrt(lay,i) + &
            sec_gr(lay) *                  predictors_ad(i)%co(10,lay)
        ENDDO
      ELSE
        sec_gr         = path(:,i)
        sec_gr_sqrt    = path_sqrt(:,i)
        ztemp_1d       = 0.4_jprb * predictors(i)%co(11,:) / path(:,i)

        DO lay = 1, nlayers
          path_ad(lay,i) = path_ad(lay,i) + &
                                           predictors_ad(i)%co(8,lay) + &
                                           predictors_ad(i)%co(10,lay) + &
            ztemp_1d(lay) *                predictors_ad(i)%co(11,lay) + &
                                           predictors_ad(i)%co(14,lay) + &
            2._jprb * path(lay,i) *       (predictors_ad(i)%co(13,lay) + &
                                           predictors_ad(i)%co(16,lay))

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                                           predictors_ad(i)%co(9,lay)
        ENDDO
        sec_gr_ad      = 0._jprb
        sec_gr_sqrt_ad = 0._jprb
      ENDIF

      DO lay = 1, nlayers
        sec_gr_ad(lay) = sec_gr_ad(lay) + &
                                           predictors_ad(i)%co(1,lay) + &
          aux%dt(lay,i) *                  predictors_ad(i)%co(3,lay) + &
          2._jprb * sec_gr(lay) *          predictors_ad(i)%co(4,lay) + &
          aux%dtabsdt(lay,i) *             predictors_ad(i)%co(7,lay)
        
        sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
                                           predictors_ad(i)%co(2,lay) + &
          aux%dt(lay,i) *                  predictors_ad(i)%co(5,lay)

        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
                                           predictors_ad(i)%co(6,lay)

        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
          sec_gr(lay) *                    predictors_ad(i)%co(3,lay) + &
          sec_gr_sqrt(lay) *               predictors_ad(i)%co(5,lay)

        aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
          sec_gr(lay) *                    predictors_ad(i)%co(7,lay)

        aux_ad%cowr_4rt(lay,i) = aux_ad%cowr_4rt(lay,i) + &
          path_4rt(lay,i) *                predictors_ad(i)%co(12,lay)

        path_4rt_ad(lay,i) = path_4rt_ad(lay,i) + &
          aux%cowr_4rt(lay,i) *            predictors_ad(i)%co(12,lay)
      ENDDO

      IF (opts%rt_all%co_data) THEN
        aux_ad%cow(:,i) = aux_ad%cow(:,i) + &
           path(:,i) * sec_gw_ad

        path_ad(:,i) = path_ad(:,i) + &
          aux%cow(:,i) * sec_gw_ad

        sec_gr_sqrt_ad = sec_gr_sqrt_ad + &
          2._jprb * sec_gr_sqrt * sec_gr_ad

        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
          aux%cor_sqrt(:,i) * sec_gr_sqrt_ad

        aux_ad%cor_sqrt(:,i) = aux_ad%cor_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gr_sqrt_ad

        path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
          aux%cor_4rt(:,i) * sec_gr_4rt_ad
                  
        aux_ad%cor_4rt(:,i) = aux_ad%cor_4rt(:,i) + &
          path_4rt(:,i) * sec_gr_4rt_ad
      ELSE
        path_ad(:,i) = path_ad(:,i) + sec_gr_ad
        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
        path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
      ENDIF

      sec_gwr_ad = &!sec_gwr_ad + &
        predictors_ad(i)%co(15,:)

      path_ad(:,i)  = path_ad(:,i) + &
        aux%cowr(:,i) * sec_gwr_ad

      aux_ad%cowr(:,i) = aux_ad%cowr(:,i) + &
         path(:,i) * sec_gwr_ad
    ENDIF

    ! ---------------------------------------------------------------------
    ! CH4
    ! ---------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%rt_all%ch4_data) THEN
        sec_gr_sqrt = path_sqrt(:,i) * aux%ch4r_sqrt(:,i)
        sec_gr      = path(:,i) * aux%ch4r(:,i)
        sec_gw      = path(:,i) * aux%ch4w(:,i)
        sec_gw_4rt  = path_4rt(:,i) * aux%ch4w_4rt(:,i)

        DO lay = 1, nlayers
          sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
            aux%ch4rw_r(lay,i) *  predictors_ad(i)%ch4(11,lay)

          aux_ad%ch4r(lay,i) = aux_ad%ch4r(lay,i) + &
            aux%dt(lay,i) *       predictors_ad(i)%ch4(5,lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
            aux%ch4r(lay,i) *     predictors_ad(i)%ch4(5,lay)

          aux_ad%ch4rw_r(lay,i) = aux_ad%ch4rw_r(lay,i) + &
            sec_gr_sqrt(lay) *    predictors_ad(i)%ch4(11,lay)
        ENDDO
      ELSE
        sec_gr_sqrt = path_sqrt(:,i)
        sec_gr      = path(:,i)
        sec_gw      = path(:,i)
        sec_gw_4rt  = path_4rt(:,i)

        DO lay = 1, nlayers
          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                                  predictors_ad(i)%ch4(11,lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                                  predictors_ad(i)%ch4(5,lay)
        ENDDO
        sec_gr_sqrt_ad = 0._jprb
      ENDIF

      DO lay = 1, nlayers
        sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
          predictors_ad(i)%ch4(1,lay) + &
          aux%dt(lay,i) * predictors_ad(i)%ch4(3,lay) + &
          2._jprb * sec_gr(lay) * predictors_ad(i)%ch4(4,lay)

        path_ad(lay,i) = path_ad(lay,i) + &
          aux%ch4wr(lay,i) *      predictors_ad(i)%ch4(7,lay)
          
        sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
                                  predictors_ad(i)%ch4(2,lay)
          
        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
                                  predictors_ad(i)%ch4(6,lay)

        sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
          2._jprb * sec_gw(lay) * predictors_ad(i)%ch4(9,lay) + &
                                  predictors_ad(i)%ch4(10,lay) + &
          sec_gw_4rt(lay) *       predictors_ad(i)%ch4(12,lay)

        aux_ad%ch4wr(lay,i) = aux_ad%ch4wr(lay,i) + &
          path(lay,i) *           predictors_ad(i)%ch4(7,lay) + &
                                  predictors_ad(i)%ch4(8,lay)

        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
          sec_gr(lay) *           predictors_ad(i)%ch4(3,lay)

        sec_gw_4rt_ad(lay) = &!sec_gw_4rt_ad(lay) + &
          sec_gw(lay) *           predictors_ad(i)%ch4(12,lay)
      ENDDO

      IF (opts%rt_all%ch4_data) THEN
! using sec_gr = path * ch4r rather than sec_g_sqrt**2
        path_ad(:,i) = path_ad(:,i) + &
          aux%ch4w(:,i) * sec_gw_ad + &
          aux%ch4r(:,i) * sec_gr_ad

        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
          aux%ch4r_sqrt(:,i) * sec_gr_sqrt_ad

        path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
          aux%ch4r_4rt(:,i) * sec_gr_4rt_ad + &
          aux%ch4w_4rt(:,i) * sec_gw_4rt_ad

        aux_ad%ch4w(:,i) = aux_ad%ch4w(:,i) + &
          path(:,i) * sec_gw_ad

        aux_ad%ch4w_4rt(:,i) = aux_ad%ch4w_4rt(:,i) + &
          path_4rt(:,i) * sec_gw_4rt_ad

        aux_ad%ch4r(:,i) = aux_ad%ch4r(:,i) + &
          path(:,i) * sec_gr_ad

        aux_ad%ch4r_sqrt(:,i) = aux_ad%ch4r_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gr_sqrt_ad

        aux_ad%ch4r_4rt(:,i) = aux_ad%ch4r_4rt(:,i) + &
          path_4rt(:,i) * sec_gr_4rt_ad
      ELSE
        path_ad(:,i) = path_ad(:,i) + sec_gr_ad + sec_gw_ad
        path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gw_4rt_ad
        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
        path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------
    ! SO2
    ! ---------------------------------------------------------------------

    IF (coef%nso2 > 0) THEN
      sec_gr_sqrt = path_sqrt(:,i) * aux%so2r_sqrt(:,i)
      sec_gr      = sec_gr_sqrt**2

      sec_gw_sqrt = path_sqrt(:,i) * aux%so2w_sqrt(:,i)
      sec_gw_4rt  = path_4rt(:,i) * aux%so2w_4rt(:,i)
      sec_gw      = sec_gw_sqrt**2

      DO lay=1, nlayers
        sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
          2._jprb * sec_gr(lay) *            predictors_ad(i)%so2(1,lay) + &
          aux%dt(lay,i) *                    predictors_ad(i)%so2(4,lay) + &
                                             predictors_ad(i)%so2(7,lay) + &
          sec_gr_sqrt(lay) *                 predictors_ad(i)%so2(10,lay) + &
          sec_gr_sqrt(lay) * aux%dt(lay,i) * predictors_ad(i)%so2(11,lay) + &
          aux%so2rw_r(lay,i) *               predictors_ad(i)%so2(14,lay)

        sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
                                             predictors_ad(i)%so2(2,lay) + &
          2._jprb * sec_gw(lay) *            predictors_ad(i)%so2(3,lay) + &
          sec_gw_sqrt(lay) *                 predictors_ad(i)%so2(9,lay) + &
          sec_gw_4rt(lay) *                  predictors_ad(i)%so2(13,lay)

        sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                             predictors_ad(i)%so2(5,lay) + &
          aux%so2rwr_r(lay,i) *              predictors_ad(i)%so2(8,lay) + &
          sec_gr(lay) *                      predictors_ad(i)%so2(10,lay) + &
          aux%dt(lay,i) * sec_gr(lay) *      predictors_ad(i)%so2(11,lay) + &
          aux%dt(lay,i) *                    predictors_ad(i)%so2(12,lay)

        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
                                             predictors_ad(i)%so2(6,lay)

        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
          sec_gr(lay) *                      predictors_ad(i)%so2(4,lay) + &
          sec_gr(lay) * sec_gr_sqrt(lay) *   predictors_ad(i)%so2(11,lay) + &
          sec_gr_sqrt(lay) *                 predictors_ad(i)%so2(12,lay)

        aux_ad%so2rwr_r(lay,i) = aux_ad%so2rwr_r(lay,i) + &
          sec_gr_sqrt(lay) *                 predictors_ad(i)%so2(8,lay)

        sec_gw_sqrt_ad(lay) = &!sec_gw_sqrt_ad(lay) + &
          sec_gw(lay) *                      predictors_ad(i)%so2(9,lay) + &
          path_sqrt(lay,i) *                 predictors_ad(i)%so2(15,lay)

        sec_gw_4rt_ad(lay) = &!sec_gw_4rt_ad(lay) + &
          sec_gw(lay) *                      predictors_ad(i)%so2(13,lay)

        aux_ad%so2rw_r(lay,i) = aux_ad%so2rw_r(lay,i) + &
          sec_gr(lay) *                      predictors_ad(i)%so2(14,lay)

        path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
          sec_gw_sqrt(lay) *                 predictors_ad(i)%so2(15,lay)
      ENDDO

      sec_gw_sqrt_ad = sec_gw_sqrt_ad + &
        2._jprb * sec_gw_sqrt * sec_gw_ad

      sec_gr_sqrt_ad = sec_gr_sqrt_ad + &
        2._jprb * sec_gr_sqrt * sec_gr_ad

      aux_ad%so2r_sqrt(:,i) = aux_ad%so2r_sqrt(:,i) + &
        path_sqrt(:,i) * sec_gr_sqrt_ad

      aux_ad%so2r_4rt(:,i) = aux_ad%so2r_4rt(:,i) + &
        path_4rt(:,i) * sec_gr_4rt_ad

      aux_ad%so2w_sqrt(:,i) = aux_ad%so2w_sqrt(:,i) + &
        path_sqrt(:,i) * sec_gw_sqrt_ad

      aux_ad%so2w_4rt(:,i) = aux_ad%so2w_4rt(:,i) + &
        path_4rt(:,i) * sec_gw_4rt_ad

      path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
        aux%so2w_4rt(:,i) * sec_gw_4rt_ad + &
        aux%so2r_4rt(:,i) * sec_gr_4rt_ad

      path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
        aux%so2w_sqrt(:,i) * sec_gw_sqrt_ad + &
        aux%so2r_sqrt(:,i) * sec_gr_sqrt_ad
    ENDIF
  ENDDO ! profiles

  ! DAR: gases should be done in reverse order according to AD coding theory but can almost be treated as being independent,
  ! making it easier to initialise arrays with calculations when the gas is mandatory - hence why we do mixedgas first but
  ! this leads to problems when considering some temperature AD variables which have dependencies on multiple gases).
  ! These need to be treated in the correct order
  path_ad = path_ad + aux%tr * sec_tr_ad
  aux_ad%tr = aux_ad%tr + path * sec_tr_ad

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_13_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_13_ad
