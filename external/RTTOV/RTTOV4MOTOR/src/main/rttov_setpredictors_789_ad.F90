! Description:
!> @file
!!   AD of predictor calculation.
!!
!> @brief
!!   AD of predictor calculation.
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
SUBROUTINE rttov_setpredictors_789_ad( &
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
  USE rttov_const, ONLY :  &
        inst_id_ssmis,     &
        inst_id_amsua
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
 
  INTEGER(jpim) :: lev, lay,i, iprof, ichan, nlayers, nprofiles

  REAL(jprb), POINTER :: path(:,:), path_sqrt(:,:), path_4rt(:,:)
  REAL(jprb), POINTER :: path_ad(:,:), path_sqrt_ad(:,:), path_4rt_ad(:,:)

  REAL(jprb) :: sec_tr(prof(1)%nlayers, SIZE(prof)), sec_twr(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r(prof(1)%nlayers), sec_wrwrtr_r(prof(1)%nlayers)
  REAL(jprb) :: sec_gr(prof(1)%nlayers), sec_gr_sqrt(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr(prof(1)%nlayers), sec_gwr_r(prof(1)%nlayers), sec_gw(prof(1)%nlayers)
  REAL(jprb) :: sec_gw_sqrt(prof(1)%nlayers), sec_gw_4rt(prof(1)%nlayers)

  REAL(jprb) :: sec_tr_ad(prof(1)%nlayers, SIZE(prof)), sec_twr_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r_ad(prof(1)%nlayers), sec_wrwrtr_r_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gr_ad(prof(1)%nlayers), sec_gr_sqrt_ad(prof(1)%nlayers), sec_gr_4rt_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gwr_ad(prof(1)%nlayers), sec_gw_ad(prof(1)%nlayers)
  REAL(jprb) :: sec_gw_sqrt_ad(prof(1)%nlayers), sec_gw_4rt_ad(prof(1)%nlayers)

  ! Pressure-modulated cell (PMC) variables
  REAL(jprb) :: Lcel_cm
  ! REAL(jprb) :: Tcel
  REAL(jprb) :: betaplus1
  REAL(jprb) :: acm
  REAL(jprb) :: Pcel_Lev
  REAL(jprb) :: Pnom_LevM1
  REAL(jprb) :: Pcel_LevM1
  REAL(jprb) :: Pnom_Lev
  REAL(jprb) :: Pnom(coef%fmv_chn)
  REAL(jprb) :: Pcel(coef%fmv_chn)
  REAL(jprb) :: acm_AD
  REAL(jprb) :: Pcel_Lev_AD
  REAL(jprb) :: Pnom_LevM1_AD
  REAL(jprb) :: Pcel_LevM1_AD
  REAL(jprb) :: Pnom_Lev_AD

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_789_AD', 0_jpim, ZHOOK_HANDLE)

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

  sec_tr = path(:, :) * aux%tr(:, :)
  
  IF (coef%fmv_model_ver <= 8) THEN

    path_4rt_ad = 0.0_jprb
    IF (coef%fmv_model_ver == 8) aux_ad%tr_sqrt = 0.0_jprb

    DO i = 1, nprofiles
      DO lay = 1, nlayers
        ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        
        path_ad(lay,i) = path_ad(lay,i) + &
                                                              predictors_ad(i)%mixedgas(1, lay) + &
                          2.0_jprb * path(lay,i) *           predictors_ad(i)%mixedgas(2, lay) + &
                          aux%tw(lay,i) *                    predictors_ad(i)%mixedgas(7, lay) + &
                          aux%tw(lay,i) * aux%tr_r(lay,i) * predictors_ad(i)%mixedgas(8, lay)

        sec_tr_ad(lay,i) = &!sec_tr_ad(lay,i) + &          
                            predictors_ad(i)%mixedgas(3, lay) + &
                            aux%tr(lay,i) * predictors_ad(i)%mixedgas(4, lay)

        aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
                            sec_tr(lay,i) * predictors_ad(i)%mixedgas(4, lay) + &
                            predictors_ad(i)%mixedgas(5, lay)  

        aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
                             predictors_ad(i)%mixedgas(6, lay)

        aux_ad%tw(lay,i) = aux_ad%tw(lay,i) + &
                            path(lay,i) * (                       predictors_ad(i)%mixedgas(7, lay) + &
                                            aux%tr_r(lay,i) *     predictors_ad(i)%mixedgas(8, lay))

        aux_ad%tr_r(lay,i) = aux_ad%tr_r(lay,i) + &
                              path(lay,i) * aux%tw(lay,i) * predictors_ad(i)%mixedgas(8, lay)

        path_sqrt_ad(lay,i) = &!path_sqrt_ad(lay,i) + &
                               predictors_ad(i)%mixedgas(9, lay) + &
                               aux%tw_4rt(lay,i) * predictors_ad(i)%mixedgas(10, lay)

        aux_ad%tw_4rt(lay,i) = aux_ad%tw_4rt(lay,i) + &
                                path_sqrt(lay,i) * predictors_ad(i)%mixedgas(10, lay)
      ENDDO
    ENDDO
  ELSE
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        path_ad(lay,i) = path_ad(lay,i) + &
                                                              predictors_ad(i)%mixedgas(1, lay) + &
                          2.0_jprb * path(lay,i) *           predictors_ad(i)%mixedgas(2, lay) + &
                          aux%tuw(lay,i) *                  (predictors_ad(i)%mixedgas(7, lay) + &
                                                              predictors_ad(i)%mixedgas(8, lay)) + &
                          1.5_jprb * aux%tr_sqrt(lay,i) * path_sqrt(lay,i) * predictors_ad(i)%mixedgas(10, lay)

        sec_tr_ad(lay,i) = &!sec_tr_ad(lay,i) + &
                            predictors_ad(i)%mixedgas(3, lay) + &
                            aux%tr(lay,i) * predictors_ad(i)%mixedgas(4, lay) + &
                            aux%tr2(lay,i) * predictors_ad(i)%mixedgas(9, lay)

        aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
                            sec_tr(lay,i) * predictors_ad(i)%mixedgas(4, lay) + &
                            predictors_ad(i)%mixedgas(5, lay)  

        aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
                             predictors_ad(i)%mixedgas(6, lay) + &
                             sec_tr(lay,i) * predictors_ad(i)%mixedgas(9, lay)

        aux_ad%tuw(lay,i) = aux_ad%tuw(lay,i) + &
                             path(lay,i) * (predictors_ad(i)%mixedgas(7, lay) + predictors_ad(i)%mixedgas(8, lay))

        aux_ad%tr_sqrt(lay,i) = aux_ad%tr_sqrt(lay,i) + &
                                 path_sqrt(lay,i) * path(lay,i) * predictors_ad(i)%mixedgas(10, lay)
      ENDDO
    ENDDO
    path_sqrt_ad = 0.0_jprb
  ENDIF

  IF (coef%id_inst == inst_id_ssmis .AND. coef%inczeeman) THEN
    aux_ad%t_layer = 0.0_jprb
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        ! SSMIS with Zeeman coefficient file
        ! geomagnetic field variables (Be, cosbk) are part of the user input

        ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        ! NB require prof(i) % Be >0. (divisor)
        
        ! X11 -> X21
        ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        predictors_ad(i)%mixedgas(20, lay) = predictors_ad(i)%mixedgas(20, lay) + &
             predictors_ad(i)%mixedgas(21, lay) * prof(i)%Be
        predictors_ad(i)%mixedgas(13, lay) = predictors_ad(i)%mixedgas(13, lay) + &
             predictors_ad(i)%mixedgas(20, lay) * prof(i)%Be

        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
             prof(i)%Be**3_jpim * predictors_ad(i)%mixedgas(19, lay)
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
             prof(i)%Be * predictors_ad(i)%mixedgas(18, lay)
        predictors_ad(i)%mixedgas(16, lay) = predictors_ad(i)%mixedgas(16, lay) + &
             predictors_ad(i)%mixedgas(17, lay) / prof(i)%Be
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
             predictors_ad(i)%mixedgas(16, lay) / prof(i)%Be
        predictors_ad(i)%mixedgas(12, lay) = predictors_ad(i)%mixedgas(12, lay) + &
             predictors_ad(i)%mixedgas(15, lay) * prof(i)%cosbk**2_jpim
        predictors_ad(i)%mixedgas(12, lay) = predictors_ad(i)%mixedgas(12, lay) + &
             predictors_ad(i)%mixedgas(14, lay) / prof(i)%Be

        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
             prof(i)%cosbk**2_jpim * predictors_ad(i)%mixedgas(13, lay)

        aux_ad%t_layer(lay,i) = aux_ad%t_layer(lay,i) - &
            (predictors(i)%mixedgas(12, lay) / &
             aux%t_layer(lay,i)) * predictors_ad(i)%mixedgas(12, lay)

        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
             (300.0_jprb / aux%t_layer(lay,i)) * predictors_ad(i)%mixedgas(12, lay)
        raytracing_ad%pathsat(lay,i) =  raytracing_ad%pathsat(lay,i) + predictors_ad(i)%mixedgas(11, lay)
      ENDDO
    ENDDO
  ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%inczeeman) THEN
    ! AMSU-A with Zeeman coefficient file
    ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
    ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
    DO i = 1, nprofiles
      DO lay = 1, nlayers
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
          prof(i)%cosbk**2_jpim * predictors_ad(i)%mixedgas(11, lay)
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) +     &
          2.0_jprb * prof(i)%Be * raytracing%pathsat(lay,i) * predictors_ad(i)%mixedgas(12, lay)
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
          prof(i)%Be**3_jpim * predictors_ad(i)%mixedgas(13, lay)
        raytracing_ad%pathsat(lay,i) = raytracing_ad%pathsat(lay,i) + &
          2.0_jprb * (prof(i)%cosbk * prof(i)%Be)**2_jpim * &
          raytracing%pathsat(lay,i) * predictors_ad(i)%mixedgas(14, lay)
      ENDDO
    ENDDO
  ENDIF

  ! ---------------------------------------------------------------------------
  ! Water vapour
  ! ---------------------------------------------------------------------------

  DO i = 1, nprofiles
    IF (coef%fmv_model_ver == 7) THEN
      sec_gr_sqrt(:) = path_sqrt(:,i) * aux%wr_sqrt(:,i)
      sec_gr(:) = sec_gr_sqrt(:)**2_jpim
      sec_gw(:) = path(:,i) * aux%ww(:,i)

      DO lay = 1, nlayers
        aux_ad%wr(lay,i) = aux_ad%wr(lay,i) + &
          sec_gr(lay) * aux%tr_r(lay,i) *         predictors_ad(i)%watervapour(14, lay) + &
          sec_gr(lay) * aux%tr_r(lay,i)**4_jpim * predictors_ad(i)%watervapour(15, lay)

        aux_ad%wrw_r(lay,i) = aux_ad%wrw_r(lay,i) + &
          sec_gr(lay) * predictors_ad(i)%watervapour(3, lay) + &
          sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(8, lay)

        aux_ad%tr_r(lay,i) = aux_ad%tr_r(lay,i) + &
          sec_gr(lay) * aux%wr(lay,i) *           predictors_ad(i)%watervapour(14, lay) + & 
          4._jprb * aux%tr_r(lay,i)**2_jpim * predictors(i)%watervapour(14, lay) * &
                                                   predictors_ad(i)%watervapour(15, lay)
        sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
          4.0_jprb * sec_gw(lay)**3_jpim * predictors_ad(i)%watervapour(12, lay) + &
          2.0_jprb * sec_gw(lay) * predictors_ad(i)%watervapour(13, lay)

        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
          sec_gr(lay) *                                 predictors_ad(i)%watervapour(4, lay)  + & !4
          sec_gr_sqrt(lay) *                            predictors_ad(i)%watervapour(6, lay) 

        aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
          sec_gr(lay) * predictors_ad(i)%watervapour(11, lay)
      ENDDO

      DO lay = 1, nlayers
        sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
                                                      predictors_ad(i)%watervapour(1, lay) + &  !  7
          aux%wrw_r(lay,i) *                         predictors_ad(i)%watervapour(3, lay) + &  ! 12
          aux%dt(lay,i) *                            predictors_ad(i)%watervapour(4, lay) + &  !  4
          2._jprb * sec_gr(lay) *                     predictors_ad(i)%watervapour(5, lay) + &  !  1
          3._jprb * sec_gr(lay)**2_jpim *             predictors_ad(i)%watervapour(9, lay) + &
          4._jprb * sec_gr(lay)**3_jpim *             predictors_ad(i)%watervapour(10, lay)+ &
          aux%dtabsdt(lay,i) *                       predictors_ad(i)%watervapour(11, lay)+ &  ! 10
          aux%wr(lay,i) * aux%tr_r(lay,i) *         predictors_ad(i)%watervapour(14, lay)+ &  ! 14
          aux%wr(lay,i) * aux%tr_r(lay,i)**4_jpim * predictors_ad(i)%watervapour(15, lay)

        sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                               predictors_ad(i)%watervapour(2, lay) + &                  !  5
          aux%dt(lay,i) *                     predictors_ad(i)%watervapour(6, lay) + &!11
          aux%wr(lay,i) * aux%ww_r(lay,i) *  predictors_ad(i)%watervapour(8, lay) !13

        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
          predictors_ad(i)%watervapour(7, lay)
      ENDDO

    ELSEIF (coef%fmv_model_ver >= 8) THEN
      sec_gr_sqrt(:) = path_sqrt(:,i) * aux%wr_sqrt(:,i)
      sec_gr(:) = sec_gr_sqrt(:)**2_jpim

      IF (coef%fmv_model_ver == 9) sec_gw_sqrt(:) = path_sqrt(:,i) * aux%ww_sqrt(:,i)
      sec_gw(:) = path(:,i) * aux%ww(:,i)

      IF (coef%nwvcont > 0) THEN
        sec_wrtr_r(:) = sec_gr(:) * aux%tr_r(:,i)
        sec_wrwrtr_r(:) = sec_wrtr_r(:) * aux%wr(:,i)

        DO lay = 1, nlayers
          sec_wrwrtr_r_ad(lay) = &!sec_wrwrtr_r_ad(lay) + &
            predictors_ad(i)%wvcont(1, lay) + &
            aux%tr_r(lay,i)**3_jpim * predictors_ad(i)%wvcont(2, lay)

          aux_ad%tr_r(lay,i) = aux_ad%tr_r(lay,i) + &
            3._jprb * aux%tr_r(lay,i)**2_jpim * sec_wrwrtr_r(lay) * predictors_ad(i)%wvcont(2, lay) + &
            sec_wrtr_r(lay) * predictors_ad(i)%wvcont(4, lay)

          sec_wrtr_r_ad(lay) = &!sec_wrtr_r_ad(lay) + &
            predictors_ad(i)%wvcont(3, lay) + &
            aux%tr_r(lay,i) * predictors_ad(i)%wvcont(4, lay)
        ENDDO

        aux_ad%wr(:,i) = aux_ad%wr(:,i) + &
          sec_wrtr_r(:) * sec_wrwrtr_r_ad(:)

        sec_wrtr_r_ad(:) = sec_wrtr_r_ad(:) + &
          aux%wr(:,i) * sec_wrwrtr_r_ad(:)

        sec_gr_ad(:) = &!sec_gr_ad(:) + &
           aux%tr_r(:,i) * sec_wrtr_r_ad(:)

        aux_ad%tr_r(:,i) = aux_ad%tr_r(:,i) + &
           sec_gr(:) * sec_wrtr_r_ad(:)
      ELSE
        sec_gr_ad(:) = 0._jprb ! need sec_gr_ad to be initialised
      ENDIF

      DO lay = 1, nlayers
        sec_gr_ad(lay) = sec_gr_ad(lay) + &
          2.0_jprb * sec_gr(lay) * predictors_ad(i)%watervapour(1, lay) + &
          aux%dt(lay,i) * predictors_ad(i)%watervapour(4, lay) + &
          predictors_ad(i)%watervapour(7, lay) + &
          3.0_jprb * sec_gr(lay)**2_jpim * predictors_ad(i)%watervapour(8, lay)

        sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
          predictors_ad(i)%watervapour(2, lay) + &
          2.0_jprb * sec_gw(lay) * predictors_ad(i)%watervapour(3, lay)

        sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
          predictors_ad(i)%watervapour(5, lay)

        sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
          predictors_ad(i)%watervapour(6, lay)
        
        aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                            sec_gr(lay) * predictors_ad(i)%watervapour(4, lay)
      ENDDO
      
      IF (coef%fmv_model_ver == 8) THEN
        DO lay = 1, nlayers
          sec_gr_ad(lay) = sec_gr_ad(lay) + &
             aux%dtabsdt(lay,i) * predictors_ad(i)%watervapour(9, lay) + &
             aux%wrwr_r(lay,i) * predictors_ad(i)%watervapour(11, lay)

          aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(9, lay)

          sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
            aux%dt(lay,i) * predictors_ad(i)%watervapour(10, lay) + &
            aux%wrwr_r(lay,i) * predictors_ad(i)%watervapour(12, lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
            sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(10, lay)

          aux_ad%wrwr_r(lay,i) = aux_ad%wrwr_r(lay,i) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(11, lay) + &
            sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(12, lay)
        ENDDO
      ELSE !v9
        sec_gw_4rt(:) = path_4rt(:,i) * aux%ww_4rt(:,i)

        DO lay = 1, nlayers
          sec_gr_ad(lay) = sec_gr_ad(lay) + &
            4.0_jprb * sec_gr(lay)**3_jpim * predictors_ad(i)%watervapour(9, lay) + &
            aux%dtabsdt(lay,i)  *             predictors_ad(i)%watervapour(10, lay) + &
            aux%wrwr_r(lay,i) *               predictors_ad(i)%watervapour(12, lay) + &
            sec_gr_sqrt(lay) *                 predictors_ad(i)%watervapour(15, lay) + &
            sec_gr_sqrt(lay) * aux%dt(lay,i) * predictors_ad(i)%watervapour(18, lay) + &
            aux%wrw_r(lay,i) *                predictors_ad(i)%watervapour(19, lay)

          sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
            aux%dt(lay,i) * predictors_ad(i)%watervapour(11, lay) + &
            aux%wrwr_r(lay,i) * predictors_ad(i)%watervapour(13, lay) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(15, lay) + &
            aux%wr(lay,i) * predictors_ad(i)%watervapour(17, lay) + &
            sec_gr(lay) * aux%dt(lay,i) * predictors_ad(i)%watervapour(18, lay)

          aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(10, lay)

          aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
            sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(11, lay) + &
            sec_gr(lay) * sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(18, lay)

          aux_ad%wrwr_r(lay,i) = aux_ad%wrwr_r(lay,i) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(12, lay) + &
            sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(13, lay)

          sec_gw_sqrt_ad(lay) = &!sec_gw_sqrt_ad(lay) + &
            sec_gw(lay) * predictors_ad(i)%watervapour(14, lay)

          sec_gw_ad(lay) =  sec_gw_ad(lay) + &
            sec_gw_sqrt(lay) * predictors_ad(i)%watervapour(14, lay) + &
            sec_gw_4rt(lay) * predictors_ad(i)%watervapour(16, lay)

          sec_gw_4rt_ad(lay) = &!sec_gw_4rt_ad(lay) + &
            sec_gw(lay) * predictors_ad(i)%watervapour(16, lay)

          aux_ad%wr(lay,i) = aux_ad%wr(lay,i) + &
             sec_gr_sqrt(lay) * predictors_ad(i)%watervapour(17, lay)

          aux_ad%wrw_r(lay,i) = aux_ad%wrw_r(lay,i) + &
            sec_gr(lay) * predictors_ad(i)%watervapour(19, lay)
        ENDDO

        path_4rt_ad(:,i) = &!path_4rt_ad(:,i) + &
          aux%ww_4rt(:,i) * sec_gw_4rt_ad(:)

        aux_ad%ww_4rt(:,i) = aux_ad%ww_4rt(:,i) + &
          path_4rt(:,i) * sec_gw_4rt_ad(:)
      ENDIF
    ENDIF

    sec_gr_sqrt_ad(:) = sec_gr_sqrt_ad(:) + &
      2._jprb * sec_gr_sqrt(:) * sec_gr_ad(:)
    
    aux_ad%wr_sqrt(:,i) = aux_ad%wr_sqrt(:,i) + &
      path_sqrt(:,i)  * sec_gr_sqrt_ad(:)
    
    path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
      aux%wr_sqrt(:,i) * sec_gr_sqrt_ad(:)
    
    aux_ad%wr_4rt(:,i) = aux_ad%wr_4rt(:,i) + &
      path_4rt(:,i) * sec_gr_4rt_ad(:)
    
    path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
      aux%wr_4rt(:,i) * sec_gr_4rt_ad(:)
    
    IF (coef%fmv_model_ver == 7) THEN
      aux_ad%ww(:,i) = aux_ad%ww(:,i) + &
        path(:,i) * sec_gw_ad(:) 
    
      path_ad(:,i) = path_ad(:,i) + &
        aux%ww(:,i) * sec_gw_ad(:)
    ELSE
      IF (coef%fmv_model_ver == 9) THEN
        aux_ad%ww_sqrt(:,i) = aux_ad%ww_sqrt(:,i) + &
                               path_sqrt(:,i)  * sec_gw_sqrt_ad(:)
    
        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
                             aux%ww_sqrt(:,i) * sec_gw_sqrt_ad(:)
      ENDIF
    
      aux_ad%ww(:,i) = aux_ad%ww(:,i) + &
        path(:,i) * sec_gw_ad(:)
    
      path_ad(:,i) = path_ad(:,i) + & 
        aux%ww(:,i) * sec_gw_ad(:)
    ENDIF

  ENDDO

  ! ---------------------------------------------------------------------------
  ! Ozone
  ! ---------------------------------------------------------------------------

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nprofiles
        sec_gr_sqrt(:)             = path_sqrt(:,i) * aux%or_sqrt(:,i)
        sec_gr(:)                  = sec_gr_sqrt(:)**2_jpim
        sec_gw_sqrt(:)             = path_sqrt(:,i) * aux%ow_sqrt(:,i)
        sec_gw(:)                  = sec_gw_sqrt(:)**2_jpim

        DO lay = 1, nlayers
          sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
                                                                  predictors_ad(i)%ozone(1, lay) + &
            aux%dto(lay,i) *                                      predictors_ad(i)%ozone(3, lay) + &
            2._jprb * sec_gr(lay) *                               predictors_ad(i)%ozone(4, lay) + &
            aux%or(lay,i) * aux%ow(lay,i) *                       predictors_ad(i)%ozone(6, lay) + &
            aux%ow(lay,i) *                                       predictors_ad(i)%ozone(8, lay) + &
            sec_gw_sqrt(lay) *                                    predictors_ad(i)%ozone(9, lay)

          sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                                predictors_ad(i)%ozone(2, lay) + &
            aux%dto(lay,i) *                   predictors_ad(i)%ozone(5, lay) + &
            aux%or(lay,i) * aux%ow_r(lay,i) * predictors_ad(i)%ozone(7, lay)

          aux_ad%or(lay,i) = aux_ad%or(lay,i) + &
            sec_gr(lay) * aux%ow(lay,i) *        predictors_ad(i)%ozone(6, lay) + &
            sec_gr_sqrt(lay) * aux%ow_r(lay,i) * predictors_ad(i)%ozone(7, lay)

          aux_ad%dto(lay,i) = aux_ad%dto(lay,i) + &
            sec_gr(lay) *      predictors_ad(i)%ozone(3, lay) + &
            sec_gr_sqrt(lay) * predictors_ad(i)%ozone(5, lay)

          sec_gw_sqrt_ad(lay) = &! sec_gw_sqrt_ad(lay) + &
            sec_gr(lay) * predictors_ad(i)%ozone(9, lay)

          aux_ad%ow(lay,i) = aux_ad%ow(lay,i) + &
             sec_gr(lay) * aux%or(lay,i) * predictors_ad(i)%ozone(6, lay) + &
             sec_gr(lay) *                  predictors_ad(i)%ozone(8, lay)

          aux_ad%ow_r(lay,i) = aux_ad%ow_r(lay,i) + &
            sec_gr_sqrt(lay) * aux%or(lay,i) * predictors_ad(i)%ozone(7, lay)

          sec_gw_ad(lay) = &! sec_gw_ad(lay) + &
            predictors_ad(i)%ozone(10, lay) + &
            2.0_jprb * sec_gw(lay) * predictors_ad(i)%ozone(11, lay)
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN       
          sec_gw_4rt(:)              = path_4rt(:,i) * aux%ow_4rt(:,i)

          DO lay = 1, nlayers
            sec_gr_ad(lay) = sec_gr_ad(lay) + &
              aux%ow_r(lay,i) * predictors_ad(i)%ozone(12, lay) 

            sec_gw_ad(lay) = sec_gw_ad(lay) + &
              1.75_jprb * sec_gw_sqrt(lay) * sec_gw_4rt(lay) * predictors_ad(i)%ozone(13, lay)

            aux_ad%ow_r(lay,i) = aux_ad%ow_r(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%ozone(12, lay)
            
            aux_ad%ow(lay,i) = aux_ad%ow(lay,i) + &
              2.0_jprb * aux%ow(lay,i) * path_sqrt(lay,i) * aux%dto(lay,i) * predictors_ad(i)%ozone(14, lay)

            aux_ad%dto(lay,i) = aux_ad%dto(lay,i) + &
              path_sqrt(lay,i) * aux%ow(lay,i)**2_jpim * predictors_ad(i)%ozone(14, lay)
            
            aux_ad%tro(lay,i) = aux_ad%tro(lay,i) + &
              3.0_jprb * path(lay,i) * aux%tro(lay,i)**2_jpim * predictors_ad(i)%ozone(15, lay)

            path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
              aux%dto(lay,i) * aux%ow(lay,i)**2_jpim * predictors_ad(i)%ozone(14, lay)

            path_ad(lay,i) = path_ad(lay,i) + &
              aux%tro(lay,i)**3_jpim * predictors_ad(i)%ozone(15, lay)
          ENDDO
        ENDIF

        sec_gr_sqrt_ad(:) = sec_gr_sqrt_ad(:) + &
          2._jprb * sec_gr_sqrt(:) * sec_gr_ad(:)
        sec_gw_sqrt_ad(:) = sec_gw_sqrt_ad(:) + &
          2._jprb * sec_gw_sqrt(:) * sec_gw_ad(:)

        aux_ad%or_sqrt(:,i) = aux_ad%or_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gr_sqrt_ad(:)
        aux_ad%ow_sqrt(:,i) = aux_ad%ow_sqrt(:,i) + &
          path_sqrt(:,i) * sec_gw_sqrt_ad(:)

        path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
          aux%ow_sqrt(:,i) * sec_gw_sqrt_ad(:) + &
          aux%or_sqrt(:,i) * sec_gr_sqrt_ad(:)

      ENDDO
    ELSE ! no user ozone data
      DO i = 1, nprofiles   
        DO lay = 1, nlayers
          path_ad(lay,i) = path_ad(lay,i) + &
            predictors_ad(i)%ozone(1, lay) + &
            aux%dto(lay,i) * predictors_ad(i)%ozone(3, lay) + &
            predictors_ad(i)%ozone(6, lay) + &
            predictors_ad(i)%ozone(8, lay) + &
            1.5_jprb * path_sqrt(lay,i) * predictors_ad(i)%ozone(9, lay) + &
            predictors_ad(i)%ozone(10, lay) + &
            2.0_jprb * path(lay,i) * (predictors_ad(i)%ozone(11, lay) + &
                                       predictors_ad(i)%ozone(4, lay))

          path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
            predictors_ad(i)%ozone(2, lay) + &
            aux%dto(lay,i) * predictors_ad(i)%ozone(5, lay) + &
            predictors_ad(i)%ozone(7, lay)

          aux_ad%dto(lay,i) = aux_ad%dto(lay,i) + &
            path(lay,i) * predictors_ad(i)%ozone(3, lay) + &
            path_sqrt(lay,i) * predictors_ad(i)%ozone(5, lay)
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN       
          DO lay = 1, nlayers
            path_ad(lay,i) = path_ad(lay,i) + &
              predictors_ad(i)%ozone(12, lay) + &
              1.75_jprb * path_sqrt(lay,i) * path_4rt(lay,i) * predictors_ad(i)%ozone(13, lay) + &
              aux%tro(lay,i)**3_jpim * predictors_ad(i)%ozone(15, lay) 

            path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
              aux%dto(lay,i) * predictors_ad(i)%ozone(14, lay)

            aux_ad%dto(lay,i) = aux_ad%dto(lay,i) + &
              path_sqrt(lay,i) * predictors_ad(i)%ozone(14, lay)

            aux_ad%tro(lay,i) = aux_ad%tro(lay,i) + &
              3.0_jprb * path(lay,i) * aux%tro(lay,i)**2_jpim * predictors_ad(i)%ozone(15, lay)
          ENDDO
        ENDIF
        
      ENDDO
    ENDIF
  ENDIF

  IF (coef%fmv_model_ver >= 8) THEN
    DO i = 1, nprofiles

      ! -----------------------------------------------------------------------
      ! CO2
      ! -----------------------------------------------------------------------

      IF (coef%nco2 > 0) THEN
        sec_twr(:) = path(:,i) * aux%twr(:,i)

        DO lay = 1, nlayers
          aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
                                                   predictors_ad(i)%co2(2, lay)

          sec_tr_ad(lay,i) = sec_tr_ad(lay,i) + &
                                                   predictors_ad(i)%co2(3, lay) + &
            aux%tr(lay,i) *                       predictors_ad(i)%co2(4, lay)

          aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
            sec_tr(lay,i) *                        predictors_ad(i)%co2(4, lay) + &
                                                   predictors_ad(i)%co2(5, lay)
          sec_twr_ad(lay) = &!sec_twr_ad(lay) + &
                                                   predictors_ad(i)%co2(7, lay) + &
            aux%tr_sqrt(lay,i) *                  predictors_ad(i)%co2(10, lay)

          aux_ad%twr(lay,i) = aux_ad%twr(lay,i) + &
            3._jprb * aux%twr(lay,i)**2_jpim *  predictors_ad(i)%co2(9, lay)
          
          aux_ad%tr_sqrt(lay,i) = aux_ad%tr_sqrt(lay,i) + &
            sec_twr(lay) *                         predictors_ad(i)%co2(10, lay)                        
            
          IF (opts%rt_all%co2_data) THEN
            path_ad(lay,i) = path_ad(lay,i) + &
              aux%co2r(lay,i) *                     predictors_ad(i)%co2(1, lay) + &
                                                     predictors_ad(i)%co2(6, lay) + &
              2.0_jprb * path(lay,i) * aux%co2w(lay,i) **2_jpim * predictors_ad(i)%co2(8, lay)

            aux_ad%co2r(lay,i) = aux_ad%co2r(lay,i) + &
              path(lay,i) * predictors_ad(i)%co2(1, lay)

            aux_ad%co2w(lay,i) = aux_ad%co2w(lay,i) + &
              2._jprb * path(lay,i)**2_jpim * aux%co2w(lay,i) * predictors_ad(i)%co2(8, lay)
          ELSE
            path_ad(lay,i) = path_ad(lay,i) + &
                                                  predictors_ad(i)%co2(1, lay) + &
                                                  predictors_ad(i)%co2(6, lay) + &
              2.0_jprb * path(lay,i) *           predictors_ad(i)%co2(8, lay)
          ENDIF
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN
          IF (opts%rt_all%co2_data) THEN
            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
              aux%co2r_sqrt(:,i) *                      predictors_ad(i)%co2(11, :)
            
            aux_ad%co2r_sqrt(:,i) = aux_ad%co2r_sqrt(:,i) + &
              path_sqrt(:,i) * predictors_ad(i)%co2(11, :)
          ELSE
            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
                                predictors_ad(i)%co2(11, :)
          ENDIF
                   
          DO lay = 1, nlayers
            path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
              aux%tr2(lay,i) * aux%twr(lay,i)**3_jpim * predictors_ad(i)%co2(14, lay)

            aux_ad%tr(lay,i) = aux_ad%tr(lay,i) + &
              3.0_jprb * aux%tr2(lay,i) * predictors_ad(i)%co2(12, lay) + &
              aux%tr(lay,i) * sec_tr(lay,i) * 2.0_jprb * predictors_ad(i)%co2(13, lay)              

            sec_tr_ad(lay,i) = sec_tr_ad(lay,i) + &
              aux%tr(lay,i)**2_jpim * predictors_ad(i)%co2(13, lay)

            aux_ad%twr(lay,i) = aux_ad%twr(lay,i) + & 
              aux%twr(lay,i) * aux%tr2(lay,i) * (3.0_jprb * path_sqrt(lay,i) * aux%twr(lay,i) * predictors_ad(i)%co2(14, lay) + &
                                                 2.0_jprb * predictors_ad(i)%co2(15, lay))

            aux_ad%tr2(lay,i) = aux_ad%tr2(lay,i) + &
              aux%twr(lay,i)** 2_jpim * (aux%twr(lay,i) * path_sqrt(lay,i) * predictors_ad(i)%co2(14, lay) + &
                                                                             predictors_ad(i)%co2(15, lay))
          ENDDO
        ENDIF
        aux_ad%twr(:,i) = aux_ad%twr(:,i) + &
          path(:,i) * sec_twr_ad(:)

        path_ad(:,i) = path_ad(:,i) + &
          aux%twr(:,i) * sec_twr_ad(:)
      ENDIF

      IF (coef%fmv_model_ver == 9) THEN

        ! ---------------------------------------------------------------------
        ! N2O
        ! ---------------------------------------------------------------------

        IF (coef%nn2o > 0) THEN
          sec_gwr(:) = path(:,i) * aux%n2owr(:,i)

          IF (opts%rt_all%n2o_data) THEN
            sec_gr_sqrt(:)    = path_sqrt(:,i) * aux%n2or_sqrt(:,i)
            sec_gr(:)         = path(:,i) * aux%n2or(:,i)

            DO lay = 1, nlayers
              aux_ad%n2or(lay,i) = aux_ad%n2or(lay,i) + &
                aux%dt(lay,i) * predictors_ad(i)%n2o(5, lay) + &
                sec_gr_sqrt(lay) * aux%n2ow_r(lay,i) * predictors_ad(i)%n2o(10, lay)

              aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                                  aux%n2or(lay,i) * predictors_ad(i)%n2o(5, lay)

              aux_ad%n2ow(lay,i) = aux_ad%n2ow(lay,i) + &
                path(lay,i) * predictors_ad(i)%n2o(7, lay)

              path_ad(lay,i) = path_ad(lay,i) + &
                aux%n2ow(lay,i) * predictors_ad(i)%n2o(7, lay)

              aux_ad%n2ow_r(lay,i) = aux_ad%n2ow_r(lay,i) + &
                sec_gr_sqrt(lay) * aux%n2or(lay,i) * predictors_ad(i)%n2o(10, lay)

              sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                                    aux%n2or(lay,i) * aux%n2ow_r(lay,i) * predictors_ad(i)%n2o(10, lay)
            ENDDO
          ELSE
            sec_gr_sqrt    = path_sqrt(:,i)
            sec_gr         = path(:,i)

            aux_ad%dt(:,i) = aux_ad%dt(:,i) + predictors_ad(i)%n2o(5, :)

            path_ad(:,i) = path_ad(:,i) + &
              predictors_ad(i)%n2o(7, :)

            sec_gr_sqrt_ad(:) = &!sec_gr_sqrt_ad(:) + &
                                  predictors_ad(i)%n2o(10, :)
          ENDIF

          DO lay = 1, nlayers
            path_ad(lay,i) = path_ad(lay,i) + &
              sec_gwr(lay) * aux%dt(lay,i) * predictors_ad(i)%n2o(13, lay)

            sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
              predictors_ad(i)%n2o(1, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%n2o(3, lay) + &
              2.0_jprb * sec_gr(lay) * predictors_ad(i)%n2o(4, lay)

            sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
              predictors_ad(i)%n2o(2, lay)

            sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
              predictors_ad(i)%n2o(6, lay)
            
            aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%n2o(3, lay) + &
              path(lay,i) * sec_gwr(lay)  * predictors_ad(i)%n2o(13, lay)

            sec_gwr_ad(lay) = &!sec_gwr_ad(lay) + &
              predictors_ad(i)%n2o(8, lay) + &
              2.0_jprb * sec_gwr(lay) * predictors_ad(i)%n2o(11, lay) + &
              3.0_jprb * sec_gwr(lay)**2_jpim * predictors_ad(i)%n2o(12, lay) + &
              path(lay,i) * aux%dt(lay,i) * predictors_ad(i)%n2o(13, lay)

            aux_ad%n2owr(lay,i) = aux_ad%n2owr(lay,i) + &
              predictors_ad(i)%n2o(9, lay)           
          ENDDO

          IF (opts%rt_all%n2o_data) THEN
            sec_gr_sqrt_ad(:) = sec_gr_sqrt_ad(:) + &
              2.0_jprb * sec_gr_sqrt(:) * sec_gr_ad(:)

            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
              aux%n2or_sqrt(:,i) * sec_gr_sqrt_ad(:)

            aux_ad%n2or_sqrt(:,i) = aux_ad%n2or_sqrt(:,i) + &
              path_sqrt(:,i) * sec_gr_sqrt_ad(:)

            path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
              aux%n2or_4rt(:,i) * sec_gr_4rt_ad(:)

            aux_ad%n2or_4rt(:,i) = aux_ad%n2or_4rt(:,i) + &
              path_4rt(:,i) * sec_gr_4rt_ad(:)
          ELSE
            path_ad(:,i) = path_ad(:,i) + sec_gr_ad
            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
            path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
          ENDIF

          path_ad(:,i) = path_ad(:,i) + &
            aux%n2owr(:,i) * sec_gwr_ad(:)

          aux_ad%n2owr(:,i) = aux_ad%n2owr(:,i) + &
            path(:,i) * sec_gwr_ad(:)
        ENDIF

        ! ---------------------------------------------------------------------
        ! CO
        ! ---------------------------------------------------------------------

        IF (coef%nco > 0) THEN
          sec_gwr(:)          = path(:,i) * aux%cowr(:,i)
          sec_gwr_r(:)        = 1.0_jprb / sec_gwr(:)

          IF (opts%rt_all%co_data) THEN
            sec_gr_sqrt(:)    = path_sqrt(:,i) * aux%cor_sqrt(:,i)
            sec_gr(:)         = sec_gr_sqrt(:)**2_jpim

            DO lay = 1, nlayers
              sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
                aux%corw_r(lay,i) * predictors_ad(i)%co(8, lay) + &
                aux%corw_rsqrt(lay,i) * predictors_ad(i)%co(10, lay) + &
                aux%corw_r4rt(lay,i) * predictors_ad(i)%co(11, lay)

              sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                aux%corw_r(lay,i) * predictors_ad(i)%co(9, lay)

              aux_ad%corw_r(lay,i) = aux_ad%corw_r(lay,i) + &
                 sec_gr(lay) * predictors_ad(i)%co(8, lay) + &
                 sec_gr_sqrt(lay) * predictors_ad(i)%co(9, lay)

              aux_ad%corw_rsqrt(lay,i) = aux_ad%corw_rsqrt(lay,i) + &
                sec_gr(lay) * predictors_ad(i)%co(10, lay)

              aux_ad%corw_r4rt(lay,i) = aux_ad%corw_r4rt(lay,i) + &
                sec_gr(lay) * predictors_ad(i)%co(11, lay)

              path_ad(lay,i) = path_ad(lay,i) + &
                aux%cow(lay,i) *  predictors_ad(i)%co(14, lay)

              aux_ad%cow(lay,i) = aux_ad%cow(lay,i) + &
                path(lay,i) * predictors_ad(i)%co(14, lay)
            ENDDO
          ELSE
            sec_gr         = path(:,i)
            sec_gr_sqrt    = path_sqrt(:,i)

            DO lay = 1, nlayers
              path_ad(lay,i) = path_ad(lay,i) + &
                predictors_ad(i)%co(8, lay) + &
                predictors_ad(i)%co(10, lay) + &
                predictors_ad(i)%co(11, lay) + &
                predictors_ad(i)%co(14, lay)

              path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                predictors_ad(i)%co(9, lay)
            ENDDO
            sec_gr_ad      = 0.0_jprb
            sec_gr_sqrt_ad = 0.0_jprb
          ENDIF

          DO lay = 1, nlayers
            sec_gr_ad(lay) = sec_gr_ad(lay) + &
              predictors_ad(i)%co(1, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%co(3, lay) + &
              2.0_jprb * sec_gr(lay) * predictors_ad(i)%co(4, lay) + &
              aux%dtabsdt(lay,i) * predictors_ad(i)%co(7, lay)
            
            sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
              predictors_ad(i)%co(2, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%co(5, lay)

            sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
              predictors_ad(i)%co(6, lay)

            aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%co(3, lay) + &
              sec_gr_sqrt(lay) * predictors_ad(i)%co(5, lay)

            aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%co(7, lay)

            aux_ad%cowr_4rt(lay,i) = aux_ad%cowr_4rt(lay,i) + &
              path_4rt(lay,i) * predictors_ad(i)%co(13, lay)

            path_4rt_ad(lay,i) = path_4rt_ad(lay,i) + &
                                 aux%cowr_4rt(lay,i) * predictors_ad(i)%co(13, lay)
          ENDDO

          IF (opts%rt_all%co_data) THEN           
            sec_gr_sqrt_ad(:) = sec_gr_sqrt_ad(:) + &
              2.0_jprb * sec_gr_sqrt(:) * sec_gr_ad(:)

            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
              aux%cor_sqrt(:,i) * sec_gr_sqrt_ad(:)

            aux_ad%cor_sqrt(:,i) = aux_ad%cor_sqrt(:,i) + &
              path_sqrt(:,i) * sec_gr_sqrt_ad(:)

            path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
              aux%cor_4rt(:,i) * sec_gr_4rt_ad(:)
                      
            aux_ad%cor_4rt(:,i) = aux_ad%cor_4rt(:,i) + &
              path_4rt(:,i) * sec_gr_4rt_ad(:)
          ELSE
            path_ad(:,i) = path_ad(:,i) + sec_gr_ad
            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
            path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
          ENDIF

          sec_gwr_ad(:) = &!sec_gwr_ad(:) + &
            0.4_jprb * predictors(i)%co(12, :) * sec_gwr_r(:) * predictors_ad(i)%co(12, :)

          path_ad(:,i)  = path_ad(:,i) + &
            aux%cowr(:,i) * sec_gwr_ad(:)

          aux_ad%cowr(:,i) = aux_ad%cowr(:,i) + &
             path(:,i) * sec_gwr_ad(:)
        ENDIF

        ! ---------------------------------------------------------------------
        ! CH4
        ! ---------------------------------------------------------------------

        IF (coef%nch4 > 0) THEN
          IF (opts%rt_all%ch4_data) THEN
            sec_gr_sqrt(:)    = path_sqrt(:,i) * aux%ch4r_sqrt(:,i)
            sec_gr(:)         = path(:,i) * aux%ch4r(:,i)
            sec_gw(:)         = path(:,i) * aux%ch4w(:,i)

            DO lay = 1, nlayers
              sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
                aux%ch4rw_r(lay,i) * predictors_ad(i)%ch4(11, lay)

              aux_ad%ch4r(lay,i) = aux_ad%ch4r(lay,i) + &
                aux%dt(lay,i) * predictors_ad(i)%ch4(5, lay)

              aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                aux%ch4r(lay,i) * predictors_ad(i)%ch4(5, lay)

              aux_ad%ch4rw_r(lay,i) = aux_ad%ch4rw_r(lay,i) + &
                sec_gr_sqrt(lay) * predictors_ad(i)%ch4(11, lay)
            ENDDO
          ELSE
            sec_gr_sqrt    = path_sqrt(:,i)
            sec_gr         = path(:,i)
            sec_gw         = path(:,i)

            DO lay = 1, nlayers
              path_sqrt_ad(lay,i) = path_sqrt_ad(lay,i) + &
                predictors_ad(i)%ch4(11, lay)

              aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
                                 predictors_ad(i)%ch4(5, lay)
            ENDDO
            sec_gr_sqrt_ad = 0.0_jprb
          ENDIF

          DO lay = 1, nlayers
            sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
              predictors_ad(i)%ch4(1, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%ch4(3, lay) + &
              2.0_jprb * sec_gr(lay) * predictors_ad(i)%ch4(4, lay)

            path_ad(lay,i) = path_ad(lay,i) + &
              aux%ch4wr(lay,i) * predictors_ad(i)%ch4(7, lay)
              
            sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
              predictors_ad(i)%ch4(2, lay)
              
            sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
              predictors_ad(i)%ch4(6, lay)

            sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
              2.0_jprb * sec_gw(lay) * predictors_ad(i)%ch4(9, lay) + &
              predictors_ad(i)%ch4(10, lay)

            aux_ad%ch4wr(lay,i) = aux_ad%ch4wr(lay,i) + &
              path(lay,i) * predictors_ad(i)%ch4(7, lay) + &
              predictors_ad(i)%ch4(8, lay)

            aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%ch4(3, lay)
          ENDDO

          IF (opts%rt_all%ch4_data) THEN
! using sec_gr = path * ch4r rather than sec_g_sqrt**2
            path_ad(:,i) = path_ad(:,i) + &
              aux%ch4w(:,i) * sec_gw_ad(:) + &
              aux%ch4r(:,i) * sec_gr_ad(:)

            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
              aux%ch4r_sqrt(:,i) * sec_gr_sqrt_ad(:)

            path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
              aux%ch4r_4rt(:,i) * sec_gr_4rt_ad(:)

            aux_ad%ch4w(:,i) = aux_ad%ch4w(:,i) + &
              path(:,i) * sec_gw_ad(:)

            aux_ad%ch4r(:,i) = aux_ad%ch4r(:,i) + &
              path(:,i) * sec_gr_ad(:)

            aux_ad%ch4r_sqrt(:,i) = aux_ad%ch4r_sqrt(:,i) + &
              path_sqrt(:,i) * sec_gr_sqrt_ad(:)

            aux_ad%ch4r_4rt(:,i) = aux_ad%ch4r_4rt(:,i) + &
              path_4rt(:,i) * sec_gr_4rt_ad(:)
          ELSE
            path_ad(:,i) = path_ad(:,i) + sec_gr_ad + sec_gw_ad
            path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + sec_gr_sqrt_ad
            path_4rt_ad(:,i) = path_4rt_ad(:,i) + sec_gr_4rt_ad
          ENDIF
        ENDIF

        ! ---------------------------------------------------------------------
        ! SO2
        ! ---------------------------------------------------------------------

        IF (coef%nso2 > 0) THEN
          sec_gr_sqrt(:) = path_sqrt(:,i) * aux%so2r_sqrt(:,i)
          sec_gr(:) = sec_gr_sqrt(:)**2_jpim

          sec_gw_sqrt(:) = path_sqrt(:,i) * aux%so2w_sqrt(:,i)
          sec_gw_4rt(:) = path_4rt(:,i) * aux%so2w_4rt(:,i)
          sec_gw(:) = sec_gw_sqrt(:)**2_jpim

          DO lay=1, nlayers
            aux_ad%so2rwr_r(lay,i) = aux_ad%so2rwr_r(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%so2(12, lay) + &
              sec_gr_sqrt(lay) * predictors_ad(i)%so2(13, lay)

            aux_ad%so2rw_r(lay,i) = aux_ad%so2rw_r(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%so2(19, lay)

            sec_gr_ad(lay) = &!sec_gr_ad(lay) + &
              aux%so2rwr_r(lay,i) * predictors_ad(i)%so2(12, lay) + &
              sec_gr_sqrt(lay) * predictors_ad(i)%so2(15, lay) + &
              aux%so2rw_r(lay,i) * predictors_ad(i)%so2(19, lay)

            sec_gr_sqrt_ad(lay) = &!sec_gr_sqrt_ad(lay) + &
              aux%so2rwr_r(lay,i) * predictors_ad(i)%so2(13, lay) + &
              sec_gr(lay) * predictors_ad(i)%so2(15, lay) + &
              aux%so2r(lay,i) * predictors_ad(i)%so2(17, lay)

            sec_gw_ad(lay) = &!sec_gw_ad(lay) + &
              sec_gw_sqrt(lay) * predictors_ad(i)%so2(14, lay) + &
              sec_gw_4rt(lay) * predictors_ad(i)%so2(16, lay)

            sec_gw_sqrt_ad(lay) = &!sec_gw_sqrt_ad(lay) + &
              sec_gw(lay) * predictors_ad(i)%so2(14, lay)

            sec_gw_4rt_ad(lay) = &!sec_gw_4rt_ad(lay) + &
              sec_gw(lay) * predictors_ad(i)%so2(16, lay)

            aux_ad%so2r(lay,i) = aux_ad%so2r(lay,i) + &
              sec_gr_sqrt(lay) * predictors_ad(i)%so2(17, lay)

            sec_gr_ad(lay) = sec_gr_ad(lay) + &
              2.0_jprb * sec_gr(lay) * predictors_ad(i)%so2(1, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%so2(4, lay) + &
              predictors_ad(i)%so2(7, lay) + &
              3.0_jprb * sec_gr(lay)**2_jpim * predictors_ad(i)%so2(8, lay) + &
              4.0_jprb * sec_gr(lay)**3_jpim * predictors_ad(i)%so2(9, lay) + &
              aux%dtabsdt(lay,i) * predictors_ad(i)%so2(10, lay) + &
              sec_gr_sqrt(lay) * aux%dt(lay,i) * predictors_ad(i)%so2(18, lay)

            sec_gw_ad(lay) = sec_gw_ad(lay) + &
              predictors_ad(i)%so2(2, lay) + &
               2.0_jprb * sec_gw(lay) * predictors_ad(i)%so2(3, lay) 
               
            aux_ad%dt(lay,i) = aux_ad%dt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%so2(4, lay) + &
              sec_gr_sqrt(lay) * predictors_ad(i)%so2(11, lay) + &
              sec_gr(lay) * sec_gr_sqrt(lay) * predictors_ad(i)%so2(18, lay)

            aux_ad%dtabsdt(lay,i) = aux_ad%dtabsdt(lay,i) + &
              sec_gr(lay) * predictors_ad(i)%so2(10, lay)

            sec_gr_sqrt_ad(lay) = sec_gr_sqrt_ad(lay) + &
              predictors_ad(i)%so2(5, lay) + &
              aux%dt(lay,i) * predictors_ad(i)%so2(11, lay) + &
              aux%dt(lay,i) * sec_gr(lay) * predictors_ad(i)%so2(18, lay)

            sec_gr_4rt_ad(lay) = &!sec_gr_4rt_ad(lay) + &
              predictors_ad(i)%so2(6, lay)
          ENDDO

          sec_gw_sqrt_ad(:) = sec_gw_sqrt_ad(:) + &
            2._jprb * sec_gw_sqrt(:) * sec_gw_ad(:)
          
          sec_gr_sqrt_ad(:) = sec_gr_sqrt_ad(:) + &
            2._jprb * sec_gr_sqrt(:) * sec_gr_ad(:)
          
          path_sqrt_ad(:,i) = path_sqrt_ad(:,i) + &
            aux%so2w_sqrt(:,i) * sec_gw_sqrt_ad(:) + &
            aux%so2r_sqrt(:,i) * sec_gr_sqrt_ad(:)
        
          path_4rt_ad(:,i) = path_4rt_ad(:,i) + &
            aux%so2w_4rt(:,i) * sec_gw_4rt_ad(:) + &
            aux%so2r_4rt(:,i) * sec_gr_4rt_ad(:)

          aux_ad%so2w_sqrt(:,i) = aux_ad%so2w_sqrt(:,i) + &
            path_sqrt(:,i) * sec_gw_sqrt_ad(:)
          
          aux_ad%so2r_sqrt(:,i) = aux_ad%so2r_sqrt(:,i) + &
            path_sqrt(:,i) * sec_gr_sqrt_ad(:)
          
          aux_ad%so2w_4rt(:,i) = aux_ad%so2w_4rt(:,i) + &
            path_4rt(:,i) * sec_gw_4rt_ad(:)
          
          aux_ad%so2r_4rt(:,i) = aux_ad%so2r_4rt(:,i) + &
            path_4rt(:,i) * sec_gr_4rt_ad(:)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  ! DAR: gases should be done in reverse order according to AD coding theory but can almost be treated as being independent,
  ! making it easier to initialise arrays with calculations when the gas is mandatory - hence why we do mixedgas first but
  ! this leads to problems when considering some temperature AD variables which have dependencies on multiple gases).
  ! These need to be treated in the correct order
  path_ad(:, :) = path_ad(:, :) + aux%tr(:,:) * sec_tr_ad(:,:)
  aux_ad%tr(:, :) = aux_ad%tr(:, :) + path(:,:) * sec_tr_ad(:,:)

  ! -------------------------------------------------------------------------
  ! Pressure-modulated cell (PMC)
  ! -------------------------------------------------------------------------

  IF (coef%pmc_shift) THEN

    ! FWD
    ! Constants
    Lcel_cm=coef%pmc_lengthcell
    ! Tcel=coef%pmc_tempcell
    betaplus1=coef%pmc_betaplus1
    DO ichan = 1, coef%fmv_chn
      Pnom(ichan)=coef%pmc_pnominal(ichan)
      Pcel(ichan)=coef%pmc_ppmc(ichan)
    ENDDO
    nlayers = coef%pmc_nlay

    ! This is on coef levels so pressure is not an active variable in the AD

    DO iprof = 1, nprofiles
      DO lay = 1, nlayers
        lev = lay + 1

        acm = aux%co2_cm(iprof)/(2._jprb*betaplus1*Lcel_cm)
        DO ichan=1, coef%fmv_chn
          Pcel_Lev= acm*raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev  )**2_jpim + Pcel(ichan)**2_jpim
          Pnom_LevM1=acm*raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev-1)**2_jpim + Pnom(ichan)**2_jpim
          Pcel_LevM1= acm*raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev-1)**2_jpim + Pcel(ichan)**2_jpim
          Pnom_Lev=acm*raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev  )**2_jpim + Pnom(ichan)**2_jpim
          ! AD
          Pcel_Lev_AD =  0._jprb
          Pnom_LevM1_AD = 0._jprb
          Pcel_LevM1_AD =  0._jprb
          Pnom_Lev_AD = 0._jprb
          acm_AD = 0._jprb
          ! for PMC predictor 2
          predictors_AD(iprof)%pmc(2,lay,ichan)=0._jprb
          ! for PMC predictor 1
          Pcel_Lev_AD = Pcel_Lev_AD &
          + predictors_AD(iprof)%pmc(1,lay,ichan) * Pnom_LevM1/(Pcel_Lev*Pnom_LevM1)
          Pnom_LevM1_AD = Pnom_LevM1_AD &
          + predictors_AD(iprof)%pmc(1,lay,ichan) * Pcel_Lev/(Pcel_Lev*Pnom_LevM1)
          Pcel_LevM1_AD = Pcel_LevM1_AD &
          - predictors_AD(iprof)%pmc(1,lay,ichan) * Pnom_Lev/(Pcel_LevM1*Pnom_Lev)
          Pnom_Lev_AD = Pnom_Lev_AD &
          - predictors_AD(iprof)%pmc(1,lay,ichan) * Pcel_LevM1/(Pcel_LevM1*Pnom_Lev)
          predictors_AD(iprof)%pmc(1,lay,ichan) = 0._jprb

          acm_AD = acm_AD &
          + Pnom_Lev_AD * raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev  )**2_jpim
          raytracing_AD%pathsat(lay,iprof) = raytracing_AD%pathsat(lay,iprof) &
          + Pnom_Lev_AD * acm * prof(iprof)%p(lev  )**2_jpim
          Pnom_Lev_AD = 0._jprb

          acm_AD = acm_AD &
          + Pcel_LevM1_AD * raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev-1)**2_jpim
          raytracing_AD%pathsat(lay,iprof) = raytracing_AD%pathsat(lay,iprof) &
          + Pcel_LevM1_AD * acm * prof(iprof)%p(lev-1)**2_jpim
          Pcel_LevM1_AD = 0._jprb

          acm_AD = acm_AD &
          + Pnom_LevM1_AD * raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev-1)**2_jpim
          raytracing_AD%pathsat(lay,iprof) = raytracing_AD%pathsat(lay,iprof) &
          + Pnom_LevM1_AD * acm * prof(iprof)%p(lev-1)**2_jpim
          Pnom_LevM1_AD = 0._jprb

          acm_AD = acm_AD &
          + Pcel_Lev_AD * raytracing%pathsat(lay,iprof)*prof(iprof)%p(lev  )**2_jpim
          raytracing_AD%pathsat(lay,iprof) = raytracing_AD%pathsat(lay,iprof) &
          + Pcel_Lev_AD * acm * prof(iprof)%p(lev  )**2_jpim
          Pcel_Lev_AD = 0._jprb

          aux_AD%co2_cm(iprof) = aux_AD%co2_cm(iprof) &
          + acm_AD/(2._jprb*betaplus1*Lcel_cm)
          acm_AD = 0._jprb
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_789_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_789_ad
