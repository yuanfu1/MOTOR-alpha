! Description:
!> @file
!!   K of predictor calculation.
!!
!> @brief
!!   K of predictor calculation.
!!
!! @details
!!   This subroutine operates on coefficient layers/levels.
!!
!! @param[in]     opts            RTTOV options
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     profiles        profiles on coefficient levels
!! @param[in]     coef            rttov_coef structure
!! @param[in]     aux             coef level auxiliary profile data structure
!! @param[in,out] aux_k           coef level auxiliary profile data perturbations
!! @param[in]     predictors      calculated predictors
!! @param[in,out] predictors_k    calculated predictor perturbations
!! @param[in]     raytracing      RTTOV raytracing structure
!! @param[in,out] raytracing_k    raytracing structure perturbations
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
SUBROUTINE rttov_setpredictors_789_k( &
             opts,         &
             chanprof,     &
             profiles,     &
             coef,         &
             aux,          &
             aux_k,        &
             predictors,   &
             predictors_k, &
             raytracing,   &
             raytracing_k, &
             raypath)

  USE rttov_types, ONLY :  &
        rttov_chanprof,         &
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
  TYPE(rttov_chanprof),         INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),          INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef),             INTENT(IN)    :: coef
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_k
  TYPE(rttov_path_pred),        INTENT(IN)    :: predictors(:)
  TYPE(rttov_path_pred),        INTENT(INOUT) :: predictors_k(:)
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_k
  INTEGER(jpim),                INTENT(IN)    :: raypath
!INTF_END
 
  INTEGER(jpim) :: lev, lay, i, j, k, ichan, prof, nlayers, nchannels
  INTEGER :: last_prof

  REAL(jprb), POINTER :: path(:,:), path_sqrt(:,:), path_4rt(:,:)
  REAL(jprb), POINTER :: path_k(:,:), path_sqrt_k(:,:), path_4rt_k(:,:)

  REAL(jprb) :: ztemp(12,profiles(1)%nlayers), ztemp_1d(profiles(1)%nlayers)

  REAL(jprb) :: sec_tr(profiles(1)%nlayers, SIZE(profiles)), sec_twr(profiles(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r(profiles(1)%nlayers), sec_wrwrtr_r(profiles(1)%nlayers)
  REAL(jprb) :: sec_gr(profiles(1)%nlayers), sec_gr_sqrt(profiles(1)%nlayers)
  REAL(jprb) :: sec_gwr(profiles(1)%nlayers), sec_gwr_r(profiles(1)%nlayers), sec_gw(profiles(1)%nlayers)
  REAL(jprb) :: sec_gw_sqrt(profiles(1)%nlayers), sec_gw_4rt(profiles(1)%nlayers)

  REAL(jprb) :: sec_tr_k(profiles(1)%nlayers, SIZE(chanprof)), sec_twr_k(profiles(1)%nlayers)
  REAL(jprb) :: sec_wrtr_r_k(profiles(1)%nlayers), sec_wrwrtr_r_k(profiles(1)%nlayers)
  REAL(jprb) :: sec_gr_k(profiles(1)%nlayers), sec_gr_sqrt_k(profiles(1)%nlayers), sec_gr_4rt_k(profiles(1)%nlayers)
  REAL(jprb) :: sec_gwr_k(profiles(1)%nlayers), sec_gw_k(profiles(1)%nlayers)
  REAL(jprb) :: sec_gw_sqrt_k(profiles(1)%nlayers), sec_gw_4rt_k(profiles(1)%nlayers)

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
  REAL(jprb) :: acm_K
  REAL(jprb) :: Pcel_Lev_K
  REAL(jprb) :: Pnom_LevM1_K
  REAL(jprb) :: Pcel_LevM1_K
  REAL(jprb) :: Pnom_Lev_K

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_789_K', 0_jpim, ZHOOK_HANDLE)

  IF (raypath == 1) THEN
    path => raytracing%pathsat
    path_sqrt => aux%pathsat_sqrt
    path_4rt => aux%pathsat_4rt
    path_k => raytracing_k%pathsat
    path_sqrt_k => aux_k%pathsat_sqrt
    path_4rt_k => aux_k%pathsat_4rt
  ELSEIF (raypath == 2) THEN
    path => raytracing%patheff
    path_sqrt => aux%patheff_sqrt
    path_4rt => aux%patheff_4rt
    path_k => raytracing_k%patheff
    path_sqrt_k => aux_k%patheff_sqrt
    path_4rt_k => aux_k%patheff_4rt
  ENDIF

  nlayers = profiles(1)%nlayers
  nchannels = SIZE(chanprof)

  ! ---------------------------------------------------------------------------
  ! Mixed gases
  ! ---------------------------------------------------------------------------

  sec_tr = path(:, :) * aux%tr(:, :)

  IF (coef%fmv_model_ver <= 8) THEN

    path_4rt_k = 0.0_jprb
    IF (coef%fmv_model_ver == 8) aux_k%tr_sqrt = 0.0_jprb
    
    last_prof = -1
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      
      DO lay = 1, nlayers
      ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22

        path_k(lay,i) = path_k(lay,i) + &
                                                                  predictors_k(i)%mixedgas(1, lay) + &
                          2.0_jprb * path(lay, prof) *            predictors_k(i)%mixedgas(2, lay) + &
                          aux%tw(lay, prof) *                     predictors_k(i)%mixedgas(7, lay) + &
                          aux%tw(lay, prof) * aux%tr_r(lay, prof) * predictors_k(i)%mixedgas(8, lay)

        sec_tr_k(lay,i) = &!sec_tr_k(lay,i) + &          
                            predictors_k(i)%mixedgas(3, lay) + &
                            aux%tr(lay, prof) * predictors_k(i)%mixedgas(4, lay)
 
        aux_k%tr(lay,i) = aux_k%tr(lay,i) + &
                            sec_tr(lay, prof) * predictors_k(i)%mixedgas(4, lay) + &
                            predictors_k(i)%mixedgas(5, lay)

        aux_k%tr2(lay,i) = aux_k%tr2(lay,i) + &
                                     predictors_k(i)%mixedgas(6, lay)

        aux_k%tw(lay,i) = aux_k%tw(lay,i) + &
                                    path(lay, prof) * (                       predictors_k(i)%mixedgas(7, lay) + &
                                                    aux%tr_r(lay, prof) *     predictors_k(i)%mixedgas(8, lay))

        aux_k%tr_r(lay,i) = aux_k%tr_r(lay,i) + &
                                      path(lay, prof) * aux%tw(lay, prof) * predictors_k(i)%mixedgas(8, lay)
  
        path_sqrt_k(lay,i) = &!path_sqrt_k(lay,i) + &
                                       predictors_k(i)%mixedgas(9, lay) + &
                                       aux%tw_4rt(lay, prof) * predictors_k(i)%mixedgas(10, lay)
        aux_k%tw_4rt(lay,i) = aux_k%tw_4rt(lay,i) + &
                                        path_sqrt(lay, prof) * predictors_k(i)%mixedgas(10, lay)
      ENDDO
    ENDDO

  ELSE !v9  
    last_prof = -1
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (last_prof .NE. prof) THEN
        DO lay = 1, nlayers
          ztemp(1,lay) = 1.5_jprb * aux%tr_sqrt(lay, prof) * path_sqrt(lay, prof)
          ztemp(2,lay) = path_sqrt(lay, prof) * path(lay, prof)
        ENDDO
        last_prof = prof
      ENDIF

      DO lay = 1, nlayers
        path_k(lay,i) = path_k(lay,i) + &
                                                           predictors_k(i)%mixedgas(1, lay) + &
                          2.0_jprb * path(lay, prof) *     predictors_k(i)%mixedgas(2, lay) + &
                          aux%tuw(lay, prof) *           ( predictors_k(i)%mixedgas(7, lay) + &
                                                           predictors_k(i)%mixedgas(8, lay)) + &
                          ztemp(1,lay) *                   predictors_k(i)%mixedgas(10, lay)

        sec_tr_k(lay,i) = &!sec_tr_k(lay,i) + &
                                                           predictors_k(i)%mixedgas(3, lay) + &
                          aux%tr(lay, prof) *              predictors_k(i)%mixedgas(4, lay) + &
                          aux%tr2(lay, prof) *             predictors_k(i)%mixedgas(9, lay)

        aux_k%tr_sqrt(lay,i) = aux_k%tr_sqrt(lay,i) + & 
                          ztemp(2,lay) *                   predictors_k(i)%mixedgas(10, lay)

        aux_k%tr(lay,i) = aux_k%tr(lay,i) + &
                          sec_tr(lay, prof) *              predictors_k(i)%mixedgas(4, lay) + &
                                                           predictors_k(i)%mixedgas(5, lay)

        aux_k%tr2(lay,i) = aux_k%tr2(lay,i) + &
                                                           predictors_k(i)%mixedgas(6, lay) + &
                          sec_tr(lay, prof) *              predictors_k(i)%mixedgas(9, lay)

        aux_k%tuw(lay,i) = aux_k%tuw(lay,i) + &
                          path(lay,prof) * (               predictors_k(i)%mixedgas(7, lay) + &
                                                           predictors_k(i)%mixedgas(8, lay))
      ENDDO
    ENDDO
  ENDIF

  IF (coef%id_inst == inst_id_ssmis .AND. coef%inczeeman) THEN

    aux_k%t_layer = 0.0_jprb

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      DO lay = 1, nlayers
        ! SSMIS with Zeeman coefficient file
        ! geomagnetic field variables (Be, cosbk) are part of the user input

        ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        ! NB require prof(i) % Be >0. (divisor)
        
        ! X11 -> X21
        ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        predictors_k(i)%mixedgas(20, lay) = predictors_k(i)%mixedgas(20, lay) + &
          predictors_k(i)%mixedgas(21, lay) * profiles(prof)%Be
        predictors_k(i)%mixedgas(13, lay) = predictors_k(i)%mixedgas(13, lay) + &
          predictors_k(i)%mixedgas(20, lay) * profiles(prof)%Be

        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
             profiles(prof)%Be**3_jpim * predictors_k(i)%mixedgas(19, lay)
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
             profiles(prof)%Be * predictors_k(i)%mixedgas(18, lay)
        predictors_k(i)%mixedgas(16, lay) = predictors_k(i)%mixedgas(16, lay) + &
             predictors_k(i)%mixedgas(17, lay) / profiles(prof)%Be
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
             predictors_k(i)%mixedgas(16, lay) / profiles(prof)%Be
        predictors_k(i)%mixedgas(12, lay) = predictors_k(i)%mixedgas(12, lay) + &
             predictors_k(i)%mixedgas(15, lay) * profiles(prof)%cosbk**2_jpim
        predictors_k(i)%mixedgas(12, lay) = predictors_k(i)%mixedgas(12, lay) + &
             predictors_k(i)%mixedgas(14, lay) / profiles(prof)%Be

        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
             profiles(prof)%cosbk**2_jpim * predictors_k(i)%mixedgas(13, lay)

        aux_k%t_layer(lay,i) = aux_k%t_layer(lay,i) - &
            (predictors(prof)%mixedgas(12, lay) / &
             aux%t_layer(lay, prof)) * predictors_k(i)%mixedgas(12, lay)

        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
             (300.0_jprb / aux%t_layer(lay, prof)) * predictors_k(i)%mixedgas(12, lay)
        raytracing_k%pathsat(lay,i) =  raytracing_k%pathsat(lay,i) + predictors_k(i)%mixedgas(11, lay)
      ENDDO
    ENDDO
  ELSEIF (coef%id_inst == inst_id_amsua .AND. coef%inczeeman) THEN
    ! AMSU-A with Zeeman coefficient file
    ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
    ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      DO lay = 1, nlayers
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
          profiles(prof)%cosbk**2_jpim * predictors_k(i)%mixedgas(11, lay)
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) +     &
          2.0_jprb * profiles(prof)%Be * raytracing%pathsat(lay,prof) * predictors_k(i)%mixedgas(12, lay)
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
          profiles(prof)%Be**3_jpim * predictors_k(i)%mixedgas(13, lay)
        raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
          2.0_jprb * (profiles(prof)%cosbk * profiles(prof)%Be)**2_jpim * &
          raytracing%pathsat(lay,prof) * predictors_k(i)%mixedgas(14, lay)
      ENDDO
    ENDDO
  ENDIF

  ! ---------------------------------------------------------------------------
  ! Water vapour
  ! ---------------------------------------------------------------------------

  last_prof = -1
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    
    IF (coef%fmv_model_ver == 7) THEN

      IF (prof .NE. last_prof) THEN
        sec_gr_sqrt(:) = path_sqrt(:, prof) * aux%wr_sqrt(:, prof)
        sec_gr(:) = sec_gr_sqrt(:)**2_jpim
        sec_gw(:) = path(:, prof) * aux%ww(:, prof)

        ztemp(1,:) = sec_gr(:) * aux%tr_r(:, prof)
        ztemp(2,:) = sec_gr(:) * aux%tr_r(:, prof)**4_jpim
        ztemp(3,:) = sec_gr(:) * aux%wr(:, prof)
        ztemp(4,:) = 4.0_jprb * aux%tr_r(:, prof)**2_jpim * predictors(prof)%watervapour(14, :)
        ztemp(5,:) = 4.0_jprb * sec_gw(:)**3_jpim
        ztemp(6,:) = 2.0_jprb * sec_gw(:)
        ztemp(7,:) = 3.0_jprb * sec_gr(:)**2_jpim
        ztemp(8,:) = 4.0_jprb * sec_gr(:)**3_jpim
        ztemp(9,:) = aux%wr(:, prof) * aux%tr_r(:, prof)
        ztemp(10,:) = aux%wr(:, prof) * aux%tr_r(:, prof)**4_jpim
        ztemp(11,:) = aux%wr(:, prof) * aux%ww_r(:, prof)

        last_prof = prof
      ENDIF

      DO lay = 1, nlayers
        aux_k%wr(lay,i) = aux_k%wr(lay,i) + &
          ztemp(1,lay) *           predictors_k(i)%watervapour(14, lay) + &
          ztemp(2,lay) *           predictors_k(i)%watervapour(15, lay)

        aux_k%wrw_r(lay,i) = aux_k%wrw_r(lay,i) + &
          sec_gr(lay) *            predictors_k(i)%watervapour(3, lay) + &
          sec_gr_sqrt(lay) *       predictors_k(i)%watervapour(8, lay)

        aux_k%tr_r(lay,i) = aux_k%tr_r(lay,i) + &
          ztemp(3,lay) *           predictors_k(i)%watervapour(14, lay) + & 
          ztemp(4,lay) *           predictors_k(i)%watervapour(15, lay)

        sec_gw_k(lay) = &!sec_gw_k(lay) + &
          ztemp(5,lay) *           predictors_k(i)%watervapour(12, lay) + &
          ztemp(6,lay) *           predictors_k(i)%watervapour(13, lay)

        aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
          sec_gr(lay) *            predictors_k(i)%watervapour(4, lay)  + & !4
          sec_gr_sqrt(lay) *       predictors_k(i)%watervapour(6, lay) 

        aux_k%dtabsdt(lay,i) = aux_k%dtabsdt(lay,i) + &
          sec_gr(lay) *            predictors_k(i)%watervapour(11, lay)
        
      ENDDO

      DO lay = 1, nlayers
        sec_gr_k(lay) = &!sec_gr_k(lay) + &
                                   predictors_k(i)%watervapour(1, lay) + &  !  7
          aux%wrw_r(lay, prof) *   predictors_k(i)%watervapour(3, lay) + &  ! 12
          aux%dt(lay, prof) *      predictors_k(i)%watervapour(4, lay) + &  !  4
          2.0_jprb * sec_gr(lay) * predictors_k(i)%watervapour(5, lay) + &  !  1
          ztemp(7,lay) *           predictors_k(i)%watervapour(9, lay) + &
          ztemp(8,lay) *           predictors_k(i)%watervapour(10, lay) + &
          aux%dtabsdt(lay, prof) * predictors_k(i)%watervapour(11, lay) + &  ! 10
          ztemp(9, lay) *          predictors_k(i)%watervapour(14, lay) + &  ! 14
          ztemp(10,lay) *          predictors_k(i)%watervapour(15, lay)

        sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                                   predictors_k(i)%watervapour(2, lay) + &                  !  5
          aux%dt(lay, prof) *      predictors_k(i)%watervapour(6, lay) + &!11
          ztemp(11,lay) *          predictors_k(i)%watervapour(8, lay) !13

        sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                   predictors_k(i)%watervapour(7, lay)
      ENDDO

    ELSEIF (coef%fmv_model_ver >= 8) THEN

      IF (prof .NE. last_prof) THEN
        sec_gr_sqrt(:) = path_sqrt(:, prof) * aux%wr_sqrt(:, prof)
        sec_gr(:) = sec_gr_sqrt(:)**2_jpim
        sec_gw(:) = path(:, prof) * aux%ww(:, prof)

        IF (coef%fmv_model_ver == 9) THEN
          sec_gw_sqrt(:) = path_sqrt(:, prof) * aux%ww_sqrt(:, prof)        
          sec_gw_4rt(:) = path_4rt(:, prof) * aux%ww_4rt(:, prof)
        ENDIF
        
        IF (coef%nwvcont > 0) THEN
          sec_wrtr_r(:) = sec_gr(:) * aux%tr_r(:, prof)
          sec_wrwrtr_r(:) = sec_wrtr_r(:) * aux%wr(:, prof)

          ztemp(1,:) = aux%tr_r(:, prof)**3_jpim
          ztemp(2,:) = 3._jprb * aux%tr_r(:, prof)**2_jpim * sec_wrwrtr_r(:)
        ENDIF
        
        ztemp(3,:) = 3.0_jprb * sec_gr(:)**2_jpim
        ztemp(4,:) = 4.0_jprb * sec_gr(:)**3_jpim
        ztemp(5,:) = sec_gr_sqrt(:) * aux%dt(:, prof)
        ztemp(6,:) = sec_gr(:) * aux%dt(:, prof) 
        ztemp(7,:) = sec_gr(:) * sec_gr_sqrt(:)

        last_prof = prof
      ENDIF

      IF (coef%nwvcont > 0) THEN
          
        DO lay = 1, nlayers
          sec_wrwrtr_r_k(lay) = &!sec_wrwrtr_r_k(lay) + &
            predictors_k(i)%wvcont(1, lay) + &
            ztemp(1,lay) * predictors_k(i)%wvcont(2, lay)

          aux_k%tr_r(lay,i) = aux_k%tr_r(lay,i) + &
            ztemp(2,lay) * predictors_k(i)%wvcont(2, lay) + &
            sec_wrtr_r(lay) * predictors_k(i)%wvcont(4, lay)

          sec_wrtr_r_k(lay) = &!sec_wrtr_r_k(lay) + &
            predictors_k(i)%wvcont(3, lay) + &
            aux%tr_r(lay, prof) * predictors_k(i)%wvcont(4, lay)
        ENDDO

        aux_k%wr(:,i) = aux_k%wr(:,i) + &
          sec_wrtr_r(:) * sec_wrwrtr_r_k(:)

        sec_wrtr_r_k(:) = sec_wrtr_r_k(:) + &
          aux%wr(:,prof) * sec_wrwrtr_r_k(:)

        sec_gr_k(:) = &!sec_gr_k(:) + &
           aux%tr_r(:,prof) * sec_wrtr_r_k(:)

        aux_k%tr_r(:,i) = aux_k%tr_r(:,i) + &
           sec_gr(:) * sec_wrtr_r_k(:)
      ELSE
        sec_gr_k(:) = 0._jprb ! need sec_gr_k to be initialised
        ! aux_k%wr(:,i) = 0._jprb
      ENDIF

      DO lay = 1, nlayers
        sec_gr_k(lay) = sec_gr_k(lay) + &
          2.0_jprb * sec_gr(lay) *        predictors_k(i)%watervapour(1, lay) + &
          aux%dt(lay, prof) *             predictors_k(i)%watervapour(4, lay) + &
                                          predictors_k(i)%watervapour(7, lay) + &
          ztemp(3,lay) *                  predictors_k(i)%watervapour(8, lay)

        sec_gw_k(lay) = &!sec_gw_k(lay) + &
                                          predictors_k(i)%watervapour(2, lay) + &
          2.0_jprb * sec_gw(lay) *        predictors_k(i)%watervapour(3, lay)

        sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                                          predictors_k(i)%watervapour(5, lay)

        sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                          predictors_k(i)%watervapour(6, lay)
        
        aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
                            sec_gr(lay) * predictors_k(i)%watervapour(4, lay)
      ENDDO
      
      IF (coef%fmv_model_ver == 8) THEN
        DO lay = 1, nlayers
          sec_gr_k(lay) = sec_gr_k(lay) + &
             aux%dtabsdt(lay, prof) *     predictors_k(i)%watervapour(9, lay) + &
             aux%wrwr_r(lay, prof) *      predictors_k(i)%watervapour(11, lay)

          aux_k%dtabsdt(lay,i) = aux_k%dtabsdt(lay,i) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(9, lay)

          sec_gr_sqrt_k(lay) = sec_gr_sqrt_k(lay) + &
            aux%dt(lay, prof) *           predictors_k(i)%watervapour(10, lay) + &
            aux%wrwr_r(lay, prof) *       predictors_k(i)%watervapour(12, lay)

          aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
            sec_gr_sqrt(lay) *            predictors_k(i)%watervapour(10, lay)

          aux_k%wrwr_r(lay,i) = aux_k%wrwr_r(lay,i) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(11, lay) + &
            sec_gr_sqrt(lay) *            predictors_k(i)%watervapour(12, lay)
        ENDDO
      ELSE !v9
        DO lay = 1, nlayers
          sec_gr_k(lay) = sec_gr_k(lay) + &
            ztemp(4,lay) *                predictors_k(i)%watervapour(9, lay) + &
            aux%dtabsdt(lay, prof)  *     predictors_k(i)%watervapour(10, lay) + &
            aux%wrwr_r(lay, prof) *       predictors_k(i)%watervapour(12, lay) + &
            sec_gr_sqrt(lay) *            predictors_k(i)%watervapour(15, lay) + &
            ztemp(5,lay) *                predictors_k(i)%watervapour(18, lay) + &
            aux%wrw_r(lay, prof) *        predictors_k(i)%watervapour(19, lay)

          sec_gr_sqrt_k(lay) = sec_gr_sqrt_k(lay) + &
            aux%dt(lay, prof) *           predictors_k(i)%watervapour(11, lay) + &
            aux%wrwr_r(lay, prof) *       predictors_k(i)%watervapour(13, lay) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(15, lay) + &
            aux%wr(lay, prof) *           predictors_k(i)%watervapour(17, lay) + &
            ztemp(6,lay) *                predictors_k(i)%watervapour(18, lay)

          aux_k%dtabsdt(lay,i) = aux_k%dtabsdt(lay,i) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(10, lay)

          aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
            sec_gr_sqrt(lay) *            predictors_k(i)%watervapour(11, lay) + &
            ztemp(7,lay) *                predictors_k(i)%watervapour(18, lay)

          aux_k%wrwr_r(lay,i) = aux_k%wrwr_r(lay,i) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(12, lay) + &
            sec_gr_sqrt(lay) *            predictors_k(i)%watervapour(13, lay)

          sec_gw_sqrt_k(lay) = &!sec_gw_sqrt_k(lay) + &
            sec_gw(lay) *                 predictors_k(i)%watervapour(14, lay)

          sec_gw_k(lay) =  sec_gw_k(lay) + &
            sec_gw_sqrt(lay) *            predictors_k(i)%watervapour(14, lay) + &
            sec_gw_4rt(lay) *             predictors_k(i)%watervapour(16, lay)

          sec_gw_4rt_k(lay) = &!sec_gw_4rt_k(lay) + &
            sec_gw(lay) *                 predictors_k(i)%watervapour(16, lay)

          aux_k%wr(lay,i) = aux_k%wr(lay,i) + &
             sec_gr_sqrt(lay) *           predictors_k(i)%watervapour(17, lay)

          aux_k%wrw_r(lay,i) = aux_k%wrw_r(lay,i) + &
            sec_gr(lay) *                 predictors_k(i)%watervapour(19, lay)
        ENDDO

!v8 and v9
        path_4rt_k(:,i) = &!path_4rt_k(:,i) + &
          aux%ww_4rt(:, prof) * sec_gw_4rt_k(:)

        aux_k%ww_4rt(:,i) = aux_k%ww_4rt(:,i) + &
          path_4rt(:, prof) * sec_gw_4rt_k(:)
      ENDIF
    ENDIF

    sec_gr_sqrt_k(:) = sec_gr_sqrt_k(:) + &
      2.0_jprb * sec_gr_sqrt(:) * sec_gr_k(:)
    
    aux_k%wr_sqrt(:,i) = aux_k%wr_sqrt(:,i) + &
      path_sqrt(:, prof)  * sec_gr_sqrt_k(:)
    
    aux_k%wr_4rt(:,i) = aux_k%wr_4rt(:,i) + &
      path_4rt(:, prof) * sec_gr_4rt_k(:)
    
    aux_k%ww(:,i) = aux_k%ww(:,i) + &
      path(:, prof) * sec_gw_k(:)
    
    path_4rt_k(:,i) = path_4rt_k(:,i) + &
      aux%wr_4rt(:, prof) * sec_gr_4rt_k(:)
    
    path_k(:,i) = path_k(:,i) + & 
      aux%ww(:, prof) * sec_gw_k(:)
    
    IF (coef%fmv_model_ver == 9) THEN
      path_sqrt_k(:,i) = &! path_sqrt_k(:,i) + &
        aux%wr_sqrt(:, prof) * sec_gr_sqrt_k(:) + &
        aux%ww_sqrt(:, prof) * sec_gw_sqrt_k(:)
    
      aux_k%ww_sqrt(:,i) = aux_k%ww_sqrt(:,i) + &
        path_sqrt(:, prof)  * sec_gw_sqrt_k(:)
    ELSE
      path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
        aux%wr_sqrt(:, prof) * sec_gr_sqrt_k(:)
    ENDIF
  ENDDO

  ! ---------------------------------------------------------------------------
  ! Ozone
  ! ---------------------------------------------------------------------------

  last_prof = -1
  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF (prof .NE. last_prof) THEN
          sec_gr_sqrt(:)             = path_sqrt(:, prof) * aux%or_sqrt(:, prof)
          sec_gr(:)                  = sec_gr_sqrt(:)**2_jpim
          sec_gw_sqrt(:)             = path_sqrt(:, prof) * aux%ow_sqrt(:, prof)
          sec_gw(:)                  = sec_gw_sqrt(:)**2_jpim

          ztemp(1,:) = aux%or(:, prof) * aux%ow(:, prof)
          ztemp(2,:) = aux%or(:, prof) * aux%ow_r(:, prof)
          ztemp(3,:) = sec_gr(:) * aux%ow(:, prof)
          ztemp(4,:) = sec_gr_sqrt(:) * aux%ow_r(:, prof)
          ztemp(5,:) = sec_gr(:) * aux%or(:, prof)
          ztemp(6,:) = sec_gr_sqrt(:) * aux%or(:, prof)

          IF (coef%fmv_model_ver == 9) THEN       
            ztemp(7,:) = 1.75_jprb * sec_gw_sqrt(:) * path_4rt(:, prof) * aux%ow_4rt(:, prof)
            ztemp(8,:) = 2.0_jprb * aux%ow(:,prof) * path_sqrt(:, prof) * aux%dto(:,prof)
            ztemp(9,:) = path_sqrt(:, prof) * aux%ow(:,prof)**2_jpim
            ztemp(10,:) = 3.0_jprb * path(:, prof) * aux%tro(:,prof)**2_jpim
            ztemp(11,:) = aux%dto(:,prof) * aux%ow(:,prof)**2_jpim
            ztemp(12,:) = aux%tro(:,prof)**3_jpim
          ENDIF
          last_prof = prof
        ENDIF
        
        DO lay = 1, nlayers
          sec_gr_k(lay) = &!sec_gr_k(lay) + &
                                     predictors_k(i)%ozone(1, lay) + &
            aux%dto(lay, prof) *     predictors_k(i)%ozone(3, lay) + &
            2._jprb * sec_gr(lay) *  predictors_k(i)%ozone(4, lay) + &
            ztemp(1,lay) *           predictors_k(i)%ozone(6, lay) + &
            aux%ow(lay, prof) *      predictors_k(i)%ozone(8, lay) + &
            sec_gw_sqrt(lay) *       predictors_k(i)%ozone(9, lay)

          sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                                     predictors_k(i)%ozone(2, lay) + &
            aux%dto(lay, prof) *     predictors_k(i)%ozone(5, lay) + &
            ztemp(2,lay) *           predictors_k(i)%ozone(7, lay)

          aux_k%or(lay,i) = aux_k%or(lay,i) + &
            ztemp(3,lay) *           predictors_k(i)%ozone(6, lay) + &
            ztemp(4,lay) *           predictors_k(i)%ozone(7, lay)

          aux_k%dto(lay,i) = aux_k%dto(lay,i) + &
            sec_gr(lay) *            predictors_k(i)%ozone(3, lay) + &
            sec_gr_sqrt(lay) *       predictors_k(i)%ozone(5, lay)

          sec_gw_sqrt_k(lay) = &! sec_gw_sqrt_k(lay) + &
            sec_gr(lay) *            predictors_k(i)%ozone(9, lay)

          aux_k%ow(lay,i) = aux_k%ow(lay,i) + &
             ztemp(5,lay) *          predictors_k(i)%ozone(6, lay) + &
             sec_gr(lay) *           predictors_k(i)%ozone(8, lay)

          aux_k%ow_r(lay,i) = aux_k%ow_r(lay,i) + &
            ztemp(6,lay) *           predictors_k(i)%ozone(7, lay)

          sec_gw_k(lay) = &! sec_gw_k(lay) + &
                                     predictors_k(i)%ozone(10, lay) + &
            2.0_jprb * sec_gw(lay) * predictors_k(i)%ozone(11, lay)
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN       
          DO lay = 1, nlayers
            sec_gr_k(lay) = sec_gr_k(lay) + &
              aux%ow_r(lay, prof) *  predictors_k(i)%ozone(12, lay) 

            sec_gw_k(lay) = sec_gw_k(lay) + &
              ztemp(7,lay) *         predictors_k(i)%ozone(13, lay)

            aux_k%ow_r(lay,i) = aux_k%ow_r(lay,i) + &
              sec_gr(lay) *          predictors_k(i)%ozone(12, lay)
            
            aux_k%ow(lay,i) = aux_k%ow(lay,i) + &
              ztemp(8,lay) *         predictors_k(i)%ozone(14, lay)

            aux_k%dto(lay,i) = aux_k%dto(lay,i) + &
              ztemp(9,lay) *        predictors_k(i)%ozone(14, lay)
            
            aux_k%tro(lay,i) = aux_k%tro(lay,i) + &
              ztemp(10,lay) *        predictors_k(i)%ozone(15, lay)

            path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
              ztemp(11,lay) *         predictors_k(i)%ozone(14, lay)

            path_k(lay,i) = path_k(lay,i) + &
              ztemp(12,lay) *        predictors_k(i)%ozone(15, lay)
          ENDDO
        ENDIF

        sec_gr_sqrt_k(:) = sec_gr_sqrt_k(:) + &
          2._jprb * sec_gr_sqrt(:) * sec_gr_k(:)
        sec_gw_sqrt_k(:) = sec_gw_sqrt_k(:) + &
          2._jprb * sec_gw_sqrt(:) * sec_gw_k(:)

        aux_k%or_sqrt(:,i) = aux_k%or_sqrt(:,i) + &
          path_sqrt(:, prof) * sec_gr_sqrt_k(:)
        aux_k%ow_sqrt(:,i) = aux_k%ow_sqrt(:,i) + &
          path_sqrt(:, prof) * sec_gw_sqrt_k(:)

        path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
          aux%ow_sqrt(:, prof) * sec_gw_sqrt_k(:) + &
          aux%or_sqrt(:, prof) * sec_gr_sqrt_k(:)
      ENDDO
    ELSE ! no user ozone data
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF (prof .NE. last_prof) THEN
          IF (coef%fmv_model_ver == 9) THEN       
            ztemp(1,:) = 1.75_jprb * path_sqrt(:,prof) * path_4rt(:,prof)
            ztemp(2,:) = aux%tro(:,prof)**3_jpim
            ztemp(3,:) = 3.0_jprb * path(:, prof) * aux%tro(:,prof)**2_jpim
          ENDIF
          last_prof = prof
        ENDIF

        DO lay = 1, nlayers
          path_k(lay,i) = path_k(lay,i) + &
                                              predictors_k(i)%ozone(1, lay) + &
            aux%dto(lay, prof) *              predictors_k(i)%ozone(3, lay) + &
                                              predictors_k(i)%ozone(6, lay) + &
                                              predictors_k(i)%ozone(8, lay) + &
            1.5_jprb * path_sqrt(lay, prof) * predictors_k(i)%ozone(9, lay) + &
                                              predictors_k(i)%ozone(10, lay) + &
            2.0_jprb * path(lay, prof) *    ( predictors_k(i)%ozone(11, lay) + &
                                              predictors_k(i)%ozone(4, lay))

          path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
                                              predictors_k(i)%ozone(2, lay) + &
            aux%dto(lay, prof) *              predictors_k(i)%ozone(5, lay) + &
                                              predictors_k(i)%ozone(7, lay)

          aux_k%dto(lay,i) = aux_k%dto(lay,i) + &
            path(lay, prof) *                 predictors_k(i)%ozone(3, lay) + &
            path_sqrt(lay, prof) *            predictors_k(i)%ozone(5, lay)
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN       
          DO lay = 1, nlayers
            path_k(lay,i) = path_k(lay,i) + &
                                              predictors_k(i)%ozone(12, lay) + &
              ztemp(1,lay) *                  predictors_k(i)%ozone(13, lay) + &
              ztemp(2,lay) *                  predictors_k(i)%ozone(15, lay) 

            path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
              aux%dto(lay,prof) *             predictors_k(i)%ozone(14, lay)

            aux_k%dto(lay,i) = aux_k%dto(lay,i) + &
              path_sqrt(lay, prof) *          predictors_k(i)%ozone(14, lay)

            aux_k%tro(lay,i) = aux_k%tro(lay,i) + &
              ztemp(3,lay) *                  predictors_k(i)%ozone(15, lay)
          ENDDO
        ENDIF
        
      ENDDO
    ENDIF
  ENDIF

  IF (coef%fmv_model_ver >= 8) THEN

    ! -------------------------------------------------------------------------
    ! CO2
    ! -------------------------------------------------------------------------

    last_prof = -1
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (prof .NE. last_prof) THEN
        ztemp(1,:) = 3.0_jprb * aux%twr(:, prof)**2_jpim

        IF (coef%nco2 > 0) THEN
          sec_twr(:) = path(:, prof) * aux%twr(:, prof)
          IF (opts%rt_all%co2_data) THEN
            ztemp(2,:) = 2.0_jprb * path(:, prof) * aux%co2w(:, prof)**2_jpim
            ztemp(3,:) = 2.0_jprb * path(:, prof)**2_jpim * aux%co2w(:, prof)
          ENDIF
        ENDIF

        IF (coef%fmv_model_ver == 9) THEN 
          ztemp(4,:) = 2.0_jprb * aux%tr(:,prof) * sec_tr(:, prof)
          ztemp(5,:) = aux%twr(:,prof)**3_jpim * path_sqrt(:,prof)
          ztemp(6,:) = aux%twr(:,prof)**2_jpim
          ztemp(7,:) = aux%twr(:,prof)**2_jpim * aux%tr2(:,prof) * 3.0_jprb * path_sqrt(:,prof)
          ztemp(8,:) = aux%twr(:,prof) * aux%tr2(:,prof) * 2.0_jprb
          ztemp(9,:) = aux%twr(:,prof)**3_jpim * aux%tr2(:,prof)
        ENDIF

        last_prof = prof
      ENDIF

      IF (coef%nco2 > 0) THEN
        DO lay = 1, nlayers
          aux_k%tr2(lay,i) = aux_k%tr2(lay,i) + &
                                              predictors_k(i)%co2(2, lay)

          sec_tr_k(lay,i) = sec_tr_k(lay,i) +     &
                                              predictors_k(i)%co2(3, lay) + &
            aux%tr(lay, prof) *               predictors_k(i)%co2(4, lay)

          aux_k%tr(lay,i) = aux_k%tr(lay,i) +   &
            sec_tr(lay,prof) *                predictors_k(i)%co2(4, lay) + &
                                              predictors_k(i)%co2(5, lay)
          sec_twr_k(lay) = &!sec_twr_k(lay) + &
                                              predictors_k(i)%co2(7, lay) + &
            aux%tr_sqrt(lay, prof) *          predictors_k(i)%co2(10, lay)

          aux_k%twr(lay,i) = aux_k%twr(lay,i) + &
            ztemp(1,lay) *                    predictors_k(i)%co2(9, lay)
          
          aux_k%tr_sqrt(lay,i) = aux_k%tr_sqrt(lay,i) + &
            sec_twr(lay) *                    predictors_k(i)%co2(10, lay)                        
            
          IF (opts%rt_all%co2_data) THEN
            path_k(lay,i) = path_k(lay,i) + &
              aux%co2r(lay, prof) *           predictors_k(i)%co2(1, lay) + &
                                              predictors_k(i)%co2(6, lay) + &
              ztemp(2,lay) *                  predictors_k(i)%co2(8, lay)

            aux_k%co2r(lay,i) = aux_k%co2r(lay,i) + &
              path(lay, prof) *               predictors_k(i)%co2(1, lay)

            aux_k%co2w(lay,i) = aux_k%co2w(lay,i) + &
              ztemp(3,lay) *                  predictors_k(i)%co2(8, lay)
          ELSE
            path_k(lay,i) = path_k(lay,i) + &
                                              predictors_k(i)%co2(1, lay) + &
                                              predictors_k(i)%co2(6, lay) + &
              2.0_jprb * path(lay, prof) *    predictors_k(i)%co2(8, lay)
          ENDIF
        ENDDO

        IF (coef%fmv_model_ver == 9) THEN
          IF (opts%rt_all%co2_data) THEN
            path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
              aux%co2r_sqrt(:,prof) *         predictors_k(i)%co2(11, :)
            
            aux_k%co2r_sqrt(:,i) = aux_k%co2r_sqrt(:,i) + &
              path_sqrt(:, prof) *            predictors_k(i)%co2(11, :)
          ELSE
            path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
                                              predictors_k(i)%co2(11, :)
          ENDIF
                   
          DO lay = 1, nlayers
            aux_k%tr(lay,i) = aux_k%tr(lay,i) + &
              3.0_jprb * aux%tr2(lay, prof) * predictors_k(i)%co2(12, lay) + &
              ztemp(4,lay) *                  predictors_k(i)%co2(13, lay)              

            aux_k%tr2(lay,i) = aux_k%tr2(lay,i) + &
             ztemp(5,lay) *                   predictors_k(i)%co2(14, lay) + &
             ztemp(6,lay) *                   predictors_k(i)%co2(15, lay)

            aux_k%twr(lay,i) = aux_k%twr(lay,i) + & 
              ztemp(7,lay) *                  predictors_k(i)%co2(14, lay) + &
              ztemp(8,lay) *                  predictors_k(i)%co2(15, lay)

            path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
              ztemp(9,lay) *                  predictors_k(i)%co2(14, lay)

            sec_tr_k(lay,i) = sec_tr_k(lay,i) + &
              aux%tr2(lay,prof) *             predictors_k(i)%co2(13, lay)
          ENDDO
        ENDIF

        aux_k%twr(:,i) = aux_k%twr(:,i) + &
          path(:, prof) * sec_twr_k(:)

        path_k(:,i) = path_k(:,i) + &
          aux%twr(:, prof) * sec_twr_k(:)
      ENDIF
    ENDDO

!the rest of the trace gases come as a package. If one is in, they're all in so we treat them together

    IF (coef%fmv_model_ver == 9) THEN

      ! -----------------------------------------------------------------------
      ! N2O
      ! -----------------------------------------------------------------------

      IF (coef%nn2o > 0) THEN

        last_prof = -1
        DO i = 1, nchannels
          prof = chanprof(i)%prof
          IF (prof .NE. last_prof) THEN
            sec_gwr(:)               = path(:, prof) * aux%n2owr(:,prof)

            IF (opts%rt_all%n2o_data) THEN
              sec_gr_sqrt(:)    = path_sqrt(:, prof) * aux%n2or_sqrt(:, prof)
              sec_gr(:)         = path(:, prof) * aux%n2or(:,prof)

              ztemp(1,:) = sec_gr_sqrt(:) * aux%n2ow_r(:,prof)
              ztemp(2,:) = sec_gr_sqrt(:) * aux%n2or(:,prof)
              ztemp(3,:) = aux%n2or(:,prof) * aux%n2ow_r(:,prof)
            ELSE
              sec_gr_sqrt    = path_sqrt(:, prof)
              sec_gr         = path(:, prof)
            ENDIF

            ztemp(4,:) = 3.0_jprb * sec_gwr(:)**2_jpim
            ztemp(5,:) = path(:, prof) * aux%dt(:,prof)
            ztemp(6,:) = sec_gwr(:) * aux%dt(:,prof)
            ztemp(7,:) = path(:, prof) * sec_gwr(:)

            last_prof = prof
          ENDIF

          IF (opts%rt_all%n2o_data) THEN

            DO lay = 1, nlayers
              sec_gr_k(lay) = &!sec_gr_k(lay) + &
                                         predictors_k(i)%n2o(1, lay) + &
                aux%dt(lay,prof) *       predictors_k(i)%n2o(3, lay) + &
                2.0_jprb * sec_gr(lay) * predictors_k(i)%n2o(4, lay)

              sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                ztemp(3,lay) *           predictors_k(i)%n2o(10, lay) + &
                                         predictors_k(i)%n2o(2, lay)

              sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                         predictors_k(i)%n2o(6, lay)

              aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
                aux%n2or(lay,prof) *     predictors_k(i)%n2o(5, lay)

              aux_k%n2or(lay,i) = aux_k%n2or(lay,i) + &
                aux%dt(lay,prof) *       predictors_k(i)%n2o(5, lay) + &
                ztemp(1,lay) *           predictors_k(i)%n2o(10, lay)

              aux_k%n2ow_r(lay,i) = aux_k%n2ow_r(lay,i) + &
                ztemp(2,lay) *           predictors_k(i)%n2o(10, lay)

              aux_k%n2ow(lay,i) = aux_k%n2ow(lay,i) + &
                path(lay, prof) *        predictors_k(i)%n2o(7, lay)

              path_k(lay,i) = path_k(lay,i) + &
                aux%n2ow(lay,prof) *     predictors_k(i)%n2o(7, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              sec_gr_k(lay) = &!sec_gr_k(lay) + &
                                         predictors_k(i)%n2o(1, lay) + &
              aux%dt(lay,prof) *         predictors_k(i)%n2o(3, lay) + &
              2.0_jprb * sec_gr(lay) *   predictors_k(i)%n2o(4, lay)

              sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                                         predictors_k(i)%n2o(2, lay) + &
                                         predictors_k(i)%n2o(10, lay)

              sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                         predictors_k(i)%n2o(6, lay)
              
              aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
                                         predictors_k(i)%n2o(5, lay)

              path_k(lay,i) = path_k(lay,i) + &
                                         predictors_k(i)%n2o(7, lay)
            ENDDO
          ENDIF

          DO lay = 1, nlayers
            sec_gwr_k(lay) = &!sec_gwr_k(lay) + &
                                         predictors_k(i)%n2o(8, lay) + &
              2.0_jprb * sec_gwr(lay) *  predictors_k(i)%n2o(11, lay) + &
              ztemp(4,lay) *             predictors_k(i)%n2o(12, lay) + &
              ztemp(5,lay) *             predictors_k(i)%n2o(13, lay)

            path_k(lay,i) = path_k(lay,i) + &
              ztemp(6,lay) *             predictors_k(i)%n2o(13, lay)
           
            aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
              sec_gr(lay) *              predictors_k(i)%n2o(3, lay) + &
              ztemp(7,lay) *             predictors_k(i)%n2o(13, lay)

            aux_k%n2owr(lay,i) = aux_k%n2owr(lay,i) + &
                                         predictors_k(i)%n2o(9, lay)           
          ENDDO

          IF (opts%rt_all%n2o_data) THEN
            sec_gr_sqrt_k(:) = sec_gr_sqrt_k(:) + &
              2.0_jprb * sec_gr_sqrt(:) * sec_gr_k(:)

            path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
              aux%n2or_sqrt(:, prof) * sec_gr_sqrt_k(:)

            path_4rt_k(:,i) = path_4rt_k(:,i) + &
              aux%n2or_4rt(:, prof) * sec_gr_4rt_k(:)

            aux_k%n2or_sqrt(:,i) = aux_k%n2or_sqrt(:,i) + &
              path_sqrt(:, prof) * sec_gr_sqrt_k(:)

            aux_k%n2or_4rt(:,i) = aux_k%n2or_4rt(:,i) + &
              path_4rt(:, prof) * sec_gr_4rt_k(:)
          ELSE
            path_k(:,i) = path_k(:,i) + sec_gr_k
            path_sqrt_k(:,i) = path_sqrt_k(:,i) + sec_gr_sqrt_k
            path_4rt_k(:,i) = path_4rt_k(:,i) + sec_gr_4rt_k
          ENDIF

          path_k(:,i) = path_k(:,i) + &
            aux%n2owr(:,prof) * sec_gwr_k(:)

          aux_k%n2owr(:,i) = aux_k%n2owr(:,i) + &
            path(:, prof) * sec_gwr_k(:)
        ENDDO
      ENDIF

      ! -----------------------------------------------------------------------
      ! CO
      ! -----------------------------------------------------------------------

      IF (coef%nco > 0) THEN
        last_prof = -1
        DO i = 1, nchannels
          prof = chanprof(i)%prof
          IF (prof .NE. last_prof) THEN
            
            sec_gwr(:)          = path(:, prof) * aux%cowr(:,prof)
            sec_gwr_r(:)        = 1.0_jprb / sec_gwr(:)
            
            IF (opts%rt_all%co_data) THEN
              sec_gr_sqrt(:)    = path_sqrt(:, prof) * aux%cor_sqrt(:, prof)
              sec_gr(:)         = sec_gr_sqrt(:)**2_jpim
            ELSE
              sec_gr_sqrt    = path_sqrt(:, prof)
              sec_gr         = path(:, prof)
            ENDIF

            ztemp_1d(:) = 0.4_jprb * predictors(prof)%co(12, :) * sec_gwr_r(:)
            last_prof = prof
          ENDIF
            
          IF (opts%rt_all%co_data) THEN
            DO lay = 1, nlayers
              sec_gr_k(lay) = &!sec_gr_k(lay) + &
                aux%corw_r(lay,prof) *     predictors_k(i)%co(8, lay) + &
                aux%corw_rsqrt(lay,prof) * predictors_k(i)%co(10, lay) + &
                aux%corw_r4rt(lay,prof) *  predictors_k(i)%co(11, lay)

              sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                aux%corw_r(lay,prof) *     predictors_k(i)%co(9, lay)

              aux_k%corw_r(lay,i) = aux_k%corw_r(lay,i) + &
                 sec_gr(lay) *             predictors_k(i)%co(8, lay) + &
                 sec_gr_sqrt(lay) *        predictors_k(i)%co(9, lay)

              aux_k%corw_rsqrt(lay,i) = aux_k%corw_rsqrt(lay,i) + &
                sec_gr(lay) *              predictors_k(i)%co(10, lay)

              aux_k%corw_r4rt(lay,i) = aux_k%corw_r4rt(lay,i) + &
                sec_gr(lay) *              predictors_k(i)%co(11, lay)

              aux_k%cow(lay,i) = aux_k%cow(lay,i) + &
                path(lay, prof) *          predictors_k(i)%co(14, lay)

              path_k(lay,i) = path_k(lay,i) + &
                aux%cow(lay, prof) *       predictors_k(i)%co(14, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              path_k(lay,i) = path_k(lay,i) + &
                                           predictors_k(i)%co(8, lay) + &
                                           predictors_k(i)%co(10, lay) + &
                                           predictors_k(i)%co(11, lay) + &
                                           predictors_k(i)%co(14, lay)

              path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
                                           predictors_k(i)%co(9, lay)
            ENDDO
            sec_gr_k      = 0.0_jprb
            sec_gr_sqrt_k = 0.0_jprb
          ENDIF

          DO lay = 1, nlayers
            sec_gr_k(lay) = sec_gr_k(lay) + &
                                           predictors_k(i)%co(1, lay) + &
              aux%dt(lay,prof) *           predictors_k(i)%co(3, lay) + &
              2.0_jprb * sec_gr(lay) *     predictors_k(i)%co(4, lay) + &
              aux%dtabsdt(lay,prof) *      predictors_k(i)%co(7, lay)
            
            sec_gr_sqrt_k(lay) = sec_gr_sqrt_k(lay) + &
                                           predictors_k(i)%co(2, lay) + &
              aux%dt(lay,prof) *           predictors_k(i)%co(5, lay)

            sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                           predictors_k(i)%co(6, lay)

            aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
              sec_gr(lay) *                predictors_k(i)%co(3, lay) + &
              sec_gr_sqrt(lay) *           predictors_k(i)%co(5, lay)

            aux_k%dtabsdt(lay,i) = aux_k%dtabsdt(lay,i) + &
              sec_gr(lay) *                predictors_k(i)%co(7, lay)

            aux_k%cowr_4rt(lay,i) = aux_k%cowr_4rt(lay,i) + &
              path_4rt(lay,prof) *         predictors_k(i)%co(13, lay)

            path_4rt_k(lay,i) = path_4rt_k(lay,i) + &
              aux%cowr_4rt(lay, prof) *    predictors_k(i)%co(13, lay)
          ENDDO

          IF (opts%rt_all%co_data) THEN           
            sec_gr_sqrt_k(:) = sec_gr_sqrt_k(:) + &
              2.0_jprb * sec_gr_sqrt(:) * sec_gr_k(:)

            path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
              aux%cor_sqrt(:, prof) * sec_gr_sqrt_k(:)

            aux_k%cor_sqrt(:,i) = aux_k%cor_sqrt(:,i) + &
              path_sqrt(:, prof) * sec_gr_sqrt_k(:)

            path_4rt_k(:,i) = path_4rt_k(:,i) + &
              aux%cor_4rt(:, prof) * sec_gr_4rt_k(:)
                      
            aux_k%cor_4rt(:,i) = aux_k%cor_4rt(:,i) + &
              path_4rt(:, prof) * sec_gr_4rt_k(:)
          ELSE
            path_k(:,i) = path_k(:,i) + sec_gr_k
            path_sqrt_k(:,i) = path_sqrt_k(:,i) + sec_gr_sqrt_k
            path_4rt_k(:,i) = path_4rt_k(:,i) + sec_gr_4rt_k
          ENDIF

          sec_gwr_k(:) = &!sec_gwr_k(:) + &
            ztemp_1d(:) * predictors_k(i)%co(12, :)
          
          path_k(:,i)  = path_k(:,i) + &
            aux%cowr(:,prof) * sec_gwr_k(:)

          aux_k%cowr(:,i) = aux_k%cowr(:,i) + &
             path(:, prof) * sec_gwr_k(:)
        ENDDO
      ENDIF

      ! -----------------------------------------------------------------------
      ! CH4
      ! -----------------------------------------------------------------------

      IF (coef%nch4 > 0) THEN

        last_prof = -1
        DO i = 1, nchannels
          prof = chanprof(i)%prof
          IF (prof .NE. last_prof) THEN
            IF (opts%rt_all%ch4_data) THEN
              sec_gr_sqrt(:)    = path_sqrt(:, prof) * aux%ch4r_sqrt(:, prof)
              sec_gr(:)         = path(:, prof) * aux%ch4r(:,prof)
              sec_gw(:)         = path(:, prof) * aux%ch4w(:,prof)       
            ELSE
              sec_gr_sqrt    = path_sqrt(:, prof)
              sec_gr         = path(:, prof)
              sec_gw         = path(:, prof)
            ENDIF
            last_prof = prof
          ENDIF

          IF (opts%rt_all%ch4_data) THEN
            DO lay = 1, nlayers
              sec_gr_k(lay) = &!sec_gr_k(lay) + &
                                         predictors_k(i)%ch4(1, lay) + &
                aux%dt(lay,prof) *       predictors_k(i)%ch4(3, lay) + &
                2.0_jprb * sec_gr(lay) * predictors_k(i)%ch4(4, lay)

              sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
                aux%ch4rw_r(lay,prof) *  predictors_k(i)%ch4(11, lay) + &
                                         predictors_k(i)%ch4(2, lay)
              
              sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
                                         predictors_k(i)%ch4(6, lay)

              sec_gw_k(lay) = &!sec_gw_k(lay) + &
                2.0_jprb * sec_gw(lay) * predictors_k(i)%ch4(9, lay) + &
                                         predictors_k(i)%ch4(10, lay)

              aux_k%ch4r(lay,i) = aux_k%ch4r(lay,i) + &
                aux%dt(lay,prof) *       predictors_k(i)%ch4(5, lay)

              aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
                aux%ch4r(lay,prof) *     predictors_k(i)%ch4(5, lay)

              aux_k%ch4rw_r(lay,i) = aux_k%ch4rw_r(lay,i) + &
                sec_gr_sqrt(lay) *       predictors_k(i)%ch4(11, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              ! copied from co2_data = T to save initialising to 0 in that case
              path_k(lay,i) = path_k(lay,i)  + &
                                         predictors_k(i)%ch4(1, lay) + &
                aux%dt(lay,prof) *       predictors_k(i)%ch4(3, lay) + &
                2.0_jprb * sec_gr(lay) * predictors_k(i)%ch4(4, lay)

              path_4rt_k(lay,i) = path_4rt_k(lay,i) + &
                                         predictors_k(i)%ch4(6, lay)

              path_k(lay,i) = path_k(lay,i)  + &
                2.0_jprb * sec_gw(lay) * predictors_k(i)%ch4(9, lay) + &
                                         predictors_k(i)%ch4(10, lay)

              aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
                                         predictors_k(i)%ch4(5, lay)

              path_sqrt_k(lay,i) = path_sqrt_k(lay,i) + &
                                         predictors_k(i)%ch4(11, lay) + &
                                         predictors_k(i)%ch4(2, lay)
            ENDDO
          ENDIF

          DO lay = 1, nlayers
            path_k(lay,i) = path_k(lay,i) + &
              aux%ch4wr(lay,prof) *      predictors_k(i)%ch4(7, lay)

            aux_k%ch4wr(lay,i) = aux_k%ch4wr(lay,i) + &
              path(lay, prof) *          predictors_k(i)%ch4(7, lay) + &
                                         predictors_k(i)%ch4(8, lay)

            aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
              sec_gr(lay) *              predictors_k(i)%ch4(3, lay)
          ENDDO

          IF (opts%rt_all%ch4_data) THEN
! using sec_gr = path * ch4r rather than sec_g_sqrt**2
            path_k(:,i) = path_k(:,i) + &
              aux%ch4w(:,prof) * sec_gw_k(:) + &
              aux%ch4r(:,prof) * sec_gr_k(:)

            path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
              aux%ch4r_sqrt(:, prof) * sec_gr_sqrt_k(:)

            path_4rt_k(:,i) = path_4rt_k(:,i) + &
              aux%ch4r_4rt(:, prof) * sec_gr_4rt_k(:)

            aux_k%ch4w(:,i) = aux_k%ch4w(:,i) + &
              path(:, prof) * sec_gw_k(:)

            aux_k%ch4r(:,i) = aux_k%ch4r(:,i) + &
              path(:, prof) * sec_gr_k(:)

            aux_k%ch4r_sqrt(:,i) = aux_k%ch4r_sqrt(:,i) + &
              path_sqrt(:, prof) * sec_gr_sqrt_k(:)

            aux_k%ch4r_4rt(:,i) = aux_k%ch4r_4rt(:,i) + &
              path_4rt(:, prof) * sec_gr_4rt_k(:)
          ENDIF
        ENDDO
      ENDIF

      ! -----------------------------------------------------------------------
      ! SO2
      ! -----------------------------------------------------------------------

      IF (coef%nso2 > 0) THEN
        last_prof = -1
        DO i = 1, nchannels
          prof = chanprof(i)%prof
          IF (prof .NE. last_prof) THEN
            sec_gr_sqrt(:) = path_sqrt(:, prof) * aux%so2r_sqrt(:, prof)
            sec_gr(:) = sec_gr_sqrt(:)**2_jpim

            sec_gw_sqrt(:) = path_sqrt(:, prof) * aux%so2w_sqrt(:, prof)
            sec_gw_4rt(:) = path_4rt(:, prof) * aux%so2w_4rt(:, prof)
            sec_gw(:) = sec_gw_sqrt(:)**2_jpim

            last_prof = prof
          ENDIF

          DO lay=1, nlayers
            aux_k%so2rwr_r(lay,i) = aux_k%so2rwr_r(lay,i) + &
              sec_gr(lay) * predictors_k(i)%so2(12, lay) + &
              sec_gr_sqrt(lay) * predictors_k(i)%so2(13, lay)
            
            aux_k%so2rw_r(lay,i) = aux_k%so2rw_r(lay,i) + &
              sec_gr(lay) * predictors_k(i)%so2(19, lay)

            sec_gr_k(lay) = &!sec_gr_k(lay) + &
              aux%so2rwr_r(lay, prof) * predictors_k(i)%so2(12, lay) + &
              sec_gr_sqrt(lay) * predictors_k(i)%so2(15, lay) + &
              aux%so2rw_r(lay, prof) * predictors_k(i)%so2(19, lay)

            sec_gr_sqrt_k(lay) = &!sec_gr_sqrt_k(lay) + &
              aux%so2rwr_r(lay, prof) * predictors_k(i)%so2(13, lay) + &
              sec_gr(lay) * predictors_k(i)%so2(15, lay) + &
              aux%so2r(lay, prof) * predictors_k(i)%so2(17, lay)

            sec_gw_k(lay) = &!sec_gw_k(lay) + &
              sec_gw_sqrt(lay) * predictors_k(i)%so2(14, lay) + &
              sec_gw_4rt(lay) * predictors_k(i)%so2(16, lay)

            sec_gw_sqrt_k(lay) = &!sec_gw_sqrt_k(lay) + &
              sec_gw(lay) * predictors_k(i)%so2(14, lay)

            sec_gw_4rt_k(lay) = &!sec_gw_4rt_k(lay) + &
              sec_gw(lay) * predictors_k(i)%so2(16, lay)

            aux_k%so2r(lay,i) = aux_k%so2r(lay,i) + &
              sec_gr_sqrt(lay) * predictors_k(i)%so2(17, lay)

            sec_gr_k(lay) = sec_gr_k(lay) + &
              2.0_jprb * sec_gr(lay) * predictors_k(i)%so2(1, lay) + &
              aux%dt(lay, prof) * predictors_k(i)%so2(4, lay) + &
              predictors_k(i)%so2(7, lay) + &
              3.0_jprb * sec_gr(lay)**2_jpim * predictors_k(i)%so2(8, lay) + &
              4.0_jprb * sec_gr(lay)**3_jpim * predictors_k(i)%so2(9, lay) + &
              aux%dtabsdt(lay, prof) * predictors_k(i)%so2(10, lay) + &
              sec_gr_sqrt(lay) * aux%dt(lay, prof) * predictors_k(i)%so2(18, lay)
                         
            sec_gw_k(lay) = sec_gw_k(lay) + &
              predictors_k(i)%so2(2, lay) + &
              2.0_jprb * sec_gw(lay) * predictors_k(i)%so2(3, lay) 
            
            aux_k%dt(lay,i) = aux_k%dt(lay,i) + &
              sec_gr(lay) * predictors_k(i)%so2(4, lay) + &
              sec_gr_sqrt(lay) * predictors_k(i)%so2(11, lay) + &
              sec_gr(lay) * sec_gr_sqrt(lay) * predictors_k(i)%so2(18, lay)
            
            aux_k%dtabsdt(lay,i) = aux_k%dtabsdt(lay,i) + &
              sec_gr(lay) * predictors_k(i)%so2(10, lay)
            
            sec_gr_sqrt_k(lay) = sec_gr_sqrt_k(lay) + &
              predictors_k(i)%so2(5, lay) + &
              aux%dt(lay, prof) * predictors_k(i)%so2(11, lay) + &
              aux%dt(lay, prof) * sec_gr(lay) * predictors_k(i)%so2(18, lay)
            
            sec_gr_4rt_k(lay) = &!sec_gr_4rt_k(lay) + &
              predictors_k(i)%so2(6, lay)
          ENDDO
          
          sec_gw_sqrt_k(:) = sec_gw_sqrt_k(:) + &
            2._jprb * sec_gw_sqrt(:) * sec_gw_k(:)
          
          sec_gr_sqrt_k(:) = sec_gr_sqrt_k(:) + &
            2._jprb * sec_gr_sqrt(:) * sec_gr_k(:)
          
          path_sqrt_k(:,i) = path_sqrt_k(:,i) + &
            aux%so2w_sqrt(:, prof) * sec_gw_sqrt_k(:) + &
            aux%so2r_sqrt(:, prof) * sec_gr_sqrt_k(:)
          
          path_4rt_k(:,i) = path_4rt_k(:,i) + &
            aux%so2w_4rt(:, prof) * sec_gw_4rt_k(:) + &
            aux%so2r_4rt(:, prof) * sec_gr_4rt_k(:)
          
          aux_k%so2w_sqrt(:,i) = aux_k%so2w_sqrt(:,i) + &
            path_sqrt(:, prof) * sec_gw_sqrt_k(:)
          
          aux_k%so2r_sqrt(:,i) = aux_k%so2r_sqrt(:,i) + &
            path_sqrt(:, prof) * sec_gr_sqrt_k(:)
          
          aux_k%so2w_4rt(:,i) = aux_k%so2w_4rt(:,i) + &
            path_4rt(:, prof) * sec_gw_4rt_k(:)
          
          aux_k%so2r_4rt(:,i) = aux_k%so2r_4rt(:,i) + &
            path_4rt(:, prof) * sec_gr_4rt_k(:)
        ENDDO
      ENDIF
    ENDIF
  ENDIF

  ! DAR: gases should be done in reverse order according to AD coding theory but can almost be treated as being independent,
  ! making it easier to initialise arrays with calculations when the gas is mandatory - hence why we do mixedgas first but
  ! this leads to problems when considering some temperature AD variables which have dependencies on multiple gases).
  ! These need to be treated in the correct order
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    path_k(:,i) = path_k(:,i) + aux%tr(:, prof) * sec_tr_k(:,i)
    aux_k%tr(:,i) = aux_k%tr(:,i) + path(:, prof) * sec_tr_k(:,i)
  ENDDO

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

    ! This is on coef levels so pressure is not an active variable in the K

    DO i = 1, nchannels
       j = chanprof(i)%prof
       k = chanprof(i)%chan

      DO lay = 1, nlayers
        lev = lay + 1

        acm = aux%co2_cm(j) / (2._jprb*betaplus1*Lcel_cm)

        Pcel_Lev = acm*raytracing%pathsat(lay,j)*profiles(j)%p(lev  )**2_jpim + Pcel(k)**2_jpim
        Pnom_LevM1 = acm*raytracing%pathsat(lay,j)*profiles(j)%p(lev-1)**2_jpim + Pnom(k)**2_jpim
        Pcel_LevM1 = acm*raytracing%pathsat(lay,j)*profiles(j)%p(lev-1)**2_jpim + Pcel(k)**2_jpim
        Pnom_Lev = acm*raytracing%pathsat(lay,j)*profiles(j)%p(lev  )**2_jpim + Pnom(k)**2_jpim
        ! AD
        Pcel_Lev_K =  0._jprb
        Pnom_LevM1_K = 0._jprb
        Pcel_LevM1_K =  0._jprb
        Pnom_Lev_K = 0._jprb
        acm_K = 0._jprb
        ! for PMC predictor 2
        predictors_K(i)%pmc(2,lay,k) = 0._jprb
        ! for PMC predictor 1
        Pcel_Lev_K = Pcel_Lev_K &
          + predictors_K(i)%pmc(1,lay,k) * Pnom_LevM1/(Pcel_Lev*Pnom_LevM1)
        Pnom_LevM1_K = Pnom_LevM1_K &
          + predictors_K(i)%pmc(1,lay,k) * Pcel_Lev/(Pcel_Lev*Pnom_LevM1)
        Pcel_LevM1_K = Pcel_LevM1_K &
          - predictors_K(i)%pmc(1,lay,k) * Pnom_Lev/(Pcel_LevM1*Pnom_Lev)
        Pnom_Lev_K = Pnom_Lev_K &
          - predictors_K(i)%pmc(1,lay,k) * Pcel_LevM1/(Pcel_LevM1*Pnom_Lev)
        predictors_K(i)%pmc(1,lay,k) = 0._jprb
        
        acm_K = acm_K &
          + Pnom_Lev_K * raytracing%pathsat(lay,j)*profiles(j)%p(lev  )**2_jpim
        raytracing_K%pathsat(lay,i) = raytracing_K%pathsat(lay,i) &
          + Pnom_Lev_K * acm * profiles(j)%p(lev  )**2_jpim
        Pnom_Lev_K = 0._jprb
        
        acm_K = acm_K &
          + Pcel_LevM1_K * raytracing%pathsat(lay,j)*profiles(j)%p(lev-1)**2_jpim
        raytracing_K%pathsat(lay,i) = raytracing_K%pathsat(lay,i) &
          + Pcel_LevM1_K * acm * profiles(j)%p(lev-1)**2_jpim
        Pcel_LevM1_K = 0._jprb
        
        acm_K = acm_K &
          + Pnom_LevM1_K * raytracing%pathsat(lay,j)*profiles(j)%p(lev-1)**2_jpim
        raytracing_K%pathsat(lay,i) = raytracing_K%pathsat(lay,i) &
          + Pnom_LevM1_K * acm * profiles(j)%p(lev-1)**2_jpim
        Pnom_LevM1_K = 0._jprb
        
        acm_K = acm_K &
          + Pcel_Lev_K * raytracing%pathsat(lay,j)*profiles(j)%p(lev  )**2_jpim
        raytracing_K%pathsat(lay,i) = raytracing_K%pathsat(lay,i) &
          + Pcel_Lev_K * acm * profiles(j)%p(lev  )**2_jpim
        Pcel_Lev_K = 0._jprb
        
        aux_k%co2_cm(i) = aux_k%co2_cm(i) &
          + acm_K/(2._jprb*betaplus1*Lcel_cm)
        acm_K = 0._jprb
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_789_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_789_k
