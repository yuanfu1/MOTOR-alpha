! Description:
!> @file
!!   AD of v7, v8 and v9 predictor pre-calculations.
!
!> @brief
!!   AD of v7, v8 and v9 predictor pre-calculations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     dosolar              flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel       flag for strict plane parallel geometry
!! @param[in]     prof                 input profiles
!! @param[in,out] prof_ad              input profile increments
!! @param[in]     coef                 optical depth coefficient structure
!! @param[in]     coef_pccomp          PC coefficients structure
!! @param[in]     aux                  coef level auxiliary profile data structure
!! @param[in,out] aux_ad               coef level auxiliary profile data increments
!! @param[in]     angles               geometry structure
!! @param[in]     raytracing           raytracing structure
!! @param[in,out] raytracing_ad        raytracing structure increments
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
SUBROUTINE rttov_predictor_precalc_789_ad( &
              opts,           &
              dosolar,        &
              plane_parallel, &
              prof,           &
              prof_ad,        &
              coef,           &
              coef_pccomp,    &
              aux,            &
              aux_ad,         &
              angles,         &
              raytracing,     &
              raytracing_ad)

  USE rttov_types, ONLY :  &
        rttov_options,          &
        rttov_coef,             &
        rttov_coef_pccomp,      &
        rttov_profile,          &
        rttov_profile_aux_coef, &
        rttov_geometry,         &
        rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : inst_id_ssmis
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod, ONLY: invsqrt_ad, reciprocal_ad, sqrt_ad
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),          INTENT(IN)    :: opts
  LOGICAL(jplm),                INTENT(IN)    :: dosolar
  LOGICAL(jplm),                INTENT(IN)    :: plane_parallel
  TYPE(rttov_profile),          INTENT(IN)    :: prof(:)
  TYPE(rttov_profile),          INTENT(INOUT) :: prof_ad(:)
  TYPE(rttov_coef),             INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp),      INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_ad
  TYPE(rttov_geometry),         INTENT(IN)    :: angles(:)
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_ad
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers, nlevels
  INTEGER(jpim) :: iprof, lay
  REAL(jprb) :: ZHOOK_HANDLE

  REAL(jprb) :: sum1_ad, sum2_ad
  INTEGER(jpim) :: iv2lay, iv2lev

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_789_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlevels

!-----------------------------------------------------------------------------------------
! AD of predictor data
!-----------------------------------------------------------------------------------------

  ! Compute the STP co2 thickness in cm
  IF (coef%pmc_shift) CALL calc_ssu_co2_thickness_ad()  

  ! t_layer shouldn't be initialised yet but the zeeman code for SSMIS does do this so don't overwrite
  IF (.NOT. (coef%id_inst == inst_id_ssmis .AND. coef%inczeeman)) aux_ad%t_layer = 0._jprb
  aux_ad%w_layer = 0.0_jprb

  ! -------------------------------------------------------------------------
  ! Calculate and store useful intermediate variables
  ! -------------------------------------------------------------------------

  IF (coef%nco > 0 .AND. opts%rt_all%co_data) THEN
    aux_ad%cor = &!aux_ad%cor + & 
                              aux%cow_r4rt *                   aux_ad%corw_r4rt + &
                              aux%cow_rsqrt * (aux%cow_rsqrt * aux_ad%corw_r + &
                                                               aux_ad%corw_rsqrt)

    aux_ad%cow_r4rt = &!aux_ad%cow_r4rt + &
                      aux%cor * aux_ad%corw_r4rt

    aux_ad%cow_rsqrt = &!aux_ad%cow_rsqrt + &
                                          aux%cor * (2.0_jprb * aux%cow_rsqrt * aux_ad%corw_r + &
                                                                                aux_ad%corw_rsqrt)

    CALL sqrt_ad(aux%cow_r4rt, aux_ad%cow_rsqrt, aux_ad%cow_r4rt, acc = .TRUE._jplm)
    CALL invsqrt_ad(aux%cow_rsqrt, aux_ad%cow, aux_ad%cow_rsqrt, acc = .TRUE._jplm)
  ENDIF

  IF (coef%nozone > 0 .and. opts%rt_all%ozone_data) THEN
      aux_ad%ow_rsqrt = &!aux_ad%ow_rsqrt + &
        aux%ow * aux_ad%ow_sqrt + &
        2._jprb * aux%ow_rsqrt * aux_ad%ow_r

      aux_ad%ow = aux_ad%ow + aux%ow_rsqrt * aux_ad%ow_sqrt

    CALL sqrt_ad(aux%or_sqrt, aux_ad%or, aux_ad%or_sqrt, acc = .TRUE._jplm)
    CALL invsqrt_ad(aux%ow_rsqrt, aux_ad%ow, aux_ad%ow_rsqrt, acc = .TRUE._jplm)
  ENDIF

  IF (coef%fmv_model_ver == 9) THEN
    CALL invsqrt_ad(aux%ww_4rt, aux_ad%ww_rsqrt, aux_ad%ww_4rt, acc = .FALSE._jplm)
    aux_ad%ww = aux_ad%ww + aux%ww_rsqrt * aux_ad%ww_sqrt
    aux_ad%ww_rsqrt = aux_ad%ww_rsqrt + aux%ww * aux_ad%ww_sqrt
    CALL invsqrt_ad(aux%ww_rsqrt, aux_ad%ww, aux_ad%ww_rsqrt, acc = .TRUE._jplm)
  ENDIF

  aux_ad%wr = aux_ad%wr + aux%wr_rsqrt * aux_ad%wr_sqrt
  aux_ad%wr_rsqrt = &!aux_ad%wr_rsqrt + 
    aux%wr * aux_ad%wr_sqrt

  CALL invsqrt_ad(aux%wr_4rt, aux_ad%wr_rsqrt, aux_ad%wr_4rt, acc = .TRUE._jplm)
  CALL invsqrt_ad(aux%wr_rsqrt, aux_ad%wr, aux_ad%wr_rsqrt, acc = .TRUE._jplm)

  IF (coef%fmv_model_ver <= 8) THEN
    CALL sqrt_ad(aux%tw_4rt, aux_ad%tw_sqrt, aux_ad%tw_4rt, acc = .FALSE._jplm)  
    CALL sqrt_ad(aux%tw_sqrt, aux_ad%tw, aux_ad%tw_sqrt, acc = .TRUE._jplm)  
  ENDIF

  IF (coef%fmv_model_ver >= 8) CALL sqrt_ad(aux%tr_sqrt, aux_ad%tr, aux_ad%tr_sqrt, acc = .TRUE._jplm)

  CALL reciprocal_ad(aux%tr_r, aux_ad%tr, aux_ad%tr_r, acc = .TRUE._jplm)

  ! -------------------------------------------------------------------------
  ! Calculate profile / reference profile sums
  ! -------------------------------------------------------------------------
   
  IF (coef%fmv_model_ver >= 8) THEN

    IF (coef%nso2 > 0) THEN
      aux_ad%so2w_sqrt = aux_ad%so2w_sqrt + &
        0.5_jprb * aux_ad%so2w_4rt / aux%so2w_4rt

      aux_ad%so2w = &!aux_ad%so2w + &
        0.5_jprb * aux_ad%so2w_sqrt / aux%so2w_sqrt

      aux_ad%so2r_sqrt = aux_ad%so2r_sqrt + &
        0.5_jprb * aux_ad%so2r_4rt / aux%so2r_4rt

      aux_ad%so2r = aux_ad%so2r + &
        0.5_jprb * aux_ad%so2r_sqrt / aux%so2r_sqrt

      aux_ad%so2wr_r = & !aux_ad%so2wr_r + &
                       aux%so2r * aux_ad%so2rwr_r
      aux_ad%so2w_r = & !aux_ad%so2w_r + &
                      aux%so2r * aux_ad%so2rw_r

      aux_ad%so2r = aux_ad%so2r + & 
                    aux%so2wr_r * aux_ad%so2rwr_r + &
                    aux%so2w_r * aux_ad%so2rw_r

      CALL reciprocal_ad(aux%so2w_r, aux_ad%so2w, aux_ad%so2w_r, acc = .TRUE._jplm)
      CALL reciprocal_ad(aux%so2wr_r, aux_ad%so2wr, aux_ad%so2wr_r, acc = .FALSE._jplm)

      DO iprof = 1, nprofiles
        sum1_ad = 0._jprb
        sum2_ad = 0._jprb

        DO lay = nlayers, 1, -1
          sum2_ad = sum2_ad + coef%so2tstar_wsum_r(lay) * aux_ad%so2wr(lay, iprof)

          aux_ad%so2_layer(lay, iprof) = & !aux_ad%so2_layer(lay, iprof) + &
            coef%dpp(lay - 1) * aux%t_layer(lay, iprof) * sum2_ad

          aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
            coef%dpp(lay - 1) * aux%so2_layer(lay, iprof) * sum2_ad
         
          sum1_ad = sum1_ad + coef%so2star_wsum_r(lay) * aux_ad%so2w(lay, iprof)
          aux_ad%so2_layer(lay, iprof) = aux_ad%so2_layer(lay, iprof) + coef%dpp(lay - 1) * sum1_ad
        ENDDO
      ENDDO
    ENDIF

    IF (coef%nch4 > 0) THEN
      IF (opts%rt_all%ch4_data) THEN
        aux_ad%ch4r = aux_ad%ch4r + aux%ch4w_r * aux_ad%ch4rw_r
        aux_ad%ch4w_r = &! aux_ad%ch4w_r + &
                        aux%ch4r * aux_ad%ch4rw_r

        CALL reciprocal_ad(aux%ch4w_r, aux_ad%ch4w, aux_ad%ch4w_r, acc = .TRUE._jplm)
      ENDIF

      DO iprof = 1, nprofiles
        sum2_ad = 0._jprb
        IF (opts%rt_all%ch4_data) THEN
          sum1_ad = 0._jprb
          DO lay = nlayers, 1, -1
            sum1_ad = sum1_ad + coef%ch4star_wsum_r(lay) * aux_ad%ch4w(lay, iprof)
            sum2_ad = sum2_ad + coef%ch4tstar_wsum_r(lay) * aux_ad%ch4wr(lay, iprof)

            aux_ad%ch4_layer(lay, iprof) = &!aux_ad%ch4_layer(lay, iprof) + 
              coef%dpp(lay - 1) * ( &
                                    sum1_ad + &
                                    aux%t_layer(lay, iprof) * sum2_ad)

            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * aux%ch4_layer(lay, iprof) * sum2_ad
          ENDDO
        ELSE
          DO lay = nlayers, 1, -1
            sum2_ad = sum2_ad + coef%ch4tstar_wsum_r(lay) * aux_ad%ch4wr(lay, iprof)
            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * coef%ch4star(lay) * sum2_ad 
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    IF (coef%nco > 0) THEN
      aux_ad%cowr = aux_ad%cowr + 0.25_jprb * (aux%cowr_4rt * aux%cowr_r) * aux_ad%cowr_4rt

      DO iprof = 1, nprofiles
        sum2_ad = 0._jprb
        IF (opts%rt_all%co_data) THEN
          sum1_ad = 0._jprb

          DO lay = nlayers, 1, -1
            sum1_ad = sum1_ad + coef%costar_wsum_r(lay) * aux_ad%cow(lay, iprof)
            sum2_ad = sum2_ad + coef%cotstar_wsum_r(lay) * aux_ad%cowr(lay, iprof)

            aux_ad%co_layer(lay, iprof) = &!aux_ad%co_layer(lay, iprof) + &!
              coef%dpp(lay - 1) * ( &
              sum1_ad + &
              aux%t_layer(lay, iprof) * sum2_ad)

            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * aux%co_layer(lay, iprof) * sum2_ad
          ENDDO
        ELSE
          DO lay = nlayers, 1, -1
            sum2_ad = sum2_ad + coef%cotstar_wsum_r(lay) * aux_ad%cowr(lay, iprof)
            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * coef%costar(lay) * sum2_ad 
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    IF (coef%nn2o > 0) THEN
      IF (opts%rt_all%n2o_data) CALL reciprocal_ad(aux%n2ow_r, aux_ad%n2ow, aux_ad%n2ow_r, acc = .TRUE._jplm)
      DO iprof = 1, nprofiles
        sum2_ad = 0._jprb
        IF (opts%rt_all%n2o_data) THEN
          sum1_ad = 0._jprb
          
          DO lay = nlayers, 1, -1
            sum1_ad = sum1_ad + coef%n2ostar_wsum_r(lay) * aux_ad%n2ow(lay, iprof)
            sum2_ad = sum2_ad + coef%n2otstar_wsum_r(lay) * aux_ad%n2owr(lay, iprof)
            
            aux_ad%n2o_layer(lay, iprof) = &!aux_ad%n2o_layer(lay, iprof) + &
              coef%dpp(lay - 1) * ( &
              sum1_ad + &
              aux%t_layer(lay, iprof) * sum2_ad)

            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * aux%n2o_layer(lay, iprof) * sum2_ad
          ENDDO
        ELSE
          DO lay = nlayers, 1, -1
            sum2_ad = sum2_ad + coef%n2otstar_wsum_r(lay) * aux_ad%n2owr(lay, iprof)
            aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
              coef%dpp(lay - 1) * coef%n2ostar(lay) * sum2_ad 
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    IF (coef%nco2 > 0) THEN
      IF (opts%rt_all%co2_data) THEN
        DO iprof = 1, nprofiles
          sum1_ad = 0._jprb
          DO lay = nlayers, 1, -1
            ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
            ! need dpp(0) to start
            sum1_ad = sum1_ad + coef%co2star_wsum_r(lay) * aux_ad%co2w(lay, iprof)
            aux_ad%co2_layer(lay, iprof) = &!aux_ad%co2_layer(lay, iprof) + &
                                           coef%dpp(lay - 1) * sum1_ad
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    IF (coef%fmv_model_ver == 9) THEN
      DO iprof = 1, nprofiles
        sum1_ad   = 0._jprb
        DO lay = nlayers, 2, -1
          sum1_ad = sum1_ad + coef%tstar_uwsum_r(lay) * aux_ad%tuw(lay, iprof)
          aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + sum1_ad
        ENDDO
      
        aux_ad%t_layer(1, iprof) = aux_ad%t_layer(1, iprof) + &
                                   coef%tstar_r(1) * aux_ad%tuw(1, iprof) + sum1_ad
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      sum2_ad   = 0._jprb
      DO lay = nlayers, 2, -1
        IF (coef%fmv_model_ver == 8) THEN
          sum2_ad = sum2_ad + aux_ad%twr(lay, iprof) * coef%tstarmod_wsum_r(lay)
          aux_ad%t_layer(lay-1, iprof) = aux_ad%t_layer(lay-1, iprof) + coef%dpp(lay - 1) * sum2_ad
        ELSE
          sum2_ad = sum2_ad + coef%tstar_wsum_r(lay) * aux_ad%twr(lay, iprof)
          aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + coef%dpp(lay - 1) * sum2_ad
        ENDIF
      ENDDO
      
      IF (coef%fmv_model_ver == 8) THEN
        aux_ad%twr(1, iprof) = 0.0_jprb
      ELSEIF (coef%fmv_model_ver == 9) THEN
        aux_ad%t_layer(1, iprof) = aux_ad%t_layer(1, iprof) + &
                                   coef%tstar_r(1) * aux_ad%twr(1, iprof) + &
                                   coef%dpp(0) * sum2_ad
      ENDIF
    ENDDO
  ENDIF

  ! All versions again

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO iprof = 1, nprofiles
        sum1_ad = 0._jprb
        DO lay = nlayers, 1, -1
          sum1_ad = sum1_ad + coef%ostar_wsum_r(lay) * aux_ad%ow(lay, iprof)
          aux_ad%o3_layer(lay, iprof) = &!aux_ad%o3_layer(lay, iprof) + &
                                        coef%dpp(lay - 1) * sum1_ad
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF (coef%fmv_model_ver >= 8) THEN
    aux_ad%wr = aux_ad%wr + aux%wwr_r * aux_ad%wrwr_r
    aux_ad%wwr_r = &!aux_ad%wwr_r + &
                   aux%wr * aux_ad%wrwr_r
    CALL reciprocal_ad(aux%wwr_r, aux_ad%wwr, aux_ad%wwr_r, acc = .FALSE._jplm)

    DO iprof = 1, nprofiles
      sum1_ad = 0._jprb
      DO lay = nlayers, 1, -1
        sum1_ad = sum1_ad + coef%wtstar_wsum_r(lay) * aux_ad%wwr(lay, iprof)

        aux_ad%w_layer(lay, iprof) = aux_ad%w_layer(lay, iprof) + &
          coef%dpp(lay - 1) * aux%t_layer(lay, iprof) * sum1_ad

        aux_ad%t_layer(lay, iprof) = aux_ad%t_layer(lay, iprof) + &
          coef%dpp(lay - 1) * aux%w_layer(lay, iprof) * sum1_ad
      ENDDO
    ENDDO
  ENDIF

  IF (.NOT. coef%fmv_model_ver == 8) THEN 
    aux_ad%wr =  aux_ad%wr + aux%ww_r * aux_ad%wrw_r
    aux_ad%ww_r = &!aux_ad%ww_r + &
                  aux%wr * aux_ad%wrw_r
    CALL reciprocal_ad(aux%ww_r, aux_ad%ww, aux_ad%ww_r, acc = .TRUE._jplm)
  ENDIF

  DO iprof = 1, nprofiles
    sum1_ad = 0._jprb
    DO lay = nlayers, 1, -1
      sum1_ad = sum1_ad + coef%wstar_wsum_r(lay) * aux_ad%ww(lay, iprof)
      aux_ad%w_layer(lay, iprof) = aux_ad%w_layer(lay, iprof) + &
                                   coef%dpp(lay - 1) * sum1_ad
    ENDDO
  ENDDO

  IF (coef%fmv_model_ver <= 8) THEN  
    DO iprof = 1, nprofiles
      DO lay = nlayers, 2, -1
        aux_ad%tw(lay - 1, iprof) = aux_ad%tw(lay - 1, iprof) + aux_ad%tw(lay, iprof)
        aux_ad%tr(lay - 1, iprof) = aux_ad%tr(lay - 1, iprof) + coef%dpp(lay - 1) * aux_ad%tw(lay, iprof) 
      ENDDO
    ENDDO
  ENDIF

  DO iprof = 1, nprofiles

    ! -----------------------------------------------------------------------
    ! Calculate (profile / reference profile) ratios
    ! -----------------------------------------------------------------------

    aux_ad%tr(:, iprof) = aux_ad%tr(:, iprof) + 2.0_jprb * aux%tr(:, iprof) * aux_ad%tr2(:, iprof)

    aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + &
                               aux_ad%tr(:, iprof) * coef%tstar_r(:)

    aux_ad%w_layer(:, iprof) = aux_ad%w_layer(:, iprof) + aux_ad%wr(:, iprof) * coef%wstar_r(:)

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux_ad%o3_layer(:, iprof) = aux_ad%o3_layer(:, iprof) + aux_ad%or(:, iprof) * coef%ostar_r(:)
    ENDIF

    IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
      IF (coef%fmv_model_ver == 9) &
        CALL sqrt_ad(aux%co2r_sqrt(:, iprof), aux_ad%co2r(:, iprof), aux_ad%co2r_sqrt(:, iprof), acc = .TRUE._jplm)
      aux_ad%co2_layer(:, iprof) = aux_ad%co2_layer(:, iprof) + aux_ad%co2r(:, iprof) * coef%co2star_r(:)
    ENDIF

    IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
      CALL sqrt_ad(aux%n2or_4rt(:, iprof), aux_ad%n2or_sqrt(:, iprof), aux_ad%n2or_4rt(:, iprof), acc = .TRUE._jplm)
      CALL sqrt_ad(aux%n2or_sqrt(:, iprof), aux_ad%n2or(:, iprof), aux_ad%n2or_sqrt(:, iprof), acc = .TRUE._jplm)
      aux_ad%n2o_layer(:, iprof) = aux_ad%n2o_layer(:, iprof) + aux_ad%n2or(:, iprof) * coef%n2ostar_r(:)
    ENDIF

    IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
      CALL sqrt_ad(aux%cor_4rt(:, iprof), aux_ad%cor_sqrt(:, iprof), aux_ad%cor_4rt(:, iprof), acc = .TRUE._jplm)
      CALL sqrt_ad(aux%cor_sqrt(:, iprof), aux_ad%cor(:, iprof), aux_ad%cor_sqrt(:, iprof), acc = .TRUE._jplm)
      aux_ad%co_layer(:, iprof) = aux_ad%co_layer(:, iprof) + aux_ad%cor(:, iprof) * coef%costar_r(:)
    ENDIF

    IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
      CALL sqrt_ad(aux%ch4r_4rt(:, iprof), aux_ad%ch4r_sqrt(:, iprof), aux_ad%ch4r_4rt(:, iprof), acc = .TRUE._jplm)
      CALL sqrt_ad(aux%ch4r_sqrt(:, iprof), aux_ad%ch4r(:, iprof), aux_ad%ch4r_sqrt(:, iprof), acc = .TRUE._jplm)
      aux_ad%ch4_layer(:, iprof) = aux_ad%ch4_layer(:, iprof) + aux_ad%ch4r(:, iprof) * coef%ch4star_r(:)
    ENDIF

    IF (coef%nso2 > 0) THEN
      aux_ad%so2_layer(:, iprof) = aux_ad%so2_layer(:, iprof) + aux_ad%so2r(:, iprof) * coef%so2star_r(:)
    ENDIF

    ! -----------------------------------------------------------------------
    ! Calculate deviations from reference profile for each layer
    ! -----------------------------------------------------------------------

    aux_ad%dt(:,iprof) = aux_ad%dt(:,iprof) + 2.0 * ABS(aux%dt(:, iprof)) * aux_ad%dtabsdt(:, iprof)
    aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + aux_ad%dt(:, iprof) 

    IF (coef%nozone > 0) THEN
      IF (coef%fmv_model_ver == 9) THEN
       aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + &
         aux_ad%tro(:, iprof) * coef%to3star_r(:)
      ENDIF
      aux_ad%t_layer(:, iprof) = aux_ad%t_layer(:, iprof) + aux_ad%dto(:, iprof) 
    ENDIF

    ! -----------------------------------------------------------------------
    ! Profile layer quantities: layer n lies between levels n and n+1
    ! -----------------------------------------------------------------------

    IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
      prof_ad(iprof)%so2(1) = prof_ad(iprof)%so2(1) + &
        0.5_jprb * aux_ad%so2_layer(1, iprof)
      prof_ad(iprof)%so2(2:nlevels-1) = prof_ad(iprof)%so2(2:nlevels-1) + &
        0.5_jprb * (aux_ad%so2_layer(1:nlevels-2, iprof) + aux_ad%so2_layer(2:nlevels-1, iprof))
      prof_ad(iprof)%so2(nlevels) = prof_ad(iprof)%so2(nlevels) + &
        0.5_jprb * aux_ad%so2_layer(nlevels-1, iprof)
    ENDIF

    IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
      IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
        prof_ad(iprof)%ch4(1) = prof_ad(iprof)%ch4(1) + &
          0.5_jprb * aux_ad%ch4_layer(1, iprof)
        prof_ad(iprof)%ch4(2:nlevels-1) = prof_ad(iprof)%ch4(2:nlevels-1) + &
          0.5_jprb * (aux_ad%ch4_layer(1:nlevels-2, iprof) + aux_ad%ch4_layer(2:nlevels-1, iprof))
        prof_ad(iprof)%ch4(nlevels) = prof_ad(iprof)%ch4(nlevels) + &
          0.5_jprb * aux_ad%ch4_layer(nlevels-1, iprof)
      ENDIF

      IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
        prof_ad(iprof)%co(1) = prof_ad(iprof)%co(1) + &
          0.5_jprb * aux_ad%co_layer(1, iprof)
        prof_ad(iprof)%co(2:nlevels-1) = prof_ad(iprof)%co(2:nlevels-1) + &
          0.5_jprb * (aux_ad%co_layer(1:nlevels-2, iprof) + aux_ad%co_layer(2:nlevels-1, iprof))
        prof_ad(iprof)%co(nlevels) = prof_ad(iprof)%co(nlevels) + &
          0.5_jprb * aux_ad%co_layer(nlevels-1, iprof)
      ENDIF
      
      IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
        prof_ad(iprof)%n2o(1) = prof_ad(iprof)%n2o(1) + &
          0.5_jprb * aux_ad%n2o_layer(1, iprof)
        prof_ad(iprof)%n2o(2:nlevels-1) = prof_ad(iprof)%n2o(2:nlevels-1) + &
          0.5_jprb * (aux_ad%n2o_layer(1:nlevels-2, iprof) + aux_ad%n2o_layer(2:nlevels-1, iprof))
        prof_ad(iprof)%n2o(nlevels) = prof_ad(iprof)%n2o(nlevels) + &
          0.5_jprb * aux_ad%n2o_layer(nlevels-1, iprof)
      ENDIF

      IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
        prof_ad(iprof)%co2(1) = prof_ad(iprof)%co2(1) + &
          0.5_jprb * aux_ad%co2_layer(1, iprof)
        prof_ad(iprof)%co2(2:nlevels-1) = prof_ad(iprof)%co2(2:nlevels-1) + &
          0.5_jprb * (aux_ad%co2_layer(1:nlevels-2, iprof) + aux_ad%co2_layer(2:nlevels-1, iprof))
        prof_ad(iprof)%co2(nlevels) = prof_ad(iprof)%co2(nlevels) + &
          0.5_jprb * aux_ad%co2_layer(nlevels-1, iprof)
      ENDIF
    ENDIF

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      prof_ad(iprof)%o3(1) = prof_ad(iprof)%o3(1) + &
        0.5_jprb * aux_ad%o3_layer(1, iprof)
      prof_ad(iprof)%o3(2:nlevels-1) = prof_ad(iprof)%o3(2:nlevels-1) + &
        0.5_jprb * (aux_ad%o3_layer(1:nlevels-2, iprof) + aux_ad%o3_layer(2:nlevels-1, iprof))
      prof_ad(iprof)%o3(nlevels) = prof_ad(iprof)%o3(nlevels) + &
        0.5_jprb * aux_ad%o3_layer(nlevels-1, iprof)
    ENDIF

    IF (opts%rt_all%use_q2m) THEN
      iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level below surface

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        prof_ad(iprof)%s2m%q =  prof_ad(iprof)%s2m%q + &
          aux_ad%w_layer(iv2lay, iprof) * 0.5_jprb

        ! This line subtracts the quantity that is added to prof_ad(iprof)%q(iv2lev) below
        ! so that it cancels out: in the case where use_q2m is true prof(iprof)%q(iv2lev) is
        ! replaced by prof(iprof)%s2m%q in the calculation of aux%w_layer(iv2lay, iprof)

        prof_ad(iprof)%q(iv2lev) = prof_ad(iprof)%q(iv2lev) - &
          aux_ad%w_layer(iv2lay, iprof) * 0.5_jprb
      ENDIF
    ENDIF

    prof_ad(iprof)%q(1) = prof_ad(iprof)%q(1) + &
      0.5_jprb * aux_ad%w_layer(1, iprof)
    prof_ad(iprof)%q(2:nlevels-1) = prof_ad(iprof)%q(2:nlevels-1) + &
      0.5_jprb * (aux_ad%w_layer(1:nlevels-2, iprof) + aux_ad%w_layer(2:nlevels-1, iprof))
    prof_ad(iprof)%q(nlevels) = prof_ad(iprof)%q(nlevels) + &
      0.5_jprb * aux_ad%w_layer(nlevels-1, iprof)

    IF (opts%rt_all%use_t2m_opdep) THEN
      iv2lev = aux%s(iprof)%nearestlev_surf    ! nearest level below surface

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        prof_ad(iprof)%s2m%t =  prof_ad(iprof)%s2m%t + &
          aux_ad%t_layer(iv2lay, iprof) * 0.5_jprb

        ! This line subtracts the quantity that is added to prof_ad(iprof)%t(iv2lev) below
        ! so that it cancels out: in the case where use_t2m_opdep is true prof(iprof)%t(iv2lev) is
        ! replaced by prof(iprof)%s2m%t in the calculation of aux%t_layer(iv2lay, iprof)

        prof_ad(iprof)%t(iv2lev) = prof_ad(iprof)%t(iv2lev) - &
          aux_ad%t_layer(iv2lay, iprof) * 0.5_jprb
      ENDIF
    ENDIF

    prof_ad(iprof)%t(1) = prof_ad(iprof)%t(1) + &
      0.5_jprb * aux_ad%t_layer(1, iprof)
    prof_ad(iprof)%t(2:nlevels-1) = prof_ad(iprof)%t(2:nlevels-1) + &
      0.5_jprb * (aux_ad%t_layer(1:nlevels-2, iprof) + aux_ad%t_layer(2:nlevels-1, iprof))
    prof_ad(iprof)%t(nlevels) = prof_ad(iprof)%t(nlevels) + &
      0.5_jprb * aux_ad%t_layer(nlevels-1, iprof)
  ENDDO

  ! -------------------------------------------------------------------------
  ! Raytracing path data required for predictor calculations
  ! -------------------------------------------------------------------------

  IF (plane_parallel) THEN
    aux_ad%pathsat_sqrt  = 0._jprb
    aux_ad%pathsat_rsqrt = 0._jprb
  ELSE
    raytracing_ad%pathsat = raytracing_ad%pathsat + &
      aux_ad%pathsat_sqrt * aux%pathsat_rsqrt

    aux_ad%pathsat_rsqrt = &!aux_ad%pathsat_rsqrt + &
      raytracing%pathsat * aux_ad%pathsat_sqrt

    CALL invsqrt_ad(aux%pathsat_4rt, aux_ad%pathsat_rsqrt, aux_ad%pathsat_4rt, acc=.TRUE._jplm)
    CALL invsqrt_ad(aux%pathsat_rsqrt, raytracing_ad%pathsat, aux_ad%pathsat_rsqrt, acc=.TRUE._jplm)

    IF (dosolar) THEN
      raytracing_ad%patheff = raytracing_ad%patheff + &
        aux_ad%patheff_sqrt * aux%patheff_rsqrt

      aux_ad%patheff_rsqrt = &!aux_ad%patheff_rsqrt + &
        raytracing%patheff * aux_ad%patheff_sqrt

      CALL invsqrt_ad(aux%patheff_4rt, aux_ad%patheff_rsqrt, aux_ad%patheff_4rt, acc=.TRUE._jplm)
      CALL invsqrt_ad(aux%patheff_rsqrt, raytracing_ad%patheff, aux_ad%patheff_rsqrt, acc=.TRUE._jplm)
    ENDIF
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_789_AD', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_ssu_co2_thickness_ad
    USE rttov_const, ONLY : p0, t0, flatt, eqrad, omega, grave
    REAL(jprb) :: rearth(SIZE(prof))
    REAL(jprb) :: gravh(SIZE(prof))
    REAL(jprb) :: gravl(SIZE(prof)) ! gravity at latitude lat [m/s^2]
    REAL(jprb) :: eta, beta ! coefficients of the international gravity formula
    REAL(jprb) :: dflat, fac
    REAL(jprb) :: pstd
    REAL(jprb) :: tstd
    REAL(jprb) :: pa(prof(1)%nlayers),    pa_ad(prof(1)%nlayers)
    REAL(jprb) :: dpa(prof(1)%nlayers),   dpa_ad(prof(1)%nlayers)
    REAL(jprb) :: rho(prof(1)%nlayers),   rho_ad(prof(1)%nlayers)
    REAL(jprb) :: temp(prof(1)%nlayers),  temp_ad(prof(1)%nlayers)
    REAL(jprb) :: co2(prof(1)%nlayers),   co2_ad(prof(1)%nlayers)
    REAL(jprb) :: dz, dz_ad
    REAL(jprb) :: dzstp_ad
    INTEGER(jpim) :: iprof, lay

    pstd = p0 * 100._jprb ! pa
    tstd = t0             ! K

    dflat = (1._jprb - flatt) ** 2_jpim
    fac = (omega ** 2_jpim * (eqrad * 1000._jprb)) / (grave)
    beta = 5._jprb * fac / 2._jprb - flatt - 17._jprb * fac * flatt / 14._jprb
    eta = flatt * (5._jprb * fac - flatt) / 8._jprb

    ! Calculate the earth's radius at latitude lat assuming the
    ! earth is an ellipsoid of revolution
    rearth(:) = SQRT(eqrad ** 2_jpim * dflat / &
                 (angles(:)%sinlat ** 2_jpim + &
                  dflat * angles(:)%coslat ** 2_jpim))

    ! The value of earth's gravity at surface at latitude lat is
    ! computed using the international gravity formula.
    gravl(:) = grave * &
      (1._jprb + beta * angles(:)%sinlat ** 2_jpim +      &
      eta * (2._jprb * angles(:)%sinlat * angles(:)%coslat) ** 2_jpim)

    gravh(:) = 2.E-3_jprb * gravl(:) / rearth(:)

    DO iprof = 1, nprofiles

      DO lay = 1, nlayers
! fwd
        pa(lay)   = 100._jprb * (prof(iprof)%p(lay) + prof(iprof)%p(lay+1)) / 2._jprb
        dpa(lay)  = 100._jprb * (prof(iprof)%p(lay) - prof(iprof)%p(lay+1))
        rho(lay)  = (raytracing%dmair(lay,iprof) + raytracing%dmair(lay+1,iprof)) / 2._jprb
        temp(lay) = (prof(iprof)%t(lay) + prof(iprof)%t(lay+1)) / 2._jprb
        co2(lay)  = (prof(iprof)%co2(lay) + prof(iprof)%co2(lay+1)) / 2._jprb
! ad
        pa_ad(lay)   = 0._jprb
        dpa_ad(lay)  = 0._jprb
        rho_ad(lay)  = 0._jprb
        temp_ad(lay) = 0._jprb
        co2_ad(lay)  = 0._jprb
      ENDDO

      DO lay = 1, nlayers
! fwd
        dz = -dpa(lay)/&
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1,iprof)) * rho(lay) )
! ad
        dzstp_ad = 0._jprb
        dzstp_ad = dzstp_ad + aux_ad%co2_cm(iprof)
        dz_ad = 0._jprb
        dz_ad =  dz_ad + &
          dzstp_ad * pa(lay)/pstd * tstd/temp(lay)* co2(lay)*1.E-6_jprb*100._jprb

        IF (opts%interpolation%lgradp) THEN
          pa_ad(lay) = pa_ad(lay) + &
            dzstp_ad * dz / pstd * tstd/temp(lay)* co2(lay)*1.E-6_jprb*100._jprb
!        ELSE
!          pa_ad(lay) = 0._jprb
        ENDIF

        temp_ad(lay) = temp_ad(lay) - &
          dzstp_ad * dz * pa(lay)/pstd * tstd/temp(lay)**2_jpim * &
          co2(lay)*1.E-6_jprb*100._jprb

        co2_ad(lay) = co2_ad(lay) + &
          dzstp_ad * dz * pa(lay)/pstd * tstd/temp(lay)*1.E-6_jprb*100._jprb
!        dzstp_ad = 0._jprb

        IF (opts%interpolation%lgradp) THEN
          dpa_ad(lay) = dpa_ad(lay) - dz_ad / ((gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1,iprof)) * rho(lay))
!        ELSE
!          dpa_ad(lay) = 0._jprb
        ENDIF

        raytracing_ad%hgpl(lay+1,iprof) = raytracing_ad%hgpl(lay+1,iprof) &
          - dz_ad * dpa(lay) * gravh(iprof)*rho(lay)/ &
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1,iprof)) * rho(lay) )**2_jpim

        rho_ad(lay) = rho_ad(lay)                                              &
          + dz_ad * dpa(lay)*(gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1,iprof)) / &
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1,iprof)) * rho(lay) )**2_jpim

!        dz_ad = 0._jprb

        prof_ad(iprof)%co2(lay) = prof_ad(iprof)%co2(lay) + &
          co2_ad(lay)/2._jprb
        prof_ad(iprof)%co2(lay+1) = prof_ad(iprof)%co2(lay+1) + &
          co2_ad(lay)/2._jprb
!        co2_ad(lay) = 0._jprb

        prof_ad(iprof)%t(lay) = prof_ad(iprof)%t(lay) + &
          temp_ad(lay)/2._jprb

        prof_ad(iprof)%t(lay+1) = prof_ad(iprof)%t(lay+1) + &
          temp_ad(lay)/2._jprb
!        temp_ad(lay) = 0._jprb

        raytracing_ad%dmair(lay,iprof) = raytracing_ad%dmair(lay,iprof) + rho_ad(lay) / 2._jprb
        raytracing_ad%dmair(lay+1,iprof) = raytracing_ad%dmair(lay+1,iprof) + rho_ad(lay) / 2._jprb

!        rho_ad(lay) = 0._jprb

        IF (opts%interpolation%lgradp) THEN
          prof_ad(iprof)%p(lay) = prof_ad(iprof)%p(lay) + &
            dpa_ad(lay)*100._jprb
          prof_ad(iprof)%p(lay+1) = prof_ad(iprof)%p(lay+1) - &
            dpa_ad(lay)*100._jprb
!          dpa_ad(lay) = 0._jprb
          prof_ad(iprof)%p(lay) = prof_ad(iprof)%p(lay) + &
            pa_ad(lay)*100._jprb/2._jprb
          prof_ad(iprof)%p(lay+1) = prof_ad(iprof)%p(lay+1) + &
            pa_ad(lay)*100._jprb/2._jprb
!          pa_ad(lay) = 0._jprb
!        ELSE
!          prof_ad(iprof)%p(lay) = 0._jprb
!          prof_ad(iprof)%p(lay+1) = 0._jprb
!          dpa_ad(lay) = 0._jprb
!          pa_ad(lay) = 0._jprb
        ENDIF
      ENDDO
    END DO

  END SUBROUTINE calc_ssu_co2_thickness_ad

END SUBROUTINE rttov_predictor_precalc_789_ad
