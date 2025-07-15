! Description:
!> @file
!!   TL of v13 predictor pre-calculations.
!
!> @brief
!!   TL of v13 predictor pre-calculations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     dosolar              flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel       flag for strict plane parallel geometry
!! @param[in]     prof                 input profiles
!! @param[in]     prof_tl              input profile perturbations
!! @param[in]     coef                 optical depth coefficient structure
!! @param[in]     coef_pccomp          PC coefficients structure
!! @param[in]     aux                  coef level auxiliary profile data structure
!! @param[in,out] aux_tl               coef level auxiliary profile data perturbations
!! @param[in]     raytracing           raytracing structure
!! @param[in,out] raytracing_tl        raytracing structure perturbations
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
SUBROUTINE rttov_predictor_precalc_13_tl( &
              opts,           &
              dosolar,        &
              plane_parallel, &
              prof,           &
              prof_tl,        &
              coef,           &
              coef_pccomp,    &
              aux,            &
              aux_tl,         &
              raytracing,     &
              raytracing_tl)

  USE rttov_types, ONLY :  &
        rttov_options,          &
        rttov_coef,             &
        rttov_coef_pccomp,      &
        rttov_profile,          &
        rttov_profile_aux_coef, &
        rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY: invsqrt_tl, reciprocal_tl, sqrt_tl
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),          INTENT(IN)    :: opts
  LOGICAL(jplm),                INTENT(IN)    :: dosolar
  LOGICAL(jplm),                INTENT(IN)    :: plane_parallel
  TYPE(rttov_profile),          INTENT(IN)    :: prof(:)
  TYPE(rttov_profile),          INTENT(IN)    :: prof_tl(:)
  TYPE(rttov_coef),             INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp),      INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_tl
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_tl
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers, nlevels
  INTEGER(jpim) :: iprof, lay
  REAL(jprb) :: ZHOOK_HANDLE

  REAL(jprb) :: sum1, sum2
  INTEGER(jpim) :: iv2lay, iv2lev, iv3lev

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_13_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlevels

!-----------------------------------------------------------------------------------------
! Calculate predictor data
!-----------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ! Raytracing path data required for predictor calculations
  ! -------------------------------------------------------------------------

  IF (plane_parallel) THEN
    aux_tl%pathsat_sqrt  = 0._jprb
    aux_tl%pathsat_rsqrt = 0._jprb
    aux_tl%pathsat_4rt   = 0._jprb

    IF (dosolar) THEN
      aux_tl%patheff_sqrt  = 0._jprb
      aux_tl%patheff_rsqrt = 0._jprb
      aux_tl%patheff_4rt   = 0._jprb
    ENDIF
  ELSE
    CALL invsqrt_tl(aux%pathsat_rsqrt, raytracing_tl%pathsat, aux_tl%pathsat_rsqrt)
    CALL invsqrt_tl(aux%pathsat_4rt, aux_tl%pathsat_rsqrt, aux_tl%pathsat_4rt)
    aux_tl%pathsat_sqrt = raytracing_tl%pathsat * aux%pathsat_rsqrt + &
                          raytracing%pathsat * aux_tl%pathsat_rsqrt
    IF (dosolar) THEN
      CALL invsqrt_tl(aux%patheff_rsqrt, raytracing_tl%patheff, aux_tl%patheff_rsqrt)
      CALL invsqrt_tl(aux%patheff_4rt, aux_tl%patheff_rsqrt, aux_tl%patheff_4rt)
      aux_tl%patheff_sqrt = raytracing_tl%patheff * aux%patheff_rsqrt + &
                            raytracing%patheff * aux_tl%patheff_rsqrt
    ENDIF
  ENDIF

  ! -------------------------------------------------------------------------
  ! Profile layer quantities: layer n lies between levels n and n+1
  ! -------------------------------------------------------------------------

  DO iprof = 1, nprofiles
    aux_tl%t_layer(1:nlayers, iprof) = &
      (prof_tl(iprof)%t(1:nlevels-1) + prof_tl(iprof)%t(2:nlevels)) * 0.5_jprb

    IF (opts%rt_all%use_t2m_opdep) THEN
      iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        aux_tl%t_layer(iv2lay, iprof) = &
          (prof_tl(iprof)%s2m%t + prof_tl(iprof)%t(iv3lev)) * 0.5_jprb
      ENDIF
    ENDIF

    aux_tl%w_layer(1:nlayers, iprof) = &
      (prof_tl(iprof)%q(1:nlevels-1) + prof_tl(iprof)%q(2:nlevels)) * 0.5_jprb

    IF (opts%rt_all%use_q2m) THEN
      iv3lev = aux%s(iprof)%nearestlev_surf - 1! nearest level above surface
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        aux_tl%w_layer(iv2lay, iprof) = &
          (prof_tl(iprof)%s2m%q + prof_tl(iprof)%q(iv3lev)) * 0.5_jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux_tl%o3_layer(1:nlayers, iprof) = &
        (prof_tl(iprof)%o3(1:nlevels-1) + prof_tl(iprof)%o3(2:nlevels)) * 0.5_jprb
    ENDIF
    
    IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux_tl%co2_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%co2(1:nlevels-1) + prof_tl(iprof)%co2(2:nlevels)) * 0.5_jprb
      ELSE
        aux_tl%co2_layer(1:nlayers, iprof) = 0._jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux_tl%n2o_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%n2o(1:nlevels-1) + prof_tl(iprof)%n2o(2:nlevels)) * 0.5_jprb
      ELSE
        aux_tl%n2o_layer(1:nlayers, iprof) = 0._jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux_tl%co_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%co(1:nlevels-1) + prof_tl(iprof)%co(2:nlevels)) * 0.5_jprb
      ELSE
        aux_tl%co_layer(1:nlayers, iprof) = 0._jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux_tl%ch4_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%ch4(1:nlevels-1) + prof_tl(iprof)%ch4(2:nlevels)) * 0.5_jprb
      ELSE
        aux_tl%ch4_layer(1:nlayers, iprof) = 0._jprb
      ENDIF
    ENDIF
    
    IF (coef%nso2 > 0) THEN
      IF (opts%rt_all%so2_data) THEN
        aux_tl%so2_layer(1:nlayers, iprof) = &
          (prof_tl(iprof)%so2(1:nlevels-1) + prof_tl(iprof)%so2(2:nlevels)) * 0.5_jprb
      ELSE
        aux_tl%so2_layer(1:nlayers, iprof) = 0._jprb
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    ! Calculate deviations from reference profile for each layer
    ! -----------------------------------------------------------------------

    aux_tl%dt(:, iprof) = aux_tl%t_layer(:, iprof)

    ! -----------------------------------------------------------------------
    ! Calculate (profile / reference profile) ratios
    ! -----------------------------------------------------------------------
    
    aux_tl%tr(:, iprof) = aux_tl%t_layer(:, iprof) * coef%tstar_r(:)
    aux_tl%tr2(:, iprof) = 2._jprb * aux%tr(:, iprof) * aux_tl%tr(:, iprof)
    aux_tl%wr(:, iprof) = aux_tl%w_layer(:, iprof) * coef%wstar_r(:)
    
    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux_tl%or(:, iprof) = aux_tl%o3_layer(:, iprof) * coef%ostar_r(:)
    ENDIF

    IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
      aux_tl%co2r(:, iprof) = aux_tl%co2_layer(:, iprof) * coef%co2star_r(:)
      CALL sqrt_tl(aux%co2r_sqrt(:, iprof), aux_tl%co2r(:, iprof), aux_tl%co2r_sqrt(:, iprof))
    ENDIF

    IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
      aux_tl%n2or(:, iprof) = aux_tl%n2o_layer(:, iprof) * coef%n2ostar_r(:)
      CALL sqrt_tl(aux%n2or_sqrt(:, iprof), aux_tl%n2or(:, iprof), aux_tl%n2or_sqrt(:, iprof))
      CALL sqrt_tl(aux%n2or_4rt(:, iprof), aux_tl%n2or_sqrt(:, iprof), aux_tl%n2or_4rt(:, iprof))
    ENDIF

    IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
      aux_tl%cor(:, iprof) = aux_tl%co_layer(:, iprof) * coef%costar_r(:)
      CALL sqrt_tl(aux%cor_sqrt(:, iprof), aux_tl%cor(:, iprof), aux_tl%cor_sqrt(:, iprof))
      CALL sqrt_tl(aux%cor_4rt(:, iprof), aux_tl%cor_sqrt(:, iprof), aux_tl%cor_4rt(:, iprof))
    ENDIF

    IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
      aux_tl%ch4r(:, iprof) = aux_tl%ch4_layer(:, iprof) * coef%ch4star_r(:)
      CALL sqrt_tl(aux%ch4r_sqrt(:, iprof), aux_tl%ch4r(:, iprof), aux_tl%ch4r_sqrt(:, iprof))
      CALL sqrt_tl(aux%ch4r_4rt(:, iprof), aux_tl%ch4r_sqrt(:, iprof), aux_tl%ch4r_4rt(:, iprof))
    ENDIF

    IF (coef%nso2 > 0) THEN
      aux_tl%so2r(:, iprof) = aux_tl%so2_layer(:, iprof) * coef%so2star_r(:)
    ENDIF
    
  ENDDO

  ! -------------------------------------------------------------------------
  ! Calculate profile / reference profile sums
  ! -------------------------------------------------------------------------

  DO iprof = 1, nprofiles
    ! cumulating column overlying layer and layer itself
    sum1 = 0._jprb
    DO lay = 1, nlayers
      ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
      ! need dpp(0) to start
      sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%w_layer(lay, iprof)
      aux_tl%ww(lay, iprof) = sum1 * coef%wstar_wsum_r(lay)
    ENDDO
  ENDDO

  CALL reciprocal_tl(aux%ww_r, aux_tl%ww, aux_tl%ww_r)
  aux_tl%wrw_r = aux%wr * aux_tl%ww_r + aux_tl%wr * aux%ww_r

  DO iprof = 1, nprofiles
    sum1 = 0._jprb
    DO lay = 1, nlayers
      sum1 = sum1 + coef%dpp(lay - 1) * (aux_tl%w_layer(lay, iprof) * aux%t_layer(lay, iprof) + &
                                         aux%w_layer(lay, iprof) * aux_tl%t_layer(lay, iprof))
      aux_tl%wwr(lay, iprof) = sum1 * coef%wtstar_wsum_r(lay)
    ENDDO
  ENDDO

  CALL reciprocal_tl(aux%wwr_r, aux_tl%wwr, aux_tl%wwr_r)
  aux_tl%wrwr_r = aux%wr * aux_tl%wwr_r + aux_tl%wr * aux%wwr_r

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO iprof = 1, nprofiles
        sum1 = 0._jprb
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%o3_layer(lay, iprof) ! OK, dpp(0) defined 
          aux_tl%ow(lay, iprof) = sum1 * coef%ostar_wsum_r(lay)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO iprof = 1, nprofiles
    aux_tl%twr(1, iprof) = aux_tl%t_layer(1, iprof) * coef%tstar_r(1)
    sum1 = coef%dpp(0) * aux_tl%t_layer(1, iprof)

    DO lay = 2, nlayers
      ! cumulate overlying layers (t, tstar relate to same layer as dpp)
      ! do not need dpp(0) to start
      sum1       = sum1 + coef%dpp(lay - 1) * aux_tl%t_layer(lay, iprof)
      aux_tl%twr(lay, iprof) = sum1 * coef%tstar_wsum_r(lay)
    ENDDO
  ENDDO

  IF (coef%nco2 > 0) THEN
    DO iprof = 1, nprofiles
      sum1 = 0._jprb
      IF (opts%rt_all%co2_data) THEN
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%co2_layer(lay, iprof) ! OK, dpp(0) defined 
          aux_tl%co2w(lay, iprof) = sum1 * coef%co2star_wsum_r(lay)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (coef%nn2o > 0) THEN
    DO iprof = 1, nprofiles         
      sum1 = 0._jprb
      sum2 = 0._jprb
      IF (opts%rt_all%n2o_data) THEN
        DO lay = 1, nlayers
          sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%n2o_layer(lay, iprof)
          aux_tl%n2ow(lay, iprof) = sum1 * coef%n2ostar_wsum_r(lay)
          sum2 = sum2 + coef%dpp(lay - 1) * (aux_tl%n2o_layer(lay, iprof) * aux%t_layer(lay, iprof) + &
                                             aux%n2o_layer(lay, iprof) * aux_tl%t_layer(lay, iprof))
          aux_tl%n2owr(lay, iprof) = sum2 * coef%n2otstar_wsum_r(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          sum2 = sum2 + coef%dpp(lay - 1) * coef%n2ostar(lay) * aux_tl%t_layer(lay, iprof)
          aux_tl%n2owr(lay, iprof) = sum2 * coef%n2otstar_wsum_r(lay)
        ENDDO
      ENDIF
    ENDDO
    IF (opts%rt_all%n2o_data) CALL reciprocal_tl(aux%n2ow_r, aux_tl%n2ow, aux_tl%n2ow_r)
  ENDIF

  IF (coef%nco > 0) THEN
    aux_tl%dtabsdt = 2._jprb * ABS(aux%dt) * aux_tl%dt

    DO iprof = 1, nprofiles
      sum1 = 0._jprb
      sum2 = 0._jprb
      IF (opts%rt_all%co_data) THEN
        DO lay = 1, nlayers
          sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%co_layer(lay, iprof)
          aux_tl%cow(lay, iprof) = sum1 * coef%costar_wsum_r(lay)
          sum2 = sum2 + coef%dpp(lay - 1) * (aux_tl%co_layer(lay, iprof) * aux%t_layer(lay, iprof) + &
                                             aux%co_layer(lay, iprof) * aux_tl%t_layer(lay, iprof))
          aux_tl%cowr(lay, iprof) = sum2 * coef%cotstar_wsum_r(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          sum2 = sum2 + coef%dpp(lay - 1) * coef%costar(lay) * aux_tl%t_layer(lay, iprof)
          aux_tl%cowr(lay, iprof) = sum2 * coef%cotstar_wsum_r(lay)
        ENDDO
      ENDIF

    ENDDO
    aux_tl%cowr_4rt = 0.25_jprb * (aux%cowr_4rt * aux%cowr_r) * aux_tl%cowr
  ENDIF

  IF (coef%nch4 > 0) THEN
    DO iprof = 1, nprofiles
      sum1 = 0._jprb
      sum2 = 0._jprb
      IF (opts%rt_all%ch4_data) THEN
        DO lay = 1, nlayers
          sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%ch4_layer(lay, iprof)
          aux_tl%ch4w(lay, iprof) = sum1 * coef%ch4star_wsum_r(lay)
          sum2 = sum2 + coef%dpp(lay - 1) * (aux_tl%ch4_layer(lay, iprof) * aux%t_layer(lay, iprof) + &
                                             aux%ch4_layer(lay, iprof) * aux_tl%t_layer(lay, iprof))
          aux_tl%ch4wr(lay, iprof) = sum2 * coef%ch4tstar_wsum_r(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          sum2 = sum2 + coef%dpp(lay - 1) * coef%ch4star(lay) * aux_tl%t_layer(lay, iprof)
          aux_tl%ch4wr(lay, iprof) = sum2 * coef%ch4tstar_wsum_r(lay)
        ENDDO
      ENDIF
    ENDDO
    IF (opts%rt_all%ch4_data) THEN
      CALL reciprocal_tl(aux%ch4w_r, aux_tl%ch4w, aux_tl%ch4w_r)
      aux_tl%ch4rw_r = aux_tl%ch4r * aux%ch4w_r + aux%ch4r * aux_tl%ch4w_r
      aux_tl%ch4w_4rt = 0.25_jprb * (aux%ch4w_4rt * aux%ch4w_r) * aux_tl%ch4w
    ENDIF
  ENDIF

  IF (coef%nso2 > 0) THEN
    DO iprof = 1, nprofiles
      sum1 = 0._jprb
      sum2 = 0._jprb
      DO lay = 1, nlayers
        ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
        ! need dpp(0) to start
        sum1 = sum1 + coef%dpp(lay - 1) * aux_tl%so2_layer(lay, iprof) ! OK, dpp(0) defined 
        aux_tl%so2w(lay, iprof) = sum1 * coef%so2star_wsum_r(lay)

        sum2 = sum2 + coef%dpp(lay - 1) * ( &
                                           aux_tl%so2_layer(lay, iprof) * aux%t_layer(lay, iprof) + &
                                           aux%so2_layer(lay, iprof) * aux_tl%t_layer(lay, iprof))
        aux_tl%so2wr(lay, iprof) = sum2 * coef%so2tstar_wsum_r(lay)
      ENDDO
    ENDDO

    CALL reciprocal_tl(aux%so2w_r, aux_tl%so2w, aux_tl%so2w_r)
    CALL reciprocal_tl(aux%so2wr_r, aux_tl%so2wr, aux_tl%so2wr_r)
    aux_tl%so2rwr_r = aux_tl%so2r * aux%so2wr_r + aux%so2r * aux_tl%so2wr_r
    aux_tl%so2rw_r = aux_tl%so2r * aux%so2w_r + aux%so2r * aux_tl%so2w_r
    aux_tl%so2r_sqrt = 0.5_jprb * aux_tl%so2r / aux%so2r_sqrt
    aux_tl%so2r_4rt = 0.5_jprb * aux_tl%so2r_sqrt / aux%so2r_4rt
    aux_tl%so2w_sqrt = 0.5_jprb * aux_tl%so2w / aux%so2w_sqrt
    aux_tl%so2w_4rt = 0.5_jprb * aux_tl%so2w_sqrt / aux%so2w_4rt
  ENDIF

  ! -------------------------------------------------------------------------
  ! Calculate and store useful intermediate variables
  ! -------------------------------------------------------------------------

  CALL reciprocal_tl(aux%tr_r, aux_tl%tr, aux_tl%tr_r)

  CALL sqrt_tl(aux%tr_sqrt, aux_tl%tr, aux_tl%tr_sqrt)

  CALL invsqrt_tl(aux%wr_rsqrt, aux_tl%wr, aux_tl%wr_rsqrt)
  CALL invsqrt_tl(aux%wr_4rt, aux_tl%wr_rsqrt, aux_tl%wr_4rt)
  aux_tl%wr_sqrt = aux_tl%wr * aux%wr_rsqrt + &
                   aux%wr * aux_tl%wr_rsqrt

  CALL invsqrt_tl(aux%ww_rsqrt, aux_tl%ww, aux_tl%ww_rsqrt)
  aux_tl%ww_sqrt = aux%ww * aux_tl%ww_rsqrt + &
                   aux_tl%ww * aux%ww_rsqrt
  CALL invsqrt_tl(aux%ww_4rt, aux_tl%ww_rsqrt, aux_tl%ww_4rt)

  IF (coef%nozone > 0 .AND. opts%rt_all%ozone_data) THEN
    CALL invsqrt_tl(aux%ow_rsqrt, aux_tl%ow, aux_tl%ow_rsqrt)
    aux_tl%ow_r = 2._jprb * aux%ow_rsqrt * aux_tl%ow_rsqrt
    CALL sqrt_tl(aux%or_sqrt, aux_tl%or, aux_tl%or_sqrt)
    aux_tl%ow_sqrt = aux%ow * aux_tl%ow_rsqrt + aux_tl%ow * aux%ow_rsqrt
  ENDIF

  IF (coef%nco > 0 .AND. opts%rt_all%co_data) THEN
    CALL invsqrt_tl(aux%cow_rsqrt, aux_tl%cow, aux_tl%cow_rsqrt)
    aux_tl%corw_rsqrt = aux%cor * aux_tl%cow_rsqrt + &
                        aux_tl%cor * aux%cow_rsqrt
    aux_tl%corw_r = aux%cow_rsqrt * (aux%cor * 2._jprb * aux_tl%cow_rsqrt + &
                                     aux_tl%cor * aux%cow_rsqrt)
  ENDIF    

  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_13_TL', 1_jpim, ZHOOK_HANDLE)

END SUBROUTINE rttov_predictor_precalc_13_tl
