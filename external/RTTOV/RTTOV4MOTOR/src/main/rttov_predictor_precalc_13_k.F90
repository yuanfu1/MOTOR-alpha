! Description:
!> @file
!!   K of v13 predictor pre-calculations.
!
!> @brief
!!   K of v13 predictor pre-calculations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     dosolar              flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel       flag for strict plane parallel geometry
!! @param[in]     chanprof             channels and profiles to simulate
!! @param[in]     profiles             input profiles
!! @param[in,out] profiles_k           input profile increments
!! @param[in]     coef                 optical depth coefficient structure
!! @param[in]     coef_pccomp          PC coefficients structure
!! @param[in]     aux                  coef level auxiliary profile data structure
!! @param[in,out] aux_k                coef level auxiliary profile data increments
!! @param[in]     raytracing           raytracing structure
!! @param[in,out] raytracing_k         raytracing structure increments
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
SUBROUTINE rttov_predictor_precalc_13_k( &
              opts,           &
              dosolar,        &
              plane_parallel, &
              chanprof,       &
              profiles,       &
              profiles_k,     &
              coef,           &
              coef_pccomp,    &
              aux,            &
              aux_k,          &
              raytracing,     &
              raytracing_k)

  USE rttov_types, ONLY :  &
        rttov_chanprof,         &
        rttov_options,          &
        rttov_coef,             &
        rttov_coef_pccomp,      &
        rttov_profile,          &
        rttov_profile_aux_coef, &
        rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_math_mod, ONLY: invsqrt_k, reciprocal_k, sqrt_k
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),          INTENT(IN)    :: opts
  LOGICAL(jplm),                INTENT(IN)    :: dosolar
  LOGICAL(jplm),                INTENT(IN)    :: plane_parallel
  TYPE(rttov_chanprof),         INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),          INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),          INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_coef),             INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp),      INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile_aux_coef), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux_coef), INTENT(INOUT) :: aux_k
  TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_k
!INTF_END

  INTEGER(jpim) :: nchannels, nlayers, nlevels
  INTEGER(jpim) :: i, prof, lay
  REAL(jprb) :: ZHOOK_HANDLE

  REAL(jprb) :: sum1_k, sum2_k
  INTEGER(jpim) :: iv2lay, iv2lev
  INTEGER(jpim) :: map(SIZE(chanprof),2), prof_stat

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_13_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = profiles(1)%nlevels

  map(1,1) = chanprof(1)%prof
  map(1,2) = chanprof(1)%chan

  prof_stat = 1 ! assume profs are contiguous and monotonic
  DO i = 2, nchannels
    map(i,1) = chanprof(i)%prof
    map(i,2) = chanprof(i)%chan

    IF (map(i,1) < map(i-1,1)) THEN ! they are not.
      prof_stat = -1
    ENDIF
  ENDDO

!-----------------------------------------------------------------------------------------
! K of predictor data
!-----------------------------------------------------------------------------------------

  aux_k%t_layer = 0._jprb

  ! -------------------------------------------------------------------------
  ! Calculate and store useful intermediate variables
  ! -------------------------------------------------------------------------

  IF (coef%nco > 0 .AND. opts%rt_all%co_data) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof

      aux_k%cor(:,i) = &!aux_k%cor + &
        aux%cow_rsqrt(:, prof) * (aux%cow_rsqrt(:, prof) * aux_k%corw_r(:,i) + &
                                         aux_k%corw_rsqrt(:,i))

      aux_k%cow_rsqrt(:,i) = &!aux_k%cow_rsqrt(:,i) + &
        aux%cor(:, prof) * (2._jprb * aux%cow_rsqrt(:, prof) * aux_k%corw_r(:,i) + &
                                              aux_k%corw_rsqrt(:,i))
    ENDDO
    CALL invsqrt_k(aux%cow_rsqrt, aux_k%cow, aux_k%cow_rsqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
  ENDIF

  IF (coef%nozone > 0 .AND. opts%rt_all%ozone_data) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      aux_k%ow_rsqrt(:,i) = &!aux_k%ow_rsqrt(:,i) + &
        aux%ow(:, prof) * aux_k%ow_sqrt(:,i) + &
        2._jprb * aux%ow_rsqrt(:, prof) * aux_k%ow_r(:,i)

      aux_k%ow(:,i) = aux_k%ow(:,i) + aux%ow_rsqrt(:, prof) * aux_k%ow_sqrt(:,i)
    END DO

    CALL sqrt_k(aux%or_sqrt, aux_k%or, aux_k%or_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

    CALL invsqrt_k(aux%ow_rsqrt, aux_k%ow, aux_k%ow_rsqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
  ENDIF

  CALL invsqrt_k(aux%ww_4rt, aux_k%ww_rsqrt, aux_k%ww_4rt, acc = .FALSE._jplm, map = map, prof_stat = prof_stat)
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    aux_k%ww(:,i)= aux_k%ww(:,i) + aux%ww_rsqrt(:, prof) * aux_k%ww_sqrt(:,i)
    aux_k%ww_rsqrt(:,i)= aux_k%ww_rsqrt(:,i) + aux%ww(:, prof) * aux_k%ww_sqrt(:,i)
  END DO
  CALL invsqrt_k(aux%ww_rsqrt, aux_k%ww, aux_k%ww_rsqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  DO i = 1, nchannels
    prof = chanprof(i)%prof

    aux_k%wr(:,i)= aux_k%wr(:,i) + aux%wr_rsqrt(:, prof) * aux_k%wr_sqrt(:,i)
    aux_k%wr_rsqrt(:,i)= &!aux_k%wr_rsqrt(:,i) + &
      aux%wr(:, prof) * aux_k%wr_sqrt(:,i)
  ENDDO

  CALL invsqrt_k(aux%wr_4rt, aux_k%wr_rsqrt, aux_k%wr_4rt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
  CALL invsqrt_k(aux%wr_rsqrt, aux_k%wr, aux_k%wr_rsqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  CALL sqrt_k(aux%tr_sqrt, aux_k%tr, aux_k%tr_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  CALL reciprocal_k(aux%tr_r, aux_k%tr, aux_k%tr_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
   
  ! -------------------------------------------------------------------------
  ! Calculate profile / reference profile sums
  ! -------------------------------------------------------------------------

  IF (coef%nso2 > 0) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof

      aux_k%so2w_sqrt(:,i) = aux_k%so2w_sqrt(:,i) + &
        0.5_jprb * aux_k%so2w_4rt(:,i) / aux%so2w_4rt(:,prof)

      aux_k%so2w(:,i) = &!aux_k%so2w(:,i) + &
        0.5_jprb * aux_k%so2w_sqrt(:,i) / aux%so2w_sqrt(:,prof)

      aux_k%so2r_sqrt(:,i) = aux_k%so2r_sqrt(:,i) + &
        0.5_jprb * aux_k%so2r_4rt(:,i) / aux%so2r_4rt(:,prof)

      aux_k%so2r(:,i) = aux_k%so2r(:,i) + &
        0.5_jprb * aux_k%so2r_sqrt(:,i) / aux%so2r_sqrt(:,prof)

      aux_k%so2wr_r(:,i) = & !aux_k%so2wr_r(:,i) + &
                       aux%so2r(:,prof) * aux_k%so2rwr_r(:,i)
      aux_k%so2w_r(:,i) = & !aux_k%so2w_r(:,i) + &
                      aux%so2r(:,prof) * aux_k%so2rw_r(:,i)

      aux_k%so2r(:,i) = aux_k%so2r(:,i) + & 
                                  aux%so2wr_r(:, prof) * aux_k%so2rwr_r(:,i) + &
                                  aux%so2w_r(:, prof) * aux_k%so2rw_r(:,i)
    ENDDO

    CALL reciprocal_k(aux%so2w_r, aux_k%so2w, aux_k%so2w_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    CALL reciprocal_k(aux%so2wr_r, aux_k%so2wr, aux_k%so2wr_r, acc = .FALSE._jplm, map = map, prof_stat = prof_stat)

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      sum1_k = 0._jprb
      sum2_k = 0._jprb

      DO lay = nlayers, 1, -1
        sum2_k = sum2_k + coef%so2tstar_wsum_r(lay) * aux_k%so2wr(lay, i)

        aux_k%so2_layer(lay, i) = & !aux_k%so2_layer(lay, i) + &
          coef%dpp(lay - 1) * aux%t_layer(lay, prof) * sum2_k

        aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
          coef%dpp(lay - 1) * aux%so2_layer(lay, prof) * sum2_k
       
        IF (opts%rt_all%so2_data) THEN
          sum1_k = sum1_k + coef%so2star_wsum_r(lay) * aux_k%so2w(lay, i)
          aux_k%so2_layer(lay, i) = aux_k%so2_layer(lay, i) + coef%dpp(lay - 1) * sum1_k
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  IF (coef%nch4 > 0) THEN
    IF (opts%rt_all%ch4_data) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof

        aux_k%ch4w(:,i)= aux_k%ch4w(:,i) + 0.25_jprb * (aux%ch4w_4rt(:, prof) * aux%ch4w_r(:, prof)) * aux_k%ch4w_4rt(:,i)
        aux_k%ch4r(:,i) = aux_k%ch4r(:,i) + aux%ch4w_r(:, prof) * aux_k%ch4rw_r(:,i)
        aux_k%ch4w_r(:,i) = &! aux_k%ch4w_r(:,i) + &
                        aux%ch4r(:, prof) * aux_k%ch4rw_r(:,i)
      ENDDO
      CALL reciprocal_k(aux%ch4w_r, aux_k%ch4w, aux_k%ch4w_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    ENDIF

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      sum2_k = 0._jprb
      IF (opts%rt_all%ch4_data) THEN
        sum1_k = 0._jprb
        DO lay = nlayers, 1, -1
          sum1_k = sum1_k + coef%ch4star_wsum_r(lay) * aux_k%ch4w(lay, i)
          sum2_k = sum2_k + coef%ch4tstar_wsum_r(lay) * aux_k%ch4wr(lay, i)

          aux_k%ch4_layer(lay, i) = &!aux_k%ch4_layer(lay, i) + 
            coef%dpp(lay - 1) * ( &
                                  sum1_k + &
                                  aux%t_layer(lay, prof) * sum2_k)

          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * aux%ch4_layer(lay, prof) * sum2_k
        ENDDO
      ELSE
        DO lay = nlayers, 1, -1
          sum2_k = sum2_k + coef%ch4tstar_wsum_r(lay) * aux_k%ch4wr(lay, i)
          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * coef%ch4star(lay) * sum2_k 
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (coef%nco > 0) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      aux_k%cowr(:,i)= aux_k%cowr(:,i) + 0.25_jprb * (aux%cowr_4rt(:, prof) * aux%cowr_r(:, prof)) * aux_k%cowr_4rt(:,i)
    ENDDO

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      sum2_k = 0._jprb
      IF (opts%rt_all%co_data) THEN
        sum1_k = 0._jprb

        DO lay = nlayers, 1, -1
          sum1_k = sum1_k + coef%costar_wsum_r(lay) * aux_k%cow(lay, i)
          sum2_k = sum2_k + coef%cotstar_wsum_r(lay) * aux_k%cowr(lay, i)

          aux_k%co_layer(lay, i) = &!aux_k%co_layer(lay, i) + &!
            coef%dpp(lay - 1) * ( &
            sum1_k + &
            aux%t_layer(lay, prof) * sum2_k)

          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * aux%co_layer(lay, prof) * sum2_k
        ENDDO
      ELSE
        DO lay = nlayers, 1, -1
          sum2_k = sum2_k + coef%cotstar_wsum_r(lay) * aux_k%cowr(lay, i)
          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * coef%costar(lay) * sum2_k 
        ENDDO
      ENDIF
    ENDDO

    aux_k%dt = aux_k%dt + 2._jprb * ABS(aux%dt(:,map(:,1))) * aux_k%dtabsdt
  ENDIF

  IF (coef%nn2o > 0) THEN
    IF (opts%rt_all%n2o_data) &
      CALL reciprocal_k(aux%n2ow_r, aux_k%n2ow, aux_k%n2ow_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      sum2_k = 0._jprb
      IF (opts%rt_all%n2o_data) THEN
        sum1_k = 0._jprb
        
        DO lay = nlayers, 1, -1
          sum1_k = sum1_k + coef%n2ostar_wsum_r(lay) * aux_k%n2ow(lay, i)
          sum2_k = sum2_k + coef%n2otstar_wsum_r(lay) * aux_k%n2owr(lay, i)
          
          aux_k%n2o_layer(lay, i) = &!aux_k%n2o_layer(lay, i) + &
            coef%dpp(lay - 1) * ( &
            sum1_k + &
            aux%t_layer(lay, prof) * sum2_k)

          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * aux%n2o_layer(lay, prof) * sum2_k
        ENDDO
      ELSE
        DO lay = nlayers, 1, -1
          sum2_k = sum2_k + coef%n2otstar_wsum_r(lay) * aux_k%n2owr(lay, i)
          aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
            coef%dpp(lay - 1) * coef%n2ostar(lay) * sum2_k 
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (coef%nco2 > 0) THEN
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (opts%rt_all%co2_data) THEN
        sum1_k = 0._jprb
        
        DO lay = nlayers, 1, -1
          sum1_k = sum1_k + coef%co2star_wsum_r(lay) * aux_k%co2w(lay, i)
          
          aux_k%co2_layer(lay, i) = &!aux_k%co2_layer(lay, i) + &
            coef%dpp(lay - 1) * sum1_k
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  DO i = 1, nchannels
    prof = chanprof(i)%prof
    sum1_k = 0._jprb
    DO lay = nlayers, 2, -1
      sum1_k = sum1_k + coef%tstar_wsum_r(lay) * aux_k%twr(lay, i)
      aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + coef%dpp(lay - 1) * sum1_k
    ENDDO

    aux_k%t_layer(1, i) = aux_k%t_layer(1, i) + &
                          coef%tstar_r(1) * aux_k%twr(1, i) + &
                          coef%dpp(0) * sum1_k          
  ENDDO

  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        sum1_k = 0._jprb
        DO lay = nlayers, 1, -1
          sum1_k = sum1_k + coef%ostar_wsum_r(lay) * aux_k%ow(lay, i)
          aux_k%o3_layer(lay, i) = &!aux_k%o3_layer(lay, i) + &
                                        coef%dpp(lay - 1) * sum1_k
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO i = 1, nchannels
    prof = chanprof(i)%prof
    aux_k%wr(:,i)= aux_k%wr(:,i) + aux%wwr_r(:, prof)* aux_k%wrwr_r(:,i)
    aux_k%wwr_r(:,i)= &!aux_k%wwr_r(:,i) + &
                   aux%wr(:, prof)* aux_k%wrwr_r(:,i)
  ENDDO
  CALL reciprocal_k(aux%wwr_r, aux_k%wwr, aux_k%wwr_r, acc = .FALSE._jplm, map = map, prof_stat = prof_stat)

  DO i = 1, nchannels
    prof = chanprof(i)%prof

    sum1_k = 0._jprb
    DO lay = nlayers, 1, -1
      sum1_k = sum1_k + coef%wtstar_wsum_r(lay) * aux_k%wwr(lay, i)

      aux_k%w_layer(lay, i) = &!aux_k%w_layer(lay, i) + &
        coef%dpp(lay - 1) * aux%t_layer(lay, prof) * sum1_k

      aux_k%t_layer(lay, i) = aux_k%t_layer(lay, i) + &
        coef%dpp(lay - 1) * aux%w_layer(lay, prof) * sum1_k
    ENDDO
  ENDDO

  DO i = 1, nchannels
    prof = chanprof(i)%prof
    aux_k%wr(:,i) =  aux_k%wr(:,i) + aux%ww_r(:, prof) * aux_k%wrw_r(:,i)
    aux_k%ww_r(:,i) = &!aux_k%ww_r(:,i) + &
                  aux%wr(:, prof) * aux_k%wrw_r(:,i)
  ENDDO
  CALL reciprocal_k(aux%ww_r, aux_k%ww, aux_k%ww_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  DO i = 1, nchannels
    prof = chanprof(i)%prof
    sum1_k = 0._jprb
    DO lay = nlayers, 1, -1
      sum1_k = sum1_k + coef%wstar_wsum_r(lay) * aux_k%ww(lay, i)
      aux_k%w_layer(lay, i) = aux_k%w_layer(lay, i) + &
                                   coef%dpp(lay - 1) * sum1_k
    ENDDO
  ENDDO

  ! -------------------------------------------------------------------------
  ! Calculate (profile / reference profile) ratios
  ! -------------------------------------------------------------------------

  DO i = 1, nchannels
    prof = chanprof(i)%prof

    aux_k%tr(:, i) = aux_k%tr(:, i) + 2._jprb * aux%tr(:, prof) * aux_k%tr2(:, i)

    aux_k%t_layer(:, i) = aux_k%t_layer(:, i) + &
                               aux_k%tr(:, i) * coef%tstar_r(:)

    aux_k%w_layer(:, i) = aux_k%w_layer(:, i) + aux_k%wr(:, i) * coef%wstar_r(:)

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux_k%o3_layer(:, i) = aux_k%o3_layer(:, i) + aux_k%or(:, i) * coef%ostar_r(:)
    ENDIF
  ENDDO

  IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
    CALL sqrt_k(aux%co2r_sqrt, aux_k%co2r, aux_k%co2r_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels; aux_k%co2_layer(:, i) = aux_k%co2_layer(:, i) + aux_k%co2r(:, i) * coef%co2star_r(:); ENDDO
  ENDIF
    
  IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
    CALL sqrt_k(aux%n2or_4rt, aux_k%n2or_sqrt, aux_k%n2or_4rt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    CALL sqrt_k(aux%n2or_sqrt, aux_k%n2or, aux_k%n2or_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels; aux_k%n2o_layer(:, i) = aux_k%n2o_layer(:, i) + aux_k%n2or(:, i) * coef%n2ostar_r(:); ENDDO
  ENDIF
    
  IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
    CALL sqrt_k(aux%cor_4rt, aux_k%cor_sqrt, aux_k%cor_4rt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    CALL sqrt_k(aux%cor_sqrt, aux_k%cor, aux_k%cor_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels; aux_k%co_layer(:, i) = aux_k%co_layer(:, i) + aux_k%cor(:, i) * coef%costar_r(:); ENDDO
  ENDIF

  IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
    CALL sqrt_k(aux%ch4r_4rt, aux_k%ch4r_sqrt, aux_k%ch4r_4rt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    CALL sqrt_k(aux%ch4r_sqrt, aux_k%ch4r, aux_k%ch4r_sqrt, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)
    DO i = 1, nchannels; aux_k%ch4_layer(:, i) = aux_k%ch4_layer(:, i) + aux_k%ch4r(:, i) * coef%ch4star_r(:); ENDDO
  ENDIF

  IF (coef%nso2 > 0) THEN
    DO i = 1, nchannels; aux_k%so2_layer(:, i) = aux_k%so2_layer(:, i) + aux_k%so2r(:, i) * coef%so2star_r(:); ENDDO
  ENDIF

  ! -------------------------------------------------------------------------
  ! Calculate deviations from reference profile for each layer
  ! -------------------------------------------------------------------------

  aux_k%t_layer(:,:) = aux_k%t_layer(:,:) + aux_k%dt(:,:)

  ! -------------------------------------------------------------------------
  ! Profile layer quantities: layer n lies between levels n and n+1
  ! -------------------------------------------------------------------------

  DO i = 1, nchannels
    prof = chanprof(i)%prof

    IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
      profiles_k(i)%so2(1) = profiles_k(i)%so2(1) + &
        0.5_jprb * aux_k%so2_layer(1, i)
      profiles_k(i)%so2(2:nlevels-1) = profiles_k(i)%so2(2:nlevels-1) + &
        0.5_jprb * (aux_k%so2_layer(1:nlevels-2, i) + aux_k%so2_layer(2:nlevels-1, i))
      profiles_k(i)%so2(nlevels) = profiles_k(i)%so2(nlevels) + &
        0.5_jprb * aux_k%so2_layer(nlevels-1, i)
    ENDIF

    IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
      IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
        profiles_k(i)%ch4(1) = profiles_k(i)%ch4(1) + &
          0.5_jprb * aux_k%ch4_layer(1, i)
        profiles_k(i)%ch4(2:nlevels-1) = profiles_k(i)%ch4(2:nlevels-1) + &
          0.5_jprb * (aux_k%ch4_layer(1:nlevels-2, i) + aux_k%ch4_layer(2:nlevels-1, i))
        profiles_k(i)%ch4(nlevels) = profiles_k(i)%ch4(nlevels) + &
          0.5_jprb * aux_k%ch4_layer(nlevels-1, i)
      ENDIF

      IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
        profiles_k(i)%co(1) = profiles_k(i)%co(1) + &
          0.5_jprb * aux_k%co_layer(1, i)
        profiles_k(i)%co(2:nlevels-1) = profiles_k(i)%co(2:nlevels-1) + &
          0.5_jprb * (aux_k%co_layer(1:nlevels-2, i) + aux_k%co_layer(2:nlevels-1, i))
        profiles_k(i)%co(nlevels) = profiles_k(i)%co(nlevels) + &
          0.5_jprb * aux_k%co_layer(nlevels-1, i)
      ENDIF
      
      IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
        profiles_k(i)%n2o(1) = profiles_k(i)%n2o(1) + &
          0.5_jprb * aux_k%n2o_layer(1, i)
        profiles_k(i)%n2o(2:nlevels-1) = profiles_k(i)%n2o(2:nlevels-1) + &
          0.5_jprb * (aux_k%n2o_layer(1:nlevels-2, i) + aux_k%n2o_layer(2:nlevels-1, i))
        profiles_k(i)%n2o(nlevels) = profiles_k(i)%n2o(nlevels) + &
          0.5_jprb * aux_k%n2o_layer(nlevels-1, i)
      ENDIF

      IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
        profiles_k(i)%co2(1) = profiles_k(i)%co2(1) + &
          0.5_jprb * aux_k%co2_layer(1, i)
        profiles_k(i)%co2(2:nlevels-1) = profiles_k(i)%co2(2:nlevels-1) + &
          0.5_jprb * (aux_k%co2_layer(1:nlevels-2, i) + aux_k%co2_layer(2:nlevels-1, i))
        profiles_k(i)%co2(nlevels) = profiles_k(i)%co2(nlevels) + &
          0.5_jprb * aux_k%co2_layer(nlevels-1, i)
      ENDIF
    ENDIF

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      profiles_k(i)%o3(1) = profiles_k(i)%o3(1) + &
        0.5_jprb * aux_k%o3_layer(1, i)
      profiles_k(i)%o3(2:nlevels-1) = profiles_k(i)%o3(2:nlevels-1) + &
        0.5_jprb * (aux_k%o3_layer(1:nlevels-2, i) + aux_k%o3_layer(2:nlevels-1, i))
      profiles_k(i)%o3(nlevels) = profiles_k(i)%o3(nlevels) + &
        0.5_jprb * aux_k%o3_layer(nlevels-1, i)
    ENDIF

    IF (opts%rt_all%use_q2m) THEN
      iv2lev = aux%s(prof)%nearestlev_surf    ! nearest level below surface

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        profiles_k(i)%s2m%q =  profiles_k(i)%s2m%q + &
          aux_k%w_layer(iv2lay, i) * 0.5_jprb

        ! This line subtracts the quantity that is added to profiles_k(i)%q(iv2lev) below
        ! so that it cancels out: in the case where use_q2m is true prof(i)%q(iv2lev) is
        ! replaced by prof(i)%s2m%q in the calculation of aux%w_layer(iv2lay, i)

        profiles_k(i)%q(iv2lev) = profiles_k(i)%q(iv2lev) - &
          aux_k%w_layer(iv2lay, i) * 0.5_jprb
      ENDIF
    ENDIF
    
    profiles_k(i)%q(1) = profiles_k(i)%q(1) + &
      0.5_jprb * aux_k%w_layer(1, i)
    profiles_k(i)%q(2:nlevels-1) = profiles_k(i)%q(2:nlevels-1) + &
      0.5_jprb * (aux_k%w_layer(1:nlevels-2, i) + aux_k%w_layer(2:nlevels-1, i))
    profiles_k(i)%q(nlevels) = profiles_k(i)%q(nlevels) + &
      0.5_jprb * aux_k%w_layer(nlevels-1, i)

    IF (opts%rt_all%use_t2m_opdep) THEN
      iv2lev = aux%s(prof)%nearestlev_surf    ! nearest level below surface

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        profiles_k(i)%s2m%t =  profiles_k(i)%s2m%t + &
          aux_k%t_layer(iv2lay, i) * 0.5_jprb

        ! This line subtracts the quantity that is added to profiles_k(i)%t(iv2lev) below
        ! so that it cancels out: in the case where use_t2m_opdep is true prof(i)%t(iv2lev) is
        ! replaced by prof(i)%s2m%t in the calculation of aux%t_layer(iv2lay, i)

        profiles_k(i)%t(iv2lev) = profiles_k(i)%t(iv2lev) - &
          aux_k%t_layer(iv2lay, i) * 0.5_jprb
      ENDIF
    ENDIF

    profiles_k(i)%t(1) = profiles_k(i)%t(1) + &
      0.5_jprb * aux_k%t_layer(1, i)
    profiles_k(i)%t(2:nlevels-1) = profiles_k(i)%t(2:nlevels-1) + &
      0.5_jprb * (aux_k%t_layer(1:nlevels-2, i) + aux_k%t_layer(2:nlevels-1, i))
    profiles_k(i)%t(nlevels) = profiles_k(i)%t(nlevels) + &
      0.5_jprb * aux_k%t_layer(nlevels-1, i)
  ENDDO

  ! -------------------------------------------------------------------------
  ! Raytracing path data required for predictor calculations
  ! -------------------------------------------------------------------------

  IF (plane_parallel) THEN
    aux_k%pathsat_rsqrt = 0._jprb
    aux_k%pathsat_sqrt  = 0._jprb
  ELSE
    DO i = 1, nchannels
      prof = chanprof(i)%prof
      raytracing_k%pathsat(:,i) = raytracing_k%pathsat(:,i) + &
        aux_k%pathsat_sqrt(:,i) * aux%pathsat_rsqrt(:,prof)

      aux_k%pathsat_rsqrt(:,i) = &!aux_k%pathsat_rsqrt + &
        raytracing%pathsat(:,prof) * aux_k%pathsat_sqrt(:,i)
    ENDDO
    CALL invsqrt_k(aux%pathsat_4rt, aux_k%pathsat_rsqrt, aux_k%pathsat_4rt, &
      acc=.TRUE._jplm, map = map, prof_stat = prof_stat)
    CALL invsqrt_k(aux%pathsat_rsqrt, raytracing_k%pathsat, aux_k%pathsat_rsqrt, &
      acc=.TRUE._jplm, map = map, prof_stat = prof_stat)

    IF (dosolar) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        raytracing_k%patheff(:,i) = raytracing_k%patheff(:,i) + &
          aux_k%patheff_sqrt(:,i) * aux%patheff_rsqrt(:,prof)

        aux_k%patheff_rsqrt(:,i) = &!aux_k%patheff_rsqrt + &
          raytracing%patheff(:,prof) * aux_k%patheff_sqrt(:,i)
      ENDDO
      CALL invsqrt_k(aux%patheff_4rt, aux_k%patheff_rsqrt, aux_k%patheff_4rt, &
        acc=.TRUE._jplm, map = map, prof_stat = prof_stat)
      CALL invsqrt_k(aux%patheff_rsqrt, raytracing_k%patheff, aux_k%patheff_rsqrt, &
        acc=.TRUE._jplm, map = map, prof_stat = prof_stat)
    ENDIF
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_13_K', 1_jpim, ZHOOK_HANDLE)

END SUBROUTINE rttov_predictor_precalc_13_k
