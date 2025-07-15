! Description:
!> @file
!!   Calculate some quantities used in v7, v8 and v9 predictor calculations.
!
!> @brief
!!   Calculate some quantities used in v7, v8 and v9 predictor calculations.
!!
!! @details
!!   All calculations here are done only on coefficient levels.
!!
!!   Most of the quantities computed relate to the temperature and gas
!!   profiles. However, a few local path quantities in the raytracing structure
!!   are calculated here because they are required only for the predictor
!!   calculations. The CO2 thickness calculation for PMC shift sensors is
!!   also performed here as this is only computed on coefficient levels.
!!
!!   The angles argument is used only for PMC shift sensors. The angles and
!!   raytracing arguments are mandatory in general, but are optional in order
!!   to simplify the coefficient generation code.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     dosolar              flag indicating whether solar computations are being performed
!! @param[in]     prof                 input profiles
!! @param[in]     coef                 optical depth coefficient structure
!! @param[in]     coef_pccomp          PC coefficients structure
!! @param[in,out] aux                  coef level auxiliary profile data structure
!! @param[in]     angles               geometry structure, optional only for LBL code
!! @param[in,out] raytracing           raytracing structure, optional only for LBL code
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
SUBROUTINE rttov_predictor_precalc_789(  &
              opts,        &
              dosolar,     &
              prof,        &
              coef,        &
              coef_pccomp, &
              aux,         &
              angles,      &
              raytracing)

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
  USE rttov_const, ONLY : gas_id_so2
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY: invsqrt, reciprocal
!INTF_ON
  USE parkind1, ONLY : jplm

  IMPLICIT NONE
  TYPE(rttov_options),          INTENT(IN)              :: opts
  LOGICAL(jplm),                INTENT(IN)              :: dosolar
  TYPE(rttov_profile),          INTENT(IN)              :: prof(:)
  TYPE(rttov_coef),             INTENT(IN)              :: coef
  TYPE(rttov_coef_pccomp),      INTENT(IN)              :: coef_pccomp
  TYPE(rttov_profile_aux_coef), INTENT(INOUT)           :: aux
  TYPE(rttov_geometry),         INTENT(IN),    OPTIONAL :: angles(:)  ! Only optional for LBL code
  TYPE(rttov_raytracing),       INTENT(INOUT), OPTIONAL :: raytracing ! Only optional for LBL code
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers, nlevels
  INTEGER(jpim) :: iprof, lay
  REAL(jprb) :: ZHOOK_HANDLE

  REAL(jprb),POINTER :: gas_profile(:)
  REAL(jprb) :: sum1, sum2
  INTEGER(jpim) :: iv2lay, iv2lev, iv3lev

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_789', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(prof)
  nlayers = prof(1)%nlayers
  nlevels = prof(1)%nlevels

!-----------------------------------------------------------------------------------------
! Calculate predictor data
!-----------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ! Raytracing path data required for predictor calculations
  ! -------------------------------------------------------------------------
  IF (PRESENT(raytracing)) THEN
    CALL invsqrt(raytracing%pathsat, aux%pathsat_rsqrt)
    CALL invsqrt(aux%pathsat_rsqrt, aux%pathsat_4rt)
    aux%pathsat_sqrt = raytracing%pathsat * aux%pathsat_rsqrt

    IF (dosolar) THEN
      CALL invsqrt(raytracing%patheff, aux%patheff_rsqrt)
      CALL invsqrt(aux%patheff_rsqrt, aux%patheff_4rt)
      aux%patheff_sqrt = raytracing%patheff * aux%patheff_rsqrt
    ENDIF
  ENDIF

  ! -------------------------------------------------------------------------
  ! Profile layer quantities: layer n lies between levels n and n+1
  ! -------------------------------------------------------------------------

  DO iprof = 1, nprofiles
    aux%t_layer(1:nlayers, iprof) = &
      (prof(iprof)%t(1:nlevels-1) + prof(iprof)%t(2:nlevels)) * 0.5_jprb

    IF (opts%rt_all%use_t2m_opdep) THEN
      iv3lev = aux%s(iprof)%nearestlev_surf - 1 ! nearest level above surface
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        aux%t_layer(iv2lay, iprof) = &
          (prof(iprof)%s2m%t + prof(iprof)%t(iv3lev)) * 0.5_jprb
      ENDIF
    ENDIF

    aux%w_layer(1:nlayers, iprof) = &
      (prof(iprof)%q(1:nlevels-1) + prof(iprof)%q(2:nlevels)) * 0.5_jprb

    IF (opts%rt_all%use_q2m) THEN
      iv3lev = aux%s(iprof)%nearestlev_surf - 1 ! nearest level above surface
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        aux%w_layer(iv2lay, iprof) = &
          (prof(iprof)%s2m%q + prof(iprof)%q(iv3lev)) * 0.5_jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux%o3_layer(1:nlayers, iprof) = &
        (prof(iprof)%o3(1:nlevels-1) + prof(iprof)%o3(2:nlevels)) * 0.5_jprb
    ENDIF
    
    IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux%co2_layer(1:nlayers, iprof) = &
          (prof(iprof)%co2(1:nlevels-1) + prof(iprof)%co2(2:nlevels)) * 0.5_jprb
      ELSE
        aux%co2_layer(1:nlayers, iprof) = &
          (coef_pccomp%co2_pc_ref(1:nlevels-1) + coef_pccomp%co2_pc_ref(2:nlevels)) * 0.5_jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux%n2o_layer(1:nlayers, iprof) = &
          (prof(iprof)%n2o(1:nlevels-1) + prof(iprof)%n2o(2:nlevels)) * 0.5_jprb
      ELSE
        aux%n2o_layer(1:nlayers, iprof) = &
          (coef_pccomp%n2o_pc_ref(1:nlevels-1) + coef_pccomp%n2o_pc_ref(2:nlevels)) * 0.5_jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux%co_layer(1:nlayers, iprof) = &
          (prof(iprof)%co(1:nlevels-1) + prof(iprof)%co(2:nlevels)) * 0.5_jprb
      ELSE
        aux%co_layer(1:nlayers, iprof) = &
          (coef_pccomp%co_pc_ref(1:nlevels-1) + coef_pccomp%co_pc_ref(2:nlevels)) * 0.5_jprb
      ENDIF
    ENDIF

    IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
      IF (.NOT. (opts%rt_ir%pc%addpc .AND. coef_pccomp%fmv_pc_comp_pc < 5)) THEN
        aux%ch4_layer(1:nlayers, iprof) = &
          (prof(iprof)%ch4(1:nlevels-1) + prof(iprof)%ch4(2:nlevels)) * 0.5_jprb
      ELSE
        aux%ch4_layer(1:nlayers, iprof) = &
          (coef_pccomp%ch4_pc_ref(1:nlevels-1) + coef_pccomp%ch4_pc_ref(2:nlevels)) * 0.5_jprb
      ENDIF
    ENDIF
    
    IF (coef%nso2 > 0) THEN
      IF (opts%rt_all%so2_data) THEN
        aux%so2_layer(1:nlayers, iprof) = &
          (prof(iprof)%so2(1:nlevels-1) + prof(iprof)%so2(2:nlevels)) * 0.5_jprb
      ELSE
        aux%so2_layer(1:nlayers, iprof) = &
          (coef%bkg_prfl_mr(1:nlevels-1,coef%fmv_gas_pos(gas_id_so2)) + &
           coef%bkg_prfl_mr(2:nlevels,  coef%fmv_gas_pos(gas_id_so2))) * 0.5_jprb
      ENDIF
    ENDIF

    ! -------------------------------------------------------------------------
    ! Calculate deviations from reference profile for each layer
    ! -------------------------------------------------------------------------

    aux%dt(:, iprof) = aux%t_layer(:, iprof) - coef%tstar(:)
    aux%dtabsdt(:, iprof) = aux%dt(:, iprof) * ABS(aux%dt(:, iprof))

    ! If no input O3 profile we still use the input temperature profile for dto
    IF (coef%nozone > 0) THEN
      aux%dto(:, iprof) = aux%t_layer(:, iprof) - coef%to3star(:)
      IF (coef%fmv_model_ver == 9) THEN
        aux%tro(:, iprof) = aux%t_layer(:, iprof) * coef%to3star_r(:)
      ENDIF
    ENDIF

    ! -------------------------------------------------------------------------
    ! Calculate (profile / reference profile) ratios
    ! -------------------------------------------------------------------------

    aux%tr(:, iprof) = aux%t_layer(:, iprof) * coef%tstar_r(:)
    aux%tr2(:, iprof) = aux%tr(:, iprof) ** 2_jpim
    aux%wr(:, iprof) = aux%w_layer(:, iprof) * coef%wstar_r(:)
    
    ! If an optional gas profile is absent then the predictor calculations
    ! in rttov_setpredictors_789 are simplified because the ratio quantities
    ! have simple fixed values (e.g. aux%or = 1.)
    ! This does not apply to SO2 because the background SO2 profile is not
    ! the same as the reference profile.

    IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
      aux%or(:, iprof) = aux%o3_layer(:, iprof) * coef%ostar_r(:)
    ENDIF

    IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
      aux%co2r(:, iprof) = aux%co2_layer(:, iprof) * coef%co2star_r(:)
      IF (coef%fmv_model_ver == 9) aux%co2r_sqrt(:, iprof) = SQRT(aux%co2r(:, iprof))
    ENDIF

    IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
      aux%n2or(:, iprof) = aux%n2o_layer(:, iprof) * coef%n2ostar_r(:)
      aux%n2or_sqrt(:, iprof) = SQRT(aux%n2or(:, iprof))
      aux%n2or_4rt(:, iprof) = SQRT(aux%n2or_sqrt(:, iprof))
    ENDIF

    IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
      aux%cor(:, iprof) = aux%co_layer(:, iprof) * coef%costar_r(:)
      aux%cor_sqrt(:, iprof) = SQRT(aux%cor(:, iprof))
      aux%cor_4rt(:, iprof) = SQRT(aux%cor_sqrt(:, iprof))
    ENDIF

    IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
      aux%ch4r(:, iprof) = aux%ch4_layer(:, iprof) * coef%ch4star_r(:)
      aux%ch4r_sqrt(:, iprof) = SQRT(aux%ch4r(:, iprof))
      aux%ch4r_4rt(:, iprof) = SQRT(aux%ch4r_sqrt(:, iprof))
    ENDIF

    ! NB Always calculate this regardless of so2_data
    IF (coef%nso2 > 0) THEN
      aux%so2r(:, iprof) = aux%so2_layer(:, iprof) * coef%so2star_r(:)
    ENDIF
    
  ENDDO

  ! -------------------------------------------------------------------------
  ! Calculate profile / reference profile sums
  ! -------------------------------------------------------------------------

  ! Quantities such as ow and co2w are not calculated if the corresponding
  ! gas profile is not input as the predictor calculations in
  ! rttov_setpredictors_789 are simplified in this case.

  IF (coef%fmv_model_ver <= 8) THEN  
    DO iprof = 1, nprofiles
      aux%tw(1, iprof) = 0._jprb

      DO lay = 2, nlayers
        ! cumulate overlying layers: weighting tr relates to same layer as dpp
        ! do not need dpp(0) to start
        aux%tw(lay, iprof) = aux%tw(lay - 1, iprof) + &
          coef%dpp(lay - 1) * aux%tr(lay - 1, iprof)
      ENDDO
    ENDDO
  ENDIF

  DO iprof = 1, nprofiles
    ! cumulating column overlying layer and layer itself
    sum1 = 0._jprb
    DO lay = 1, nlayers
      ! cumulate overlying layers: weighting w or wstar relates to layer below dpp
      ! need dpp(0) to start
      sum1 = sum1 + coef%dpp(lay - 1) * aux%w_layer(lay, iprof)
      aux%ww(lay, iprof) = sum1 * coef%wstar_wsum_r(lay)
    ENDDO
  ENDDO

  IF (.NOT. coef%fmv_model_ver == 8) THEN
    CALL reciprocal(aux%ww, aux%ww_r)
    aux%wrw_r = aux%wr * aux%ww_r
  ENDIF

  IF (coef%fmv_model_ver >= 8) THEN
    DO iprof = 1, nprofiles
      sum1 = 0._jprb
      DO lay = 1, nlayers
        sum1 = sum1 + coef%dpp(lay - 1) * aux%w_layer(lay, iprof) * aux%t_layer(lay, iprof)
        aux%wwr(lay, iprof) = sum1 * coef%wtstar_wsum_r(lay)
      ENDDO
    ENDDO

    CALL reciprocal(aux%wwr, aux%wwr_r)
    aux%wrwr_r = aux%wr * aux%wwr_r
  ENDIF
 
  IF (coef%nozone > 0) THEN
    IF (opts%rt_all%ozone_data) THEN
      DO iprof = 1, nprofiles
        sum1 = 0._jprb
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          sum1 = sum1 + coef%dpp(lay - 1) * aux%o3_layer(lay, iprof) ! OK, dpp(0) defined 
          aux%ow(lay, iprof) = sum1 * coef%ostar_wsum_r(lay)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF (coef%fmv_model_ver >= 8) THEN
    ! twr calculation still uses the input temperature profile
    DO iprof = 1, nprofiles
      sum2   = 0._jprb
      IF (coef%fmv_model_ver == 8) THEN
        aux%twr(1, iprof) = 1.0_jprb
      ELSEIF (coef%fmv_model_ver == 9) THEN
        aux%twr(1, iprof) = aux%t_layer(1, iprof) * coef%tstar_r(1)
        sum2 = coef%dpp(0) * aux%t_layer(1, iprof)
      ENDIF

      DO lay = 2, nlayers
        ! cumulate overlying layers (t, tstar relate to same layer as dpp)
        ! do not need dpp(0) to start
        IF (coef%fmv_model_ver == 8) THEN
          sum2       = sum2 + coef%dpp(lay - 1) * aux%t_layer(lay-1, iprof)
          aux%twr(lay, iprof) = sum2 * coef%tstarmod_wsum_r(lay)
        ELSE
          sum2       = sum2 + coef%dpp(lay - 1) * aux%t_layer(lay, iprof)
          aux%twr(lay, iprof) = sum2 * coef%tstar_wsum_r(lay)
        ENDIF
      ENDDO
    ENDDO

    IF (coef%fmv_model_ver == 9) THEN
      DO iprof = 1, nprofiles
        aux%tuw(1, iprof) = aux%t_layer(1, iprof) * coef%tstar_r(1)
        sum1 = aux%t_layer(1, iprof)

        DO lay = 2, nlayers
          ! cumulate overlying layers (t, tstar relate to same layer as dpp)
          ! do not need dpp(0) to start
          sum1       = sum1 + aux%t_layer(lay, iprof) 
          aux%tuw(lay, iprof) = sum1 * coef%tstar_uwsum_r(lay)
        ENDDO
      ENDDO
    ENDIF

    IF (coef%nco2 > 0) THEN
      IF (opts%rt_all%co2_data) THEN
        DO iprof = 1, nprofiles
          sum1 = 0._jprb
          DO lay = 1, nlayers
            ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
            ! need dpp(0) to start
            sum1 = sum1 + coef%dpp(lay - 1) * aux%co2_layer(lay, iprof) ! OK, dpp(0) defined 
            aux%co2w(lay, iprof) = sum1 * coef%co2star_wsum_r(lay)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    IF (coef%nn2o > 0) THEN
      DO iprof = 1, nprofiles
        IF (opts%rt_all%n2o_data) THEN
          gas_profile => aux%n2o_layer(:, iprof)
        ELSE
          gas_profile => coef%n2ostar
        ENDIF
        
        sum1 = 0._jprb
        sum2 = 0._jprb
        DO lay = 1, nlayers
          IF (opts%rt_all%n2o_data) THEN
            sum1 = sum1 + coef%dpp(lay - 1) * gas_profile(lay)
            aux%n2ow(lay, iprof) = sum1 * coef%n2ostar_wsum_r(lay)
          ENDIF
          sum2 = sum2 + coef%dpp(lay - 1) * gas_profile(lay) * aux%t_layer(lay, iprof)
          aux%n2owr(lay, iprof) = sum2 * coef%n2otstar_wsum_r(lay)
        ENDDO
      ENDDO
      IF (opts%rt_all%n2o_data) CALL reciprocal(aux%n2ow, aux%n2ow_r)
    ENDIF

    IF (coef%nco > 0) THEN
      DO iprof = 1, nprofiles
        IF (opts%rt_all%co_data) THEN
          gas_profile => aux%co_layer(:, iprof)
        ELSE
          gas_profile => coef%costar
        ENDIF
        
        sum1 = 0._jprb
        sum2 = 0._jprb
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          IF (opts%rt_all%co_data) THEN
            sum1 = sum1 + coef%dpp(lay - 1) * gas_profile(lay)
            aux%cow(lay, iprof) = sum1 * coef%costar_wsum_r(lay)
          ENDIF
          sum2 = sum2 + coef%dpp(lay - 1) * gas_profile(lay) * aux%t_layer(lay, iprof)
          aux%cowr(lay, iprof) = sum2 * coef%cotstar_wsum_r(lay)
        ENDDO
      ENDDO
      CALL reciprocal(aux%cowr, aux%cowr_r)
      aux%cowr_4rt = SQRT(SQRT(aux%cowr))
    ENDIF

    IF (coef%nch4 > 0) THEN
      DO iprof = 1, nprofiles
        IF (opts%rt_all%ch4_data) THEN
          gas_profile => aux%ch4_layer(:, iprof)
        ELSE
          gas_profile => coef%ch4star
        ENDIF
        
        sum1 = 0._jprb
        sum2 = 0._jprb
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          IF (opts%rt_all%ch4_data) THEN
            sum1 = sum1 + coef%dpp(lay - 1) * gas_profile(lay) ! OK, dpp(0) defined 
            aux%ch4w(lay, iprof) = sum1 * coef%ch4star_wsum_r(lay)
          ENDIF
          sum2 = sum2 + coef%dpp(lay - 1) * gas_profile(lay) * aux%t_layer(lay, iprof)
          aux%ch4wr(lay, iprof) = sum2 * coef%ch4tstar_wsum_r(lay)
        ENDDO
      ENDDO
      IF (opts%rt_all%ch4_data) THEN
        CALL reciprocal(aux%ch4w, aux%ch4w_r)
        aux%ch4rw_r = aux%ch4r * aux%ch4w_r
      ENDIF
    ENDIF
    
    IF (coef%nso2 > 0) THEN
      DO iprof = 1, nprofiles
        sum1 = 0._jprb
        sum2 = 0._jprb
        DO lay = 1, nlayers
          ! cumulate overlying layers: weighting o or ostar relates to layer below dpp
          ! need dpp(0) to start
          sum1 = sum1 + coef%dpp(lay - 1) * aux%so2_layer(lay,iprof) ! OK, dpp(0) defined 
          aux%so2w(lay,iprof) = sum1 * coef%so2star_wsum_r(lay)
          sum2 = sum2 + coef%dpp(lay - 1) * aux%so2_layer(lay,iprof) * aux%t_layer(lay,iprof)
          aux%so2wr(lay,iprof) = sum2 * coef%so2tstar_wsum_r(lay)
        ENDDO
      ENDDO
      
      CALL reciprocal(aux%so2w, aux%so2w_r)
      CALL reciprocal(aux%so2wr, aux%so2wr_r)
      aux%so2rwr_r = aux%so2r * aux%so2wr_r
      aux%so2rw_r = aux%so2r * aux%so2w_r
      aux%so2r_sqrt = SQRT(aux%so2r)
      aux%so2r_4rt = SQRT(aux%so2r_sqrt)
      aux%so2w_sqrt = SQRT(aux%so2w)
      aux%so2w_4rt = SQRT(aux%so2w_sqrt)
    ENDIF
  ENDIF

  ! -------------------------------------------------------------------------
  ! Calculate and store useful intermediate variables
  ! -------------------------------------------------------------------------

  CALL reciprocal(aux%tr, aux%tr_r)

  IF (coef%fmv_model_ver >= 8) THEN
    aux%tr_sqrt = SQRT(aux%tr)
  ENDIF

  IF (coef%fmv_model_ver <= 8) THEN
    aux%tw_sqrt = SQRT(aux%tw)
    aux%tw_sqrt(1,:) = 1e-100_jprb ! DAR used for division later...
    aux%tw_4rt = SQRT(aux%tw_sqrt)
  ENDIF

  CALL invsqrt(aux%wr, aux%wr_rsqrt)
  CALL invsqrt(aux%wr_rsqrt, aux%wr_4rt)
  aux%wr_sqrt = aux%wr * aux%wr_rsqrt

  IF (coef%fmv_model_ver >= 9) THEN
    CALL invsqrt(aux%ww, aux%ww_rsqrt)
    CALL invsqrt(aux%ww_rsqrt, aux%ww_4rt)
    aux%ww_sqrt = aux%ww * aux%ww_rsqrt
  ENDIF

  IF (coef%nozone > 0 .and. opts%rt_all%ozone_data) THEN
    CALL invsqrt(aux%ow, aux%ow_rsqrt)
    IF (coef%fmv_model_ver == 9) CALL invsqrt(aux%ow_rsqrt, aux%ow_4rt)
    aux%ow_r = aux%ow_rsqrt**2_jpim
    aux%or_sqrt = SQRT(aux%or)
    aux%ow_sqrt = aux%ow * aux%ow_rsqrt
  ENDIF

  IF (coef%nco > 0 .AND. opts%rt_all%co_data) THEN
    CALL invsqrt(aux%cow, aux%cow_rsqrt)
    aux%cow_r4rt = SQRT(aux%cow_rsqrt)

    aux%corw_rsqrt = aux%cor * aux%cow_rsqrt
    aux%corw_r = aux%cor * aux%cow_rsqrt**2_jpim
    aux%corw_r4rt = aux%cor * aux%cow_r4rt
  ENDIF

  ! Compute the STP co2 thickness in cm
  IF (coef%pmc_shift) CALL calc_ssu_co2_thickness()  
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_PREDICTOR_PRECALC_789', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_ssu_co2_thickness
    USE rttov_const, ONLY : p0, t0, flatt, eqrad, omega, grave
    REAL(jprb) :: rearth(SIZE(prof))
    REAL(jprb) :: gravh(SIZE(prof))
    REAL(jprb) :: gravl(SIZE(prof)) ! gravity at latitude lat [m/s^2]
    REAL(jprb) :: eta, beta ! coefficients of the international gravity formula
    REAL(jprb) :: dflat, fac
    REAL(jprb) :: pstd
    REAL(jprb) :: tstd
    REAL(jprb) :: pa(prof(1)%nlayers)
    REAL(jprb) :: dpa(prof(1)%nlayers)
    REAL(jprb) :: rho(prof(1)%nlayers)
    REAL(jprb) :: temp(prof(1)%nlayers)
    REAL(jprb) :: co2(prof(1)%nlayers)
    REAL(jprb) :: dz
    REAL(jprb) :: dzstp
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
        pa(lay)   = (prof(iprof)%p(lay) + prof(iprof)%p(lay+1)) * 100._jprb / 2._jprb ! pa
        dpa(lay)  = (prof(iprof)%p(lay) - prof(iprof)%p(lay+1)) * 100._jprb           ! pa
        rho(lay)  = (raytracing%dmair(lay,iprof) + raytracing%dmair(lay+1,iprof)) / 2._jprb   ! kg/m3
        temp(lay) = (prof(iprof)%t(lay) + prof(iprof)%t(lay+1)) / 2._jprb             ! K
        co2(lay)  = (prof(iprof)%co2(lay) + prof(iprof)%co2(lay+1)) / 2._jprb ! ppmv dry
      ENDDO

      ! Calculation of co2 thickness at stp
      aux%co2_cm(iprof) = 0._jprb

      DO lay = 1, nlayers
        dz = -dpa(lay) / ((gravl(iprof)-gravh(iprof) * &
          raytracing%hgpl(lay+1,iprof)) * rho(lay))    ! m
        ! pv=rt so v=rt/p so h=rt/p when v is height h of unit column
        dzstp = dz * pa(lay) / pstd * tstd / temp(lay) * co2(lay) * 1.E-6_jprb * 100._jprb   ! cm
        ! NB equal v contains equal numbers of molecules (e.g. co2) - result in cm

        aux%co2_cm(iprof) = aux%co2_cm(iprof) + dzstp
      ENDDO
    ENDDO
  END SUBROUTINE calc_ssu_co2_thickness

END SUBROUTINE rttov_predictor_precalc_789
