SUBROUTINE rttov_intavg_prof_k( &
              opts,        &
              chanprof,    &
              kni,         &
              kno,         &
              ProfIn,      &
              ProfIn_k,    &
              ProfGasIn,   &
              ProfGasIn_k, &
              ProfOut,     &
              ProfOut_k,   &
              coef,        &
              coef_pccomp)
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!-----------------------------------------------------------------------
! usage of arguments
!              set of profile elements
!              -----------------------
! nprofiles     number of profiles
! kni           number of source levs
! kno           number of target levs
! ProfIn        source for interpolator
! ProfIn_k      source for interpolator
! ProfOut       target for interpolator
! ProfOut_k     target for interpolator
!
!
!--------------------------------------------------------------------------
! application within RTTOV (for NWP SAF)
!   - profiles on USER levels -> prediction of gas/atm opdeps on COEF levels
!   - gas/atm opdeps on COEF levels -> transmittances on USER levels for RTE
!---------------------------------------------------------------------------
!
!     History:
!
!     Version   Date      Comment
!     -------   ----      -------
!
!     1         12/2006   main code intavg restructured Peter Rayer
!
!     2         10/2007   Optimisation Deborah Salmond and Pascal Brunel
!
!     3         11/2007   layeravg now deals better with output levels near
!                         and beyond the upper and lower boundaries of the
!                         input levels, so there is now no need for a special
!                         case here - A. Geer
!
!     4         01/2008   lgradp option for TL/AD for user pressure levels,
!                         Niels Bormann
!     r1428     09/2013   Made coef optional for use with HTFRTC
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
  USE rttov_types, ONLY : rttov_options, rttov_coef, rttov_coef_pccomp, rttov_chanprof, rttov_profile
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb, jplm
  USE rttov_const, ONLY : &
    realtol,              &
    interp_rochon,        &
    interp_loglinear,     &
    gas_id_watervapour,   &
    gas_id_ozone,         &
    gas_id_co2,           &
    gas_id_n2o,           &
    gas_id_co,            &
    gas_id_ch4,           &
    gas_id_so2
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim),      INTENT(IN)    :: kni, kno                    ! number of levels
  TYPE(rttov_profile),     INTENT(IN)    :: ProfIn(:)                   ! atmospheric profiles
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfIn_k(SIZE(chanprof))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfGasIn(SIZE(ProfIn))     ! atmospheric gas profiles (ppmv dry)
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfGasIn_k(SIZE(chanprof))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfOut(SIZE(ProfIn))       ! atmospheric profiles
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfOut_k(SIZE(chanprof))
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
#include "rttov_layeravg.interface"
#include "rttov_layeravg_k.interface"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, jotop, kn, istart(kno, SIZE(ProfIn)), iend(kno, SIZE(ProfIn))
  INTEGER(KIND=jpim) :: prof, ig
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: pvlev(kni, SIZE(ProfIn))
  REAL   (KIND=jprb) :: pvlev_k(kni)
  REAL   (KIND=jprb) :: ppo(kno, SIZE(ProfIn))
  REAL   (KIND=jprb) :: zlnpi(kni, SIZE(ProfIn)), zlnpi_k(kni), zpz(kni, kno, SIZE(ProfIn)), zpz_k(kni, kno)
  REAL   (KIND=jprb) :: zlnpo(kno, SIZE(ProfIn)), zlnpo_k(kno)
  INTEGER(KIND=jpim) :: interp_mode
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_K', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = SIZE(ProfIn)
  nchannels         = SIZE(chanprof)
  IF (opts%interpolation%lgradp) THEN
    llevels_different = .TRUE.
    zpz_k             = 0._jprb
  ELSE
    DO prof = 1, nprofiles - 1
      DO jo = 1, kni
        IF (ABS(ProfIn(prof)%p(jo) - ProfIn(prof + 1)%p(jo)) > realtol) THEN
          llevels_different = .TRUE.
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! Set interpolation mode: default to Rochon for invalid mode
  ! Only use loglinear for profile interpolation if user asked for loglinear-only

  interp_mode = interp_rochon
  IF (opts%interpolation%interp_mode == interp_loglinear) THEN
    interp_mode = interp_loglinear
  ENDIF

  DO prof = 1, nprofiles
    IF (llevels_different .OR. prof == 1) THEN
! source levels
      pvlev(1:kni, prof) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni, prof) = LOG(pvlev(1:kni, prof))
! target levels
      ppo(1:kno, prof)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno, prof) = LOG(ppo(1:kno, prof))
!--------------------------------------------------------------------------
!1. profile interpolation user -> coef
!--------------------------------------------------------------------------
! replaced by profile assignments
!1.2 - interpolation procedure for profile elements (user -> coef)
!------------------------------------------------------------------
! for rttov same set of weights for all elements, so one call suffices
      CALL rttov_layeravg( &
              zlnpo(:, prof),  &
              zlnpi(:, prof),  &
              kno,             &
              kni,             &
              zpz(:, :, prof), &
              istart(:, prof), &
              iend(:, prof),   &
              interp_mode)
    ENDIF
  ENDDO

  DO kn = 1, nchannels
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF

    IF (opts%interpolation%reg_limit_extrap .AND. ppo(2, prof) < pvlev(1, prof)) THEN
      DO jotop = 2, kno-1
        IF (ppo(jotop+1, prof) > pvlev(1, prof)) EXIT
      ENDDO

      IF (.NOT. opts%rt_ir%pc%addpc) THEN

        CALL reg_lim_extrap_k(coef%lim_prfl_tmin(:), coef%lim_prfl_tmax(:), &
                              jotop, ProfOut_k(kn)%t(:))

        ig = coef % fmv_gas_pos(gas_id_watervapour)
        CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                              jotop, ProfOut_k(kn)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ozone)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%o3(:))
        ENDIF
        IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co2)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%co2(:))
        ENDIF
        IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_n2o)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%n2o(:))
        ENDIF
        IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%co(:))
        ENDIF
        IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ch4)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%ch4(:))
        ENDIF
        IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_so2)
          CALL reg_lim_extrap_k(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                jotop, ProfOut_k(kn)%so2(:))
        ENDIF

      ELSE ! addpc

        CALL reg_lim_extrap_k(coef_pccomp%lim_pc_prfl_tmin(:), coef_pccomp%lim_pc_prfl_tmax(:), &
                              jotop, ProfOut_k(kn)%t(:))

        CALL reg_lim_extrap_k(coef_pccomp%lim_pc_prfl_qmin(:), coef_pccomp%lim_pc_prfl_qmax(:), &
                              jotop, ProfOut_k(kn)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          CALL reg_lim_extrap_k(coef_pccomp%lim_pc_prfl_ozmin(:), coef_pccomp%lim_pc_prfl_ozmax(:), &
                                jotop, ProfOut_k(kn)%o3(:))
        ENDIF

        IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

          IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
            CALL reg_lim_extrap_k(coef_pccomp%co2_pc_min, coef_pccomp%co2_pc_max, &
                                  jotop, ProfOut_k(kn)%co2(:))
          ENDIF
          IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
            CALL reg_lim_extrap_k(coef_pccomp%n2o_pc_min, coef_pccomp%n2o_pc_max, &
                                  jotop, ProfOut_k(kn)%n2o(:))
          ENDIF
          IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
            CALL reg_lim_extrap_k(coef_pccomp%co_pc_min, coef_pccomp%co_pc_max, &
                                  jotop, ProfOut_k(kn)%co(:))
          ENDIF
          IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
            CALL reg_lim_extrap_k(coef_pccomp%ch4_pc_min, coef_pccomp%ch4_pc_max, &
                                  jotop, ProfOut_k(kn)%ch4(:))
          ENDIF

        ENDIF

      ENDIF ! addpc

    ENDIF ! reg_limit_extrap

    IF (opts%interpolation%lgradp) THEN
      DO jo = 1, kno
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = 0.0_jprb
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            ProfOut_k(kn)%t(jo) * ProfIn(prof)%t(istart(jo, prof):iend(jo, prof))
        zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
            ProfOut_k(kn)%q(jo) * ProfGasIn(prof)%q(istart(jo, prof):iend(jo, prof))
        IF (opts%rt_all%ozone_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%o3(jo) * ProfGasIn(prof)%o3(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%rt_all%co2_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%co2(jo) * ProfGasIn(prof)%co2(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%rt_all%n2o_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%n2o(jo) * ProfGasIn(prof)%n2o(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%rt_all%co_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%co(jo) * ProfGasIn(prof)%co(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%rt_all%ch4_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%ch4(jo) * ProfGasIn(prof)%ch4(istart(jo, prof):iend(jo, prof))
        ENDIF
        IF (opts%rt_all%so2_data) THEN
          zpz_k(istart(jo, prof):iend(jo, prof), jo) = zpz_k(istart(jo, prof):iend(jo, prof), jo) +      &
              ProfOut_k(kn)%so2(jo) * ProfGasIn(prof)%so2(istart(jo, prof):iend(jo, prof))
        ENDIF
      ENDDO
    ENDIF
    DO jo = 1, kno
      ProfIn_k(kn)%t(istart(jo, prof):iend(jo, prof)) = ProfIn_k(kn)%t(istart(jo, prof):iend(jo, prof)) +      &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%t(jo)
      ProfGasIn_k(kn)%q(istart(jo, prof):iend(jo, prof)) = ProfGasIn_k(kn)%q(istart(jo, prof):iend(jo, prof)) +      &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%q(jo)
      IF (opts%rt_all%ozone_data) ProfGasIn_k(kn)%o3(istart(jo, prof):iend(jo, prof))  =      &
          ProfGasIn_k(kn)%o3(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%o3(jo)
      IF (opts%rt_all%co2_data  ) ProfGasIn_k(kn)%co2(istart(jo, prof):iend(jo, prof)) =      &
          ProfGasIn_k(kn)%co2(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%co2(jo)
      IF (opts%rt_all%n2o_data  ) ProfGasIn_k(kn)%n2o(istart(jo, prof):iend(jo, prof)) =      &
          ProfGasIn_k(kn)%n2o(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%n2o(jo)
      IF (opts%rt_all%co_data   ) ProfGasIn_k(kn)%co(istart(jo, prof):iend(jo, prof))  =      &
          ProfGasIn_k(kn)%co(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%co(jo)
      IF (opts%rt_all%ch4_data  ) ProfGasIn_k(kn)%ch4(istart(jo, prof):iend(jo, prof)) =      &
          ProfGasIn_k(kn)%ch4(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%ch4(jo)
      IF (opts%rt_all%so2_data  ) ProfGasIn_k(kn)%so2(istart(jo, prof):iend(jo, prof)) =      &
          ProfGasIn_k(kn)%so2(istart(jo, prof):iend(jo, prof)) +    &
          zpz(istart(jo, prof):iend(jo, prof), jo, prof) * ProfOut_k(kn)%so2(jo)
    ENDDO ! target levels
    IF (opts%interpolation%lgradp) THEN
      zlnpi_k = 0.0_jprb
      zlnpo_k = 0.0_jprb
      CALL rttov_layeravg_k( &
              zlnpo(:, prof),  &
              zlnpo_k,         &
              zlnpi(:, prof),  &
              zlnpi_k,         &
              kno,             &
              kni,             &
              zpz(:, :, prof), &
              zpz_k,           &
              istart(:, prof), &
              iend(:, prof),   &
              interp_mode)
      pvlev_k(1:kni)        = 1 / pvlev(1:kni, prof) * zlnpi_k(1:kni)
      ProfIn_k(kn)%p(1:kni) = ProfIn_k(kn)%p(1:kni) + pvlev_k(1:kni)
    ENDIF
  ENDDO ! chan loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_K', 1_jpim, ZHOOK_HANDLE)
!------------------------------------------------------------
CONTAINS
  SUBROUTINE reg_lim_extrap_k(lim_min, lim_max, jotop, prof_k)

    REAL(KIND=jprb)   , INTENT(IN)    :: lim_min(:), lim_max(:)
    INTEGER(KIND=jpim), INTENT(IN)    :: jotop
    REAL(KIND=jprb)   , INTENT(INOUT) :: prof_k(:)

    INTEGER(KIND=jpim) :: jo
    REAL(KIND=jprb)    :: factor_k

    factor_k = 0._jprb
    DO jo = 1, jotop - 1
      factor_k = factor_k + prof_k(jo) * (lim_max(jo) - lim_min(jo))
      prof_k(jo) = 0._jprb
    ENDDO
    prof_k(jotop) = prof_k(jotop) + factor_k / (lim_max(jotop) - lim_min(jotop))

  END SUBROUTINE reg_lim_extrap_k
END SUBROUTINE rttov_intavg_prof_k
