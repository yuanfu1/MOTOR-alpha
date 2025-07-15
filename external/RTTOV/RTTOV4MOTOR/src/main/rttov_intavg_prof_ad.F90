SUBROUTINE rttov_intavg_prof_ad( &
              opts,         &
              kni,          &
              kno,          &
              ProfIn,       &
              ProfIn_ad,    &
              ProfGasIn,    &
              ProfGasIn_ad, &
              ProfOut,      &
              ProfOut_ad,   &
              coef,         &
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
! ProfIn_ad     source for interpolator
! ProfOut       target for interpolator
! ProfOut_ad    target for interpolator
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
!
!     r1428     09/2013   Made coef optional for use with HTFRTC
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
  USE rttov_types, ONLY : rttov_options, rttov_coef, rttov_coef_pccomp, rttov_profile
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
  INTEGER(KIND=jpim),      INTENT(IN)    :: kni, kno                   ! number of levels
  TYPE(rttov_profile),     INTENT(IN)    :: ProfIn(:)                  ! atmospheric profiles
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfIn_ad(SIZE(ProfIn))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfGasIn(SIZE(ProfIn))    ! atmospheric gas profiles (ppmv dry)
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfGasIn_ad(SIZE(ProfIn))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfOut(SIZE(ProfIn))      ! atmospheric profiles
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfOut_ad(SIZE(ProfIn))
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
#include "rttov_layeravg.interface"
#include "rttov_layeravg_ad.interface"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, jotop, istart(kno), iend(kno)
  INTEGER(KIND=jpim) :: prof, ig
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: pvlev(kni)
  REAL   (KIND=jprb) :: pvlev_ad(kni)
  REAL   (KIND=jprb) :: ppo(kno)
  REAL   (KIND=jprb) :: zlnpi(kni), zlnpi_ad(kni), zpz(kni, kno), zpz_ad(kni, kno)
  REAL   (KIND=jprb) :: zlnpo(kno), zlnpo_ad(kno)
  INTEGER(KIND=jpim) :: interp_mode
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_AD', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = SIZE(ProfIn)
  IF (opts%interpolation%lgradp) THEN
    llevels_different = .TRUE.
    zpz_ad            = 0._jprb
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
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      zlnpi(1:kni) = LOG(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = LOG(ppo(1:kno))
!--------------------------------------------------------------------------
!1. profile interpolation user -> coef
!--------------------------------------------------------------------------
! replaced by profile assignments
!1.2 - interpolation procedure for profile elements (user -> coef)
!------------------------------------------------------------------
! for rttov same set of weights for all elements, so one call suffices
      CALL rttov_layeravg( &
              zlnpo,  &
              zlnpi,  &
              kno,    &
              kni,    &
              zpz,    &
              istart, &
              iend,   &
              interp_mode)
    ENDIF

    IF (opts%interpolation%reg_limit_extrap .AND. ppo(2) < pvlev(1)) THEN
      DO jotop = 2, kno-1
        IF (ppo(jotop+1) > pvlev(1)) EXIT
      ENDDO

      IF (.NOT. opts%rt_ir%pc%addpc) THEN

        CALL reg_lim_extrap_ad(coef%lim_prfl_tmin(:), coef%lim_prfl_tmax(:), &
                               jotop, ProfOut_ad(prof)%t(:))

        ig = coef % fmv_gas_pos(gas_id_watervapour)
        CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                               jotop, ProfOut_ad(prof)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ozone)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%o3(:))
        ENDIF
        IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co2)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%co2(:))
        ENDIF
        IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_n2o)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%n2o(:))
        ENDIF
        IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%co(:))
        ENDIF
        IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ch4)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%ch4(:))
        ENDIF
        IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_so2)
          CALL reg_lim_extrap_ad(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_ad(prof)%so2(:))
        ENDIF

      ELSE ! addpc

        CALL reg_lim_extrap_ad(coef_pccomp%lim_pc_prfl_tmin(:), coef_pccomp%lim_pc_prfl_tmax(:), &
                               jotop, ProfOut_ad(prof)%t(:))

        CALL reg_lim_extrap_ad(coef_pccomp%lim_pc_prfl_qmin(:), coef_pccomp%lim_pc_prfl_qmax(:), &
                               jotop, ProfOut_ad(prof)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          CALL reg_lim_extrap_ad(coef_pccomp%lim_pc_prfl_ozmin(:), coef_pccomp%lim_pc_prfl_ozmax(:), &
                                 jotop, ProfOut_ad(prof)%o3(:))
        ENDIF

        IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

          IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
            CALL reg_lim_extrap_ad(coef_pccomp%co2_pc_min, coef_pccomp%co2_pc_max, &
                                   jotop, ProfOut_ad(prof)%co2(:))
          ENDIF
          IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
            CALL reg_lim_extrap_ad(coef_pccomp%n2o_pc_min, coef_pccomp%n2o_pc_max, &
                                   jotop, ProfOut_ad(prof)%n2o(:))
          ENDIF
          IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
            CALL reg_lim_extrap_ad(coef_pccomp%co_pc_min, coef_pccomp%co_pc_max, &
                                   jotop, ProfOut_ad(prof)%co(:))
          ENDIF
          IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
            CALL reg_lim_extrap_ad(coef_pccomp%ch4_pc_min, coef_pccomp%ch4_pc_max, &
                                   jotop, ProfOut_ad(prof)%ch4(:))
          ENDIF

        ENDIF

      ENDIF ! addpc

    ENDIF ! reg_limit_extrap

    IF (opts%interpolation%lgradp) THEN
      DO jo = 1, kno
        zpz_ad(istart(jo):iend(jo), jo) = 0.0_jprb
        zpz_ad(istart(jo):iend(jo), jo) =      &
            zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%t(jo) * ProfIn(prof)%t(istart(jo):iend(jo))
        zpz_ad(istart(jo):iend(jo), jo) =      &
            zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%q(jo) * ProfGasIn(prof)%q(istart(jo):iend(jo))
        IF (opts%rt_all%ozone_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%o3(jo) * ProfGasIn(prof)%o3(istart(jo):iend(jo))
        ENDIF
        IF (opts%rt_all%co2_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%co2(jo) * ProfGasIn(prof)%co2(istart(jo):iend(jo))
        ENDIF
        IF (opts%rt_all%n2o_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%n2o(jo) * ProfGasIn(prof)%n2o(istart(jo):iend(jo))
        ENDIF
        IF (opts%rt_all%co_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%co(jo) * ProfGasIn(prof)%co(istart(jo):iend(jo))
        ENDIF
        IF (opts%rt_all%ch4_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%ch4(jo) * ProfGasIn(prof)%ch4(istart(jo):iend(jo))
        ENDIF
        IF (opts%rt_all%so2_data) THEN
          zpz_ad(istart(jo):iend(jo), jo) =      &
              zpz_ad(istart(jo):iend(jo), jo) + ProfOut_ad(prof)%so2(jo) * ProfGasIn(prof)%so2(istart(jo):iend(jo))
        ENDIF
      ENDDO
    ENDIF
    DO jo = 1, kno
      ProfIn_ad(prof)%t(istart(jo):iend(jo)) =      &
          ProfIn_ad(prof)%t(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%t(jo)
      ProfGasIn_ad(prof)%q(istart(jo):iend(jo)) =      &
          ProfGasIn_ad(prof)%q(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%q(jo)
      IF (opts%rt_all%ozone_data) ProfGasIn_ad(prof)%o3(istart(jo):iend(jo))  =      &
          ProfGasIn_ad(prof)%o3(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%o3(jo)
      IF (opts%rt_all%co2_data  ) ProfGasIn_ad(prof)%co2(istart(jo):iend(jo)) =      &
          ProfGasIn_ad(prof)%co2(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%co2(jo)
      IF (opts%rt_all%n2o_data  ) ProfGasIn_ad(prof)%n2o(istart(jo):iend(jo)) =      &
          ProfGasIn_ad(prof)%n2o(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%n2o(jo)
      IF (opts%rt_all%co_data   ) ProfGasIn_ad(prof)%co(istart(jo):iend(jo))  =      &
          ProfGasIn_ad(prof)%co(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%co(jo)
      IF (opts%rt_all%ch4_data  ) ProfGasIn_ad(prof)%ch4(istart(jo):iend(jo)) =      &
          ProfGasIn_ad(prof)%ch4(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%ch4(jo)
      IF (opts%rt_all%so2_data  ) ProfGasIn_ad(prof)%so2(istart(jo):iend(jo)) =      &
          ProfGasIn_ad(prof)%so2(istart(jo):iend(jo)) + zpz(istart(jo):iend(jo), jo) * ProfOut_ad(prof)%so2(jo)
    ENDDO ! target levels
    IF (opts%interpolation%lgradp) THEN
      zlnpi_ad = 0.0_jprb
      zlnpo_ad = 0.0_jprb
      CALL rttov_layeravg_ad( &
              zlnpo,    &
              zlnpo_ad, &
              zlnpi,    &
              zlnpi_ad, &
              kno,      &
              kni,      &
              zpz,      &
              zpz_ad,   &
              istart,   &
              iend,     &
              interp_mode)
      pvlev_ad(1:kni)          = 1 / pvlev(1:kni) * zlnpi_ad(1:kni)
      ProfIn_ad(prof)%p(1:kni) = ProfIn_ad(prof)%p(1:kni) + pvlev_ad(1:kni)
    ENDIF
  ENDDO ! profiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_AD', 1_jpim, ZHOOK_HANDLE)
!------------------------------------------------------------
CONTAINS
  SUBROUTINE reg_lim_extrap_ad(lim_min, lim_max, jotop, prof_ad)

    REAL(KIND=jprb)   , INTENT(IN)    :: lim_min(:), lim_max(:)
    INTEGER(KIND=jpim), INTENT(IN)    :: jotop
    REAL(KIND=jprb)   , INTENT(INOUT) :: prof_ad(:)

    INTEGER(KIND=jpim) :: jo
    REAL(KIND=jprb)    :: factor_ad

    factor_ad = 0._jprb
    DO jo = 1, jotop - 1
      factor_ad = factor_ad + prof_ad(jo) * (lim_max(jo) - lim_min(jo))
      prof_ad(jo) = 0._jprb
    ENDDO
    prof_ad(jotop) = prof_ad(jotop) + factor_ad / (lim_max(jotop) - lim_min(jotop))

  END SUBROUTINE reg_lim_extrap_ad
END SUBROUTINE rttov_intavg_prof_ad
