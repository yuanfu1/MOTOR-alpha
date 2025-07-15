SUBROUTINE rttov_intavg_prof_tl( &
              opts,         &
              kni,          &
              kno,          &
              ProfIn,       &
              ProfIn_tl,    &
              ProfGasIn,    &
              ProfGasIn_tl, &
              ProfOut,      &
              ProfOut_tl,   &
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
! ProfIn_tl     source for interpolator
! ProfOut       target for interpolator
! ProfOut_tl    target for interpolator
!
!--------------------------------------------------------------------------
! application within RTTOV (for NWP SAF)
!   - profiles on USER levels -> prediction of gas/atm opdeps on COEF levels
!   - gas/atm opdeps on COEF levels -> transmittances on USER levels for RTE
!---------------------------------------------------------------------------
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
  TYPE(rttov_profile),     INTENT(IN)    :: ProfIn_tl(SIZE(ProfIn))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfGasIn(SIZE(ProfIn))    ! atmospheric gas profiles (ppmv dry)
  TYPE(rttov_profile),     INTENT(IN)    :: ProfGasIn_tl(SIZE(ProfIn))
  TYPE(rttov_profile),     INTENT(IN)    :: ProfOut(SIZE(ProfIn))      ! atmospheric profiles
  TYPE(rttov_profile),     INTENT(INOUT) :: ProfOut_tl(SIZE(ProfIn))
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
#include "rttov_layeravg.interface"
#include "rttov_layeravg_tl.interface"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: jo, jotop, istart(kno), iend(kno)
  INTEGER(KIND=jpim) :: prof, ig
  LOGICAL(KIND=jplm) :: llevels_different
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: pvlev(kni)
  REAL   (KIND=jprb) :: pvlev_tl(kni)
  REAL   (KIND=jprb) :: ppo(kno)
  REAL   (KIND=jprb) :: zlnpi(kni), zlnpi_tl(kni), zpz(kni, kno), zpz_tl(kni, kno)
  REAL   (KIND=jprb) :: zlnpo(kno), zlnpo_tl(kno)
  INTEGER(KIND=jpim) :: interp_mode
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
!0. initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_TL', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  nprofiles         = SIZE(ProfIn)
  IF (opts%interpolation%lgradp) THEN
    llevels_different = .TRUE.
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
      IF (opts%interpolation%lgradp) THEN
        pvlev_tl(1:kni) = ProfIn_tl(prof)%p(1:kni)
        zlnpi_tl(1:kni) = pvlev_tl(1:kni) / pvlev(1:kni)
        zlnpo_tl(1:kno) = 0.0_jprb
        CALL rttov_layeravg_tl( &
                zlnpo,    &
                zlnpo_tl, &
                zlnpi,    &
                zlnpi_tl, &
                kno,      &
                kni,      &
                zpz,      &
                zpz_tl,   &
                istart,   &
                iend,     &
                interp_mode)
      ELSE
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
    ENDIF
    DO jo = 1, kno
! TL
      ProfOut_tl(prof)%t(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfIn_tl(prof)%t(istart(jo):iend(jo))))
      ProfOut_tl(prof)%q(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%q(istart(jo):iend(jo))))
      IF (opts%rt_all%ozone_data) THEN
        ProfOut_tl(prof)%o3(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%o3(istart(jo):iend(jo))))
      ENDIF
      IF (opts%rt_all%co2_data) THEN
        ProfOut_tl(prof)%co2(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%co2(istart(jo):iend(jo))))
      ENDIF
      IF (opts%rt_all%n2o_data) THEN
        ProfOut_tl(prof)%n2o(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%n2o(istart(jo):iend(jo))))
      ENDIF
      IF (opts%rt_all%co_data) THEN
        ProfOut_tl(prof)%co(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%co(istart(jo):iend(jo))))
      ENDIF
      IF (opts%rt_all%ch4_data) THEN
        ProfOut_tl(prof)%ch4(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%ch4(istart(jo):iend(jo))))
      ENDIF
      IF (opts%rt_all%so2_data) THEN
        ProfOut_tl(prof)%so2(jo) = SUM(zpz(istart(jo):iend(jo), jo) * (ProfGasIn_tl(prof)%so2(istart(jo):iend(jo))))
      ENDIF
      IF (opts%interpolation%lgradp) THEN
        ProfOut_tl(prof)%t(jo) =      &
            ProfOut_tl(prof)%t(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfIn(prof)%t(istart(jo):iend(jo))))
        ProfOut_tl(prof)%q(jo) =      &
            ProfOut_tl(prof)%q(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%q(istart(jo):iend(jo))))
        IF (opts%rt_all%ozone_data) THEN
          ProfOut_tl(prof)%o3(jo) =      &
              ProfOut_tl(prof)%o3(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%o3(istart(jo):iend(jo))))
        ENDIF
        IF (opts%rt_all%co2_data) THEN
          ProfOut_tl(prof)%co2(jo) =      &
              ProfOut_tl(prof)%co2(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%co2(istart(jo):iend(jo))))
        ENDIF
        IF (opts%rt_all%n2o_data) THEN
          ProfOut_tl(prof)%n2o(jo) =      &
              ProfOut_tl(prof)%n2o(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%n2o(istart(jo):iend(jo))))
        ENDIF
        IF (opts%rt_all%co_data) THEN
          ProfOut_tl(prof)%co(jo) =      &
              ProfOut_tl(prof)%co(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%co(istart(jo):iend(jo))))
        ENDIF
        IF (opts%rt_all%ch4_data) THEN
          ProfOut_tl(prof)%ch4(jo) =      &
              ProfOut_tl(prof)%ch4(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%ch4(istart(jo):iend(jo))))
        ENDIF
        IF (opts%rt_all%so2_data) THEN
          ProfOut_tl(prof)%so2(jo) =      &
              ProfOut_tl(prof)%so2(jo) + SUM(zpz_tl(istart(jo):iend(jo), jo) * (ProfGasIn(prof)%so2(istart(jo):iend(jo))))
        ENDIF
      ENDIF
    ENDDO ! target levels

    IF (opts%interpolation%reg_limit_extrap .AND. ppo(2) < pvlev(1)) THEN
      DO jotop = 2, kno-1
        IF (ppo(jotop+1) > pvlev(1)) EXIT
      ENDDO

      IF (.NOT. opts%rt_ir%pc%addpc) THEN

        CALL reg_lim_extrap_tl(coef%lim_prfl_tmin(:), coef%lim_prfl_tmax(:), &
                               jotop, ProfOut_tl(prof)%t(:))

        ig = coef % fmv_gas_pos(gas_id_watervapour)
        CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                               jotop, ProfOut_tl(prof)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ozone)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%o3(:))
        ENDIF
        IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co2)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%co2(:))
        ENDIF
        IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_n2o)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%n2o(:))
        ENDIF
        IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_co)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%co(:))
        ENDIF
        IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_ch4)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%ch4(:))
        ENDIF
        IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
          ig = coef % fmv_gas_pos(gas_id_so2)
          CALL reg_lim_extrap_tl(coef%lim_prfl_gmin(:, ig), coef%lim_prfl_gmax(:, ig), &
                                 jotop, ProfOut_tl(prof)%so2(:))
        ENDIF

      ELSE ! addpc

        CALL reg_lim_extrap_tl(coef_pccomp%lim_pc_prfl_tmin(:), coef_pccomp%lim_pc_prfl_tmax(:), &
                               jotop, ProfOut_tl(prof)%t(:))

        CALL reg_lim_extrap_tl(coef_pccomp%lim_pc_prfl_qmin(:), coef_pccomp%lim_pc_prfl_qmax(:), &
                               jotop, ProfOut_tl(prof)%q(:))

        IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
          CALL reg_lim_extrap_tl(coef_pccomp%lim_pc_prfl_ozmin(:), coef_pccomp%lim_pc_prfl_ozmax(:), &
                                 jotop, ProfOut_tl(prof)%o3(:))
        ENDIF

        IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

          IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
            CALL reg_lim_extrap_tl(coef_pccomp%co2_pc_min, coef_pccomp%co2_pc_max, &
                                   jotop, ProfOut_tl(prof)%co2(:))
          ENDIF
          IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
            CALL reg_lim_extrap_tl(coef_pccomp%n2o_pc_min, coef_pccomp%n2o_pc_max, &
                                   jotop, ProfOut_tl(prof)%n2o(:))
          ENDIF
          IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
            CALL reg_lim_extrap_tl(coef_pccomp%co_pc_min, coef_pccomp%co_pc_max, &
                                   jotop, ProfOut_tl(prof)%co(:))
          ENDIF
          IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
            CALL reg_lim_extrap_tl(coef_pccomp%ch4_pc_min, coef_pccomp%ch4_pc_max, &
                                   jotop, ProfOut_tl(prof)%ch4(:))
          ENDIF

        ENDIF

      ENDIF ! addpc

    ENDIF ! reg_limit_extrap

  ENDDO ! profiles
!------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_PROF_TL', 1_jpim, ZHOOK_HANDLE)
!
CONTAINS
  SUBROUTINE reg_lim_extrap_tl(lim_min, lim_max, jotop, prof_tl)

    REAL(KIND=jprb)   , INTENT(IN)    :: lim_min(:), lim_max(:)
    INTEGER(KIND=jpim), INTENT(IN)    :: jotop
    REAL(KIND=jprb)   , INTENT(INOUT) :: prof_tl(:)

    INTEGER(KIND=jpim) :: jo
    REAL(KIND=jprb)    :: factor_tl

    factor_tl = prof_tl(jotop) / (lim_max(jotop) - lim_min(jotop))
    DO jo = 1, jotop - 1
      prof_tl(jo) = factor_tl * (lim_max(jo) - lim_min(jo))
    ENDDO

  END SUBROUTINE reg_lim_extrap_tl
END SUBROUTINE rttov_intavg_prof_tl
