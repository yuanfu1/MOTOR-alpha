! Description:
!> @file
!!   Clip profile variables to regression limits TL.
!
!> @brief
!!   Clip profile variables to regression limits TL.
!!
!! @details
!!   This subroutine is called for profiles on coefficient levels if the
!!   apply_reg_limits option is set or if the interpolator is active.
!!
!!   In the former case values are checked on all coefficient levels at
!!   least as far down as the first user pressure level below the surface.
!!   This ensures all values which contribute to the optical depth
!!   regression are covered.
!!
!!   In the latter case values are checked on all coefficient levels which
!!   contain extrapolated values from the user profile at the top of the
!!   atmosphere. This ensures that any profile extrapolation does not
!!   result in values beyond the range of the training data (and it also
!!   avoids warnings about such values).
!!
!!   For standard RTTOV the temperature and gas profiles are checked. For
!!   PC-RTTOV some surface variables are also checked.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     prof_user      profiles structure on user levels
!! @param[in]     prof           profiles structure on coefficient levels
!! @param[in,out] prof_tl        profile perturbations on coefficient levels
!! @param[in]     coef           rttov_coef coefficient structure
!! @param[in]     coef_pccomp    PC-RTTOV coefficient structure
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_apply_reg_limits_tl( &
        opts,        &
        prof_user,   &
        prof,        &
        prof_tl,     &
        coef,        &
        coef_pccomp)

  USE rttov_types, ONLY : &
    rttov_coef,        &
    rttov_options,     &
    rttov_coef_pccomp, &
    rttov_profile

!INTF_OFF
  USE parkind1, ONLY : jprb, jpim

  USE rttov_const, ONLY : &
    gas_id_watervapour, &
    gas_id_ozone,       &
    gas_id_co2,         &
    gas_id_co,          &
    gas_id_n2o,         &
    gas_id_ch4,         &
    gas_id_so2

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: prof_user(:) ! Profiles on user levels (only p(:) and 2m p used)
  TYPE(rttov_profile),     INTENT(IN)    :: prof(:)      ! Profiles on coef levels (ppmv dry)
  TYPE(rttov_profile),     INTENT(INOUT) :: prof_tl(:)   ! Profiles_tl on coef levels (ppmv dry)
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END

  REAL(KIND=jprb)    :: wind
  INTEGER(KIND=jpim) :: firstlevel, firstuserlevel, ilev
  INTEGER(KIND=jpim) :: ig
  INTEGER(KIND=jpim) :: nprofiles, iprof
  REAL(KIND=jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_REG_LIMITS_TL',0_jpim,ZHOOK_HANDLE)

  nprofiles = SIZE(prof)

  DO iprof = 1, nprofiles

    ! If apply_reg_limits is true then check all profile values from TOA to surface.
    ! Otherwise if addinterp is true then check only the values above the top of
    !   the input profile (i.e. extrapolated values).

    IF (opts%config%apply_reg_limits) THEN

      ! Find first user level at or below surface
      DO firstuserlevel = prof_user(iprof)%nlevels, 2, -1
        IF (prof_user(iprof)%p(firstuserlevel-1) < prof_user(iprof)%s2m%p) EXIT
      ENDDO

      ! Find first coef level at or below firstuserlevel
      DO firstlevel = prof(iprof)%nlevels, 2, -1
        IF (prof(iprof)%p(firstlevel-1) < prof_user(iprof)%p(firstuserlevel)) EXIT
      ENDDO

    ELSE IF (opts%interpolation%addinterp) THEN

      IF (prof(iprof)%p(2) < prof_user(iprof)%p(1)) THEN

        ! Determine first coef level above user input profile top level
        DO firstlevel = 2, prof(iprof)%nlevels - 1
          IF (prof(iprof)%p(firstlevel+1) > prof_user(iprof)%p(1)) EXIT
        ENDDO

      ELSE
        CYCLE
      ENDIF

    ELSE

      ! Don't do anything if apply_reg_limits and addinterp are both false
      EXIT

    ENDIF

    IF (.NOT. opts%rt_ir%pc%addpc) THEN

      DO ilev = 1, firstlevel
        IF (prof(iprof)%t(ilev) > coef%lim_prfl_tmax(ilev)) THEN
          prof_tl(iprof)%t(ilev) = 0._jprb
        ELSE IF (prof(iprof)%t(ilev) < coef%lim_prfl_tmin(ilev)) THEN
          prof_tl(iprof)%t(ilev) = 0._jprb
        ENDIF
      ENDDO

      ig = coef%fmv_gas_pos(gas_id_watervapour)
      DO ilev = 1, firstlevel
        IF (prof(iprof)%q(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
          prof_tl(iprof)%q(ilev) = 0._jprb
        ELSE IF (prof(iprof)%q(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
          prof_tl(iprof)%q(ilev) = 0._jprb
        ENDIF
      ENDDO

      IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_ozone)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%o3(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%o3(ilev) = 0._jprb
          ELSE IF (prof(iprof)%o3(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%o3(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_co2)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%co2(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%co2(ilev) = 0._jprb
          ELSE IF (prof(iprof)%co2(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%co2(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_co)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%co(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%co(ilev) = 0._jprb
          ELSE IF (prof(iprof)%co(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%co(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_n2o)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%n2o(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%n2o(ilev) = 0._jprb
          ELSE IF (prof(iprof)%n2o(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%n2o(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_ch4)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%ch4(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%ch4(ilev) = 0._jprb
          ELSE IF (prof(iprof)%ch4(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%ch4(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
        ig = coef%fmv_gas_pos(gas_id_so2)
        DO ilev = 1, firstlevel
          IF (prof(iprof)%so2(ilev) > coef%lim_prfl_gmax(ilev, ig)) THEN
            prof_tl(iprof)%so2(ilev) = 0._jprb
          ELSE IF (prof(iprof)%so2(ilev) < coef%lim_prfl_gmin(ilev, ig)) THEN
            prof_tl(iprof)%so2(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

    ELSE ! addpc

      DO ilev = 1, firstlevel
        IF (prof(iprof)%t(ilev) > coef_pccomp%lim_pc_prfl_tmax(ilev)) THEN
          prof_tl(iprof)%t(ilev) = 0._jprb
        ELSE IF (prof(iprof)%t(ilev) < coef_pccomp%lim_pc_prfl_tmin(ilev)) THEN
          prof_tl(iprof)%t(ilev) = 0._jprb
        ENDIF
      ENDDO

      DO ilev = 1, firstlevel
        IF (prof(iprof)%q(ilev) > coef_pccomp%lim_pc_prfl_qmax(ilev)) THEN
          prof_tl(iprof)%q(ilev) = 0._jprb
        ELSE IF (prof(iprof)%q(ilev) < coef_pccomp%lim_pc_prfl_qmin(ilev)) THEN
          prof_tl(iprof)%q(ilev) = 0._jprb
        ENDIF
      ENDDO

      IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
        DO ilev = 1, firstlevel
          IF (prof(iprof)%o3(ilev) > coef_pccomp%lim_pc_prfl_ozmax(ilev)) THEN
            prof_tl(iprof)%o3(ilev) = 0._jprb
          ELSE IF (prof(iprof)%o3(ilev) < coef_pccomp%lim_pc_prfl_ozmin(ilev)) THEN
            prof_tl(iprof)%o3(ilev) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

        IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
          DO ilev = 1, firstlevel
            IF (prof(iprof)%co2(ilev) > coef_pccomp%co2_pc_max(ilev)) THEN
              prof_tl(iprof)%co2(ilev) = 0._jprb
            ELSE IF (prof(iprof)%co2(ilev) < coef_pccomp%co2_pc_min(ilev)) THEN
              prof_tl(iprof)%co2(ilev) = 0._jprb
            ENDIF
          ENDDO
        ENDIF

        IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
          DO ilev = 1, firstlevel
            IF (prof(iprof)%n2o(ilev) > coef_pccomp%n2o_pc_max(ilev)) THEN
              prof_tl(iprof)%n2o(ilev) = 0._jprb
            ELSE IF (prof(iprof)%n2o(ilev) < coef_pccomp%n2o_pc_min(ilev)) THEN
              prof_tl(iprof)%n2o(ilev) = 0._jprb
            ENDIF
          ENDDO
        ENDIF

        IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
          DO ilev = 1, firstlevel
            IF (prof(iprof)%co(ilev) > coef_pccomp%co_pc_max(ilev)) THEN
              prof_tl(iprof)%co(ilev) = 0._jprb
            ELSE IF (prof(iprof)%co(ilev) < coef_pccomp%co_pc_min(ilev)) THEN
              prof_tl(iprof)%co(ilev) = 0._jprb
            ENDIF
          ENDDO
        ENDIF

        IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
          DO ilev = 1, firstlevel
            IF (prof(iprof)%ch4(ilev) > coef_pccomp%ch4_pc_max(ilev)) THEN
              prof_tl(iprof)%ch4(ilev) = 0._jprb
            ELSE IF (prof(iprof)%ch4(ilev) < coef_pccomp%ch4_pc_min(ilev)) THEN
              prof_tl(iprof)%ch4(ilev) = 0._jprb
            ENDIF
          ENDDO
        ENDIF

      ENDIF

      IF (prof(iprof)%s2m%p < coef_pccomp%lim_pc_prfl_pmin) THEN
        prof_tl(iprof)%s2m%p = 0._jprb
      ELSE IF (prof(iprof)%s2m%p > coef_pccomp%lim_pc_prfl_pmax) THEN
        prof_tl(iprof)%s2m%p = 0._jprb
      ENDIF

      IF (prof(iprof)%s2m%t < coef_pccomp%lim_pc_prfl_tsmin) THEN
        prof_tl(iprof)%s2m%t = 0._jprb
      ELSE IF (prof(iprof)%s2m%t > coef_pccomp%lim_pc_prfl_tsmax) THEN
        prof_tl(iprof)%s2m%t = 0._jprb
      ENDIF

      IF (prof(iprof)%skin%t < coef_pccomp%lim_pc_prfl_skmin) THEN
        prof_tl(iprof)%skin%t = 0._jprb
      ELSE IF (prof(iprof)%skin%t > coef_pccomp%lim_pc_prfl_skmax) THEN
        prof_tl(iprof)%skin%t = 0._jprb
      ENDIF

      wind = SQRT(prof(iprof)%s2m%u * prof(iprof)%s2m%u + &
                  prof(iprof)%s2m%v * prof(iprof)%s2m%v)

      IF (wind < coef_pccomp%lim_pc_prfl_wsmin) THEN
        prof_tl(iprof)%s2m%u = 0._jprb
        prof_tl(iprof)%s2m%v = 0._jprb
      ELSE IF (wind > coef_pccomp%lim_pc_prfl_wsmax) THEN
        prof_tl(iprof)%s2m%u = 0._jprb
        prof_tl(iprof)%s2m%v = 0._jprb
      ENDIF

    ENDIF ! addpc

  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_REG_LIMITS_TL',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_apply_reg_limits_tl
