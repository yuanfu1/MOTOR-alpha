! Description:
!> @file
!!   Check input profile variables on coefficient levels for against the
!!   regression limits.
!
!> @brief
!!   Check input profile variables on coefficient levels for against the
!!   regression limits.
!!
!! @details
!!   This subroutine is called internally within RTTOV for profiles on
!!   coefficient levels. Note that unphysical profile variable values are
!!   checked by the rttov_check_profiles subroutine.
!!
!!   Profiles are compared to regression limits (for RTTOV this involves the
!!   temperature and gas profiles, while for PC-RTTOV this also includes some
!!   surface variables). If any regression limit is exceeded the
!!   qflag_reg_limits bit is set in the output radiance%quality for the
!!   channels associated with that profile. If the verbose flag is true
!!   then a warning is printed out. Values are not checked on any levels which
!!   have been extrapolated at the top of the input profile as these are always
!!   clipped to the regression limits.
!!
!!   Gas profiles are always in units of ppmv over dry air on coefficient
!!   levels (and the regression limits are in units of ppmv over dry air)
!!   so all comparisons are carried out in these units and warnings are
!!   reported in these units.
!!
!!   Users can call rttov_user_profile_checkinput before calling RTTOV to
!!   check all profile variables for problematic values and to check them
!!   against the regression limits. In this case the call to this subroutine
!!   is not required and it can be disabled by setting the do_checkinput
!!   option to FALSE. Note then that the output radiance%quality will not
!!   flag profiles which exceed the regression limits.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           options to configure the simulations
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     thermal        flags for active thermal channels
!! @param[in]     solar          flags for active solar channels
!! @param[in]     prof           profiles structure on coefficient levels
!! @param[in]     prof_user      profiles structure on user levels
!! @param[in]     coef           rttov_coef coefficient structure
!! @param[in]     coef_pccomp    PC-RTTOV coefficient structure
!! @param[in,out] radiance       radiance structure
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
SUBROUTINE rttov_check_reg_limits( &
         err,         &
         opts,        &
         chanprof,    &
         thermal,     &
         solar,       &
         prof,        &
         prof_user,   &
         coef,        &
         coef_pccomp, &
         radiance)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : &
      rttov_options,     &
      rttov_chanprof,    &
      rttov_coef,        &
      rttov_coef_pccomp, &
      rttov_profile,     &
      rttov_radiance

  USE parkind1, ONLY : jpim, jplm

!INTF_OFF
  USE rttov_const, ONLY : &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_co,          &
      gas_id_ch4,         &
      gas_id_so2,         &
      gas_id_n2o,         &
      qflag_reg_limits
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)   :: err
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),      INTENT(IN)    :: thermal(:)
  LOGICAL(KIND=jplm),      INTENT(IN)    :: solar(:)
  TYPE(rttov_profile),     INTENT(IN)    :: prof(:)
  TYPE(rttov_profile),     INTENT(IN)    :: prof_user(:)
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_radiance),    INTENT(INOUT) :: radiance
!INTF_END

#include "rttov_errorreport.interface"

  REAL(KIND=jprb)    :: wind
  REAL(KIND=jprb)    :: dp(coef%nlevels)
  INTEGER(KIND=jpim) :: firstlevel, firstuserlevel, toplevel
  INTEGER(KIND=jpim) :: ig
  INTEGER(KIND=jpim) :: nprofiles, iprof
  CHARACTER(32)      :: sprof
  LOGICAL(KIND=jplm) :: ltest(SIZE(prof(1)%p)), reg_lim_flag, lreg_lim_verbose
  REAL(KIND=jprb)    :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
  TRY

  !-------------
  ! Initialize
  !-------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_REG_LIMITS',0_jpim,ZHOOK_HANDLE)

  nprofiles = SIZE(prof)

  DO iprof = 1, nprofiles
  ! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
    WRITE(sprof,'(" (profile number = ",I8,")")') iprof
!$OMP END CRITICAL

    ! Compare profile levels and model levels
    IF (prof(iprof)%nlevels /= coef%nlevels) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid profile number of levels"//sprof)
    ENDIF

    dp(:) = ABS(coef%ref_prfl_p(:) - prof(iprof)%p(:)) / coef%ref_prfl_p(:)
    IF (ANY( dp > 0.01_jprb )) THEN
      err = errorstatus_fatal
      THROWM(err .NE. 0, "invalid profile pressure levels"//sprof)
    ENDIF


    ! We should check all levels which can contribute to output radiances.
    ! Exactly which coef levels contribute depends on the user and coef
    ! levels, the surface pressure, and the interpolation mode used.
    ! We would like to avoid applying checks to levels which may contain
    ! extrapolated values below the surface.
    ! In general we should check down at least as far as the *input*
    ! pressure level which lies on or below the surface pressure.

    ! Find first user level at or below surface
    DO firstuserlevel = prof_user(iprof)%nlevels, 2, -1
      IF (prof_user(iprof)%p(firstuserlevel-1) < prof_user(iprof)%s2m%p) EXIT
    ENDDO

    ! Find first coef level at or below firstuserlevel
    DO firstlevel = prof(iprof)%nlevels, 2, -1
      IF (prof(iprof)%p(firstlevel-1) < prof_user(iprof)%p(firstuserlevel)) EXIT
    ENDDO


    ! When the interpolator is active any extrapolated profile values (where coef
    ! levels are above the top of the user profile) should be ignored because these
    ! will be clipped to the regression limits by rttov_apply_reg_limits

    IF (opts%interpolation%addinterp .AND. &
        prof(iprof)%p(2) < prof_user(iprof)%p(1)) THEN

      ! Determine first coef level below user input profile top level
      DO toplevel = 3, prof(iprof)%nlevels
        IF (prof(iprof)%p(toplevel) > prof_user(iprof)%p(1)) EXIT
      ENDDO
    ELSE
      toplevel = 1
    ENDIF


    !---------------------------------
    ! Check against regression limits
    !---------------------------------

    reg_lim_flag = .FALSE.
    ltest = .FALSE.
    lreg_lim_verbose = opts%config%verbose .AND. .NOT. opts%config%apply_reg_limits

    IF (.NOT. opts%rt_ir%pc%addpc) THEN

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) > coef%lim_prfl_tmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("Input temperature profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_tmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) < coef%lim_prfl_tmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("Input temperature profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_tmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ig = coef%fmv_gas_pos( gas_id_watervapour )
      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("Input water vapour profile exceeds upper coef limit"//sprof, &
            PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("Input water vapour profile exceeds lower coef limit"//sprof, &
            PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ozone )
        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input ozone profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input ozone profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_co2 )
        ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CO2 profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CO2 profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_co )
        ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CO profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CO profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_n2o )
        ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input N2O profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input N2O profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ch4 )
        ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CH4 profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input CH4 profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_all%so2_data .AND. coef%nso2 > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_so2 )
        ltest(toplevel:firstlevel) = (prof(iprof)%so2(toplevel:firstlevel) > coef%lim_prfl_gmax(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input SO2 profile exceeds upper coef limit"//sprof, &
              PACK(coef%lim_prfl_gmax(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%so2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%so2(toplevel:firstlevel) < coef%lim_prfl_gmin(toplevel:firstlevel, ig))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("Input SO2 profile exceeds lower coef limit"//sprof, &
              PACK(coef%lim_prfl_gmin(toplevel:firstlevel, ig), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%so2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

    ELSE ! addpc

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_tmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("PC-RTTOV: Input temperature profile exceeds upper coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_tmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%t(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_tmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("PC-RTTOV: Input temperature profile exceeds lower coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_tmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%t(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ig = coef%fmv_gas_pos( gas_id_watervapour )
      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_qmax(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("PC-RTTOV: Input water vapour profile exceeds upper coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_qmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      ltest(toplevel:firstlevel) = (prof(iprof)%q(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_qmin(toplevel:firstlevel))
      IF (ANY(ltest)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          CALL print_info("PC-RTTOV: Input water vapour profile exceeds lower coef limit"//sprof, &
            PACK(coef_pccomp%lim_pc_prfl_qmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
            PACK(prof(iprof)%q(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
        ENDIF
      ENDIF

      IF (opts%rt_all%ozone_data .AND. coef%nozone > 0) THEN
        ig = coef%fmv_gas_pos( gas_id_ozone )
        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) > coef_pccomp%lim_pc_prfl_ozmax(toplevel:firstlevel))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("PC-RTTOV: Input ozone profile exceeds upper coef limit"//sprof, &
              PACK(coef_pccomp%lim_pc_prfl_ozmax(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF

        ltest(toplevel:firstlevel) = (prof(iprof)%o3(toplevel:firstlevel) < coef_pccomp%lim_pc_prfl_ozmin(toplevel:firstlevel))
        IF (ANY(ltest)) THEN
          reg_lim_flag = .TRUE.
          IF (lreg_lim_verbose) THEN
            CALL print_info("PC-RTTOV: Input ozone profile exceeds lower coef limit"//sprof, &
              PACK(coef_pccomp%lim_pc_prfl_ozmin(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
              PACK(prof(iprof)%o3(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
          ENDIF
        ENDIF
      ENDIF

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

        IF (opts%rt_all%co2_data .AND. coef%nco2 > 0) THEN
          ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) > coef_pccomp%co2_pc_max(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CO2 profile exceeds upper coef limit"//sprof, &
                PACK(coef_pccomp%co2_pc_max(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF

          ltest(toplevel:firstlevel) = (prof(iprof)%co2(toplevel:firstlevel) < coef_pccomp%co2_pc_min(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CO2 profile exceeds lower coef limit"//sprof, &
                PACK(coef_pccomp%co2_pc_min(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%co2(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_all%n2o_data .AND. coef%nn2o > 0) THEN
          ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) > coef_pccomp%n2o_pc_max(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input N2O profile exceeds upper coef limit"//sprof, &
                PACK(coef_pccomp%n2o_pc_max(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF

          ltest(toplevel:firstlevel) = (prof(iprof)%n2o(toplevel:firstlevel) < coef_pccomp%n2o_pc_min(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input N2O profile exceeds lower coef limit"//sprof, &
                PACK(coef_pccomp%n2o_pc_min(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%n2o(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_all%co_data .AND. coef%nco > 0) THEN
          ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) > coef_pccomp%co_pc_max(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CO profile exceeds upper coef limit"//sprof, &
                PACK(coef_pccomp%co_pc_max(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF

          ltest(toplevel:firstlevel) = (prof(iprof)%co(toplevel:firstlevel) < coef_pccomp%co_pc_min(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CO profile exceeds lower coef limit"//sprof, &
                PACK(coef_pccomp%co_pc_min(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%co(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF
        ENDIF

        IF (opts%rt_all%ch4_data .AND. coef%nch4 > 0) THEN
          ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) > coef_pccomp%ch4_pc_max(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CH4 profile exceeds upper coef limit"//sprof, &
                PACK(coef_pccomp%ch4_pc_max(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF

          ltest(toplevel:firstlevel) = (prof(iprof)%ch4(toplevel:firstlevel) < coef_pccomp%ch4_pc_min(toplevel:firstlevel))
          IF (ANY(ltest)) THEN
            reg_lim_flag = .TRUE.
            IF (lreg_lim_verbose) THEN
              CALL print_info("PC-RTTOV: Input CH4 profile exceeds lower coef limit"//sprof, &
                PACK(coef_pccomp%ch4_pc_min(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%p(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)), &
                PACK(prof(iprof)%ch4(toplevel:firstlevel), mask=ltest(toplevel:firstlevel)))
            ENDIF
          ENDIF
        ENDIF

      ENDIF

      IF ((prof(iprof)%s2m%p < coef_pccomp%lim_pc_prfl_pmin) .OR. &
          (prof(iprof)%s2m%p > coef_pccomp%lim_pc_prfl_pmax)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          WARN("PC-RTTOV: surface pressure outside limits"//sprof)
        ENDIF
      ENDIF

      IF ((prof(iprof)%s2m%t < coef_pccomp%lim_pc_prfl_tsmin) .OR. &
          (prof(iprof)%s2m%t > coef_pccomp%lim_pc_prfl_tsmax)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          WARN("PC-RTTOV: surface temperature outside limits"//sprof)
        ENDIF
      ENDIF

      IF ((prof(iprof)%skin%t < coef_pccomp%lim_pc_prfl_skmin) .OR.   &
          (prof(iprof)%skin%t > coef_pccomp%lim_pc_prfl_skmax)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          WARN("PC-RTTOV: skin temperature outside limits"//sprof)
        ENDIF
      ENDIF

      wind = SQRT(prof(iprof)%s2m%u * prof(iprof)%s2m%u + &
                  prof(iprof)%s2m%v * prof(iprof)%s2m%v)
      IF ((wind < coef_pccomp%lim_pc_prfl_wsmin) .OR.  &
          (wind > coef_pccomp%lim_pc_prfl_wsmax)) THEN
        reg_lim_flag = .TRUE.
        IF (lreg_lim_verbose) THEN
          WARN("PC-RTTOV: 10m wind speed outside limits"//sprof)
        ENDIF
      ENDIF

    ENDIF ! addpc

    IF (reg_lim_flag) THEN
      ! Check thermal or solar to ensure consistent output from parallel interface
      WHERE ((chanprof%prof == iprof) .AND. (thermal .OR. solar))
        radiance%quality(1:SIZE(chanprof)) = IBSET(radiance%quality(1:SIZE(chanprof)), qflag_reg_limits)
      ENDWHERE
    ENDIF

  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_REG_LIMITS',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_REG_LIMITS',1_jpim,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE print_info(msg1, limits, levels, values)
    CHARACTER(LEN=*), INTENT(IN) :: msg1
    REAL(KIND=jprb),  INTENT(IN) :: limits(:)
    REAL(KIND=jprb),  INTENT(IN) :: levels(:)
    REAL(KIND=jprb),  INTENT(IN) :: values(:)

    CHARACTER(LEN=256) :: msg2
    INTEGER(KIND=jpim) :: imax, lmax

    imax = MIN(10, SIZE(levels))
    lmax = MIN(imax, SIZE(limits))

! Replace warn/info macros from throw.h with in-line code because
! NAG v5.3 complains otherwise.

! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
!     WARN(msg1)
    CALL rttov_errorreport(errorstatus_success, TRIM(msg1), 'rttov_check_reg_limits.F90')
    WRITE(msg2, '(a,10f10.4)') 'Limit   = ',limits(1:lmax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'p (hPa) = ',levels(1:imax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
    WRITE(msg2, '(a,10f10.4)') 'Value   = ',values(1:imax)
!     INFO(TRIM(msg2))
    CALL rttov_errorreport(errorstatus_success, TRIM(msg2))
!$OMP END CRITICAL

  END SUBROUTINE print_info

END SUBROUTINE rttov_check_reg_limits
