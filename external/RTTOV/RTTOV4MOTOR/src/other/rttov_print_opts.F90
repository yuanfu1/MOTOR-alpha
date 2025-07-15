! Description:
!> @file
!!   Print out the contents of the rttov_options structure.
!
!> @brief
!!   Print out the contents of the rttov_options structure.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   opts      options to configure the simulations
!! @param[in]   lu        logical unit for output, optional
!! @param[in]   text      additional text to print, optional
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
SUBROUTINE rttov_print_opts(opts, lu, text)

  USE rttov_types, ONLY : rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
  USE rttov_const, ONLY : cloud_overlap_simple
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options), INTENT(IN)           :: opts  ! options
  INTEGER(KIND=jpim),  INTENT(IN), OPTIONAL :: lu    ! logical unit for print
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: text  ! additional text to print
!INTF_END

  INTEGER(KIND=jpim)  :: iu       ! logical unit for print
  CHARACTER(LEN=20)   :: tmp_text ! temporary string for formatting output

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV options structure: ", TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV options structure"
  END IF

  WRITE(iu,'(a)')           "General config options"
  WRITE(iu,'(2x,a,l1)')     "apply_reg_limits        ", opts%config%apply_reg_limits
  WRITE(iu,'(2x,a,l1)')     "verbose                 ", opts%config%verbose
  WRITE(iu,'(2x,a,l1)')     "do_checkinput           ", opts%config%do_checkinput
  WRITE(iu,'(2x,a,l1)')     "fix_hgpl                ", opts%config%fix_hgpl
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "General RT options"
  WRITE(iu,'(2x,a,l1)')     "ozone_data              ", opts%rt_all%ozone_data
  WRITE(iu,'(2x,a,l1)')     "co2_data                ", opts%rt_all%co2_data
  WRITE(iu,'(2x,a,l1)')     "n2o_data                ", opts%rt_all%n2o_data
  WRITE(iu,'(2x,a,l1)')     "co_data                 ", opts%rt_all%co_data
  WRITE(iu,'(2x,a,l1)')     "ch4_data                ", opts%rt_all%ch4_data
  WRITE(iu,'(2x,a,l1)')     "so2_data                ", opts%rt_all%so2_data
  WRITE(iu,'(2x,a,l1)')     "addrefrac               ", opts%rt_all%addrefrac
  WRITE(iu,'(2x,a,l1)')     "plane_parallel          ", opts%rt_all%plane_parallel
  WRITE(iu,'(2x,a,l1)')     "switchrad               ", opts%rt_all%switchrad
  WRITE(iu,'(2x,a,l1)')     "use_t2m_opdep           ", opts%rt_all%use_t2m_opdep
  WRITE(iu,'(2x,a,l1)')     "use_q2m                 ", opts%rt_all%use_q2m
  WRITE(iu,'(2x,a,l1)')     "do_lambertian           ", opts%rt_all%do_lambertian
  IF (opts%rt_all%do_lambertian) THEN
    WRITE(iu,'(2x,a,l1)')     "lambertian_fixed_angle  ", opts%rt_all%lambertian_fixed_angle
  ENDIF
  WRITE(iu,'(2x,a,l1)')     "rad_down_lin_tau        ", opts%rt_all%rad_down_lin_tau
  WRITE(iu,'(2x,a,l1)')     "dtau_test               ", opts%rt_all%dtau_test
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "VIS/IR-only RT options  "
  WRITE(tmp_text,'(i3)') opts%rt_ir%ir_sea_emis_model
  WRITE(iu,'(2x,a,a)')      "ir_sea_emis_model       ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "addsolar                ", opts%rt_ir%addsolar
  IF (opts%rt_ir%addsolar) THEN
    WRITE(tmp_text,'(i3)') opts%rt_ir%solar_sea_brdf_model
    WRITE(iu,'(2x,a,a)')      "solar_sea_brdf_model    ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(f10.4)') opts%rt_ir%rayleigh_max_wavelength
    WRITE(iu,'(2x,a,a)')      "rayleigh_max_wavelength ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(f10.4)') opts%rt_ir%rayleigh_min_pressure
    WRITE(iu,'(2x,a,a)')      "rayleigh_min_pressure   ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)')     "rayleigh_single_scatt   ", opts%rt_ir%rayleigh_single_scatt
  ENDIF
  WRITE(iu,'(2x,a,l1)')     "do_nlte_correction      ", opts%rt_ir%do_nlte_correction
  WRITE(iu,'(2x,a,l1)')     "addaerosl               ", opts%rt_ir%addaerosl
  IF (opts%rt_ir%addaerosl) THEN
    WRITE(iu,'(2x,a,l1)')     "user_aer_opt_param      ", opts%rt_ir%user_aer_opt_param
  ENDIF
  WRITE(iu,'(2x,a,l1)')     "addclouds               ", opts%rt_ir%addclouds
  IF (opts%rt_ir%addclouds) THEN
    WRITE(iu,'(2x,a,l1)')     "user_cld_opt_param      ", opts%rt_ir%user_cld_opt_param
    WRITE(iu,'(2x,a,l1)')     "grid_box_avg_cloud      ", opts%rt_ir%grid_box_avg_cloud
    WRITE(tmp_text,'(g12.4)') opts%rt_ir%cldcol_threshold
    WRITE(iu,'(2x,a,a)')      "cldcol_threshold        ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') opts%rt_ir%cloud_overlap
    WRITE(iu,'(2x,a,a)')      "cloud_overlap           ", TRIM(ADJUSTL(tmp_text))
    IF (opts%rt_ir%cloud_overlap == cloud_overlap_simple) THEN
      WRITE(tmp_text,'(f10.4)') opts%rt_ir%cc_low_cloud_top
      WRITE(iu,'(2x,a,a)')      "cc_low_cloud_top        ", TRIM(ADJUSTL(tmp_text))
    ENDIF
  ENDIF
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    WRITE(tmp_text,'(i3)') opts%rt_ir%ir_scatt_model
    WRITE(iu,'(2x,a,a)')      "ir_scatt_model          ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') opts%rt_ir%vis_scatt_model
    WRITE(iu,'(2x,a,a)')      "vis_scatt_model         ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i5)') opts%rt_ir%dom_nstreams
    WRITE(iu,'(2x,a,a)')      "dom_nstreams            ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(f10.4)') opts%rt_ir%dom_accuracy
    WRITE(iu,'(2x,a,a)')      "dom_accuracy            ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(f10.4)') opts%rt_ir%dom_opdep_threshold
    WRITE(iu,'(2x,a,a)')      "dom_opdep_threshold     ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)')     "dom_rayleigh            ", opts%rt_ir%dom_rayleigh
  ENDIF
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "PC-RTTOV options"
  WRITE(iu,'(2x,a,l1)')     "addpc                   ", opts%rt_ir%pc%addpc
  IF (opts%rt_ir%pc%addpc) THEN
    WRITE(tmp_text,'(i3)') opts%rt_ir%pc%ipcbnd
    WRITE(iu,'(2x,a,a)')      "ipcbnd                  ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') opts%rt_ir%pc%ipcreg
    WRITE(iu,'(2x,a,a)')      "ipcreg                  ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') opts%rt_ir%pc%npcscores
    WRITE(iu,'(2x,a,a)')      "npcscores               ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)')     "addradrec               ", opts%rt_ir%pc%addradrec
  ENDIF
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "HTFRTC options"
  WRITE(iu,'(2x,a,l1)')     "htfrtc                  ", opts%htfrtc_opts%htfrtc
  IF (opts%htfrtc_opts%htfrtc) THEN
    WRITE(tmp_text,'(i3)') opts%htfrtc_opts%n_pc_in
    WRITE(iu,'(2x,a,a)')      "n_pc_in                 ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)')     "reconstruct             ", opts%htfrtc_opts%reconstruct
    WRITE(iu,'(2x,a,l1)')     "simple_cloud            ", opts%htfrtc_opts%simple_cloud
    WRITE(iu,'(2x,a,l1)')     "overcast                ", opts%htfrtc_opts%overcast
  ENDIF
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "MW-only RT options"
  WRITE(iu,'(2x,a,l1)')     "clw_data                ", opts%rt_mw%clw_data
  IF (opts%rt_mw%clw_data) THEN
    WRITE(tmp_text,'(i3)') opts%rt_mw%clw_scheme
    WRITE(iu,'(2x,a,a)')      "clw_scheme              ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(f10.4)') opts%rt_mw%clw_cloud_top
    WRITE(iu,'(2x,a,a)')      "clw_cloud_top           ", TRIM(ADJUSTL(tmp_text))
  ENDIF
  WRITE(tmp_text,'(i3)') opts%rt_mw%fastem_version
  WRITE(iu,'(2x,a,a)')      "fastem_version          ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "supply_foam_fraction    ", opts%rt_mw%supply_foam_fraction
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "Interpolation and vertical grid options"
  WRITE(iu,'(2x,a,l1)')     "addinterp               ", opts%interpolation%addinterp
  IF (opts%interpolation%addinterp) THEN
    WRITE(tmp_text,'(i3)') opts%interpolation%interp_mode
    WRITE(iu,'(2x,a,a)')      "interp_mode             ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)')     "lgradp                  ", opts%interpolation%lgradp
    WRITE(iu,'(2x,a,l1)')     "reg_limit_extrap        ", opts%interpolation%reg_limit_extrap
  ENDIF
  WRITE(iu,'(2x,a,l1)')     "spacetop                ", opts%interpolation%spacetop
  WRITE(iu,'(a)',advance='yes')

END SUBROUTINE rttov_print_opts
