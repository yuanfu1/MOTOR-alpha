! Description:
!> @file
!!   Check equality of two RTTOV options structures.
!
!> @brief
!!   Check equality of two RTTOV options structures. Only options which
!!   affect calculation results are compared (e.g. not the verbose flag).
!!
!! @param[in]     opts1                first options structure to compare
!! @param[in]     opts2                second options structure to compare
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
FUNCTION rttov_opts_eq(opts1, opts2)

  USE parkind1, ONLY : jplm
  USE rttov_types, ONLY : rttov_options

  IMPLICIT NONE

  LOGICAL(KIND=jplm) :: rttov_opts_eq

  TYPE(rttov_options), INTENT(IN) :: opts1
  TYPE(rttov_options), INTENT(IN) :: opts2
!INTF_END

! Only options which affect calculation results are compared (e.g. not verbose)
  rttov_opts_eq = &
   ( opts1%config%apply_reg_limits        .EQV. opts2%config%apply_reg_limits        ) .AND. &
   ( opts1%config%do_checkinput           .EQV. opts2%config%do_checkinput           ) .AND. &
   ( opts1%config%fix_hgpl                .EQV. opts2%config%fix_hgpl                ) .AND. &
   ( opts1%rt_all%switchrad               .EQV. opts2%rt_all%switchrad               ) .AND. &
   ( opts1%rt_all%addrefrac               .EQV. opts2%rt_all%addrefrac               ) .AND. &
   ( opts1%rt_all%use_t2m_opdep           .EQV. opts2%rt_all%use_t2m_opdep           ) .AND. &
   ( opts1%rt_all%use_q2m                 .EQV. opts2%rt_all%use_q2m                 ) .AND. &
   ( opts1%rt_all%do_lambertian           .EQV. opts2%rt_all%do_lambertian           ) .AND. &
   ( opts1%rt_all%lambertian_fixed_angle  .EQV. opts2%rt_all%lambertian_fixed_angle  ) .AND. &
   ( opts1%rt_all%plane_parallel          .EQV. opts2%rt_all%plane_parallel          ) .AND. &
   ( opts1%rt_all%rad_down_lin_tau        .EQV. opts2%rt_all%rad_down_lin_tau        ) .AND. &
   ( opts1%rt_all%dtau_test               .EQV. opts2%rt_all%dtau_test               ) .AND. &
   ( opts1%rt_ir%solar_sea_brdf_model      ==   opts2%rt_ir%solar_sea_brdf_model     ) .AND. &
   ( opts1%rt_ir%ir_sea_emis_model         ==   opts2%rt_ir%ir_sea_emis_model        ) .AND. &
   ( opts1%rt_ir%addsolar                 .EQV. opts2%rt_ir%addsolar                 ) .AND. &
   ( opts1%rt_ir%rayleigh_max_wavelength   ==   opts2%rt_ir%rayleigh_max_wavelength  ) .AND. &
   ( opts1%rt_ir%rayleigh_min_pressure     ==   opts2%rt_ir%rayleigh_min_pressure    ) .AND. &
   ( opts1%rt_ir%rayleigh_single_scatt    .EQV. opts2%rt_ir%rayleigh_single_scatt    ) .AND. &
   ( opts1%rt_ir%do_nlte_correction       .EQV. opts2%rt_ir%do_nlte_correction       ) .AND. &
   ( opts1%rt_ir%addaerosl                .EQV. opts2%rt_ir%addaerosl                ) .AND. &
   ( opts1%rt_ir%user_aer_opt_param       .EQV. opts2%rt_ir%user_aer_opt_param       ) .AND. &
   ( opts1%rt_ir%addclouds                .EQV. opts2%rt_ir%addclouds                ) .AND. &
   ( opts1%rt_ir%user_cld_opt_param       .EQV. opts2%rt_ir%user_cld_opt_param       ) .AND. &
   ( opts1%rt_ir%grid_box_avg_cloud       .EQV. opts2%rt_ir%grid_box_avg_cloud       ) .AND. &
   ( opts1%rt_ir%cldcol_threshold          ==   opts2%rt_ir%cldcol_threshold         ) .AND. &
   ( opts1%rt_ir%cloud_overlap             ==   opts2%rt_ir%cloud_overlap            ) .AND. &
   ( opts1%rt_ir%cc_low_cloud_top          ==   opts2%rt_ir%cc_low_cloud_top         ) .AND. &
   ( opts1%rt_ir%vis_scatt_model           ==   opts2%rt_ir%vis_scatt_model          ) .AND. &
   ( opts1%rt_ir%ir_scatt_model            ==   opts2%rt_ir%ir_scatt_model           ) .AND. &
   ( opts1%rt_ir%dom_nstreams              ==   opts2%rt_ir%dom_nstreams             ) .AND. &
   ( opts1%rt_ir%dom_accuracy              ==   opts2%rt_ir%dom_accuracy             ) .AND. &
   ( opts1%rt_ir%dom_opdep_threshold       ==   opts2%rt_ir%dom_opdep_threshold      ) .AND. &
   ( opts1%rt_ir%dom_rayleigh             .EQV. opts2%rt_ir%dom_rayleigh             ) .AND. &
   ( opts1%rt_all%ozone_data              .EQV. opts2%rt_all%ozone_data              ) .AND. &
   ( opts1%rt_all%co2_data                .EQV. opts2%rt_all%co2_data                ) .AND. &
   ( opts1%rt_all%n2o_data                .EQV. opts2%rt_all%n2o_data                ) .AND. &
   ( opts1%rt_all%co_data                 .EQV. opts2%rt_all%co_data                 ) .AND. &
   ( opts1%rt_all%ch4_data                .EQV. opts2%rt_all%ch4_data                ) .AND. &
   ( opts1%rt_all%so2_data                .EQV. opts2%rt_all%so2_data                ) .AND. &
   ( opts1%rt_ir%pc%addpc                 .EQV. opts2%rt_ir%pc%addpc                 ) .AND. &
   ( opts1%rt_ir%pc%ipcbnd                 ==   opts2%rt_ir%pc%ipcbnd                ) .AND. &
   ( opts1%rt_ir%pc%ipcreg                 ==   opts2%rt_ir%pc%ipcreg                ) .AND. &
   ( opts1%rt_ir%pc%npcscores              ==   opts2%rt_ir%pc%npcscores             ) .AND. &
   ( opts1%rt_ir%pc%addradrec             .EQV. opts2%rt_ir%pc%addradrec             ) .AND. &
   ( opts1%rt_mw%fastem_version            ==   opts2%rt_mw%fastem_version           ) .AND. &
   ( opts1%rt_mw%clw_data                 .EQV. opts2%rt_mw%clw_data                 ) .AND. &
   ( opts1%rt_mw%clw_scheme                ==   opts2%rt_mw%clw_scheme               ) .AND. &
   ( opts1%rt_mw%clw_cloud_top             ==   opts2%rt_mw%clw_cloud_top            ) .AND. &
   ( opts1%rt_mw%supply_foam_fraction     .EQV. opts2%rt_mw%supply_foam_fraction     ) .AND. &
   ( opts1%interpolation%addinterp        .EQV. opts2%interpolation%addinterp        ) .AND. &
   ( opts1%interpolation%interp_mode       ==   opts2%interpolation%interp_mode      ) .AND. &
   ( opts1%interpolation%reg_limit_extrap .EQV. opts2%interpolation%reg_limit_extrap ) .AND. &
   ( opts1%interpolation%spacetop         .EQV. opts2%interpolation%spacetop         ) .AND. &
   ( opts1%interpolation%lgradp           .EQV. opts2%interpolation%lgradp           )

END FUNCTION
