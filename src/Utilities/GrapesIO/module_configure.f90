!!--------------------------------------------------------------------------------------------------
! PROJECT           : GRAPES IO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Zhao Liu (liuzhao@nsccsz.cn), 2021/3/18, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # GRAPES IO Module
!!
!!  *This module defines data structures for GRAPES input namelist*
!! @author Sanshan Tu
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
MODULE module_configure
  USE kinds_m
  USE NMLRead_m
  TYPE grid_config_rec_type
    INTEGER(i_kind) :: s_we
    INTEGER(i_kind) :: e_we
    INTEGER(i_kind) :: s_sn
    INTEGER(i_kind) :: e_sn
    INTEGER(i_kind) :: s_vert
    INTEGER(i_kind) :: e_vert
    LOGICAL :: global_opt
    !  integer(i_kind) :: idep
    !  integer(i_kind) :: jdep
    !  real(r_kind)    :: polar_line
    !  integer(i_kind) :: max_dom
    integer(i_kind) :: spec_bdy_width
    !  integer(i_kind) :: time_step_max
    !  integer(i_kind) :: time_section_of_output(5)
    !  integer(i_kind) :: output_freq_every_section(5)
    !  integer(i_kind) :: time_step_count_output
    !  integer(i_kind) :: step_rain_output
    !  integer(i_kind) :: step_phytend_output
    !  integer(i_kind) :: step_diagnosis
    !  integer(i_kind) :: step_begin_output
    !  integer(i_kind) :: check_dt
    !  integer(i_kind) :: maxit
    !  integer(i_kind) :: interp
    !  integer(i_kind) :: vinterp_method
    INTEGER(i_kind) :: nh
    INTEGER(i_kind) :: iotype
    !  logical :: writetofile
    !  logical :: write_modelvar
    !  logical :: write_postvar
    !  logical :: split_output
    REAL(r_kind) :: xs_we
    REAL(r_kind) :: ys_sn
    !  real(r_kind) :: cen_lat
    REAL(r_kind) :: xd
    REAL(r_kind) :: yd
    !  real(r_kind) :: dt
    !  real(r_kind) :: epson_depart
    !  real(r_kind)    :: threshold
    !  logical :: real_data
    !  logical :: balance_flow
    !  logical :: cross_pole
    !  logical :: phys_opt
    !  integer(i_kind) :: i_parent_start
    !   integer(i_kind) :: j_parent_start
    !  integer(i_kind) :: shw
    !  integer(i_kind) :: parent_grid_ratio
    !  integer(i_kind) :: tile_sz_x
    ! integer(i_kind) :: tile_sz_y
    ! integer(i_kind) :: numtiles
    !  integer(i_kind) :: nproc_x
    !  integer(i_kind) :: nproc_y
    !  logical :: reorder_mesh
    !  integer(i_kind) :: dyn_opt
    !  integer(i_kind) :: coor_opt
    !  integer(i_kind) :: diff_opt
    !  integer(i_kind) :: km_opt
    !  integer(i_kind) :: damp_opt
    !  integer(i_kind) :: spec_zone
    !  integer(i_kind) :: relax_zone
    !  integer(i_kind) :: isfflx
    !  integer(i_kind) :: ifsnow
    !  integer(i_kind) :: icloud
    !  integer(i_kind) :: if_ozone_data
    !  integer(i_kind) :: number_o3
    !  logical :: prm
    INTEGER(i_kind) :: num_soil_layers
    !  integer(i_kind) :: num_coeffs_a_b
    !  integer(i_kind) :: julyr
    !  integer(i_kind) :: julday
    !  real(r_kind) :: gmt
    !  real(r_kind) :: radt
    !  real(r_kind) :: bldt
    !  real(r_kind) :: cudt
    !  logical :: top_nudging
    !  logical :: if_smooth
    !  integer(i_kind) :: time_step_begin_restart
    !  integer(i_kind) :: io_form_restart
    !  logical :: self_test_domain
    !  logical :: mp_lc_logical
    !  integer(i_kind) :: mp_physics
    !  logical :: gwd_mb
    !  logical :: if_tofd
    !  logical :: nmrf
    !  logical :: nsas
    ! integer(i_kind) :: lc_physics
    ! integer(i_kind) :: ra_lw_physics
    !  integer(i_kind) :: ra_sw_physics
    !  integer(i_kind) :: bl_sfclay_physics
    !  integer(i_kind) :: bl_surface_physics
    !  integer(i_kind) :: bl_pbl_physics
    !  integer(i_kind) :: cu_physics
    !  logical :: tlad_gwd
    !  integer(i_kind) :: tlad_pbl
    !  integer(i_kind) :: tlad_lc
    !  integer(i_kind) :: tlad_cu
    INTEGER(i_kind) :: start_year
    INTEGER(i_kind) :: start_month
    INTEGER(i_kind) :: start_day
    INTEGER(i_kind) :: start_hour
    INTEGER(i_kind) :: start_minute
    INTEGER(i_kind) :: start_second
    !  integer(i_kind) :: end_year
    !  integer(i_kind) :: end_month
    !  integer(i_kind) :: end_day
    !  integer(i_kind) :: end_hour
    !  integer(i_kind) :: end_minute
    !  integer(i_kind) :: end_second
    !  integer(i_kind) :: interval_seconds
    !  integer(i_kind) :: real_data_init_type
    ! real(r_kind) :: bdyfrq
    ! logical :: nested
    ! logical :: specified
    ! integer(i_kind) :: iswater
    ! logical :: do_df
    ! real(r_kind) :: df_period
    ! pengf - 1) decalare parmaters in "namelist_ens" for GEPS in "TYPE_grid_config"
    !  logical :: ens_conf
    !  logical :: sppt_opt
    !  logical :: skeb_opt
    !  integer(i_kind) :: ens_sppt_memb
    !  integer(i_kind) :: ens_sppt_ncha2
    !  integer(i_kind) :: ens_sppt_trnl2
    !  integer(i_kind) :: ens_sppt_trnh2
    !  integer(i_kind) :: ens_skeb_memb
    !  integer(i_kind) :: ens_skeb_ncha2
    !  integer(i_kind) :: ens_skeb_trnl2
    !  integer(i_kind) :: ens_skeb_trnh2
    !  real(r_kind)    :: ens_sppt_min2
    !  real(r_kind)    :: ens_sppt_max2
    !  real(r_kind)    :: ens_sppt_mean2
    !  real(r_kind)    :: ens_sppt_std2
    !  real(r_kind)    :: ens_sppt_tau2
    !  real(r_kind)    :: ens_sppt_str2
    !  real(r_kind)    :: ens_skeb_min2
    !  real(r_kind)    :: ens_skeb_max2
    !  real(r_kind)    :: ens_skeb_mean2
    !  real(r_kind)    :: ens_skeb_std2
    !  real(r_kind)    :: ens_skeb_tau2
    !  real(r_kind)    :: ens_skeb_str2
    !  real(r_kind)    :: ens_skeb_alpha
    ! pengf - end
  END TYPE grid_config_rec_type

CONTAINS

  SUBROUTINE initial_config(nmlst_fn, config_flags)

    IMPLICIT NONE
    TYPE(grid_config_rec_type) :: config_flags
    CHARACTER(*)        :: nmlst_fn

    call namelist_read(nmlst_fn, 's_we', config_flags%s_we)
    call namelist_read(nmlst_fn, 'e_we', config_flags%e_we)
    call namelist_read(nmlst_fn, 's_sn', config_flags%s_sn)
    call namelist_read(nmlst_fn, 'e_sn', config_flags%e_sn)
    call namelist_read(nmlst_fn, 's_vert', config_flags%s_vert)
    call namelist_read(nmlst_fn, 'e_vert', config_flags%e_vert)
    call namelist_read(nmlst_fn, 'spec_bdy_width', config_flags%spec_bdy_width)
    ! call namelist_read(nmlst_fn, 'global_opt', config_flags%global_opt)
    ! CALL namelist_read(nmlst_fn, 'iotype', config_flags%iotype)
    CALL namelist_read(nmlst_fn, 'xd', config_flags%xd)
    CALL namelist_read(nmlst_fn, 'yd', config_flags%yd)
    CALL namelist_read(nmlst_fn, 'start_year', config_flags%start_year)
    CALL namelist_read(nmlst_fn, 'start_month', config_flags%start_month)
    CALL namelist_read(nmlst_fn, 'start_day', config_flags%start_day)
    CALL namelist_read(nmlst_fn, 'start_hour', config_flags%start_hour)
    CALL namelist_read(nmlst_fn, 'num_soil_layers', config_flags%num_soil_layers)
    CALL namelist_read(nmlst_fn, 'xs_we', config_flags%xs_we)
    CALL namelist_read(nmlst_fn, 'ys_sn', config_flags%ys_sn)
    CALL namelist_read(nmlst_fn, 'nh', config_flags%nh)

  END SUBROUTINE initial_config

END MODULE module_configure
