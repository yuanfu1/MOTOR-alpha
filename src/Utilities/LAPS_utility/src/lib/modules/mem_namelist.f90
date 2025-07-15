MODULE mem_namelist

  INCLUDE 'lapsparms.for'

!       Globally used variables that are independent of the namelists
!       character(len=200) :: generic_data_root, cstaticdir, grid_fnam_common

!       Declarations for namelist variables
  INTEGER iflag_lapsparms
  REAL max_radar_files_nl   ! max_radar_files is in lapsparms.for
  REAL PRESSURE_INTERVAL_L
  REAL PRESSURE_0_L
  INTEGER nk_laps
  REAL standard_latitude
  REAL standard_latitude2
  REAL standard_longitude
  INTEGER NX_L
  INTEGER NY_L
  INTEGER I_PERIMETER
  REAL grid_spacing_m
  REAL grid_cen_lat
  REAL grid_cen_lon
  REAL earth_radius
  INTEGER laps_cycle_time
  INTEGER model_cycle_time, model_fcst_intvl, model_fcst_len
  REAL purge_time
  INTEGER i2_missing_data
  INTEGER iverbose
  REAL r_missing_data
  INTEGER MAX_RADARS
  REAL aod, aod_bin(3), aod_asy(3, 3), fcterm, aod_ha, ht_ha(4), ssa(3)
  REAL aero_scaleht
  REAL ref_base
  REAL ref_base_useable
  REAL r_hybrid_first_gate
  REAL aircraft_time_window
  REAL hydrometeor_scale_pcp
  REAL hydrometeor_scale_cld
  INTEGER maxstns
  INTEGER N_PIREP
  INTEGER max_snd_grid
  INTEGER max_snd_levels
  REAL redp_lvl
  REAL prtop
  INTEGER vert_rad_meso
  INTEGER vert_rad_sao
  INTEGER vert_rad_pirep
  INTEGER vert_rad_prof
  REAL silavwt_parm
  REAL toptwvl_parm
  INTEGER iwrite_output
  INTEGER i_offset_radar

  CHARACTER * 40 vertical_grid
  CHARACTER * 10 lvl_coord_cdf
  CHARACTER * 50 c50_lowres_directory
  CHARACTER * 6 c6_maproj
  CHARACTER * 8 radarext_3d
  CHARACTER * 8 radarext_3d_accum
  CHARACTER * 200 path_to_raw_pirep
  CHARACTER * 200 path_to_raw_rass
  CHARACTER * 200 path_to_raw_profiler
  CHARACTER * 200 path_to_raw_blprass
  CHARACTER * 200 path_to_raw_blpprofiler
  CHARACTER * 200 path_to_wsi_2d_radar
  CHARACTER * 200 path_to_wsi_3d_radar
  CHARACTER * 200 path_to_qc_acars
  CHARACTER * 200 path_to_wisdom
  CHARACTER * 30 fdda_model_source(maxbgmodels)
  CHARACTER * 8 c8_project
  CHARACTER * 8 c8_blpfmt
  CHARACTER * 3 c_raddat_type
  CHARACTER * 80 c80_description
  CHARACTER * 200 path_to_topt30s
  CHARACTER * 200 path_to_topt10m
  CHARACTER * 200 path_to_pctl10m
  CHARACTER * 200 path_to_soil2m
  CHARACTER * 200 path_to_landuse30s
  CHARACTER * 200 path_to_soiltype_top30s
  CHARACTER * 200 path_to_soiltype_bot30s
  CHARACTER * 200 path_to_greenfrac
  CHARACTER * 200 path_to_soiltemp1deg
  CHARACTER * 200 path_to_albedo
  CHARACTER * 200 path_to_maxsnoalb
  CHARACTER * 200 path_to_islope
  CHARACTER * 200 path_to_sst

  LOGICAL * 1 l_compress_radar, l_dual_pol, l_use_tamdar, l_3dvar, l_pad1
  LOGICAL * 1 l_accum_fg, l_accum_radar, l_accum_gauge
  LOGICAL * 1 l_superob_barnes, l_mosaic_sat, l_fsf_gridgen

  COMMON / lapsparms / iflag_lapsparms &
    , max_radar_files_nl, PRESSURE_INTERVAL_L, PRESSURE_0_L &
    , nk_laps, standard_latitude, standard_latitude2 &
    , standard_longitude, NX_L, NY_L, I_PERIMETER &
    , grid_spacing_m, grid_cen_lat, grid_cen_lon &
    , earth_radius &
    , laps_cycle_time &
    , i2_missing_data, r_missing_data, MAX_RADARS &
    , ref_base, ref_base_useable, r_hybrid_first_gate &
    , maxstns, N_PIREP &
    , max_snd_grid, max_snd_levels &
    , redp_lvl, prtop &
    , vert_rad_meso, vert_rad_sao &
    , vert_rad_pirep, vert_rad_prof &
    , silavwt_parm, toptwvl_parm &
    , iwrite_output &
    , aircraft_time_window &
    , vertical_grid, lvl_coord_cdf &
    , c50_lowres_directory, c6_maproj &
    , radarext_3d, radarext_3d_accum &
    , path_to_raw_pirep &
    , path_to_raw_rass, path_to_raw_profiler &
    , path_to_raw_blprass, path_to_raw_blpprofiler &
    , path_to_wsi_2d_radar, path_to_wsi_3d_radar &
    , path_to_qc_acars, path_to_wisdom &
    , c8_project, c8_blpfmt &
    , c_raddat_type, c80_description &
    , path_to_topt30s &
    , path_to_topt10m, path_to_pctl10m, path_to_soil2m &
    , path_to_landuse30s, path_to_soiltype_top30s &
    , path_to_soiltype_bot30s, path_to_greenfrac &
    , path_to_soiltemp1deg, path_to_albedo, path_to_sst &
    , path_to_maxsnoalb, path_to_islope &
    , fdda_model_source &
    , l_compress_radar, l_dual_pol, l_use_tamdar, l_3dvar &
    , l_pad1, l_accum_fg, l_accum_radar, l_accum_gauge

! wind_nl variables
  LOGICAL :: l_use_raob, l_use_cdw, l_use_radial_vel
  REAL    :: weight_bkg_const_wind &
             , weight_radar &
             , rms_thresh_wind &
             , stdev_thresh_radial &
             , r0_barnes_max_m &
             , brns_conv_rate_wind &
             , qc_thresh_wind_def &
             , qc_thresh_wind_pin &
             , qc_thresh_wind_cdw &
             , qc_thresh_wind_pro
  INTEGER :: thresh_2_radarobs_lvl_unfltrd &
             , thresh_4_radarobs_lvl_unfltrd &
             , thresh_9_radarobs_lvl_unfltrd &
             , thresh_25_radarobs_lvl_unfltrd
  INTEGER :: max_pr, max_pr_levels, max_wind_obs

! pressures_nl variables
  INTEGER, PARAMETER   :: max_p = 150
  REAL                 :: pressures(max_p)
  INTEGER              :: nplevs

! surface_analysis variables
  INTEGER  ::  use_lso_qc, skip_internal_qc, itheta
  LOGICAL  ::  l_require_lso, l_dev_ck
  REAL     ::  del, gam, ak &
              , bad_t, bad_td, bad_u, bad_v, bad_p &
              , bad_mp, bad_th, bad_the &
              , bad_tgd_land, bad_tgd_water, bad_vis, bad_tb8 &
              , thresh_t, thresh_td, thresh_mslp &
              , rms_wind, rms_temp, rms_dewpoint, rms_pres

  !  redp_lvl utilized in background code also

! temp_nl variables
  LOGICAL  :: l_read_raob_t, l_use_raob_t
  REAL     :: weight_bkg_const_temp, pres_mix_thresh, rms_thresh_temp, radiometer_ht_temp
  INTEGER  :: max_obs, mode_adjust_heights

! cloud_nl variables
  LOGICAL  :: l_use_vis, l_use_vis_add, l_use_vis_partial, l_use_s8a, l_use_39 &
              , l_use_pireps, l_use_metars, l_use_radar, l_corr_parallax
  INTEGER  :: latency_co2, i4_sat_window, i4_sat_window_offset, i_varadj
  REAL     :: pct_req_lvd_s8a, echotop_thr_a(3), cld_weight_modelfg

! moisture_switch_nl variables
  INTEGER :: covar_switch, print_switch, raob_switch, radiometer_switch, raob_lookback, endian, goes_switch &
             , cloud_switch, cloud_d, sounder_switch, tiros_switch, sat_skip &
             , gvap_switch, ihop_flag, time_diff, gps_switch, sfc_mix &
             , mod_4dda_1
  REAL    :: raob_radius, mod_4dda_factor, t_ref
  REAL    :: max_cdelrh_nl
  REAL    :: cf_set_nl
  REAL    :: cloud_weight_nl
  REAL    :: radio_wt_nl
  CHARACTER(len=256) ::  path_to_gvap12, path_to_gvap10, path_to_gps, path2covar

! lapsprep_nl variables (not yet used)
  LOGICAL ::  hotstart, balance, make_sfc_uv
  CHARACTER(len=16) :: output_format(10)
  REAL ::  snow_thresh, lwc2vapor_thresh, rai_frac, sno_frac
  ! RAMS processing variables
  REAL ::  sfcinf
  CHARACTER(len=256) :: var_prefix

CONTAINS

!-------------------------------------------------------------------

  SUBROUTINE read_namelist_laps(namelist_name, filename)

!use control_coms
!use mem_grid
!use lapsparms_com
!use isan_coms

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: namelist_name, filename

    NAMELIST /lapsparms_NL/ iflag_lapsparms &
      , max_radar_files_nl, PRESSURE_INTERVAL_L &
      , nk_laps, standard_latitude, standard_latitude2 &
      , standard_longitude, NX_L, NY_L, I_PERIMETER &
      , l_compress_radar, l_dual_pol, l_use_tamdar, l_3dvar &
      , l_accum_fg, l_accum_radar, l_accum_gauge &
      , l_superob_barnes, l_mosaic_sat &
      , grid_spacing_m, grid_cen_lat, grid_cen_lon &
      , earth_radius &
      , laps_cycle_time &
      , model_cycle_time, model_fcst_intvl, model_fcst_len &
      , purge_time &
      , i2_missing_data, iverbose, r_missing_data, MAX_RADARS, i_offset_radar &
      , aod, aod_bin, aod_asy, fcterm, aod_ha, ht_ha, ssa, aero_scaleht &
      , ref_base, ref_base_useable, r_hybrid_first_gate &
      , maxstns, N_PIREP &
      , max_snd_grid, max_snd_levels, redp_lvl, prtop &
      , hydrometeor_scale_pcp, hydrometeor_scale_cld &
      , vert_rad_meso, vert_rad_sao &
      , vert_rad_pirep, vert_rad_prof &
      , silavwt_parm, toptwvl_parm &
      , iwrite_output &
      , aircraft_time_window &
      , vertical_grid, lvl_coord_cdf &
      , c50_lowres_directory, c6_maproj &
      , radarext_3d, radarext_3d_accum &
      , path_to_raw_pirep &
      , path_to_raw_rass, path_to_raw_profiler &
      , path_to_raw_blprass, path_to_raw_blpprofiler &
      , path_to_wsi_2d_radar, path_to_wsi_3d_radar &
      , path_to_qc_acars, path_to_wisdom &
      , c8_project, c8_blpfmt &
      , c_raddat_type, c80_description &
      , l_fsf_gridgen &
      , path_to_topt30s, path_to_topt10m &
      , path_to_soiltype_top30s, path_to_soiltype_bot30s &
      , path_to_landuse30s, path_to_greenfrac &
      , path_to_soiltemp1deg, path_to_albedo, path_to_maxsnoalb &
      , path_to_islope, path_to_sst, fdda_model_source

    NAMELIST /pressures_nl/ pressures

    NAMELIST /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel &
      , thresh_2_radarobs_lvl_unfltrd &
      , thresh_4_radarobs_lvl_unfltrd &
      , thresh_9_radarobs_lvl_unfltrd &
      , thresh_25_radarobs_lvl_unfltrd &
      , stdev_thresh_radial &
      , weight_bkg_const_wind &
      , weight_radar &
      , rms_thresh_wind &
      , max_pr, max_pr_levels, max_wind_obs &
      , r0_barnes_max_m &
      , brns_conv_rate_wind &
      , qc_thresh_wind_def &
      , qc_thresh_wind_pin &
      , qc_thresh_wind_cdw &
      , qc_thresh_wind_pro

    NAMELIST /surface_analysis/ &
      use_lso_qc, skip_internal_qc, itheta &
      , l_require_lso, l_dev_ck &
      , del, gam, ak &
      , bad_t, bad_td, bad_u, bad_v, bad_p &
      , bad_mp, bad_th, bad_the &
      , bad_tgd_land, bad_tgd_water, bad_vis, bad_tb8 &
      , thresh_t, thresh_td, thresh_mslp &
      , rms_wind, rms_temp, rms_dewpoint, rms_pres

    NAMELIST /temp_nl/ l_read_raob_t, l_use_raob_t, mode_adjust_heights &
      , weight_bkg_const_temp, pres_mix_thresh, rms_thresh_temp &
      , max_obs, radiometer_ht_temp

    NAMELIST /cloud_nl/ l_use_vis, l_use_vis_add, l_use_vis_partial &
      , l_use_s8a, l_use_pireps, l_use_39, l_use_metars, l_use_radar, l_corr_parallax &
      , latency_co2 &
      , pct_req_lvd_s8a, echotop_thr_a &
      , cld_weight_modelfg &
      , i4_sat_window, i4_sat_window_offset &
      , i_varadj

    NAMELIST /moisture_switch_nl/ &
      covar_switch, print_switch, raob_switch, radiometer_switch, raob_lookback, endian, raob_radius &
      , goes_switch, cloud_switch, cloud_d, sounder_switch &
      , tiros_switch, sat_skip, gvap_switch &
      , max_cdelrh_nl, cf_set_nl, cloud_weight_nl, radio_wt_nl &
      , ihop_flag, time_diff, gps_switch, sfc_mix, mod_4dda_1 &
      , raob_radius, mod_4dda_factor, t_ref &
      , path_to_gvap12, path_to_gvap10, path_to_gps, path2covar

    NAMELIST /lapsprep_nl/ var_prefix, sfcinf, hotstart, balance &
      , output_format, snow_thresh, lwc2vapor_thresh &
      , make_sfc_uv, rai_frac, sno_frac

    PRINT *
    PRINT *, '======> Read_namelist_laps: ', TRIM(namelist_name), nk_laps
    PRINT *, '======>          File_name: ', TRIM(filename)
    PRINT *

! open the namelist file name

    OPEN (12, file=filename, status='old')

! read the requested namelist
    IF (namelist_name == 'ilaps_control') THEN
      ! Set some default values
!  proc_grids = 0

      ! Read ILAPS options information
!  read(12,ilaps_control)

    ELSEIF (namelist_name == 'RAMS') THEN
      ! Read RAMS grid point information
!  open(1,status='OLD',file=rams_grid_input_file)
!  read(1,model_grids)
!  close(1)

      ! Since there is more stuff in the LAPS nest7parm file, read this also.
      !  Then overwrite LAPS grid params with RAMS grid params later.
      !open(1,status='OLD',file=laps_grid_input_file)
      !open(1,status='OLD',file=filename)
      PRINT *, '======> Reading namelist: lapsparms_nl', nk_laps, nx_l, ny_l
      READ (12, lapsparms_nl)
      PRINT *, '======> Reading namelist: lapsparms_nl', nk_laps, nx_l, ny_l
      !close(1)

    ELSEIF (namelist_name == 'LAPS') THEN
      ! Set some default values
!  proc_grids = 0 ; proc_grids(1) = 1

      ! Read LAPS grid point information
      !open(1,status='OLD',file=laps_grid_input_file)
      READ (12, lapsparms_nl)
      !close(1)

    ELSEIF (namelist_name == 'pressures') THEN
      READ (12, pressures_nl)

    ELSEIF (namelist_name == 'lapsparms') THEN

!  default values
      earth_radius = 6370000. ! WRF value
      lvl_coord_cdf = 'HPA'
      prtop = 10000.          ! 100mb for SIGMA_P grid
      l_superob_barnes = .FALSE.
      l_mosaic_sat = .FALSE.
      l_dual_pol = .FALSE.
      l_fsf_gridgen = .FALSE.
      iverbose = 0
      i_offset_radar = -1
      aod = 0.05              ! default column aerosol optical depth
      aero_scaleht = 1500.    ! default aerosol scale height (m)
      fcterm = 0.05           ! range from 0.00 to 0.09 for large aerosols
      ! corresponding phase function peak from 20-110
      aod_ha = .015           ! high altitude aerosol optical depth
      ht_ha(1) = 13000.; ht_ha(2) = 13000.; ht_ha(3) = 25000.; ht_ha(4) = 31000.

!  Single scattering albedo for aerosols
      ssa(1) = .90; ssa(2) = .90; ssa(3) = .90 ! non-dust
!  ssa(1) = .92 ; ssa(2) = .92 ; ssa(3) = .92 ! urban
!  ssa(1) = .99 ; ssa(2) = .99 ; ssa(3) = .99 ! non-absorbing (sea salt)
!  ssa(1) = .95 ; ssa(2) = .85 ; ssa(3) = .75 ! hematite dust
!  ssa(1) = .95 ; ssa(2) = .85 ; ssa(3) = .75 ! smoke

!  fraction of aerosols in each bin & asymmetry factor
!  Factor of 2 back scatter increase from minimum, peak of 50
      aod_bin(1) = 0.000
      aod_bin(2) = 0.987
      aod_bin(3) = 0.013

!  Coarse mode aerosols (rgb)
      aod_asy(1, 1) = +0.957; aod_asy(1, 2) = +0.962; aod_asy(1, 3) = +0.967

!  Fine mode aerosols (rgb)
      aod_asy(2, 1) = +0.44; aod_asy(2, 2) = +0.45; aod_asy(2, 3) = +0.46

!  Backscattering aerosols
      aod_asy(3, :) = +0.55

!  Factor of 10 back scatter increase from minimum
!  aod_bin(1) = 0.70   ;   aod_asy(1) = +0.65
!  aod_bin(2) = 0.12   ;   aod_asy(2) = +0.95
!  aod_bin(3) = 0.18   ;   aod_asy(3) = -0.65

!  Factor of ~3 back scatter increase from minimum
!  aod_bin(1) = 0.80   ;   aod_asy(1) = +0.65
!  aod_bin(2) = 0.14   ;   aod_asy(2) = +0.95
!  aod_bin(3) = 0.06   ;   aod_asy(3) = -0.65

!  Original
!  aod_bin(1) = 0.85   ;   aod_asy(1) = +0.65
!  aod_bin(2) = 0.15   ;   aod_asy(2) = +0.95
!  aod_bin(3) = 0.00   ;   aod_asy(3) = -0.65

      READ (12, lapsparms_nl, err=901)

      ! QC the input variables if desired
      !  .
      !  .
      !  .
      !  .

    ELSEIF (namelist_name == 'wind') THEN

      thresh_25_radarobs_lvl_unfltrd = 450 ! default used mainly for transition

      READ (12, wind_nl, err=905)

      ! QC the input variables if desired
      IF (r0_barnes_max_m .LE. 0) THEN
        WRITE (6, *) ' Error in value of r0_barnes_max_m ', r0_barnes_max_m
        GOTO 905
      END IF

    ELSEIF (namelist_name == 'sfc_anal') THEN

      READ (12, surface_analysis, err=906)

      ! QC the input variables if desired
      !  .
      !  .
      !  .
      !  .

    ELSEIF (namelist_name == 'temp_anal') THEN

      radiometer_ht_temp = 2000.

      READ (12, temp_nl, err=907)

      ! QC the input variables if desired
      !  .
      !  .
      !  .
      !  .

    ELSEIF (namelist_name == 'cloud_anal') THEN

      i_varadj = 1                    ! keep this on as a default
      l_corr_parallax = .TRUE.
      l_use_s8a = .TRUE.
      l_use_pireps = .TRUE.

      READ (12, cloud_nl, err=908)

      ! QC the input variables if desired
      !  .
      !  .
      !  .
      !  .

    ELSEIF (namelist_name == 'moisture_anal') THEN

      READ (12, moisture_switch_nl)

      ! QC the input variables if desired
      !  .
      !  .
      !  .
      !  .

    ELSEIF (namelist_name == 'lapsprep') THEN
      ! Set some default values
      output_format = ' '

      ! Read LAPSPREP info for RAMS-vfile processing
      READ (12, lapsprep_nl)

    ELSE
      PRINT *, 'Illegal namelist_name in read_namelist_laps:', TRIM(namelist_name)
      STOP 'read_namelist_laps: illegal namelist_name'
    END IF

    CLOSE (12)

    RETURN

901 PRINT *, 'error reading lapsparms_nl'
    WRITE (*, lapsparms_nl)
    STOP

905 PRINT *, 'error reading wind_nl'
    WRITE (*, wind_nl)
    STOP

906 PRINT *, 'error reading surface_analysis'
    WRITE (*, surface_analysis)
    STOP

907 PRINT *, 'error reading temp_anal'
    WRITE (*, temp_nl)
    STOP

908 PRINT *, 'error reading cloud_anal'
    WRITE (*, cloud_nl)
    STOP

909 PRINT *, 'error reading moisture_anal'
    WRITE (*, moisture_switch_nl)
    STOP

  END SUBROUTINE

!---------------------------------------------------------------------

END MODULE
