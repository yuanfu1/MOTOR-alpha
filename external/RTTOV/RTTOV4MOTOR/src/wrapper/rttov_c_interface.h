extern void rttov_load_inst_(
    int* inst_id,
    char* opts_str,
    int* nchannels,
    int channels[],
    int l);

extern void rttov_call_direct_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                                  [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                                [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                                  [nprofiles][3]
    int surftype[],             // surftype, watertype                                                  [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5               [nprofiles][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch                          [nprofiles][6]
    double simplecloud[],       // ctp, cfraction                                                       [nprofiles][2]
    int clwscheme[],            // clw_scheme, clwde_param                                              [nprofiles][2]
    int icecloud[],             // ice_scheme, icede_param                                              [nprofiles][2]
    double zeeman[],            // Be, cosbk                                                            [nprofiles][2]
    double p[],                 // pressure                                                             [nprofiles][nlevels]
    double t[],                 // temperature                                                          [nprofiles][nlevels]
    int* gas_units,             // units for gas profiles
    int* mmr_cldaer,            // units for cloud/aerosol profiles
    int gas_id[],               // gas ID list                                                          [ngases]
    double gases[],             // gas profiles                                                         [ngases][nprofiles][nlevels]
    double surfemisrefl[],      // input/output surface emissivities, reflectances, specularities       [4][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)                           [nprofiles][nchannels]
    double rads[],              // output radiances                                                     [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_visir_scatt_call_direct_( //
    int* err,                   //
    int* inst_id,               //
    int channel_list[],         //
    int datetimes[],            // profile dates/times                                                  [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                                [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                                  [nprofiles][3]
    int surftype[],             // surftype, watertype                                                  [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5               [nprofiles][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch                          [nprofiles][6]
    int clwscheme[],            // clw_scheme, clwde_param                                              [nprofiles][2]
    int icecloud[],             // ice_scheme, icede_param                                              [nprofiles][2]
    double p[],                 // pressure                                                             [nprofiles][nlevels]
    double t[],                 // temperature                                                          [nprofiles][nlevels]
    int* gas_units,             // units for gas profiles
    int* mmr_cldaer,            // units for cloud/aerosol profiles
    int gas_id[],               // gas ID list                                                          [ngases]
    double gases[],             // gas profiles                                                         [ngases][nprofiles][nlevels]
    double aer_phangle[],       // aerosol phase fn angles                                              [aer_nphangle]
    double aer_asb[],           // aerosol abs, sca, bpr parameters                                     [3,nprofiles,nchannels,nlayers]
    double aer_legcoef[],       // aerosol phase fn Legendre coefficients                               [nprofiles,nchannels,nlayers,aer_nmom+1]
    double aer_pha[],           // aerosol phase fns                                                    [nprofiles,nchannels,nlayers,aer_nphangle]
    double cld_phangle[],       // cloud phase fn angles                                                [cld_nphangle]
    double cld_asb[],           // cloud abs, sca, bpr parameters                                       [3,nprofiles,nchannels,nlayers]
    double cld_legcoef[],       // cloud phase fn Legendre coefficients                                 [nprofiles,nchannels,nlayers,cld_nmom+1]
    double cld_pha[],           // cloud phase fns                                                      [nprofiles,nchannels,nlayers,cld_nphangle]
    double surfemisrefl[],      // input/output surface emissivities, reflectances, specularities       [4][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)                           [nprofiles][nchannels]
    double rads[],              // output radiances                                                     [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles,
    int* aer_nphangle, int* aer_nmom, int* cld_nphangle, int* cld_nmom);

extern void rttov_scatt_call_direct_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                     [nprofiles][6]
    double angles[],            // satzen, satazi angles                                   [nprofiles][2]
    double surfgeom[],          // lat, lon, elevation                                     [nprofiles][3]
    int surftype[],             // surftype                                                [nprofiles]
    double skin[],              // skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
    double zeeman[],            // Be, cosbk                                               [nprofiles][2]
    double p[],                 // pressure                                                [nprofiles][nlevels]
    double t[],                 // temperature                                             [nprofiles][nlevels]
    int* gas_units,             // units for gas profiles
    int gas_id[],               // gas ID list                                             [ngases]
    double gases[],             // gas profiles                                            [ngases][nprofiles][nlevels]
    double ph[],                // pressure half-levels                                    [nprofiles][nlevels+1]
    double cfrac[],             // user cloud fraction                                     [nprofiles]
    int* multi_hydro_frac,      // false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
    int* calc_zef,              // enable/disable radar reflectivity calculations
    double surfemis[],          // input/output surface emissivities                       [nprofiles][nchannels]
    double bt[],                // output BTs                                              [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_call_k_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                                  [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                                [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                                  [nprofiles][3]
    int surftype[],             // surftype, watertype                                                  [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5               [nprofiles][9]
    double skin_k[],            // output skin K                                                        [nprofiles][nchannels][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch                          [nprofiles][6]
    double s2m_k[],             // output 2m K                                                          [nprofiles][nchannels][6]
    double simplecloud[],       // ctp, cfraction                                                       [nprofiles][2]
    double simplecloud_k[],     // output ctp, cfraction K                                              [nprofiles][nchannels][2]
    int clwscheme[],            // clw_scheme, clwde_param                                              [nprofiles][2]
    int icecloud[],             // ice_scheme, icede_param                                              [nprofiles][2]
    double zeeman[],            // Be, cosbk                                                            [nprofiles][2]
    double p[],                 // pressure                                                             [nprofiles][nlevels]
    double p_k[],               // output pressure K                                                    [nprofiles][nchannels][nlevels]
    double t[],                 // temperature                                                          [nprofiles][nlevels]
    double t_k[],               // output temperature K                                                 [nprofiles][nchannels][nlevels]
    int* gas_units,             // units for gas profiles
    int* mmr_cldaer,            // units for cloud/aerosol profiles
    int gas_id[],               // gas ID list                                                          [ngases]
    double gases[],             // gas profiles                                                         [ngases][nprofiles][nlevels]
    double gases_k[],           // output gas profiles K                                                [ngases][nprofiles][nchannels][nlevels]
    double surfemisrefl[],      // input/output surface emissivities, reflectances, specularities       [4][nprofiles][nchannels]
    double surfemisrefl_k[],    // output surface emis, refl, spec K                                    [4][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)                           [nprofiles][nchannels]
    double rads[],              // output radiances                                                     [nprofiles][nchannels]
    double btrefl_k[],          // input BT perturbations                                               [nprofiles][nchannels]
    double rads_k[],            // input radiance perturbations                                         [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_visir_scatt_call_k_( //
    int* err,                   //
    int* inst_id,               //
    int channel_list[],         //
    int datetimes[],            // profile dates/times                                                  [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                                [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                                  [nprofiles][3]
    int surftype[],             // surftype, watertype                                                  [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5               [nprofiles][9]
    double skin_k[],            // output skin K                                                        [nprofiles][nchannels][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch                          [nprofiles][6]
    double s2m_k[],             // output 2m K                                                          [nprofiles][nchannels][6]
    int clwscheme[],            // clw_scheme, clwde_param                                              [nprofiles][2]
    int icecloud[],             // ice_scheme, icede_param                                              [nprofiles][2]
    double p[],                 // pressure                                                             [nprofiles][nlevels]
    double p_k[],               // output pressure K                                                    [nprofiles][nchannels][nlevels]
    double t[],                 // temperature                                                          [nprofiles][nlevels]
    double t_k[],               // output temperature K                                                 [nprofiles][nchannels][nlevels]
    int* gas_units,             // units for gas profiles
    int* mmr_cldaer,            // units for cloud/aerosol profiles
    int gas_id[],               // gas ID list                                                          [ngases]
    double gases[],             // gas profiles                                                         [ngases][nprofiles][nlevels]
    double gases_k[],           // output gas profiles K                                                [ngases][nprofiles][nchannels][nlevels]
    double aer_phangle[],       // aerosol phase fn angles                                              [aer_nphangle]
    double aer_asb[],           // aerosol abs, sca, bpr parameters                                     [3,nprofiles,nchannels,nlayers]
    double aer_legcoef[],       // aerosol phase fn Legendre coefficients                               [nprofiles,nchannels,nlayers,aer_nmom+1]
    double aer_pha[],           // aerosol phase fns                                                    [nprofiles,nchannels,nlayers,aer_nphangle]
    double cld_phangle[],       // cloud phase fn angles                                                [cld_nphangle]
    double cld_asb[],           // cloud abs, sca, bpr parameters                                       [3,nprofiles,nchannels,nlayers]
    double cld_legcoef[],       // cloud phase fn Legendre coefficients                                 [nprofiles,nchannels,nlayers,cld_nmom+1]
    double cld_pha[],           // cloud phase fns                                                      [nprofiles,nchannels,nlayers,cld_nphangle]
    double surfemisrefl[],      // input/output surface emissivities, reflectances, specularities       [4][nprofiles][nchannels]
    double surfemisrefl_k[],    // output surface emis, refl, spec K                                    [4][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)                           [nprofiles][nchannels]
    double rads[],              // output radiances                                                     [nprofiles][nchannels]
    double btrefl_k[],          // input BT perturbations                                               [nprofiles][nchannels]
    double rads_k[],            // input radiance perturbations                                         [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles,
    int* aer_nphangle, int* aer_nmom, int* cld_nphangle, int* cld_nmom);

extern void rttov_scatt_call_k_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                     [nprofiles][6]
    double angles[],            // satzen, satazi angles                                   [nprofiles][2]
    double surfgeom[],          // lat, lon, elevation                                     [nprofiles][3]
    int surftype[],             // surftype                                                [nprofiles]
    double skin[],              // skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
    double skin_k[],            // output skin K                                           [nprofiles][8]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
    double s2m_k[],             // output 2m K                                             [nprofiles][5]
    double zeeman[],            // Be, cosbk                                               [nprofiles][2]
    double p[],                 // pressure                                                [nprofiles][nlevels]
    double p_k[],               // output pressure K                                       [nprofiles][nchannels][nlevels]
    double t[],                 // temperature                                             [nprofiles][nlevels]
    double t_k[],               // output temperature K                                    [nprofiles][nchannels][nlevels]
    int* gas_units,             // units for gas profiles
    int gas_id[],               // gas ID list                                             [ngases]
    double gases[],             // gas profiles                                            [ngases][nprofiles][nlevels]
    double gases_k[],           // output gas profiles K                                   [ngases][nprofiles][nchannels][nlevels]
    double ph[],                // pressure half-levels                                    [nprofiles][nlevels+1]
    double ph_k[],              // output pressure half-levels K                           [nprofiles][nchannels][nlevels+1]
    double cfrac[],             // user cloud fraction                                     [nprofiles]
    double cfrac_k[],           // user cloud fraction K                                   [nprofiles][nchannels]
    int* multi_hydro_frac,      // false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
    int* calc_zef,              // enable/disable radar reflectivity calculations
    double surfemis[],          // input/output surface emissivities                       [nprofiles][nchannels]
    double surfemis_k[],        // input/output surface emissivities K                     [nprofiles][nchannels]
    double bt[],                // output BTs                                              [nprofiles][nchannels]
    double bt_k[],              // input BT perturbations                                  [nprofiles][nchannels]
    double zef_k[],             // input radar reflectivity perturbations                  [nprofiles][nchannels][nlevels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_drop_inst_(int* err, int* inst_id);

extern void rttov_drop_all_(int* err);

extern void rttov_set_options_(int* err, int* inst_id, char* opts_str, int l);

extern void rttov_print_options_(int* err, int* inst_id);


extern void rttov_load_ir_emis_atlas_(
    int* atlas_wrap_id,
    char* path,
    int* month,
    int* atlas_id,
    int* inst_id,
    int* ang_corr,
    int l);

extern void rttov_load_mw_emis_atlas_(
    int* atlas_wrap_id,
    char* path,
    int* month,
    int* atlas_id,
    int* inst_id,
    int* year,
    int l);

extern void rttov_load_brdf_atlas_(
    int* atlas_wrap_id,
    char* path,
    int* month,
    int* atlas_id,
    int* inst_id,
    int l);

extern void rttov_get_emisbrdf_(
    int* err,
    int* atlas_wrap_id,
    double latitude[],         // [nprofiles]
    double longitude[],        // [nprofiles]
    int surftype[],            // [nprofiles]
    int watertype[],           // [nprofiles]
    double zenangle[],         // [nprofiles]
    double azangle[],          // [nprofiles]
    double sunzenangle[],      // [nprofiles]
    double sunazangle[],       // [nprofiles]
    double snow_fraction[],    // [nprofiles]
    int* inst_id,
    int channel_list[],        // [nchannels]
    double emisbrdf[],         // [nprofiles][nchannels]
    int* nchannels,
    int* nprofiles);

extern void rttov_drop_atlas_(int* err, int* atlas_wrap_id);


extern void rttov_bpr_(int* err,
                       double phangle[],     // [nphangle]
                       double pha[],         // [nphangle]
                       double* bpr,
                       int* nthreads,
                       int* nphangle);

extern void rttov_legcoef_(int* err,
                           double phangle[], // [nphangle]
                           double pha[],     // [nphangle]
                           double legcoef[], // [nmom]
                           int* ngauss,
                           int* nphangle,
                           int* nmom);


extern void rttov_get_rad_clear_(int* err, int* inst_id, double rad_clear[], int* nchanprof);

extern void rttov_get_rad_total_(int* err, int* inst_id, double rad_total[], int* nchanprof);

extern void rttov_get_bt_clear_(int* err, int* inst_id, double bt_clear[], int* nchanprof);

extern void rttov_get_bt_(int* err, int* inst_id, double bt[], int* nchanprof);

extern void rttov_get_refl_clear_(int* err, int* inst_id, double refl_clear[], int* nchanprof);

extern void rttov_get_refl_(int* err, int* inst_id, double refl[], int* nchanprof);

extern void rttov_get_rad_cloudy_(int* err, int* inst_id, double rad_cloudy[], int* nchanprof);

extern void rttov_get_overcast_(int* err, int* inst_id, double overcast[], int* nchanprof, int* nlayers);

extern void rttov_get_rad_quality_(int* err, int* inst_id, int rad_quality[], int* nchanprof);

extern void rttov_get_plane_parallel_(int* err, int* inst_id, int* plane_parallel);

extern void rttov_get_geometric_height_(int* err, int* inst_id, double geometric_height[], int* nchanprof, int* nlevels);


extern void rttov_get_rad2_upclear_(int* err, int* inst_id, double rad2_upclear[], int* nchanprof);

extern void rttov_get_rad2_dnclear_(int* err, int* inst_id, double rad2_dnclear[], int* nchanprof);

extern void rttov_get_rad2_refldnclear_(int* err, int* inst_id, double rad2_refldnclear[], int* nchanprof);

extern void rttov_get_rad2_up_(int* err, int* inst_id, double rad2_up[], int* nchanprof, int* nlayers);

extern void rttov_get_rad2_down_(int* err, int* inst_id, double rad2_down[], int* nchanprof, int* nlayers);

extern void rttov_get_rad2_surf_(int* err, int* inst_id, double rad2_surf[], int* nchanprof, int* nlayers);


extern void rttov_get_tau_total_(int* err, int* inst_id, double tau_total[], int* nchanprof);

extern void rttov_get_tau_levels_(int* err, int* inst_id, double tau_levels[], int* nchanprof, int* nlevels);

extern void rttov_get_tausun_total_path2_(int* err, int* inst_id, double tau_total_path2[], int* nchanprof);

extern void rttov_get_tausun_levels_path2_(int* err, int* inst_id, double tau_levels_path2[], int* nchanprof, int* nlevels);

extern void rttov_get_tausun_total_path1_(int* err, int* inst_id, double tau_total_path1[], int* nchanprof);

extern void rttov_get_tausun_levels_path1_(int* err, int* inst_id, double tau_levels_path1[], int* nchanprof, int* nlevels);

extern void rttov_get_tau_total_cld_(int* err, int* inst_id, double tau_total_cld[], int* nchanprof);

extern void rttov_get_tau_levels_cld_(int* err, int* inst_id, double tau_levels_cld[], int* nchanprof, int* nlevels);


extern void rttov_get_emis_terms_cfrac_(int* err, int* inst_id, double cfrac[], int* nchanprof);

extern void rttov_get_emis_terms_bsfc_(int* err, int* inst_id, double bsfc[], int* nchanprof);

extern void rttov_get_emis_terms_tau_cld_(int* err, int* inst_id, double tau_cld[], int* nchanprof);

extern void rttov_get_emis_terms_up_cld_(int* err, int* inst_id, double up_cld[], int* nchanprof);

extern void rttov_get_emis_terms_down_cld_(int* err, int* inst_id, double down_cld[], int* nchanprof);

extern void rttov_get_emis_terms_tau_clr_(int* err, int* inst_id, double tau_clr[], int* nchanprof);

extern void rttov_get_emis_terms_up_clr_(int* err, int* inst_id, double up_clr[], int* nchanprof);

extern void rttov_get_emis_terms_down_clr_(int* err, int* inst_id, double down_clr[], int* nchanprof);


extern void rttov_get_zef_(int* err, int* inst_id, double zef[], int* nchanprof, int* nlevels);

extern void rttov_get_azef_(int* err, int* inst_id, double azef[], int* nchanprof, int* nlevels);


// The following are used by the Rttov and RttovSafe classes to interrogate the RTTOV coefficients structure
extern void rttov_get_coef_val_i0_(int* err, int* inst_id, char* varch, int* i0, int l);

extern void rttov_get_coef_val_r0_(int* err, int* inst_id, char* varch, double* r0, int l);

extern void rttov_get_coef_val_c0_(int* err, int* inst_id, char* varch, char* c0, int l);

extern void rttov_get_coef_val_i1_(int* err, int* inst_id, char* varch, int* m, int* i1, int l);

extern void rttov_get_coef_val_i2_(int* err, int* inst_id, char* varch, int* m, int* n, int* i2, int l);

extern void rttov_get_coef_val_r1_(int* err, int* inst_id, char* varch, int* m, double* r1, int l);




