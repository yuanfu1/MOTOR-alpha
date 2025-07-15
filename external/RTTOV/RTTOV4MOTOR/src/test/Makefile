
TOPDIR=../..
SUBDIR=test

FFLAGS_ARCH_example_aer_file_fwd               = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_aer_param_fwd              = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_atlas_fwd                  = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_cld_file_fwd               = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_cld_mfasis_fwd             = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_cld_param_fwd              = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_fwd                        = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_htfrtc_fwd                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_k                          = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_pc_fwd                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_example_rttovscatt_fwd             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_chain                        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_k_ad                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_k_bf                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_k_tl                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_lun                          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_cld_profile_inc         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_emisrefl_inc            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_opt_param               = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_opt_param_inc           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_pccomp_inc              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_profile_inc             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_radiance_inc            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_make_reflectivity_inc        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_cld_profile_inc        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_opt_param_inc          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_pccomp_inc             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_profile_inc            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_radiance_inc           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scale_reflectivity_inc       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_taylor_test                  = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_test                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_test_k_mod                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_test_mod                     = $(FFLAGS_ARCH)


F90_INTF=rttov_k_ad.F90 rttov_k_bf.F90 rttov_k_tl.F90 rttov_make_cld_profile_inc.F90 rttov_make_emisrefl_inc.F90 rttov_make_opt_param_inc.F90 rttov_make_pccomp_inc.F90 rttov_make_profile_inc.F90 rttov_make_radiance_inc.F90 rttov_make_reflectivity_inc.F90 rttov_scale_cld_profile_inc.F90 rttov_scale_opt_param_inc.F90 rttov_scale_pccomp_inc.F90 rttov_scale_profile_inc.F90 rttov_scale_radiance_inc.F90 rttov_scale_reflectivity_inc.F90 rttov_taylor_test.F90

INTF=$(INCDIR)/rttov_k_ad.interface $(INCDIR)/rttov_k_bf.interface $(INCDIR)/rttov_k_tl.interface $(INCDIR)/rttov_make_cld_profile_inc.interface $(INCDIR)/rttov_make_emisrefl_inc.interface $(INCDIR)/rttov_make_opt_param_inc.interface $(INCDIR)/rttov_make_pccomp_inc.interface $(INCDIR)/rttov_make_profile_inc.interface $(INCDIR)/rttov_make_radiance_inc.interface $(INCDIR)/rttov_make_reflectivity_inc.interface $(INCDIR)/rttov_scale_cld_profile_inc.interface $(INCDIR)/rttov_scale_opt_param_inc.interface $(INCDIR)/rttov_scale_pccomp_inc.interface $(INCDIR)/rttov_scale_profile_inc.interface $(INCDIR)/rttov_scale_radiance_inc.interface $(INCDIR)/rttov_scale_reflectivity_inc.interface $(INCDIR)/rttov_taylor_test.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_k_ad$(EXT_OBJ) $(OBJDIR)/rttov_k_bf$(EXT_OBJ) $(OBJDIR)/rttov_k_tl$(EXT_OBJ) $(OBJDIR)/rttov_make_cld_profile_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_emisrefl_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_opt_param_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_pccomp_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_profile_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_radiance_inc$(EXT_OBJ) $(OBJDIR)/rttov_make_reflectivity_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_cld_profile_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_opt_param_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_pccomp_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_profile_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_radiance_inc$(EXT_OBJ) $(OBJDIR)/rttov_scale_reflectivity_inc$(EXT_OBJ) $(OBJDIR)/rttov_taylor_test$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/rttov_chain$(EXT_OBJ) $(OBJDIR)/rttov_lun$(EXT_OBJ) $(OBJDIR)/rttov_test_k_mod$(EXT_OBJ) $(OBJDIR)/rttov_test_mod$(EXT_OBJ)

OBJBIN=$(OBJDIR)/example_aer_file_fwd$(EXT_OBJ) $(OBJDIR)/example_aer_param_fwd$(EXT_OBJ) $(OBJDIR)/example_atlas_fwd$(EXT_OBJ) $(OBJDIR)/example_cld_file_fwd$(EXT_OBJ) $(OBJDIR)/example_cld_mfasis_fwd$(EXT_OBJ) $(OBJDIR)/example_cld_param_fwd$(EXT_OBJ) $(OBJDIR)/example_fwd$(EXT_OBJ) $(OBJDIR)/example_htfrtc_fwd$(EXT_OBJ) $(OBJDIR)/example_k$(EXT_OBJ) $(OBJDIR)/example_pc_fwd$(EXT_OBJ) $(OBJDIR)/example_rttovscatt_fwd$(EXT_OBJ) $(OBJDIR)/rttov_make_opt_param$(EXT_OBJ) $(OBJDIR)/rttov_test$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_test$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/example_aer_file_fwd$(EXT_BIN) $(BINDIR)/example_aer_param_fwd$(EXT_BIN) $(BINDIR)/example_atlas_fwd$(EXT_BIN) $(BINDIR)/example_cld_file_fwd$(EXT_BIN) $(BINDIR)/example_cld_mfasis_fwd$(EXT_BIN) $(BINDIR)/example_cld_param_fwd$(EXT_BIN) $(BINDIR)/example_fwd$(EXT_BIN) $(BINDIR)/example_htfrtc_fwd$(EXT_BIN) $(BINDIR)/example_k$(EXT_BIN) $(BINDIR)/example_pc_fwd$(EXT_BIN) $(BINDIR)/example_rttovscatt_fwd$(EXT_BIN) $(BINDIR)/rttov_make_opt_param$(EXT_BIN) $(BINDIR)/rttov_test$(EXT_BIN)

MODS=$(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_test_mod$(EXT_MOD)

include $(TOPDIR)/build/Makefile.inc

lib: dir intf $(LIB)

bin: lib $(BINS)

dir: 
	@mkdir -p $(TMPDIR)/$(SUBDIR) $(BINDIR) $(OBJDIR) $(INCDIR) $(MODDIR) $(LIBDIR)

intf: 
	@$(MKINTF) $(F90_INTF) $(INTF)
	@$(CPINCH) $(INCDIR)

$(LIB): $(OBJSUB) $(OBJMOD)
	\rm -f $(LIB) $(LIBNOVER)
	$(AR) $(LIB) $(OBJSUB) $(OBJMOD)
	\ln -s $(LIBNAME) $(LIBNOVER)

clean:
	\rm -f $(INTF) $(LIB) $(LIBNOVER) $(BINS) $(OBJSUB) $(OBJMOD) $(OBJBIN) $(MODS)


$(BINDIR)/example_aer_file_fwd$(EXT_BIN): $(OBJDIR)/example_aer_file_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_aer_file_fwd$(EXT_BIN) $(OBJDIR)/example_aer_file_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_aer_param_fwd$(EXT_BIN): $(OBJDIR)/example_aer_param_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_aer_param_fwd$(EXT_BIN) $(OBJDIR)/example_aer_param_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_atlas_fwd$(EXT_BIN): $(OBJDIR)/example_atlas_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_atlas_fwd$(EXT_BIN) $(OBJDIR)/example_atlas_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_cld_file_fwd$(EXT_BIN): $(OBJDIR)/example_cld_file_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_cld_file_fwd$(EXT_BIN) $(OBJDIR)/example_cld_file_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_cld_mfasis_fwd$(EXT_BIN): $(OBJDIR)/example_cld_mfasis_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_cld_mfasis_fwd$(EXT_BIN) $(OBJDIR)/example_cld_mfasis_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_cld_param_fwd$(EXT_BIN): $(OBJDIR)/example_cld_param_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_cld_param_fwd$(EXT_BIN) $(OBJDIR)/example_cld_param_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_fwd$(EXT_BIN): $(OBJDIR)/example_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_fwd$(EXT_BIN) $(OBJDIR)/example_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_htfrtc_fwd$(EXT_BIN): $(OBJDIR)/example_htfrtc_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_htfrtc_fwd$(EXT_BIN) $(OBJDIR)/example_htfrtc_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_k$(EXT_BIN): $(OBJDIR)/example_k$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_k$(EXT_BIN) $(OBJDIR)/example_k$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_pc_fwd$(EXT_BIN): $(OBJDIR)/example_pc_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_pc_fwd$(EXT_BIN) $(OBJDIR)/example_pc_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/example_rttovscatt_fwd$(EXT_BIN): $(OBJDIR)/example_rttovscatt_fwd$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/example_rttovscatt_fwd$(EXT_BIN) $(OBJDIR)/example_rttovscatt_fwd$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_make_opt_param$(EXT_BIN): $(OBJDIR)/rttov_make_opt_param$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_make_opt_param$(EXT_BIN) $(OBJDIR)/rttov_make_opt_param$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_test$(EXT_BIN): $(OBJDIR)/rttov_test$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_test$(EXT_BIN) $(OBJDIR)/rttov_test$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_test -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov_k_ad$(EXT_OBJ): rttov_k_ad.F90 $(INCDIR)/rttov_ad.interface $(INCDIR)/rttov_parallel_ad.interface $(INCDIR)/rttov_scatt_ad.interface $(INCDIR)/rttov_parallel_scatt_ad.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_reflectivity.interface $(INCDIR)/rttov_alloc_pccomp.interface $(INCDIR)/rttov_alloc_traj.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_init_scatt_prof.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_init_pccomp.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_transmission.interface $(INCDIR)/rttov_init_reflectivity.interface $(INCDIR)/rttov_copy_prof.interface $(INCDIR)/rttov_copy_scatt_prof.interface $(INCDIR)/rttov_copy_opt_param.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_k_ad) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_k_ad.F90 -o $(OBJDIR)/rttov_k_ad$(EXT_OBJ)


$(OBJDIR)/rttov_k_bf$(EXT_OBJ): rttov_k_bf.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_parallel_scatt.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_alloc_pccomp.interface $(INCDIR)/rttov_alloc_traj.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_reflectivity.interface $(INCDIR)/rttov_alloc_pc_dimensions.interface $(INCDIR)/rttov_copy_prof.interface $(INCDIR)/rttov_copy_scatt_prof.interface $(INCDIR)/rttov_copy_opt_param.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_make_profile_inc.interface $(INCDIR)/rttov_make_cld_profile_inc.interface $(INCDIR)/rttov_make_opt_param_inc.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_k_bf) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_k_bf.F90 -o $(OBJDIR)/rttov_k_bf$(EXT_OBJ)


$(OBJDIR)/rttov_k_tl$(EXT_OBJ): rttov_k_tl.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_tl.interface $(INCDIR)/rttov_parallel_tl.interface $(INCDIR)/rttov_scatt_tl.interface $(INCDIR)/rttov_parallel_scatt_tl.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_alloc_pccomp.interface $(INCDIR)/rttov_alloc_traj.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_reflectivity.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_reflectivity.interface $(INCDIR)/rttov_init_pccomp.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_k_tl) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_k_tl.F90 -o $(OBJDIR)/rttov_k_tl$(EXT_OBJ)


$(OBJDIR)/rttov_make_cld_profile_inc$(EXT_OBJ): rttov_make_cld_profile_inc.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_cld_profile_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_cld_profile_inc.F90 -o $(OBJDIR)/rttov_make_cld_profile_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_emisrefl_inc$(EXT_OBJ): rttov_make_emisrefl_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_emisrefl_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_emisrefl_inc.F90 -o $(OBJDIR)/rttov_make_emisrefl_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_opt_param_inc$(EXT_OBJ): rttov_make_opt_param_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_opt_param_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_opt_param_inc.F90 -o $(OBJDIR)/rttov_make_opt_param_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_pccomp_inc$(EXT_OBJ): rttov_make_pccomp_inc.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_pccomp_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_pccomp_inc.F90 -o $(OBJDIR)/rttov_make_pccomp_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_profile_inc$(EXT_OBJ): rttov_make_profile_inc.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_profile_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_profile_inc.F90 -o $(OBJDIR)/rttov_make_profile_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_radiance_inc$(EXT_OBJ): rttov_make_radiance_inc.F90  $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_math_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_radiance_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_radiance_inc.F90 -o $(OBJDIR)/rttov_make_radiance_inc$(EXT_OBJ)


$(OBJDIR)/rttov_make_reflectivity_inc$(EXT_OBJ): rttov_make_reflectivity_inc.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_reflectivity_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_reflectivity_inc.F90 -o $(OBJDIR)/rttov_make_reflectivity_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_cld_profile_inc$(EXT_OBJ): rttov_scale_cld_profile_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_cld_profile_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_cld_profile_inc.F90 -o $(OBJDIR)/rttov_scale_cld_profile_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_opt_param_inc$(EXT_OBJ): rttov_scale_opt_param_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_opt_param_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_opt_param_inc.F90 -o $(OBJDIR)/rttov_scale_opt_param_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_pccomp_inc$(EXT_OBJ): rttov_scale_pccomp_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_pccomp_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_pccomp_inc.F90 -o $(OBJDIR)/rttov_scale_pccomp_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_profile_inc$(EXT_OBJ): rttov_scale_profile_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_profile_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_profile_inc.F90 -o $(OBJDIR)/rttov_scale_profile_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_radiance_inc$(EXT_OBJ): rttov_scale_radiance_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_radiance_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_radiance_inc.F90 -o $(OBJDIR)/rttov_scale_radiance_inc$(EXT_OBJ)


$(OBJDIR)/rttov_scale_reflectivity_inc$(EXT_OBJ): rttov_scale_reflectivity_inc.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scale_reflectivity_inc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scale_reflectivity_inc.F90 -o $(OBJDIR)/rttov_scale_reflectivity_inc$(EXT_OBJ)


$(OBJDIR)/rttov_taylor_test$(EXT_OBJ): rttov_taylor_test.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_scale_profile_inc.interface $(INCDIR)/rttov_scale_cld_profile_inc.interface $(INCDIR)/rttov_scale_opt_param_inc.interface $(INCDIR)/rttov_make_profile_inc.interface $(INCDIR)/rttov_make_cld_profile_inc.interface $(INCDIR)/rttov_make_opt_param_inc.interface $(INCDIR)/rttov_make_emisrefl_inc.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_alloc_tl.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_copy_prof.interface $(INCDIR)/rttov_copy_scatt_prof.interface $(INCDIR)/rttov_copy_opt_param.interface $(INCDIR)/rttov_add_prof.interface $(INCDIR)/rttov_add_scatt_prof.interface $(INCDIR)/rttov_add_opt_param.interface $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_tl.interface $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_scatt_tl.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_taylor_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_taylor_test.F90 -o $(OBJDIR)/rttov_taylor_test$(EXT_OBJ)


$(OBJDIR)/example_aer_file_fwd$(EXT_OBJ): example_aer_file_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_aer_clim_prof.interface $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_aer_file_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_aer_file_fwd.F90 -o $(OBJDIR)/example_aer_file_fwd$(EXT_OBJ)


$(OBJDIR)/example_aer_param_fwd$(EXT_OBJ): example_aer_param_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_legcoef_calc.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_aer_param_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_aer_param_fwd.F90 -o $(OBJDIR)/example_aer_param_fwd$(EXT_OBJ)


$(OBJDIR)/example_atlas_fwd$(EXT_OBJ): example_atlas_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_setup_brdf_atlas.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_camel_clim_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_atlas_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_atlas_fwd.F90 -o $(OBJDIR)/example_atlas_fwd$(EXT_OBJ)


$(OBJDIR)/example_cld_file_fwd$(EXT_OBJ): example_cld_file_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_cld_file_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_cld_file_fwd.F90 -o $(OBJDIR)/example_cld_file_fwd$(EXT_OBJ)


$(OBJDIR)/example_cld_mfasis_fwd$(EXT_OBJ): example_cld_mfasis_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_print_radiance_quality.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_cld_mfasis_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_cld_mfasis_fwd.F90 -o $(OBJDIR)/example_cld_mfasis_fwd$(EXT_OBJ)


$(OBJDIR)/example_cld_param_fwd$(EXT_OBJ): example_cld_param_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_legcoef_calc.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_cld_param_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_cld_param_fwd.F90 -o $(OBJDIR)/example_cld_param_fwd$(EXT_OBJ)


$(OBJDIR)/example_fwd$(EXT_OBJ): example_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_fwd.F90 -o $(OBJDIR)/example_fwd$(EXT_OBJ)


$(OBJDIR)/example_htfrtc_fwd$(EXT_OBJ): example_htfrtc_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs_htfrtc.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_camel_clim_atlas$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_htfrtc_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_htfrtc_fwd.F90 -o $(OBJDIR)/example_htfrtc_fwd$(EXT_OBJ)


$(OBJDIR)/example_k$(EXT_OBJ): example_k.F90 $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_k.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_transmission.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_k) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_k.F90 -o $(OBJDIR)/example_k$(EXT_OBJ)


$(OBJDIR)/example_pc_fwd$(EXT_OBJ): example_pc_fwd.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_get_pc_predictindex.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_pc_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_pc_fwd.F90 -o $(OBJDIR)/example_pc_fwd$(EXT_OBJ)


$(OBJDIR)/example_rttovscatt_fwd$(EXT_OBJ): example_rttovscatt_fwd.F90 $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_parallel_scatt.interface $(INCDIR)/rttov_read_scattcoeffs.interface $(INCDIR)/rttov_dealloc_scattcoeffs.interface $(INCDIR)/rttov_scatt_setupindex.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_print_opts_scatt.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_print_cld_profile.interface $(INCDIR)/rttov_skipcommentline.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_example_rttovscatt_fwd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/example_rttovscatt_fwd.F90 -o $(OBJDIR)/example_rttovscatt_fwd$(EXT_OBJ)


$(OBJDIR)/rttov_make_opt_param$(EXT_OBJ): rttov_make_opt_param.F90 $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_aux_prof.interface $(INCDIR)/rttov_baran2014_calc_optpar.interface $(INCDIR)/rttov_baran2018_calc_optpar.interface $(INCDIR)/rttov_baran_calc_phase.interface $(INCDIR)/rttov_convert_profile_units.interface $(INCDIR)/rttov_profaux_cldaer.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_rttov_baran2018_icldata$(EXT_MOD) $(MODDIR)/mod_rttov_baran2014_icldata$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_make_opt_param) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_make_opt_param.F90 -o $(OBJDIR)/rttov_make_opt_param$(EXT_OBJ)


$(OBJDIR)/rttov_test$(EXT_OBJ): rttov_test.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_tl.interface $(INCDIR)/rttov_ad.interface $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_k_tl.interface $(INCDIR)/rttov_k_ad.interface $(INCDIR)/rttov_k_bf.interface $(INCDIR)/rttov_taylor_test.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_parallel_tl.interface $(INCDIR)/rttov_parallel_ad.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_scatt_tl.interface $(INCDIR)/rttov_scatt_ad.interface $(INCDIR)/rttov_parallel_scatt.interface $(INCDIR)/rttov_parallel_scatt_tl.interface $(INCDIR)/rttov_parallel_scatt_ad.interface $(INCDIR)/rttov_make_profile_inc.interface $(INCDIR)/rttov_scale_profile_inc.interface $(INCDIR)/rttov_make_cld_profile_inc.interface $(INCDIR)/rttov_scale_cld_profile_inc.interface $(INCDIR)/rttov_make_opt_param_inc.interface $(INCDIR)/rttov_scale_opt_param_inc.interface $(INCDIR)/rttov_make_radiance_inc.interface $(INCDIR)/rttov_scale_radiance_inc.interface $(INCDIR)/rttov_make_pccomp_inc.interface $(INCDIR)/rttov_scale_pccomp_inc.interface $(INCDIR)/rttov_make_emisrefl_inc.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_read_coefs_htfrtc.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_read_scattcoeffs.interface $(INCDIR)/rttov_dealloc_scattcoeffs.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_reflectivity.interface $(INCDIR)/rttov_alloc_emis_ret_terms.interface $(INCDIR)/rttov_alloc_pccomp.interface $(INCDIR)/rttov_init_emis_refl.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_init_scatt_prof.interface $(INCDIR)/rttov_copy_prof.interface $(INCDIR)/rttov_copy_scatt_prof.interface $(INCDIR)/rttov_copy_pccomp.interface $(INCDIR)/rttov_alloc_traj.interface $(INCDIR)/rttov_alloc_opt_param.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_scatt_setupindex.interface $(INCDIR)/rttov_scale_ref_gas_prof.interface $(INCDIR)/rttov_print_radiance_quality.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_opts_scatt.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_print_cld_profile.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_user_profile_checkinput.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/mod_rttov_baran2014_icldata$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_rttov_baran2018_icldata$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_test_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_math_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_test.F90 -o $(OBJDIR)/rttov_test$(EXT_OBJ)


$(OBJDIR)/rttov_chain$(EXT_OBJ): rttov_chain.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_chain) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_chain.F90 -o $(OBJDIR)/rttov_chain$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_chain$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD)

$(MODDIR)/rttov_chain$(EXT_MOD): rttov_chain.F90  $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_chain) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_chain.F90 -o $(OBJDIR)/rttov_chain$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_chain$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD)


$(OBJDIR)/rttov_lun$(EXT_OBJ): rttov_lun.F90  $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_lun) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_lun.F90 -o $(OBJDIR)/rttov_lun$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_lun$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD)

$(MODDIR)/rttov_lun$(EXT_MOD): rttov_lun.F90  $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_lun) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_lun.F90 -o $(OBJDIR)/rttov_lun$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_lun$(EXT_MOD) $(MODDIR)/rttov_lun$(EXT_MOD)


$(OBJDIR)/rttov_test_k_mod$(EXT_OBJ): rttov_test_k_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_test_k_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_test_k_mod.F90 -o $(OBJDIR)/rttov_test_k_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD)

$(MODDIR)/rttov_test_k_mod$(EXT_MOD): rttov_test_k_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_test_k_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_test_k_mod.F90 -o $(OBJDIR)/rttov_test_k_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_test_k_mod$(EXT_MOD) $(MODDIR)/rttov_test_k_mod$(EXT_MOD)


$(OBJDIR)/rttov_test_mod$(EXT_OBJ): rttov_test_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_lun$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_test_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_test_mod.F90 -o $(OBJDIR)/rttov_test_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_test_mod$(EXT_MOD) $(MODDIR)/rttov_test_mod$(EXT_MOD)

$(MODDIR)/rttov_test_mod$(EXT_MOD): rttov_test_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_lun$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_chain$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_test_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_test_mod.F90 -o $(OBJDIR)/rttov_test_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_test_mod$(EXT_MOD) $(MODDIR)/rttov_test_mod$(EXT_MOD)

