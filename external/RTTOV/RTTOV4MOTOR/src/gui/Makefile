
TOPDIR=../..
SUBDIR=gui

FFLAGS_ARCH_rttov_gui_context        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_gui_f2py           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_gui_handle         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_gui_test_run       = $(FFLAGS_ARCH)


F90_INTF=rttov_gui_f2py.F90

INTF=$(INCDIR)/rttov_gui_f2py.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_gui_f2py$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/rttov_gui_context$(EXT_OBJ) $(OBJDIR)/rttov_gui_handle$(EXT_OBJ)

OBJBIN=$(OBJDIR)/rttov_gui_test_run$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_gui$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov_gui_test_run$(EXT_BIN)

MODS=$(MODDIR)/rttov_gui_context$(EXT_MOD) $(MODDIR)/rttov_gui_handle$(EXT_MOD)

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


$(BINDIR)/rttov_gui_test_run$(EXT_BIN): $(OBJDIR)/rttov_gui_test_run$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_gui_test_run$(EXT_BIN) $(OBJDIR)/rttov_gui_test_run$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_gui -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov_gui_f2py$(EXT_OBJ): rttov_gui_f2py.F90 $(INCDIR)/rttov_print_info.interface $(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_setup_brdf_atlas.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_alloc_pccomp.interface $(INCDIR)/rttov_aer_clim_prof.interface $(MODDIR)/rttov_gui_handle$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_gui_context$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/mod_camel_clim_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_f2py) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_f2py.F90 -o $(OBJDIR)/rttov_gui_f2py$(EXT_OBJ)


$(OBJDIR)/rttov_gui_test_run$(EXT_OBJ): rttov_gui_test_run.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_test_run) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_test_run.F90 -o $(OBJDIR)/rttov_gui_test_run$(EXT_OBJ)


$(OBJDIR)/rttov_gui_context$(EXT_OBJ): rttov_gui_context.F90 $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_gui_handle$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_context) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_context.F90 -o $(OBJDIR)/rttov_gui_context$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_gui_context$(EXT_MOD) $(MODDIR)/rttov_gui_context$(EXT_MOD)

$(MODDIR)/rttov_gui_context$(EXT_MOD): rttov_gui_context.F90 $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_gui_handle$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_context) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_context.F90 -o $(OBJDIR)/rttov_gui_context$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_gui_context$(EXT_MOD) $(MODDIR)/rttov_gui_context$(EXT_MOD)


$(OBJDIR)/rttov_gui_handle$(EXT_OBJ): rttov_gui_handle.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_dealloc_coefs.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_camel_clim_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_handle) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_handle.F90 -o $(OBJDIR)/rttov_gui_handle$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_gui_handle$(EXT_MOD) $(MODDIR)/rttov_gui_handle$(EXT_MOD)

$(MODDIR)/rttov_gui_handle$(EXT_MOD): rttov_gui_handle.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_dealloc_coefs.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_camel_clim_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_gui_handle) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_gui_handle.F90 -o $(OBJDIR)/rttov_gui_handle$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_gui_handle$(EXT_MOD) $(MODDIR)/rttov_gui_handle$(EXT_MOD)

