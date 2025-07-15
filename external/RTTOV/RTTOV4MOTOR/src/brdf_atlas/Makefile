
TOPDIR=../..
SUBDIR=brdf_atlas

FFLAGS_ARCH_mod_brdf_atlas                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_rttov_brdf_atlas              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_brdf_atlas_test             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_deallocate_brdf_atlas       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_get_brdf                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_setup_brdf_atlas            = $(FFLAGS_ARCH)


F90_INTF=rttov_deallocate_brdf_atlas.F90 rttov_get_brdf.F90 rttov_setup_brdf_atlas.F90

INTF=$(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_setup_brdf_atlas.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_deallocate_brdf_atlas$(EXT_OBJ) $(OBJDIR)/rttov_get_brdf$(EXT_OBJ) $(OBJDIR)/rttov_setup_brdf_atlas$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/mod_brdf_atlas$(EXT_OBJ) $(OBJDIR)/mod_rttov_brdf_atlas$(EXT_OBJ)

OBJBIN=$(OBJDIR)/rttov_brdf_atlas_test$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov_brdf_atlas_test$(EXT_BIN)

MODS=$(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD)

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


$(BINDIR)/rttov_brdf_atlas_test$(EXT_BIN): $(OBJDIR)/rttov_brdf_atlas_test$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_brdf_atlas_test$(EXT_BIN) $(OBJDIR)/rttov_brdf_atlas_test$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov_deallocate_brdf_atlas$(EXT_OBJ): rttov_deallocate_brdf_atlas.F90  $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_deallocate_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_deallocate_brdf_atlas.F90 -o $(OBJDIR)/rttov_deallocate_brdf_atlas$(EXT_OBJ)


$(OBJDIR)/rttov_get_brdf$(EXT_OBJ): rttov_get_brdf.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_get_brdf) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_get_brdf.F90 -o $(OBJDIR)/rttov_get_brdf$(EXT_OBJ)


$(OBJDIR)/rttov_setup_brdf_atlas$(EXT_OBJ): rttov_setup_brdf_atlas.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_setup_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_setup_brdf_atlas.F90 -o $(OBJDIR)/rttov_setup_brdf_atlas$(EXT_OBJ)


$(OBJDIR)/rttov_brdf_atlas_test$(EXT_OBJ): rttov_brdf_atlas_test.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_setup_brdf_atlas.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_brdf_atlas_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_brdf_atlas_test.F90 -o $(OBJDIR)/rttov_brdf_atlas_test$(EXT_OBJ)


$(OBJDIR)/mod_brdf_atlas$(EXT_OBJ): mod_brdf_atlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_hdf_load.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_lapack_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_solar_refl_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_brdf_atlas.F90 -o $(OBJDIR)/mod_brdf_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD)

$(MODDIR)/mod_brdf_atlas$(EXT_MOD): mod_brdf_atlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_hdf_load.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_lapack_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_solar_refl_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_brdf_atlas.F90 -o $(OBJDIR)/mod_brdf_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD)


$(OBJDIR)/mod_rttov_brdf_atlas$(EXT_OBJ): mod_rttov_brdf_atlas.F90  $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_rttov_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_rttov_brdf_atlas.F90 -o $(OBJDIR)/mod_rttov_brdf_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD)

$(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD): mod_rttov_brdf_atlas.F90  $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_rttov_brdf_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_rttov_brdf_atlas.F90 -o $(OBJDIR)/mod_rttov_brdf_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD)

