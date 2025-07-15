
TOPDIR=../..
SUBDIR=hdf

FFLAGS_ARCH_rttov_hdf_chanprof_io                  = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_coefs                        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_emissivity_io                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_load                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_mod                          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_opt_param_io                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_config_io            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_htfrtc_opts_io       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_interp_io            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_io                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_pc_io                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_rt_all_io            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_rt_ir_io             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_options_rt_mw_io             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_pccomp_io                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_profile_io                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_profiles                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_radiance2_io                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_radiance_io                  = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_reflectance_io               = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_coef_io                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc1_io           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc2_io           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc_io            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_fast_coef_io           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_rttov_nlte_coef_io           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_s2m_io                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_save                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_skin_io                      = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_hdf_transmission_io              = $(FFLAGS_ARCH)


F90_INTF=rttov_hdf_load.F90 rttov_hdf_save.F90

INTF=$(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_hdf_load$(EXT_OBJ) $(OBJDIR)/rttov_hdf_save$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/rttov_hdf_chanprof_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_coefs$(EXT_OBJ) $(OBJDIR)/rttov_hdf_emissivity_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_mod$(EXT_OBJ) $(OBJDIR)/rttov_hdf_opt_param_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_config_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_interp_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_pc_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_rt_all_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_rt_ir_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_options_rt_mw_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_pccomp_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_profile_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_profiles$(EXT_OBJ) $(OBJDIR)/rttov_hdf_radiance2_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_radiance_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_reflectance_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_coef_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_s2m_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_skin_io$(EXT_OBJ) $(OBJDIR)/rttov_hdf_transmission_io$(EXT_OBJ)

OBJBIN=

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=

MODS=$(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD)

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


$(OBJDIR)/rttov_hdf_load$(EXT_OBJ): rttov_hdf_load.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_load) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_load.F90 -o $(OBJDIR)/rttov_hdf_load$(EXT_OBJ)


$(OBJDIR)/rttov_hdf_save$(EXT_OBJ): rttov_hdf_save.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_save) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_save.F90 -o $(OBJDIR)/rttov_hdf_save$(EXT_OBJ)


$(OBJDIR)/rttov_hdf_chanprof_io$(EXT_OBJ): rttov_hdf_chanprof_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_chanprof_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_chanprof_io.F90 -o $(OBJDIR)/rttov_hdf_chanprof_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD): rttov_hdf_chanprof_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_chanprof_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_chanprof_io.F90 -o $(OBJDIR)/rttov_hdf_chanprof_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_coefs$(EXT_OBJ): rttov_hdf_coefs.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_alloc_phfn_int.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_coefs) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_coefs.F90 -o $(OBJDIR)/rttov_hdf_coefs$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD)

$(MODDIR)/rttov_hdf_coefs$(EXT_MOD): rttov_hdf_coefs.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_alloc_phfn_int.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_coefs) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_coefs.F90 -o $(OBJDIR)/rttov_hdf_coefs$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD)


$(OBJDIR)/rttov_hdf_emissivity_io$(EXT_OBJ): rttov_hdf_emissivity_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_emissivity_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_emissivity_io.F90 -o $(OBJDIR)/rttov_hdf_emissivity_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD): rttov_hdf_emissivity_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_emissivity_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_emissivity_io.F90 -o $(OBJDIR)/rttov_hdf_emissivity_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_mod$(EXT_OBJ): rttov_hdf_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_mod.F90 -o $(OBJDIR)/rttov_hdf_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)

$(MODDIR)/rttov_hdf_mod$(EXT_MOD): rttov_hdf_mod.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_mod.F90 -o $(OBJDIR)/rttov_hdf_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)


$(OBJDIR)/rttov_hdf_opt_param_io$(EXT_OBJ): rttov_hdf_opt_param_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_opt_param_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_opt_param_io.F90 -o $(OBJDIR)/rttov_hdf_opt_param_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD): rttov_hdf_opt_param_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_opt_param_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_opt_param_io.F90 -o $(OBJDIR)/rttov_hdf_opt_param_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_config_io$(EXT_OBJ): rttov_hdf_options_config_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_config_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_config_io.F90 -o $(OBJDIR)/rttov_hdf_options_config_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD): rttov_hdf_options_config_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_config_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_config_io.F90 -o $(OBJDIR)/rttov_hdf_options_config_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_OBJ): rttov_hdf_options_htfrtc_opts_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_htfrtc_opts_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_htfrtc_opts_io.F90 -o $(OBJDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_htfrtc_opts_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_MOD): rttov_hdf_options_htfrtc_opts_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_htfrtc_opts_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_htfrtc_opts_io.F90 -o $(OBJDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_htfrtc_opts_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_htfrtc_opts_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_interp_io$(EXT_OBJ): rttov_hdf_options_interp_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_interp_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_interp_io.F90 -o $(OBJDIR)/rttov_hdf_options_interp_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD): rttov_hdf_options_interp_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_interp_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_interp_io.F90 -o $(OBJDIR)/rttov_hdf_options_interp_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_io$(EXT_OBJ): rttov_hdf_options_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_io.F90 -o $(OBJDIR)/rttov_hdf_options_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_io$(EXT_MOD): rttov_hdf_options_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_io.F90 -o $(OBJDIR)/rttov_hdf_options_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_pc_io$(EXT_OBJ): rttov_hdf_options_pc_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_pc_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_pc_io.F90 -o $(OBJDIR)/rttov_hdf_options_pc_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD): rttov_hdf_options_pc_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_pc_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_pc_io.F90 -o $(OBJDIR)/rttov_hdf_options_pc_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_rt_all_io$(EXT_OBJ): rttov_hdf_options_rt_all_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_all_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_all_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_all_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD): rttov_hdf_options_rt_all_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_all_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_all_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_all_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_rt_ir_io$(EXT_OBJ): rttov_hdf_options_rt_ir_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_ir_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_ir_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_ir_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD): rttov_hdf_options_rt_ir_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_ir_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_ir_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_ir_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_options_rt_mw_io$(EXT_OBJ): rttov_hdf_options_rt_mw_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_mw_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_mw_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_mw_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD): rttov_hdf_options_rt_mw_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_options_rt_mw_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_options_rt_mw_io.F90 -o $(OBJDIR)/rttov_hdf_options_rt_mw_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_pccomp_io$(EXT_OBJ): rttov_hdf_pccomp_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_pccomp_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_pccomp_io.F90 -o $(OBJDIR)/rttov_hdf_pccomp_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD): rttov_hdf_pccomp_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_pccomp_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_pccomp_io.F90 -o $(OBJDIR)/rttov_hdf_pccomp_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_profile_io$(EXT_OBJ): rttov_hdf_profile_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_profile_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_profile_io.F90 -o $(OBJDIR)/rttov_hdf_profile_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_profile_io$(EXT_MOD): rttov_hdf_profile_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_profile_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_profile_io.F90 -o $(OBJDIR)/rttov_hdf_profile_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_profiles$(EXT_OBJ): rttov_hdf_profiles.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_profiles) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_profiles.F90 -o $(OBJDIR)/rttov_hdf_profiles$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD)

$(MODDIR)/rttov_hdf_profiles$(EXT_MOD): rttov_hdf_profiles.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_profiles) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_profiles.F90 -o $(OBJDIR)/rttov_hdf_profiles$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD)


$(OBJDIR)/rttov_hdf_radiance2_io$(EXT_OBJ): rttov_hdf_radiance2_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_radiance2_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_radiance2_io.F90 -o $(OBJDIR)/rttov_hdf_radiance2_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD): rttov_hdf_radiance2_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_radiance2_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_radiance2_io.F90 -o $(OBJDIR)/rttov_hdf_radiance2_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_radiance_io$(EXT_OBJ): rttov_hdf_radiance_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_radiance_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_radiance_io.F90 -o $(OBJDIR)/rttov_hdf_radiance_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD): rttov_hdf_radiance_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_radiance_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_radiance_io.F90 -o $(OBJDIR)/rttov_hdf_radiance_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_reflectance_io$(EXT_OBJ): rttov_hdf_reflectance_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_reflectance_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_reflectance_io.F90 -o $(OBJDIR)/rttov_hdf_reflectance_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD): rttov_hdf_reflectance_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_reflectance_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_reflectance_io.F90 -o $(OBJDIR)/rttov_hdf_reflectance_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_coef_io$(EXT_OBJ): rttov_hdf_rttov_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD): rttov_hdf_rttov_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_OBJ): rttov_hdf_rttov_coef_pcc1_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc1_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc1_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD): rttov_hdf_rttov_coef_pcc1_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc1_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc1_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc1_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_OBJ): rttov_hdf_rttov_coef_pcc2_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc2_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc2_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD): rttov_hdf_rttov_coef_pcc2_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc2_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc2_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc2_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_OBJ): rttov_hdf_rttov_coef_pcc_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD): rttov_hdf_rttov_coef_pcc_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_coef_pcc_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_coef_pcc_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_coef_pcc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_pcc_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_OBJ): rttov_hdf_rttov_fast_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_fast_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_fast_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD): rttov_hdf_rttov_fast_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_fast_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_fast_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_OBJ): rttov_hdf_rttov_nlte_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_nlte_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_nlte_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD): rttov_hdf_rttov_nlte_coef_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_rttov_nlte_coef_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_rttov_nlte_coef_io.F90 -o $(OBJDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_s2m_io$(EXT_OBJ): rttov_hdf_s2m_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_s2m_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_s2m_io.F90 -o $(OBJDIR)/rttov_hdf_s2m_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD): rttov_hdf_s2m_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_s2m_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_s2m_io.F90 -o $(OBJDIR)/rttov_hdf_s2m_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_skin_io$(EXT_OBJ): rttov_hdf_skin_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_skin_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_skin_io.F90 -o $(OBJDIR)/rttov_hdf_skin_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_skin_io$(EXT_MOD): rttov_hdf_skin_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_skin_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_skin_io.F90 -o $(OBJDIR)/rttov_hdf_skin_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD)


$(OBJDIR)/rttov_hdf_transmission_io$(EXT_OBJ): rttov_hdf_transmission_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_transmission_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_transmission_io.F90 -o $(OBJDIR)/rttov_hdf_transmission_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD)

$(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD): rttov_hdf_transmission_io.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_hdf_transmission_io) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_hdf_transmission_io.F90 -o $(OBJDIR)/rttov_hdf_transmission_io$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD)

