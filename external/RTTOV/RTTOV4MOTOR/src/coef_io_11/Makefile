
TOPDIR=../..
SUBDIR=coef_io_11

FFLAGS_ARCH_rttov11_conv_coef_11to12       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov11_conv_coef_12to11       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov11_read_ascii_coef        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov11_read_hdf5_coef         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov11_write_ascii_coef       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov11_write_hdf5_coef        = $(FFLAGS_ARCH)


F90_INTF=rttov11_read_ascii_coef.F90 rttov11_read_hdf5_coef.F90 rttov11_write_ascii_coef.F90 rttov11_write_hdf5_coef.F90

INTF=$(INCDIR)/rttov11_read_ascii_coef.interface $(INCDIR)/rttov11_read_hdf5_coef.interface $(INCDIR)/rttov11_write_ascii_coef.interface $(INCDIR)/rttov11_write_hdf5_coef.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov11_read_ascii_coef$(EXT_OBJ) $(OBJDIR)/rttov11_read_hdf5_coef$(EXT_OBJ) $(OBJDIR)/rttov11_write_ascii_coef$(EXT_OBJ) $(OBJDIR)/rttov11_write_hdf5_coef$(EXT_OBJ) 

OBJMOD=

OBJBIN=$(OBJDIR)/rttov11_conv_coef_11to12$(EXT_OBJ) $(OBJDIR)/rttov11_conv_coef_12to11$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io_11$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov11_conv_coef_11to12$(EXT_BIN) $(BINDIR)/rttov11_conv_coef_12to11$(EXT_BIN)

MODS=

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


$(BINDIR)/rttov11_conv_coef_11to12$(EXT_BIN): $(OBJDIR)/rttov11_conv_coef_11to12$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov11_conv_coef_11to12$(EXT_BIN) $(OBJDIR)/rttov11_conv_coef_11to12$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_coef_io_11 -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov11_conv_coef_12to11$(EXT_BIN): $(OBJDIR)/rttov11_conv_coef_12to11$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov11_conv_coef_12to11$(EXT_BIN) $(OBJDIR)/rttov11_conv_coef_12to11$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_coef_io_11 -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov11_read_ascii_coef$(EXT_OBJ): rttov11_read_ascii_coef.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_deletecomment.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov_findnextsection.interface $(INCDIR)/rttov_nullify_coef.interface $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_read_ascii_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_read_ascii_coef.F90 -o $(OBJDIR)/rttov11_read_ascii_coef$(EXT_OBJ)


$(OBJDIR)/rttov11_read_hdf5_coef$(EXT_OBJ): rttov11_read_hdf5_coef.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_nullify_coef.interface $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_read_hdf5_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_read_hdf5_coef.F90 -o $(OBJDIR)/rttov11_read_hdf5_coef$(EXT_OBJ)


$(OBJDIR)/rttov11_write_ascii_coef$(EXT_OBJ): rttov11_write_ascii_coef.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_write_ascii_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_write_ascii_coef.F90 -o $(OBJDIR)/rttov11_write_ascii_coef$(EXT_OBJ)


$(OBJDIR)/rttov11_write_hdf5_coef$(EXT_OBJ): rttov11_write_hdf5_coef.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_write_hdf5_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_write_hdf5_coef.F90 -o $(OBJDIR)/rttov11_write_hdf5_coef$(EXT_OBJ)


$(OBJDIR)/rttov11_conv_coef_11to12$(EXT_OBJ): rttov11_conv_coef_11to12.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov11_read_ascii_coef.interface $(INCDIR)/rttov_write_ascii_coef.interface $(INCDIR)/rttov11_read_hdf5_coef.interface $(INCDIR)/rttov_init_coef.interface $(INCDIR)/rttov_nullify_coef.interface $(INCDIR)/rttov_dealloc_coef.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_solar_refl_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_conv_coef_11to12) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_conv_coef_11to12.F90 -o $(OBJDIR)/rttov11_conv_coef_11to12$(EXT_OBJ)


$(OBJDIR)/rttov11_conv_coef_12to11$(EXT_OBJ): rttov11_conv_coef_12to11.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov_read_ascii_coef.interface $(INCDIR)/rttov11_write_ascii_coef.interface $(INCDIR)/rttov11_write_hdf5_coef.interface $(INCDIR)/rttov_init_coef.interface $(INCDIR)/rttov_nullify_coef.interface $(INCDIR)/rttov_dealloc_coef.interface $(INCDIR)/rttov_hdf_load.interface $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_nlte_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_fast_coef_utils_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_skin_io$(EXT_MOD) $(MODDIR)/rttov_solar_refl_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov11_conv_coef_12to11) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov11_conv_coef_12to11.F90 -o $(OBJDIR)/rttov11_conv_coef_12to11$(EXT_OBJ)

