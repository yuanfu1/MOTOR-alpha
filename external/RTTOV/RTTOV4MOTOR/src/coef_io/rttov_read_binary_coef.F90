! Description:
!> @file
!!   Read a binary coefficient file, optionally extracting a subset of channels.
!
!> @brief
!!   Read a binary coefficient file, optionally extracting a subset of channels.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!!   Note that after reading a subset of channels RTTOV will identify them by
!!   indexes 1...SIZE(channels), not by the original channel numbers.
!!
!! @param[out]    err             status on exit
!! @param[in,out] coef            RTTOV optical depth coefficient structure
!! @param[in]     file_id         logical unit for input rtcoef file
!! @param[in]     channels        list of channels to read, optional
!! @param[in]     lbl             flag to indicate this was called from the LBL code, optional
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_read_binary_coef(err, coef, file_id, channels, lbl)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
        rttov_magic_string,     &
        rttov_magic_number,     &
        version_compatible_min, &
        version_compatible_max, &
        gas_id_mixed,           &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_wvcont,          &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        gas_id_so2,             &
        ngases_max,             &
        sensor_id_mw,           &
        sensor_id_po
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)            :: err
  TYPE(rttov_coef),   INTENT(INOUT)          :: coef
  INTEGER(KIND=jpim), INTENT(IN)             :: file_id
  INTEGER(KIND=jpim), INTENT(IN),   OPTIONAL :: channels(:)
  LOGICAL(KIND=jplm), INTENT(IN),   OPTIONAL :: lbl
!INTF_END

#include "rttov_errorreport.interface"

  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: n
  INTEGER(KIND=jpim) :: chn
  INTEGER(KIND=jpim) :: i, inczeeman
  INTEGER(KIND=jpim) :: nlte_count, nlte_start, nlte_file_nchan
  INTEGER(KIND=jpim), ALLOCATABLE :: nlte_chans(:)

  COMPLEX(KIND=jprb), POINTER :: values_c_0(:)
  COMPLEX(KIND=jprb), POINTER :: values_c_1(:)
  REAL   (KIND=jprb), POINTER :: values0(:)
  REAL   (KIND=jprb), POINTER :: values1(:)
  REAL   (KIND=jprb), POINTER :: values2(:)
  REAL   (KIND=jprb), POINTER :: values3(:)
  REAL   (KIND=jprb), POINTER :: values4(:)
  REAL   (KIND=jprb), POINTER :: iremis_values(:,:)
  REAL   (KIND=jprb), POINTER :: nlte_values(:,:,:,:)
  INTEGER(KIND=jpim), POINTER :: ivalues0(:)
  INTEGER(KIND=jpim), POINTER :: ivalues1(:)
  CHARACTER(LEN=16)  :: bin_check_string
  REAL(KIND=jprb)    :: bin_check_number
  REAL(KIND=jprb)    :: bin_check_value
  CHARACTER(LEN=80)  :: errMessage
  LOGICAL(KIND=jplm) :: section_present
  LOGICAL(KIND=jplm) :: lbl1
!- End of header --------------------------------------------------------

  TRY

  lbl1 = .FALSE.
  IF(PRESENT(lbl)) lbl1 = lbl

  all_channels = .NOT. PRESENT(channels)

  READ (file_id, iostat=err)bin_check_string, bin_check_number
  THROWM(err.NE.0, 'io status while reading header')

! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal
  THROWM(err.NE.0,'Wrong header string in file')

! Verification of single/double precision using a 5 digit number
! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._jprb - ABS(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_jprb .OR. bin_check_value < 0.99_jprb) err = errorstatus_fatal
  THROWM(err.NE.0,'File created with a different real precision (R4<->R8)')


  errMessage = 'io status while reading IDENTIFICATION'
  READ (file_id, iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, coef%id_sensor
  THROWM(err.NE.0,errMessage)

  READ (file_id, iostat=err)coef%id_comp_lvl, coef%id_creation_date, coef%id_creation, coef%id_common_name
  THROWM(err.NE.0,errMessage)

! Error if the compatibility version of the coefficient file
! is not in the range defined by the constant module

  IF (coef%id_comp_lvl < version_compatible_min .OR. &
      coef%id_comp_lvl > version_compatible_max) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,"Version of coefficient file is incompatible with RTTOV library")
  ENDIF

  errMessage = 'io status while reading FAST_MODEL_VARIABLES'

  READ (file_id, iostat=err)coef%fmv_model_def, coef%fmv_model_ver, coef%fmv_ori_nchn, coef%fmv_gas
  THROWM(err.NE.0,errMessage)

  IF (coef%fmv_model_ver > 9) THEN
    READ (file_id, iostat=err)coef%nlevels
    THROWM(err.NE.0,errMessage)
  ENDIF

  READ (file_id, iostat=err)coef%id_comp_pc
  THROWM(err.NE.0,errMessage)

  READ (file_id, iostat=err)inczeeman
  THROWM(err.NE.0,errMessage)
  coef%inczeeman = (inczeeman == 1)

! Take care of the user list of channels
! coef%fmv_ori_nchn store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests
  IF (all_channels) THEN
    coef%fmv_chn = coef%fmv_ori_nchn
  ELSE
    IF (MAXVAL(channels) > coef%fmv_ori_nchn) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Channel index out of range for coefficient file")
    ENDIF
    coef%fmv_chn = SIZE(channels)
  ENDIF

  ALLOCATE (coef%fmv_gas_id(coef%fmv_gas), &
            coef%fmv_gas_pos(ngases_max), &
            coef%fmv_var(coef%fmv_gas), &
            coef%fmv_coe(coef%fmv_gas), &
            coef%fmv_ncorr(coef%fmv_gas), &
            coef%fmv_lvl(coef%fmv_gas), &
            coef%ff_ori_chn(coef%fmv_chn), &
            coef%ff_val_chn(coef%fmv_chn), &
            coef%ff_cwn(coef%fmv_chn), &
            coef%ff_bco(coef%fmv_chn), &
            coef%ff_bcs(coef%fmv_chn), &
            coef%ff_gam(coef%fmv_chn), STAT = err)
  THROWM(err.NE.0, "allocation of coef arrays")

  coef%fmv_gas_id(:)  = 0_jpim
  coef%fmv_gas_pos(:) = 0_jpim
  coef%fmv_var(:)     = 0_jpim
  coef%fmv_lvl(:)     = 0_jpim
  coef%fmv_coe(:)     = 0_jpim
  coef%fmv_ncorr(:)   = 0_jpim

  IF (coef%fmv_model_ver <= 9) THEN
    READ (file_id, iostat=err)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var, coef%fmv_lvl
    THROWM(err.NE.0,errMessage)

    READ (file_id, iostat=err)coef%fmv_coe
    THROWM(err.NE.0,errMessage)

    coef%nlevels = coef%fmv_lvl(1)
  ELSE
    READ (file_id, iostat=err)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var
    THROWM(err.NE.0,errMessage)

    READ (file_id, iostat=err)coef%fmv_coe, coef%fmv_ncorr
    THROWM(err.NE.0,errMessage)

    coef%fmv_lvl = coef%nlevels
  ENDIF

  DO n = 1, coef%fmv_gas
    SELECT CASE (coef%fmv_gas_id(n))
    CASE (gas_id_mixed)
      coef%nmixed   = coef%fmv_var(n)
      coef%ncmixed  = coef%fmv_coe(n)
      coef%nccmixed = coef%fmv_ncorr(n)
    CASE (gas_id_watervapour)
      coef%nwater   = coef%fmv_var(n)
      coef%ncwater  = coef%fmv_coe(n)
      coef%nccwater = coef%fmv_ncorr(n)
    CASE (gas_id_ozone)
      coef%nozone   = coef%fmv_var(n)
      coef%ncozone  = coef%fmv_coe(n)
      coef%nccozone = coef%fmv_ncorr(n)
    CASE (gas_id_wvcont)
      coef%nwvcont   = coef%fmv_var(n)
      coef%ncwvcont  = coef%fmv_coe(n)
      coef%nccwvcont = coef%fmv_ncorr(n)
    CASE (gas_id_co2)
      coef%nco2   = coef%fmv_var(n)
      coef%ncco2  = coef%fmv_coe(n)
      coef%nccco2 = coef%fmv_ncorr(n)
    CASE (gas_id_n2o)
      coef%nn2o   = coef%fmv_var(n)
      coef%ncn2o  = coef%fmv_coe(n)
      coef%nccn2o = coef%fmv_ncorr(n)
    CASE (gas_id_co)
      coef%nco   = coef%fmv_var(n)
      coef%ncco  = coef%fmv_coe(n)
      coef%nccco = coef%fmv_ncorr(n)
    CASE (gas_id_ch4)
      coef%nch4   = coef%fmv_var(n)
      coef%ncch4  = coef%fmv_coe(n)
      coef%nccch4 = coef%fmv_ncorr(n)
    CASE (gas_id_so2)
      coef%nso2   = coef%fmv_var(n)
      coef%ncso2  = coef%fmv_coe(n)
      coef%nccso2 = coef%fmv_ncorr(n)
      END SELECT
  ENDDO

  coef%nlayers = coef%nlevels - 1


  READ (file_id)section_present
  errMessage = 'io status while reading TRANSMITTANCE_TRESHOLD'

  IF (section_present) THEN
    ALLOCATE (coef%tt_val_chn(coef%fmv_chn), coef%tt_a0(coef%fmv_chn), &
              coef%tt_a1(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of transmittance threshold arrays")

    IF (all_channels) THEN
      READ (file_id, iostat=err)coef%tt_val_chn, coef%tt_a0, coef%tt_a1
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (ivalues1(coef%fmv_ori_nchn), values1(coef%fmv_ori_nchn), &
                values2(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of temporary arrays")

      READ (file_id, iostat=err)ivalues1, values1, values2
      THROWM(err.NE.0,errMessage)

      coef%tt_val_chn(:) = ivalues1(channels(:))
      coef%tt_a0(:)      = values1(channels(:))
      coef%tt_a1(:)      = values2(channels(:))

      DEALLOCATE (ivalues1, values1, values2, STAT = err)
      THROWM(err.NE.0, "deallocation of temporary arrays")
    ENDIF
  ENDIF


  READ (file_id)section_present
  errMessage = 'io status while reading SOLAR_SPECTRUM'

  IF (section_present) THEN
    ALLOCATE (coef%ss_val_chn(coef%fmv_chn), &
              coef%ss_solar_spectrum(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of solar spectrum arrays")

    IF (coef%fmv_model_ver > 9) THEN
      ALLOCATE (coef%ss_rayleigh_ext(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ss_rayleigh_ext")
    ENDIF

    IF (all_channels) THEN
      IF (coef%fmv_model_ver > 9) THEN
        READ (file_id, iostat=err)coef%ss_val_chn, coef%ss_solar_spectrum, coef%ss_rayleigh_ext
      ELSE
        READ (file_id, iostat=err)coef%ss_val_chn, coef%ss_solar_spectrum
      ENDIF
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (ivalues1(coef%fmv_ori_nchn), &
                values1(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of temporary arrays")

      IF (coef%fmv_model_ver > 9) THEN
        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values2")

        READ (file_id, iostat=err)ivalues1, values1, values2
        THROWM(err.NE.0,errMessage)

        coef%ss_rayleigh_ext(:) = values2(channels(:))

        DEALLOCATE (values2, STAT = err)
        THROWM(err.NE.0, "deallocation of values2")
      ELSE
        READ (file_id, iostat=err)ivalues1, values1
        THROWM(err.NE.0,errMessage)
      ENDIF

      coef%ss_val_chn(:)        = ivalues1(channels(:))
      coef%ss_solar_spectrum(:) = values1(channels(:))

      DEALLOCATE (ivalues1, values1, STAT = err)
      THROWM(err.NE.0, "deallocation of temporary arrays")
    ENDIF
  ENDIF


  READ (file_id)section_present
  errMessage = 'io status while reading WATER_OPTICAL_CONSTANT'

  IF (section_present) THEN
    ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), &
              coef%woc_waopc_fw(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of water optical constant arrays")

    IF (all_channels) THEN
      READ (file_id, iostat=err)coef%woc_waopc_ow, coef%woc_waopc_fw
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (values_c_0(coef%fmv_ori_nchn), &
                values_c_1(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of temporary arrays")

      READ (file_id, iostat=err)values_c_0, values_c_1
      THROWM(err.NE.0,errMessage)

      coef%woc_waopc_ow(:) = values_c_0(channels(:))
      coef%woc_waopc_fw(:) = values_c_1(channels(:))

      DEALLOCATE (values_c_0, values_c_1, STAT = err)
      THROWM(err.NE.0, "deallocation of temporary arrays")
    ENDIF
  ENDIF


  READ (file_id)section_present
  errMessage = 'io status while reading WAVE_SPECTRUM'

  IF (section_present) THEN
    READ (file_id, iostat=err)coef%ws_nomega
    THROWM(err.NE.0,errMessage)

    ALLOCATE (coef%ws_k_omega(coef%ws_nomega), &
              coef%ws_npoint(coef%ws_nomega), STAT = err)
    THROWM(err.NE.0, "allocation of wave spectrum arrays")

    READ (file_id, iostat=err)coef%ws_k_omega, coef%ws_npoint
    THROWM(err.NE.0,errMessage)
  ENDIF


  errMessage = 'io status while reading FILTER_FUNCTIONS'

  IF (all_channels) THEN
    READ (file_id, iostat=err)coef%ff_ori_chn, coef%ff_val_chn, coef%ff_cwn, coef%ff_bco, coef%ff_bcs, coef%ff_gam
    THROWM(err.NE.0,errMessage)
  ELSE
    ALLOCATE (ivalues0(coef%fmv_ori_nchn), ivalues1(coef%fmv_ori_nchn), &
              values0(coef%fmv_ori_nchn), values1(coef%fmv_ori_nchn), &
              values2(coef%fmv_ori_nchn), values3(coef%fmv_ori_nchn), STAT = err)
    THROWM(err.NE.0, "allocation of temporary arrays")

    READ (file_id, iostat=err)ivalues0, ivalues1, values0, values1, values2, values3
    THROWM(err.NE.0,errMessage)

    coef%ff_ori_chn(:) = ivalues0(channels(:))
    coef%ff_val_chn(:) = ivalues1(channels(:))
    coef%ff_cwn(:)     = values0(channels(:))
    coef%ff_bco(:)     = values1(channels(:))
    coef%ff_bcs(:)     = values2(channels(:))
    coef%ff_gam(:)     = values3(channels(:))

    DEALLOCATE (ivalues0, ivalues1, values0, values1, values2, values3, STAT = err)
    THROWM(err.NE.0, "deallocation of temporary arrays")
  ENDIF


  errMessage = 'io status while reading FUNDAMENTAL_CONSTANTS'
  READ (file_id, iostat=err)coef%fc_planck_c1, coef%fc_planck_c2, coef%fc_sat_height
  THROWM(err.NE.0,errMessage)


  errMessage = 'io status while reading FASTEM'

  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
    ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of fastem_polar")

    IF (all_channels) THEN
      READ (file_id, iostat=err)coef%fastem_polar
      THROWM(err.NE.0,errMessage)

      IF (ANY(coef%fastem_polar == 7)) THEN
        ALLOCATE (coef%pol_phi(coef%fmv_chn), STAT = err)
        THROWM(err.NE.0, "allocation of pol_phi")

        READ (file_id, iostat=err)coef%pol_phi
        THROWM(err.NE.0,errMessage)
      ENDIF
    ELSE
      ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of ivalues0")

      READ (file_id, iostat=err)ivalues0
      THROWM(err.NE.0,errMessage)

      coef%fastem_polar(:) = ivalues0(channels(:))
      DEALLOCATE (ivalues0, STAT = err)
      THROWM(err.NE.0, "deallocation of ivalues0")

      IF (ANY(coef%fastem_polar == 7)) THEN
        ALLOCATE (coef%pol_phi(coef%fmv_chn), STAT = err)
        THROWM(err.NE.0, "allocation of pol_phi")
        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values0")

        READ (file_id, iostat=err)values0
        THROWM(err.NE.0,errMessage)
        coef%pol_phi(:) = values0(channels(:))
        DEALLOCATE (values0, STAT = err)
        THROWM(err.NE.0, "deallocation of values0")
      ENDIF
    ENDIF
  ENDIF


  errMessage = 'io status while reading IR emissivity model versions'
  READ (file_id, iostat=err)coef%ssirem_ver, coef%iremis_version
  THROWM(err.NE.0,errMessage)


  errMessage = 'io status while reading SSIREM'

  IF (coef%ssirem_ver >= 1) THEN
    ALLOCATE (coef%ssirem_a0(coef%fmv_chn), coef%ssirem_a1(coef%fmv_chn), &
              coef%ssirem_a2(coef%fmv_chn), coef%ssirem_xzn1(coef%fmv_chn), &
              coef%ssirem_xzn2(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of ssirem arrays")

    IF (all_channels) THEN
      READ (file_id, iostat=err)     &
          coef%ssirem_a0, coef%ssirem_a1, coef%ssirem_a2, coef%ssirem_xzn1, coef%ssirem_xzn2
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (values0(coef%fmv_ori_nchn), values1(coef%fmv_ori_nchn), &
                values2(coef%fmv_ori_nchn), values3(coef%fmv_ori_nchn), &
                values4(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of temporary arrays")

      READ (file_id, iostat=err)values0, values1, values2, values3, values4
      THROWM(err.NE.0,errMessage)

      coef%ssirem_a0(:)   = values0(channels(:))
      coef%ssirem_a1(:)   = values1(channels(:))
      coef%ssirem_a2(:)   = values2(channels(:))
      coef%ssirem_xzn1(:) = values3(channels(:))
      coef%ssirem_xzn2(:) = values4(channels(:))

      DEALLOCATE (values0, values1, values2, values3, values4, STAT = err)
      THROWM(err.NE.0, "deallocation of temporary arrays")
    ENDIF
  ENDIF


  errMessage = 'io status while reading IR_SEA_EMIS'

  IF (coef%iremis_version >= 1) THEN
    READ (file_id, iostat=err) coef%iremis_ncoef
    READ (file_id, iostat=err) coef%iremis_angle0, coef%iremis_tskin0

    ALLOCATE (coef%iremis_coef(coef%iremis_ncoef,coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of iremis_coef")

    IF (all_channels) THEN
      READ (file_id, iostat=err) coef%iremis_coef
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (iremis_values(coef%iremis_ncoef,coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of iremis_values")

      READ (file_id, iostat=err)iremis_values
      THROWM(err.NE.0,errMessage)

      coef%iremis_coef(:,:) = iremis_values(:,channels(:))
      DEALLOCATE (iremis_values, STAT = err)
      THROWM(err.NE.0, "deallocation of iremis_values")
    ENDIF
  ENDIF


  ALLOCATE (coef%ref_prfl_p(coef%fmv_lvl(gas_id_mixed)), &
            coef%ref_prfl_t(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), &
            coef%ref_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), &
            coef%bkg_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), &
            coef%lim_prfl_p(coef%fmv_lvl(gas_id_mixed)), &
            coef%env_prfl_tmax(coef%fmv_lvl(gas_id_mixed)), &
            coef%env_prfl_tmin(coef%fmv_lvl(gas_id_mixed)), &
            coef%env_prfl_gmin(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), &
            coef%env_prfl_gmax(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
  THROWM(err.NE.0, "allocation of reference and envelope profiles")

  errMessage = 'io status while reading REFERENCE PROFILE'
  READ (file_id, iostat=err)coef%ref_prfl_p, coef%ref_prfl_t, coef%ref_prfl_mr, coef%bkg_prfl_mr
  THROWM(err.NE.0,errMessage)

  errMessage = 'io status while reading PROFILE ENVELOPE'
  READ (file_id, iostat=err)     &
      coef%lim_prfl_p, coef%env_prfl_tmax, coef%env_prfl_tmin, coef%env_prfl_gmax, coef%env_prfl_gmin
  THROWM(err.NE.0,errMessage)


! FAST COEFFICIENT section

  ALLOCATE (coef%thermal(coef%fmv_chn), STAT = err)
  THROWM(err.NE.0, "allocation of thermal fast coefs")
  IF (coef%fmv_model_ver > 9) THEN
    ALLOCATE (coef%thermal_corr(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of thermal_corr fast coefs")
  ENDIF
  DO i = 1, coef%fmv_chn
    ALLOCATE (coef%thermal(i)%gasarray(coef%fmv_gas), STAT = err)
    THROWM(err.NE.0, "allocation of gasarray")
    CALL nullify_gas_coef_pointers(coef%thermal(i))
    IF (coef%fmv_model_ver > 9) THEN
      ALLOCATE (coef%thermal_corr(i)%gasarray(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of gasarray")
      CALL nullify_gas_coef_pointers(coef%thermal_corr(i))
    ENDIF
  ENDDO

  DO n = 1, coef%fmv_gas
    ! test existence of gas
    IF(coef%fmv_coe(n) > 0_jpim) THEN
      CALL read_binary_fast_coef(err, coef%thermal, n, coef%fmv_gas_id(n), &
                                 coef%fmv_coe(n), coef%fmv_model_ver, lbl1)
      THROWM(err.NE.0, "error occurred inside read_binary_fast_coef reading gas")

      IF (coef%fmv_ncorr(n) > 0_jpim) THEN
        CALL read_binary_fast_coef(err, coef%thermal_corr, n, coef%fmv_gas_id(n), &
                                   coef%fmv_ncorr(n), coef%fmv_model_ver, lbl1)
        THROWM(err.NE.0,"error occurred inside read_binary_fast_coef reading gas")
      ENDIF
    ENDIF
  ENDDO


! SOLAR_FAST_COEFFICIENTS SECTION
! We need to know if they are present, and if so read them just as for FAST_COEFS
  READ (file_id, iostat=err) coef%solarcoef
  THROWM(err.NE.0,'io status while reading solar coef section')

  IF (coef%solarcoef) THEN

    ALLOCATE (coef%solar(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of solar fast coefs")
    IF (coef%fmv_model_ver > 9) THEN
      ALLOCATE (coef%solar_corr(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of solar_corr fast coefs")
    ENDIF
    DO i = 1, coef%fmv_chn
      ALLOCATE (coef%solar(i)%gasarray(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of gasarray")
      CALL nullify_gas_coef_pointers(coef%solar(i))
      IF (coef%fmv_model_ver > 9) THEN
        ALLOCATE (coef%solar_corr(i)%gasarray(coef%fmv_gas), STAT = err)
        THROWM(err.NE.0, "allocation of gasarray")
        CALL nullify_gas_coef_pointers(coef%solar_corr(i))
      ENDIF
    ENDDO

    DO n = 1, coef%fmv_gas
      ! test existence of gas
      IF(coef%fmv_coe(n) > 0_jpim) THEN
        CALL read_binary_fast_coef(err, coef%solar, n, coef%fmv_gas_id(n), &
                                   coef%fmv_coe(n), coef%fmv_model_ver, lbl1)
        THROWM(err.NE.0, "error occurred inside read_binary_fast_coef")

        IF (coef%fmv_ncorr(n) > 0_jpim) THEN
          CALL read_binary_fast_coef(err, coef%solar_corr, n, coef%fmv_gas_id(n), &
                                     coef%fmv_ncorr(n), coef%fmv_model_ver, lbl1)
          THROWM(err.NE.0,"error occurred inside read_binary_fast_coef reading gas")
        ENDIF
      ENDIF
    ENDDO

  ENDIF


  READ (file_id) section_present
  errMessage = 'io status while reading PLANCK_WEIGHTED'

  IF (section_present) THEN
    ALLOCATE (coef%pw_val_chn(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of pw_val_chn")

    IF (all_channels) THEN
      READ (file_id, iostat=err)coef%pw_val_chn
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of ivalues1")

      READ (file_id, iostat=err)ivalues1
      THROWM(err.NE.0,errMessage)

      coef%pw_val_chn(:) = ivalues1(channels(:))

      DEALLOCATE (ivalues1, STAT = err)
      THROWM(err.NE.0, "deallocation of ivalues1")
    ENDIF
  ENDIF


  READ (file_id)section_present
  errMessage = 'io status while reading README_SPECTRAL_RESPONSE_FUNCTION'

  IF (section_present) THEN
    READ (file_id, iostat=err)coef%readme_srf
    THROWM(err.NE.0,errMessage)
  ENDIF


  READ (file_id) section_present
  errMessage = 'io status while reading NLTE_RADIANCE_COEFS'

  IF (section_present) THEN
    coef%nltecoef = .TRUE.
    ALLOCATE(coef%nlte_coef)

    READ(file_id, iostat=err) coef%nlte_coef%ncoef, &
    coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
    coef%nlte_coef%start_chan
    THROWM(err.NE.0,errMessage)

    NULLIFY (coef%nlte_coef%coef)
    NULLIFY (coef%nlte_coef%sol_zen_angle, coef%nlte_coef%cos_sol)
    NULLIFY (coef%nlte_coef%sat_zen_angle, coef%nlte_coef%sec_sat)

    ! For any monotonic channel selection we must find those selected channels
    ! which lie within the range of NLTE channels in the coef file. This
    ! constitutes another contiguous block of channels in the coef structure.
    IF (.NOT. all_channels) THEN
      ALLOCATE(nlte_chans(SIZE(channels)))
      nlte_count = 0
      nlte_start = 0
      DO i = 1, SIZE(channels)
        IF (i > 1) THEN
          IF (channels(i) < channels(i-1)) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0, "non-monotonic channel selection incompatible with NLTE coefficients")
          ENDIF
        ENDIF
        IF (channels(i) >= coef%nlte_coef%start_chan .AND. &
            channels(i) < coef%nlte_coef%start_chan + coef%nlte_coef%nchan) THEN
          nlte_count = nlte_count + 1
          nlte_chans(nlte_count) = channels(i) - coef%nlte_coef%start_chan + 1
          IF (nlte_count == 1) nlte_start = i
        ENDIF
      ENDDO

      IF (nlte_count > 0) THEN
        nlte_file_nchan           = coef%nlte_coef%nchan
        coef%nlte_coef%start_chan = nlte_start
        coef%nlte_coef%nchan      = nlte_count
      ELSE
        coef%nlte_coef%start_chan = 0_jpim
        coef%nlte_coef%nchan = 0_jpim
        DEALLOCATE(coef%nlte_coef, nlte_chans)
        coef%nltecoef = .FALSE.
      ENDIF
    ENDIF

    IF (coef%nltecoef) THEN

      ALLOCATE(coef%nlte_coef%sol_zen_angle(coef%nlte_coef%nsol), &
               coef%nlte_coef%sec_sat(coef%nlte_coef%nsat), &
               coef%nlte_coef%coef(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
               coef%nlte_coef%nsol, coef%nlte_coef%nchan), STAT = err)
      THROWM(err.NE.0, "allocation of nlte_coef")

      READ(file_id, iostat=err) coef%nlte_coef%sol_zen_angle
      THROWM(err.NE.0,errMessage)
      READ(file_id, iostat=err) coef%nlte_coef%sec_sat
      THROWM(err.NE.0,errMessage)

      IF (all_channels) THEN
        READ(file_id, iostat=err) coef%nlte_coef%coef
        THROWM(err.NE.0,errMessage)
      ELSE
        ALLOCATE(nlte_values(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
                 coef%nlte_coef%nsol, nlte_file_nchan), STAT = err)
        THROWM(err.NE.0, "allocation of nlte_values array")
        READ(file_id, iostat=err) nlte_values
        THROWM(err.NE.0,errMessage)
        coef%nlte_coef%coef(:,:,:,:) = nlte_values(:,:,:,nlte_chans(1:nlte_count))
        DEALLOCATE(nlte_chans, nlte_values)
      ENDIF

    ENDIF
  ENDIF


  READ (file_id) section_present
  errMessage = 'io status while reading PRESSURE_MODULATED_CELL'

  IF (section_present) THEN
    coef%pmc_shift = .TRUE.

    READ (file_id, iostat=err) coef%pmc_lengthcell
    THROWM(err.NE.0,errMessage)

    ALLOCATE (coef%pmc_pnominal(coef%fmv_chn), STAT = err)
    THROWM(err.NE.0, "allocation of coef%pmc_pnominal array")
    IF (all_channels) THEN
      READ (file_id, iostat=err) coef%pmc_pnominal
      THROWM(err.NE.0,errMessage)
    ELSE
      ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of values1")
      READ (file_id, iostat=err) values1
      THROWM(err.NE.0,errMessage)
      coef%pmc_pnominal(:) = values1(channels(:))
      DEALLOCATE (values1, STAT = err)
      THROWM(err.NE.0, "deallocation of values1")
    ENDIF

    READ (file_id, iostat=err) coef%pmc_tempcell, coef%pmc_betaplus1, &
        coef%pmc_nlay, coef%pmc_nvar
    THROWM(err.NE.0,errMessage)

    ALLOCATE (coef%pmc_coef(coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar), STAT = err)
    THROWM(err.NE.0, "allocation of coef%pmc_coef array")

    DO chn = 1, coef%fmv_ori_nchn
      IF (all_channels) THEN
        READ (file_id, iostat=err) coef%pmc_coef(:, chn, :)
      ELSE
        DO i = 1, coef%fmv_chn
          IF (chn == channels(i)) EXIT
        END DO
        IF (i > coef%fmv_chn) THEN
          READ (file_id, iostat=err)
        ELSE
          READ (file_id, iostat=err)coef%pmc_coef(:, i, :)
        END IF
      END IF
      THROWM(err.NE.0,errMessage)
    ENDDO
  ENDIF


    READ (file_id)section_present
    errMessage = 'io status while reading LINE_BY_LINE'

    IF (section_present) THEN
      READ (file_id, iostat=err)coef%line_by_line
      THROWM(err.NE.0,errMessage)
    ENDIF

!
! Here add reading of new sections for binary format in order to keep compatibility with
! previous versions
!
  CATCH

CONTAINS

  SUBROUTINE read_binary_fast_coef(err, fast_coef, gas_pos, gas_id, ncoef, version, lbl_mode)

    USE rttov_types, ONLY : rttov_fast_coef
    USE rttov_fast_coef_utils_mod, ONLY : set_pointers

    INTEGER(jpim),         INTENT(OUT)          :: err
    TYPE(rttov_fast_coef), INTENT(INOUT)        :: fast_coef(:)
    INTEGER(jpim),         INTENT(IN)           :: gas_pos
    INTEGER(jpim),         INTENT(IN)           :: gas_id
    INTEGER(jpim),         INTENT(IN)           :: ncoef
    INTEGER(jpim),         INTENT(IN)           :: version
    LOGICAL(jplm),         INTENT(IN)           :: lbl_mode

    REAL(jprb), ALLOCATABLE :: coef_temp(:,:)
    INTEGER(KIND=jpim) :: i, chn
    LOGICAL(jplm) :: data_present

    TRY
    IF (version <= 9) THEN
      !This is transposed from pre-RTTOV 11.3 order
      !Coefs will be transposed for more efficent access in memory
      !but remain same on disk (for pre-RTTOV 11.3 coefs)
      ALLOCATE(coef_temp(coef%nlayers,ncoef), stat=err)
      THROW(err.NE.0)
    ENDIF

    chn = 1
    DO i = 1, coef%fmv_chn
! Reading a subset of channels from coefficient file. Discarding unwanted
      IF (.NOT. all_channels) THEN
        DO WHILE (chn /= channels(i))
          READ (file_id, iostat=err) data_present
          THROW(err.NE.0)
          IF (data_present) THEN
            READ (file_id, iostat=err) ! dump contents
            THROW(err.NE.0)
          ENDIF
          chn = chn + 1
        ENDDO ! now chn == channels(i)
      ENDIF

      READ (file_id, iostat=err) data_present
      THROW(err.NE.0)

      chn = chn + 1

      IF (lbl_mode .AND. .NOT. data_present) THEN
        ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef,coef%nlayers), stat=err)
        THROW(err.NE.0)
        fast_coef(i)%gasarray(gas_pos)%coef = 0._jprb
        CALL set_pointers(fast_coef(i), gas_pos, gas_id)
      ELSE
        IF (data_present) THEN
          ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef,coef%nlayers), stat=err)
          THROW(err.NE.0)
          IF (version > 9) THEN
            READ (file_id, iostat=err) fast_coef(i)%gasarray(gas_pos)%coef
            THROW(err.NE.0)
          ELSE
            READ (file_id, iostat=err) coef_temp(:,:)
            THROW (err.NE.0)! RETURN
            fast_coef(i)%gasarray(gas_pos)%coef = TRANSPOSE(coef_temp)
          ENDIF
          CALL set_pointers(fast_coef(i), gas_pos, gas_id)
        ENDIF
      ENDIF
    ENDDO

    IF (version <= 9) THEN
      DEALLOCATE(coef_temp, stat = err)
      THROW(err.NE.0)
    ENDIF
    DO WHILE (chn <= coef%fmv_ori_nchn)
      READ (file_id, iostat=err) data_present
      THROW(err.NE.0)
      IF (data_present) THEN
        READ (file_id, iostat=err) ! dump contents to end of record
        THROW(err.NE.0)
      ENDIF
      chn = chn + 1
    ENDDO
    CATCH
  END SUBROUTINE read_binary_fast_coef

  SUBROUTINE nullify_gas_coef_pointers(fast_coef)
    USE rttov_types, ONLY : rttov_fast_coef
    TYPE(rttov_fast_coef), INTENT(INOUT) :: fast_coef
    NULLIFY (fast_coef%mixedgas,    &
             fast_coef%watervapour, &
             fast_coef%ozone,       &
             fast_coef%wvcont,      &
             fast_coef%co2,         &
             fast_coef%n2o,         &
             fast_coef%co,          &
             fast_coef%ch4,         &
             fast_coef%so2)
  END SUBROUTINE nullify_gas_coef_pointers

END SUBROUTINE rttov_read_binary_coef
