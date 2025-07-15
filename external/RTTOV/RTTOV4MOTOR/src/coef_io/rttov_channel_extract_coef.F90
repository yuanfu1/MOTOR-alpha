! Description:
!> @file
!!   Extract data for given channel list from an optical depth coefficients
!!   structure.
!
!> @brief
!!   Extract data for given channel list from an optical depth coefficients
!!   structure.
!!
!! @details
!!   This is used by HDF5 I/O code to read in a subset of channels from a
!!   coefficient file. The first coef argument contains the coefficients
!!   from the file. The second argument is an uninitialised structure
!!   which contains the extracted coefficients on exit.
!!
!! @param[out]     err       status on exit
!! @param[in]      coef1     input coefficients read from file
!! @param[in,out]  coef2     output coefficients, uninitialised on entry
!! @param[in]      channels  list of channels to extract
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
SUBROUTINE rttov_channel_extract_coef(err, coef1, coef2, channels)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef
!INTF_OFF
  USE parkind1, ONLY : jplm
  USE rttov_types, ONLY : rttov_fast_coef
  USE rttov_fast_coef_utils_mod, ONLY : set_pointers
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  TYPE(rttov_coef), INTENT(IN)    :: coef1
  TYPE(rttov_coef), INTENT(INOUT) :: coef2
  INTEGER(jpim),    INTENT(IN)    :: channels(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim)              :: i
  INTEGER(jpim)              :: nlte_start, nlte_count
  INTEGER(jpim), ALLOCATABLE :: nlte_chans(:)
! ----------------------------------------------------------------------------

TRY

  IF (MAXVAL(channels) > coef1%fmv_ori_nchn) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Channel index out of range for coefficient file')
  ENDIF
  ! Scalar variables

  ! The number of extracted channels is different...
  coef2%fmv_chn          = SIZE(channels)

  ! ... but everything else is the same
  coef2%id_platform      = coef1%id_platform
  coef2%id_sat           = coef1%id_sat
  coef2%id_inst          = coef1%id_inst
  coef2%id_sensor        = coef1%id_sensor
  coef2%id_comp_lvl      = coef1%id_comp_lvl
  coef2%id_creation_date = coef1%id_creation_date
  coef2%line_by_line     = coef1%line_by_line
  coef2%readme_srf       = coef1%readme_srf
  coef2%id_creation      = coef1%id_creation
  coef2%id_common_name   = coef1%id_common_name
  coef2%fmv_model_def    = coef1%fmv_model_def
  coef2%fmv_model_ver    = coef1%fmv_model_ver
  coef2%id_comp_pc       = coef1%id_comp_pc
  coef2%inczeeman        = coef1%inczeeman
  coef2%fc_planck_c1     = coef1%fc_planck_c1
  coef2%fc_planck_c2     = coef1%fc_planck_c2
  coef2%fc_sat_height    = coef1%fc_sat_height
  coef2%fmv_ori_nchn     = coef1%fmv_ori_nchn
  coef2%fmv_gas          = coef1%fmv_gas
  coef2%nmixed           = coef1%nmixed
  coef2%nwater           = coef1%nwater
  coef2%nozone           = coef1%nozone
  coef2%nwvcont          = coef1%nwvcont
  coef2%nco2             = coef1%nco2
  coef2%nn2o             = coef1%nn2o
  coef2%nco              = coef1%nco
  coef2%nch4             = coef1%nch4
  coef2%nso2             = coef1%nso2
  coef2%nlevels          = coef1%nlevels
  coef2%nlayers          = coef1%nlayers
  coef2%ncmixed          = coef1%ncmixed
  coef2%ncwater          = coef1%ncwater
  coef2%ncozone          = coef1%ncozone
  coef2%ncwvcont         = coef1%ncwvcont
  coef2%ncco2            = coef1%ncco2
  coef2%ncn2o            = coef1%ncn2o
  coef2%ncco             = coef1%ncco
  coef2%ncch4            = coef1%ncch4
  coef2%ncso2            = coef1%ncso2
  coef2%nccmixed         = coef1%nccmixed
  coef2%nccwater         = coef1%nccwater
  coef2%nccozone         = coef1%nccozone
  coef2%nccwvcont        = coef1%nccwvcont
  coef2%nccco2           = coef1%nccco2
  coef2%nccn2o           = coef1%nccn2o
  coef2%nccco            = coef1%nccco
  coef2%nccch4           = coef1%nccch4
  coef2%nccso2           = coef1%nccso2
  coef2%ws_nomega        = coef1%ws_nomega
  coef2%ssirem_ver       = coef1%ssirem_ver
  coef2%iremis_version   = coef1%iremis_version
  coef2%iremis_ncoef     = coef1%iremis_ncoef
  coef2%iremis_angle0    = coef1%iremis_angle0
  coef2%iremis_tskin0    = coef1%iremis_tskin0
  coef2%solarcoef        = coef1%solarcoef
  coef2%nltecoef         = coef1%nltecoef
  coef2%pmc_shift        = coef1%pmc_shift
  coef2%pmc_lengthcell   = coef1%pmc_lengthcell
  coef2%pmc_tempcell     = coef1%pmc_tempcell
  coef2%pmc_betaplus1    = coef1%pmc_betaplus1
  coef2%pmc_nlay         = coef1%pmc_nlay
  coef2%pmc_nvar         = coef1%pmc_nvar


  ! Fast model variables

  IF (ASSOCIATED(coef1%fmv_gas_id)) THEN

    ALLOCATE(coef2%fmv_gas_id(coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_chn')
    coef2%fmv_gas_id = coef1%fmv_gas_id

    ALLOCATE(coef2%fmv_gas_pos(SIZE(coef1%fmv_gas_pos)), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_gas_pos')
    coef2%fmv_gas_pos = coef1%fmv_gas_pos

    ALLOCATE(coef2%fmv_var(coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_var')
    coef2%fmv_var = coef1%fmv_var

    ALLOCATE(coef2%fmv_lvl(coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_lvl')
    coef2%fmv_lvl = coef1%fmv_lvl

    ALLOCATE(coef2%fmv_coe(coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_coe')
    coef2%fmv_coe = coef1%fmv_coe

    ALLOCATE(coef2%fmv_ncorr(coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of fmv_ncorr')
    coef2%fmv_ncorr = coef1%fmv_ncorr

  ENDIF


  ! Filter

  IF (ASSOCIATED(coef1%ff_ori_chn)) THEN

    ALLOCATE(coef2%ff_ori_chn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_ori_chn')
    coef2%ff_ori_chn = coef1%ff_ori_chn(channels)

    ALLOCATE(coef2%ff_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_val_chn')
    coef2%ff_val_chn = coef1%ff_val_chn(channels)

    ALLOCATE(coef2%ff_cwn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_cwn')
    coef2%ff_cwn = coef1%ff_cwn(channels)

    ALLOCATE(coef2%ff_bco(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_bco')
    coef2%ff_bco = coef1%ff_bco(channels)

    ALLOCATE(coef2%ff_bcs(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_bcs')
    coef2%ff_bcs = coef1%ff_bcs(channels)

    ALLOCATE(coef2%ff_gam(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ff_gam')
    coef2%ff_gam = coef1%ff_gam(channels)

  ENDIF


  ! Transmittance threshold

  IF (ASSOCIATED(coef1%tt_a0)) THEN

    ALLOCATE(coef2%tt_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of tt_val_chn')
    coef2%tt_val_chn = coef1%tt_val_chn(channels)

    ALLOCATE(coef2%tt_a0(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of tt_a0')
    coef2%tt_a0 = coef1%tt_a0(channels)

    ALLOCATE(coef2%tt_a1(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of tt_a1')
    coef2%tt_a1 = coef1%tt_a1(channels)

  ENDIF


  ! Planck-weighted

  IF (ASSOCIATED(coef1%pw_val_chn)) THEN

    ALLOCATE(coef2%pw_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of pw_val_chn')
    coef2%pw_val_chn = coef1%pw_val_chn(channels)

  ENDIF


  ! Solar spectrum

  IF (ASSOCIATED(coef1%ss_solar_spectrum)) THEN

    ALLOCATE(coef2%ss_val_chn(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ss_val_chn')
    coef2%ss_val_chn = coef1%ss_val_chn(channels)

    ALLOCATE(coef2%ss_solar_spectrum(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ss_solar_spectrum')
    coef2%ss_solar_spectrum = coef1%ss_solar_spectrum(channels)

    IF (ASSOCIATED(coef1%ss_rayleigh_ext)) THEN
      ALLOCATE(coef2%ss_rayleigh_ext(coef2%fmv_chn), stat=err)
      THROWM(err.NE.0, 'allocation of ss_rayleigh_ext')
      coef2%ss_rayleigh_ext = coef1%ss_rayleigh_ext(channels)
    ENDIF

  ENDIF


  ! Water optical constants

  IF (ASSOCIATED(coef1%woc_waopc_ow)) THEN

    ALLOCATE(coef2%woc_waopc_ow(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of woc_waopc_ow')
    coef2%woc_waopc_ow = coef1%woc_waopc_ow(channels)

    ALLOCATE(coef2%woc_waopc_fw(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of woc_waopc_fw')
    coef2%woc_waopc_fw = coef1%woc_waopc_fw(channels)

  ENDIF


  ! Wave spectrum

  IF (ASSOCIATED(coef1%ws_npoint)) THEN

    ALLOCATE(coef2%ws_npoint(coef2%ws_nomega), stat=err)
    THROWM(err.NE.0, 'allocation of ws_npoint')
    coef2%ws_npoint = coef1%ws_npoint

    ALLOCATE(coef2%ws_k_omega(coef2%ws_nomega), stat=err)
    THROWM(err.NE.0, 'allocation of ws_k_omega')
    coef2%ws_k_omega = coef1%ws_k_omega

  ENDIF


  ! FASTEM

  IF (ASSOCIATED(coef1%fastem_polar)) THEN

    ALLOCATE(coef2%fastem_polar(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of fastem_polar')
    coef2%fastem_polar = coef1%fastem_polar(channels)

    IF (ASSOCIATED(coef1%pol_phi)) THEN

      ALLOCATE(coef2%pol_phi(coef2%fmv_chn), stat=err)
      THROWM(err.NE.0, 'allocation of pol_phi')
      coef2%pol_phi = coef1%pol_phi(channels)

    ENDIF

  ENDIF


  ! SSIREM

  IF (ASSOCIATED(coef1%ssirem_a0)) THEN

    ALLOCATE(coef2%ssirem_a0(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ssirem_a0')
    coef2%ssirem_a0 = coef1%ssirem_a0(channels)

    ALLOCATE(coef2%ssirem_a1(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ssirem_a1')
    coef2%ssirem_a1 = coef1%ssirem_a1(channels)

    ALLOCATE(coef2%ssirem_a2(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ssirem_a2')
    coef2%ssirem_a2 = coef1%ssirem_a2(channels)

    ALLOCATE(coef2%ssirem_xzn1(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ssirem_xzn1')
    coef2%ssirem_xzn1 = coef1%ssirem_xzn1(channels)

    ALLOCATE(coef2%ssirem_xzn2(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of ssirem_xzn2')
    coef2%ssirem_xzn2 = coef1%ssirem_xzn2(channels)

  ENDIF


  ! IR SEA EMIS

  IF (ASSOCIATED(coef1%iremis_coef)) THEN

    ALLOCATE(coef2%iremis_coef(coef2%iremis_ncoef,coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of iremis_coef')
    coef2%iremis_coef = coef1%iremis_coef(:,channels)

  ENDIF


  ! Reference profile

  IF (ASSOCIATED(coef1%ref_prfl_p)) THEN

    ALLOCATE(coef2%ref_prfl_p(coef2%nlevels), stat=err)
    THROWM(err.NE.0, 'allocation of ref_prfl_p')
    coef2%ref_prfl_p = coef1%ref_prfl_p

    ALLOCATE(coef2%ref_prfl_t(coef2%nlevels, coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of ref_prfl_t')
    coef2%ref_prfl_t = coef1%ref_prfl_t

    ALLOCATE(coef2%ref_prfl_mr(coef2%nlevels, coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of ref_prfl_mr')
    coef2%ref_prfl_mr = coef1%ref_prfl_mr

    ALLOCATE(coef2%bkg_prfl_mr(coef2%nlevels, coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of bkg_prfl_mr')
    coef2%bkg_prfl_mr = coef1%bkg_prfl_mr

  ENDIF


  ! Profile envelope

  IF (ASSOCIATED(coef1%lim_prfl_p)) THEN

    ALLOCATE(coef2%lim_prfl_p(coef2%nlevels), stat=err)
    THROWM(err.NE.0, 'allocation of lim_prfl_p')
    coef2%lim_prfl_p = coef1%lim_prfl_p

    ALLOCATE(coef2%env_prfl_tmax(coef2%nlevels), stat=err)
    THROWM(err.NE.0, 'allocation of env_prfl_tmax')
    coef2%env_prfl_tmax = coef1%env_prfl_tmax

    ALLOCATE(coef2%env_prfl_tmin(coef2%nlevels), stat=err)
    THROWM(err.NE.0, 'allocation of env_prfl_tmin')
    coef2%env_prfl_tmin = coef1%env_prfl_tmin

    ALLOCATE(coef2%env_prfl_gmax(coef2%nlevels, coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of env_prfl_gmax')
    coef2%env_prfl_gmax = coef1%env_prfl_gmax

    ALLOCATE(coef2%env_prfl_gmin(coef2%nlevels, coef2%fmv_gas), stat=err)
    THROWM(err.NE.0, 'allocation of env_prfl_gmin')
    coef2%env_prfl_gmin = coef1%env_prfl_gmin

  ENDIF


  ! Thermal fast coefs

  CALL extract_fast_coef(err, coef1%thermal, coef2%thermal, .FALSE._jplm)
  THROW(err.NE.0)

  CALL extract_fast_coef(err, coef1%thermal_corr, coef2%thermal_corr, .TRUE._jplm)
  THROW(err.NE.0)


  ! Solar fast coefs

  CALL extract_fast_coef(err, coef1%solar, coef2%solar, .FALSE._jplm)
  THROW(err.NE.0)

  CALL extract_fast_coef(err, coef1%solar_corr, coef2%solar_corr, .TRUE._jplm)
  THROW(err.NE.0)


  ! NLTE coefs

  IF (ASSOCIATED(coef1%nlte_coef)) THEN

    ! For any monotonic channel selection we must find those selected channels
    ! which lie within the range of NLTE channels in the coef file. This
    ! constitutes another contiguous block of channels in the coef structure.

    ALLOCATE(nlte_chans(SIZE(channels))) ! Index of selected channels in nlte_coefs array in the file
    nlte_count = 0  ! Number of NLTE channels being read in
    nlte_start = 0  ! Index (in input channel list) of first NLTE channel being read in
    DO i = 1, SIZE(channels)
      IF (i > 1) THEN
        IF (channels(i) < channels(i-1)) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "non-monotonic channel selection incompatible with NLTE coefficients")
        ENDIF
      ENDIF
      IF (channels(i) >= coef1%nlte_coef%start_chan .AND. &
          channels(i) < coef1%nlte_coef%start_chan + coef1%nlte_coef%nchan) THEN
        nlte_count = nlte_count + 1
        nlte_chans(nlte_count) = channels(i) - coef1%nlte_coef%start_chan + 1
        IF (nlte_count == 1) nlte_start = i
      ENDIF
    ENDDO

    IF (nlte_count > 0) THEN
      ! We have some NLTE channels in the selection so go ahead and allocate space

      ALLOCATE(coef2%nlte_coef)

      coef2%nlte_coef%ncoef = coef1%nlte_coef%ncoef
      coef2%nlte_coef%nsol  = coef1%nlte_coef%nsol
      coef2%nlte_coef%nsat  = coef1%nlte_coef%nsat

      NULLIFY (coef2%nlte_coef%coef)
      NULLIFY (coef2%nlte_coef%sol_zen_angle, coef2%nlte_coef%cos_sol)
      NULLIFY (coef2%nlte_coef%sat_zen_angle, coef2%nlte_coef%sec_sat)

      coef2%nlte_coef%start_chan = nlte_start
      coef2%nlte_coef%nchan      = nlte_count
    ELSE

      ! No NLTE channels selected so no need for NLTE coef section
      coef2%nltecoef = .FALSE.

    ENDIF

    IF (coef2%nltecoef) THEN

      ALLOCATE(coef2%nlte_coef%sol_zen_angle(coef2%nlte_coef%nsol), stat=err)
      THROWM(err.NE.0, "allocation of NLTE solar zenith angle array")
      coef2%nlte_coef%sol_zen_angle = coef1%nlte_coef%sol_zen_angle

      ALLOCATE(coef2%nlte_coef%sec_sat(coef2%nlte_coef%nsat), stat=err)
      THROWM(err.NE.0, "allocation of NLTE satellite zenith angle array")
      coef2%nlte_coef%sec_sat = coef1%nlte_coef%sec_sat

      ALLOCATE(coef2%nlte_coef%coef(coef2%nlte_coef%ncoef, coef2%nlte_coef%nsat, &
               coef2%nlte_coef%nsol, coef2%nlte_coef%nchan), stat=err)
      THROWM(err.NE.0, "allocation of NLTE coef array")
      coef2%nlte_coef%coef(:,:,:,:) = coef1%nlte_coef%coef(:,:,:,nlte_chans(1:nlte_count))

    ENDIF

    DEALLOCATE(nlte_chans)

  ENDIF


  ! PMC coefs

  IF (ASSOCIATED(coef1%pmc_coef)) THEN

    ALLOCATE(coef2%pmc_coef(coef2%pmc_nlay, coef2%fmv_chn, coef2%pmc_nvar), stat=err)
    THROWM(err.NE.0, 'allocation of pmc_coef')
    coef2%pmc_coef(:,:,:) = coef1%pmc_coef(:,channels,:)

    ALLOCATE(coef2%pmc_pnominal(coef2%fmv_chn), stat=err)
    THROWM(err.NE.0, 'allocation of pmc_pnominal')
    coef2%pmc_pnominal = coef1%pmc_pnominal(channels)

  ENDIF

CATCH

CONTAINS

  SUBROUTINE nullify_gas_coef_pointers(fast_coef)
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

  SUBROUTINE extract_fast_coef(err, fast_coef1, fast_coef2, corr)
    INTEGER(jpim),                  INTENT(OUT)   :: err
    TYPE(rttov_fast_coef), POINTER, INTENT(IN)    :: fast_coef1(:)
    TYPE(rttov_fast_coef), POINTER, INTENT(INOUT) :: fast_coef2(:)
    LOGICAL(jplm),                  INTENT(IN)    :: corr

    INTEGER(jpim) :: ichan, igas, gas_id, gas_pos, ncoef

    TRY
    IF (ASSOCIATED(fast_coef1)) THEN
      ALLOCATE(fast_coef2(coef2%fmv_chn), stat=err)
      THROWM(err.NE.0, 'allocation of fast_coef')
      DO ichan = 1, coef2%fmv_chn
        ALLOCATE (fast_coef2(ichan)%gasarray(coef2%fmv_gas), stat=err)
        THROWM(err.NE.0, 'allocation of fast_coef(:)%gasarray')
        CALL nullify_gas_coef_pointers(fast_coef2(ichan))

        DO igas = 1, coef1%fmv_gas
          gas_id = coef1%fmv_gas_id(igas)
          gas_pos = coef1%fmv_gas_pos(gas_id)
          IF (corr) THEN
            ncoef = coef1%fmv_ncorr(gas_pos)
          ELSE
            ncoef = coef1%fmv_coe(gas_pos)
          ENDIF

          IF (ASSOCIATED(fast_coef1(channels(ichan))%gasarray(gas_pos)%coef)) THEN
            ALLOCATE(fast_coef2(ichan)%gasarray(gas_pos)%coef(ncoef, coef2%nlayers), stat=err)
            THROWM(err.NE.0, 'allocation of fast_coef(:)%gasarray(:)%coef')
            fast_coef2(ichan)%gasarray(gas_pos)%coef = &
              fast_coef1(channels(ichan))%gasarray(gas_pos)%coef

            CALL set_pointers(fast_coef2(ichan), gas_pos, gas_id)
          ENDIF
        ENDDO

      ENDDO
    ENDIF
    CATCH
  END SUBROUTINE

END SUBROUTINE rttov_channel_extract_coef
