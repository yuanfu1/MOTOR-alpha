! Description:
!> @file
!!   Write a binary coefficient file.
!
!> @brief
!!   Write a binary coefficient file.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in]     file_id         logical unit for output rtcoef file
!! @param[in]     verbose         flag to switch verbose output on/off (default TRUE), optional
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
SUBROUTINE rttov_write_binary_coef(err, coef, file_id, verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
         rttov_magic_string, &
         rttov_magic_number, &
         sensor_id_mw,       &
         sensor_id_po
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)           :: err
  TYPE(rttov_coef),   INTENT(IN)            :: coef
  INTEGER(KIND=jpim), INTENT(IN)            :: file_id
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL  :: verbose
!INTF_END

#include "rttov_errorreport.interface"

  LOGICAL(KIND=jplm) :: section_present, data_present
  LOGICAL(KIND=jplm) :: lverbose
  INTEGER(KIND=jpim) :: i, k
  CHARACTER(LEN=80)  :: errMessage
  INTEGER(KIND=jpim) :: inczeeman

!- End of header --------------------------------------------------------
  TRY

  lverbose = .TRUE._jplm
  IF (PRESENT(verbose)) lverbose = verbose

  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")')file_id
    INFO(errMessage)
  END IF

! Write a string that could be displayed
! Write a real number to be able to check single/double precision
  WRITE (file_id, iostat=err)rttov_magic_string, rttov_magic_number
  THROW(err.NE.0)


  WRITE (file_id, iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, coef%id_sensor
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef%id_comp_lvl, coef%id_creation_date, coef%id_creation, coef%id_common_name
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef%fmv_model_def, coef%fmv_model_ver, coef%fmv_chn, coef%fmv_gas
  THROW(err.NE.0)

  IF (coef%fmv_model_ver > 9) THEN
    WRITE (file_id, iostat=err)coef%nlevels
    THROW(err.NE.0)
  ENDIF

  WRITE (file_id, iostat=err)coef%id_comp_pc
  THROW(err.NE.0)

  inczeeman = 0_jpim
  IF (coef%inczeeman) inczeeman = 1_jpim
  WRITE (file_id, iostat=err) inczeeman
  THROW(err.NE.0)

  IF (coef%fmv_model_ver <= 9) THEN
    WRITE (file_id, iostat=err)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var, coef%fmv_lvl
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef%fmv_coe
    THROW(err.NE.0)
  ELSE
    WRITE (file_id, iostat=err)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef%fmv_coe, coef%fmv_ncorr
    THROW(err.NE.0)
  ENDIF

!TRANSMITTANCE_TRESHOLD
  section_present = .FALSE.
  IF (ASSOCIATED(coef%tt_val_chn)) section_present = ANY(coef%tt_val_chn > 0)

  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%tt_val_chn, coef%tt_a0, coef%tt_a1
    THROW(err.NE.0)
  ENDIF


!SOLAR_SPECTRUM
  section_present = .FALSE.
  IF (ASSOCIATED(coef%ss_val_chn)) section_present = ANY(coef%ss_val_chn > 0)

  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    IF (coef%fmv_model_ver > 9) THEN
      WRITE (file_id, iostat=err)coef%ss_val_chn, coef%ss_solar_spectrum, coef%ss_rayleigh_ext
    ELSE
      WRITE (file_id, iostat=err)coef%ss_val_chn, coef%ss_solar_spectrum
    ENDIF
    THROW(err.NE.0)
  ENDIF


!WATER_OPTICAL_CONSTANT
  section_present = ASSOCIATED(coef%woc_waopc_ow)
  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%woc_waopc_ow, coef%woc_waopc_fw
    THROW(err.NE.0)
  ENDIF


!WAVE_SPECTRUM
  section_present = ASSOCIATED(coef%ws_k_omega)
  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%ws_nomega
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef%ws_k_omega, coef%ws_npoint
    THROW(err.NE.0)
  ENDIF


  WRITE (file_id, iostat=err)coef%ff_ori_chn, coef%ff_val_chn, coef%ff_cwn, coef%ff_bco, coef%ff_bcs, coef%ff_gam
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef%fc_planck_c1, coef%fc_planck_c2, coef%fc_sat_height
  THROW(err.NE.0)


  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
    WRITE (file_id, iostat=err)coef%fastem_polar
    THROW(err.NE.0)

    IF (ANY(coef%fastem_polar == 7)) THEN
      WRITE (file_id, iostat=err)coef%pol_phi
      THROW(err.NE.0)
    ENDIF
  ENDIF


  WRITE (file_id, iostat=err)coef%ssirem_ver, coef%iremis_version
  THROW(err.NE.0)

  IF (coef%ssirem_ver >= 1) THEN
    WRITE (file_id, iostat=err)     &
        coef%ssirem_a0, coef%ssirem_a1, coef%ssirem_a2, coef%ssirem_xzn1, coef%ssirem_xzn2
    THROW(err.NE.0)
  ENDIF

  IF (coef%iremis_version >= 1) THEN
    WRITE (file_id, iostat=err) coef%iremis_ncoef
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%iremis_angle0, coef%iremis_tskin0
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%iremis_coef
    THROW(err.NE.0)
  ENDIF


  WRITE (file_id, iostat=err)coef%ref_prfl_p, coef%ref_prfl_t, coef%ref_prfl_mr, coef%bkg_prfl_mr
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)     &
      coef%lim_prfl_p, coef%env_prfl_tmax, coef%env_prfl_tmin, coef%env_prfl_gmax, coef%env_prfl_gmin
  THROW(err.NE.0)


  DO k = 1, coef%fmv_gas
    IF (coef%fmv_model_ver > 9) THEN
      IF (coef%fmv_coe(k) > 0_jpim) THEN
        DO i = 1, coef%fmv_chn
          IF (ASSOCIATED(coef%thermal(i)%gasarray(k)%coef)) THEN
            data_present = .TRUE._jplm
            WRITE (file_id, iostat=err) data_present
            THROW(err.NE.0)
            WRITE (file_id, iostat=err) coef%thermal(i)%gasarray(k)%coef
            THROW(err.NE.0)
          ELSE
            data_present = .FALSE._jplm
            WRITE (file_id, iostat=err) data_present
            THROW(err.NE.0)
          ENDIF
        ENDDO

        DO i = 1, coef%fmv_chn
          ! Gas correction coefficients
          IF (coef%fmv_ncorr(k) > 0_jpim) THEN
            IF (ASSOCIATED(coef%thermal_corr(i)%gasarray(k)%coef)) THEN
              data_present = .TRUE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
              WRITE (file_id, iostat=err) coef%thermal_corr(i)%gasarray(k)%coef
              THROW(err.NE.0)
            ELSE
              data_present = .FALSE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ELSE
      IF (coef%fmv_var(k) > 0_jpim) THEN
        DO i = 1, coef%fmv_chn
          IF (ASSOCIATED(coef%thermal(i)%gasarray(k)%coef)) THEN
            data_present = .TRUE._jplm
            WRITE (file_id, iostat=err) data_present
            THROW(err.NE.0)
            WRITE (file_id, iostat=err) TRANSPOSE(coef%thermal(i)%gasarray(k)%coef) ! transpose because of rttov11 coefficient format
            THROW(err.NE.0)
          ELSE
            data_present = .FALSE._jplm
            WRITE (file_id, iostat=err) data_present
            THROW(err.NE.0)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO

  ! SOLAR_FAST_COEFFICIENTS

  WRITE (file_id, iostat=err) coef%solarcoef
  THROW(err.NE.0)
  IF (coef%solarcoef) THEN
    DO k = 1, coef%fmv_gas
      IF (coef%fmv_model_ver > 9) THEN
        IF (coef%fmv_coe(k) > 0_jpim) THEN
          DO i = 1, coef%fmv_chn
            IF (ASSOCIATED(coef%solar(i)%gasarray(k)%coef)) THEN
              data_present = .TRUE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
              WRITE (file_id, iostat=err) coef%solar(i)%gasarray(k)%coef
              THROW(err.NE.0)
            ELSE
              data_present = .FALSE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
            ENDIF
          ENDDO

          DO i = 1, coef%fmv_chn
            ! Gas correction coefficients
            IF (coef%fmv_ncorr(k) > 0_jpim) THEN
              IF (ASSOCIATED(coef%solar_corr(i)%gasarray(k)%coef)) THEN
                data_present = .TRUE._jplm
                WRITE (file_id, iostat=err) data_present
                THROW(err.NE.0)
                WRITE (file_id, iostat=err) coef%solar_corr(i)%gasarray(k)%coef
                THROW(err.NE.0)
              ELSE
                data_present = .FALSE._jplm
                WRITE (file_id, iostat=err) data_present
                THROW(err.NE.0)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ELSE
        IF (coef%fmv_var(k) > 0_jpim) THEN
          DO i = 1, coef%fmv_chn
            IF (ASSOCIATED(coef%solar(i)%gasarray(k)%coef)) THEN
              data_present = .TRUE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
              WRITE (file_id, iostat=err) TRANSPOSE(coef%solar(i)%gasarray(k)%coef) ! transpose because of rttov11 coefficient format
              THROW(err.NE.0)
            ELSE
              data_present = .FALSE._jplm
              WRITE (file_id, iostat=err) data_present
              THROW(err.NE.0)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF


!PLANCK_WEIGHTED
  section_present = ANY(coef%pw_val_chn > 0)
  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%pw_val_chn
    THROW(err.NE.0)
  ENDIF


!README_SRF
  section_present = coef%readme_srf(1)(1:4) .NE. "xxxx"
  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%readme_srf
    THROW(err.NE.0)
  ENDIF


!NLTE_RADIANCE_COEFS
  section_present = coef%nltecoef
  WRITE (file_id, iostat=err) section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err) coef%nlte_coef%ncoef, &
    coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
    coef%nlte_coef%start_chan
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%nlte_coef%sol_zen_angle
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%nlte_coef%sec_sat
    THROW(err.NE.0)
    WRITE(file_id, iostat=ERR) coef%nlte_coef%coef
    THROW(err.NE.0)
  ENDIF

!PRESSURE_MODULATED_CELL
  section_present = coef%pmc_shift
  WRITE (file_id, iostat=err) section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err) coef%pmc_lengthcell
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%pmc_pnominal
    THROW(err.NE.0)
    WRITE (file_id, iostat=err) coef%pmc_tempcell, coef%pmc_betaplus1, &
        coef%pmc_nlay, coef%pmc_nvar
    THROW(err.NE.0)
    DO i = 1, coef%fmv_chn
      WRITE (file_id, iostat=err)coef%pmc_coef(:, i, :)
      THROW(err.NE.0)
    ENDDO
  ENDIF

!LINE_BY_LINE
  section_present = coef%line_by_line(1)(1:4) .NE. "xxxx"
  WRITE (file_id, iostat=err)section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    WRITE (file_id, iostat=err)coef%line_by_line
    THROW(err.NE.0)
  ENDIF

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE rttov_write_binary_coef
