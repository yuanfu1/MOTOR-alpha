! Description:
!> @file
!!   Print information about an RTTOV optical depth coefficient file.
!
!> @brief
!!   Print some basic information about RTTOV and, optionally
!!   also print information about the supplied coefs structure.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!!   If supplied the verbose flag prints out additional
!!   information about every channel in the coefficient file.
!!
!! @param[in]   coefs     coefficients structure, optional
!! @param[in]   lu        logical unit for output, optional
!! @param[in]   text      additional text to print, optional
!! @param[in]   verbose   enable verbose output, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_print_info(coefs, lu, text, verbose)

  USE rttov_types, ONLY : rttov_coefs
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_global, ONLY : error_unit
  USE rttov_const, ONLY : &
          version,        &
          release,        &
          minor_version,  &
          platform_name,  &
          inst_name,      &
          sensor_name,    &
          sensor_id_ir,   &
          sensor_id_hi,   &
          gas_name
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_coefs),  INTENT(IN), OPTIONAL :: coefs  ! coefs structure
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lu     ! logical unit for print
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: text   ! additional text to print
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: verbose
!INTF_END

  INTEGER(KIND=jpim)  :: iu        ! logical unit for print
  CHARACTER(LEN=128)  :: tmp_text  ! temporary string for formatting output
  INTEGER(KIND=jpim)  :: int_kind  ! Standard integer
  REAL(KIND=jprb)     :: real_kind ! Standard real
  INTEGER(KIND=jpim)  :: c, l, g, gas_id ! Indexes
  LOGICAL(KIND=jplm)  :: lverbose

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  lverbose = .TRUE.
  IF (PRESENT(verbose)) lverbose = verbose

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV information: ",TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV information"
  ENDIF

  WRITE(tmp_text,'(i2,".",i1,".",i1)') version, release, minor_version
  WRITE(iu,'(2x,a,a)')    "RTTOV library version        ", TRIM(tmp_text)
  WRITE(tmp_text,'(es19.9)') HUGE(real_kind)
  WRITE(iu,'(2x,a,a)')    "Maximum real size            ", TRIM(ADJUSTL(tmp_text))
  WRITE(tmp_text,'(i25)') HUGE(int_kind)
  WRITE(iu,'(2x,a,a)')    "Maximum integer size         ", TRIM(ADJUSTL(tmp_text))

  IF (PRESENT(coefs)) THEN
    WRITE(iu,'(/,a)') "Coefficient file"
    WRITE(tmp_text,'(i4.4,"/",i2.2,"/",i2.2)') coefs%coef%id_creation_date(1), &
                                               coefs%coef%id_creation_date(2), &
                                               coefs%coef%id_creation_date(3)
    WRITE(iu,'(2x,a,a)')  "Creation date                ", TRIM(tmp_text)
    WRITE(tmp_text,'(i2)') coefs%coef%id_comp_lvl
    WRITE(iu,'(2x,a,a)')  "Coefficient file version     ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(a)',advance='yes')

    WRITE(tmp_text,'(a," - ",i2)') TRIM(platform_name(coefs%coef%id_platform)), coefs%coef%id_sat
    WRITE(iu,'(2x,a,a)')  "Platform                     ", TRIM(tmp_text)
    WRITE(iu,'(2x,a,a)')  "Instrument                   ", TRIM(inst_name(coefs%coef%id_inst))
    WRITE(iu,'(2x,a,a)')  "Sensor type                  ", TRIM(sensor_name(coefs%coef%id_sensor))
    WRITE(tmp_text,'(f9.2)') coefs%coef%fc_sat_height
    WRITE(iu,'(2x,a,a)')  "Nominal satellite height     ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i5)') coefs%coef%fmv_ori_nchn
    WRITE(iu,'(2x,a,a)')  "Number of channels in file   ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i5)') coefs%coef%fmv_chn
    WRITE(iu,'(2x,a,a)')  "Number of channels extracted ", TRIM(ADJUSTL(tmp_text))
    WRITE(iu,'(2x,a,l1)') "Zeeman coefficient file?     ", coefs%coef%inczeeman
    IF (coefs%coef%id_comp_pc == 0) THEN
      WRITE(iu,'(2x,a,l1)') "PC compatibility?            ", .FALSE._jplm
    ELSE
      WRITE(iu,'(2x,a,l1)') "PC compatibility?            ", .TRUE._jplm
      WRITE(tmp_text,'(i3)') coefs%coef%id_comp_pc
      WRITE(iu,'(2x,a,a)')  "PC compatibility version     ", TRIM(ADJUSTL(tmp_text))
    ENDIF
    IF (coefs%coef%fmv_model_ver >= 9) THEN
      WRITE(iu,'(2x,a,l1)')   "Solar compatibility?         ", .TRUE._jplm
      WRITE(iu,'(2x,a,l1)')   "Separate solar fast coefs?   ", coefs%coef%solarcoef
    ELSE
      WRITE(iu,'(2x,a,l1)')   "Solar compatibility?         ", .FALSE._jplm
    ENDIF
    WRITE(iu,'(a)',advance='yes')

    WRITE(tmp_text,'(i2)') coefs%coef%fmv_model_ver
    WRITE(iu,'(2x,a,a)')  "Fast model version           ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') coefs%coef%nlevels
    WRITE(iu,'(2x,a,a)')  "Number of levels             ", TRIM(ADJUSTL(tmp_text))
    WRITE(tmp_text,'(i3)') coefs%coef%fmv_gas
    WRITE(iu,'(2x,a,a)')  "Number of gases              ", TRIM(ADJUSTL(tmp_text))
    IF (coefs%coef%fmv_model_ver > 9) THEN
      WRITE(iu,'(2x,3a12,a19)')  "Gas name    ", "# predictors"," # gas coefs"," # correction coefs"
      DO g = 1, coefs%coef%fmv_gas
        gas_id = coefs%coef%fmv_gas_id(g)
        WRITE(iu,'(2x,a12,i8,4x,i8,7x,i8)') &
          gas_name(gas_id), coefs%coef%fmv_var(g), coefs%coef%fmv_coe(g), coefs%coef%fmv_ncorr(g)
      END DO
    ELSE
      WRITE(iu,'(2x,3a12)')  "Gas name    ", "# predictors","  # coeffs  "
      DO g = 1, coefs%coef%fmv_gas
        gas_id = coefs%coef%fmv_gas_id(g)
        WRITE(iu,'(2x,a12,i8,4x,i8)') &
          gas_name(gas_id), coefs%coef%fmv_var(g), coefs%coef%fmv_coe(g)
      END DO
    ENDIF
    WRITE(iu,'(a)',advance='yes')

    WRITE(iu,'(2x,a)')    "Line-By-Line information"
    DO l = 1, SIZE(coefs%coef%line_by_line)
      IF (coefs%coef%line_by_line(l) .EQ. 'xxxx') EXIT
      IF (coefs%coef%line_by_line(l) .EQ. '') EXIT
      WRITE(iu,'(4x,a)') TRIM(coefs%coef%line_by_line(l))
    ENDDO
    WRITE(iu,'(a)',advance='yes')

    IF (coefs%coef%readme_srf(1)(1:4) .NE. "xxxx") THEN
      WRITE(iu,'(2x,a)')    "Spectral Response Function information"
      DO l = 1, SIZE(coefs%coef%readme_srf)
        IF (coefs%coef%readme_srf(l) .EQ. 'xxxx') EXIT
        WRITE (iu, '(4x,a)') TRIM(coefs%coef%readme_srf(l))
      ENDDO
    ELSE
      WRITE (iu, '(2x,a)') "No README_SPECTRAL_RESPONSE_FUNCTION section"
    ENDIF
    WRITE(iu,'(a)',advance='yes')

    IF (coefs%coef%id_sensor == sensor_id_ir .OR. &
        coefs%coef%id_sensor == sensor_id_hi) THEN
      IF (coefs%coef%ssirem_ver > 0) THEN
        WRITE(tmp_text,'(2x,i3)') coefs%coef%ssirem_ver
        WRITE(iu,'(2x,a,a)')  "SSIREM model version ",TRIM(ADJUSTL(tmp_text))
      ENDIF
      IF (coefs%coef%iremis_version > 0) THEN
        WRITE(tmp_text,'(2x,i3)') coefs%coef%iremis_version
        WRITE(iu,'(2x,a,a)')  "IREMIS model version ",TRIM(ADJUSTL(tmp_text))
      ENDIF
    ELSE
      WRITE(tmp_text,'(a,i3.3,a)') "(2x,",coefs%coef%fmv_chn,"i3)"
      WRITE(iu,'(2x,a)')    "Channel polarisations"
      WRITE(iu,tmp_text) coefs%coef%fastem_polar(:)
    ENDIF
    WRITE(iu,'(a)',advance='yes')

    IF (lverbose) THEN
      IF (ASSOCIATED(coefs%coef%ss_val_chn)) THEN
        WRITE(iu,'(2x,a)')    "Thermal/Solar and Planck-weighted flags"
        DO c = 1, coefs%coef%fmv_chn
          IF (coefs%coef%ss_val_chn(c) .EQ. 0_jpim) THEN
            WRITE(iu,'(2x,i5,f10.3," cm-1    th  ")',advance='no')  c,coefs%coef%ff_cwn(c)
          ELSE IF (coefs%coef%ss_val_chn(c) .EQ. 1_jpim) THEN
            WRITE(iu,'(2x,i5,f10.3," cm-1  th+sol")',advance='no')  c,coefs%coef%ff_cwn(c)
          ELSE IF (coefs%coef%ss_val_chn(c) .EQ. 2_jpim) THEN
            WRITE(iu,'(2x,i5,f10.3," cm-1    sol ")',advance='no')  c,coefs%coef%ff_cwn(c)
          ENDIF
          IF (coefs%coef%pw_val_chn(c) > 0_jpim) THEN
            WRITE(iu,'(" pw ")',advance='no')
          ENDIF
          WRITE(iu,'(a)',advance='yes')
        ENDDO
      ENDIF
      WRITE(iu,'(a)',advance='yes')
    ENDIF

    IF (coefs%coef%nltecoef .AND. coefs%coef%id_sensor == sensor_id_hi) THEN
      WRITE(iu,'(2x,a)') "NLTE aware coefficient file"
      WRITE(iu,'(4x,a,i5,a,i5)') "Correction applied between channels ", &
        coefs%coef%nlte_coef%start_chan, " and ", &
        coefs%coef%nlte_coef%start_chan + coefs%coef%nlte_coef%nchan

      WRITE(iu,'(a)',advance='yes')
    ENDIF

    IF (coefs%coef%pmc_shift) THEN
      WRITE(iu,'(2x,a)') "PMC shift aware Coefficient File"

      WRITE(iu,'(4x,a,f7.3)') "Cell length (cm):     ", coefs%coef%pmc_lengthcell
      WRITE(iu,'(4x,a,f7.3)') "Cell temperature (K): ", coefs%coef%pmc_tempcell
      WRITE(iu,'(4x,a)') "Nominal cell pressures"
      DO c = 1, coefs%coef%fmv_chn
        WRITE(iu,'(4x,i5,f7.3," hPa")') c,coefs%coef%pmc_pnominal(c)
      ENDDO

      WRITE(iu,'(a)',advance='yes')
    ENDIF

  ENDIF

END SUBROUTINE rttov_print_info
