! Description:
!> @file
!!   Write a v10/v11 ASCII coefficient file from the v12 rttov_coef structure.
!
!> @brief
!!   Write a v10/v11 ASCII coefficient file from the v12 rttov_coef structure.
!!
!! @details
!!   Where necessary this subroutine converts the v12 coefficient structure to
!!   v10/v11 format. Suitable defaults are chosen for values which no longer
!!   exist in the v12 coefficient structure.
!!
!!   NB This subroutine should not be called for variable SO2 coefficients.
!!
!! @param[out]    err             status on exit
!! @param[in]     file_id         logical unit for output rtcoef file
!! @param[in]     coef            RTTOV optical depth coefficient structure
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
SUBROUTINE rttov11_write_ascii_coef(err, file_id, coef)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :    &
         version,            &
         release,            &
         minor_version,      &
         gas_id_mixed,       &
         gas_name,           &
         tscale_v11,         &
         gscale_v11,         &
         lensection,         &
         speedl,             &
         sensor_id_ir,       &
         sensor_id_mw,       &
         sensor_id_hi,       &
         sensor_id_po
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err
  INTEGER(KIND=jpim), INTENT(IN)  :: file_id
  TYPE(rttov_coef)  , INTENT(IN)  :: coef
!INTF_END
#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim)  :: i, j, l, k
!   INTEGER(KIND=jpim)  :: isat, isol, ichan
  INTEGER(KIND=jpim)  :: IncZeeman
  CHARACTER(LEN=2)          :: sensor
  CHARACTER(LEN=lensection) :: section
  CHARACTER(LEN=80)         :: errMessage
  CHARACTER(LEN=80)         :: version_name
!   CHARACTER(LEN=20)         :: FMT_sol, FMT_sat,
  CHARACTER(LEN=20)         :: FMT_coef, FMT_pnom

  INTEGER(KIND=jpim), ALLOCATABLE :: tt_val_chn(:)
  REAL(KIND=jprb), ALLOCATABLE :: coef_temp(:,:,:), tt_a0(:), tt_a1(:)
  INTEGER(KIND=jpim) :: ncoef, gas_pos, gas_id
  CHARACTER(LEN=*), PARAMETER :: routinename = 'rttov11_write_ascii_coef'
!- End of header --------------------------------------------------------
  TRY

  IF (coef%fmv_gas > 8) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'Too many gases: cannot write variable SO2 files in v10/v11 format')
  ENDIF

  WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
  INFO(errMessage)

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV coefficient file '//TRIM(coef%id_Common_name)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Automatic creation by subroutine '//routinename
  THROW(err.NE.0)

  IF (release < 10 .and. minor_version < 10 ) THEN
    WRITE (version_name, '(I2.2,".",i1,".",i1)', iostat=err) version, release, minor_version
  ELSE
    WRITE (version_name, '(I2.2,".",i2.2,".",i2.2)', iostat=err) version, release, minor_version
  ENDIF

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV library version '//TRIM(version_name)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  section = 'IDENTIFICATION'

  SELECT CASE (coef%id_sensor)
  CASE (sensor_id_ir)
    sensor = 'ir'
  CASE (sensor_id_mw)
    sensor = 'mw'
  CASE (sensor_id_hi)
    sensor = 'hi'
  CASE (sensor_id_po)
    sensor = 'po'
  END SELECT

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(3i4,T20,a)', iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, '! Platform  sat_id  instrument'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_Common_name)
  THROW(err.NE.0)

  WRITE (file_id, '(1x,a,T20,a)', iostat=err)sensor, '! Sensor type [ir,mw,hi,po]'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i2,T20,a)', iostat=err)10, '! RTTOV coefficient file version number'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_creation)
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i4,1x,i2.2,1x,i2.2,t20,a)', iostat=err)coef%id_creation_date, '! Creation date'
  THROW(err.NE.0)


  IF (coef%line_by_line(1) .NE. 'xxxx') THEN
    section = 'LINE-BY-LINE'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Line-by-line and other information'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    DO i = 1, SIZE(coef%line_by_line)
      IF (coef%line_by_line(i) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)TRIM(coef%line_by_line(i))
      THROW(err.NE.0)
    ENDDO
  ENDIF


  section = 'FAST_MODEL_VARIABLES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' !'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,a,t20,a)', iostat=err)coef%fmv_model_def, '! Fast model name'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_model_ver, '! Fast model version compatibility level'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i6,t20,a)', iostat=err)coef%fmv_chn, '! Number of channels described in the coef file'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_gas, '! Number of gases described in the coef file'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,i4,t20,a)' , iostat=ERR)coef%id_comp_pc, '! PC compatibility level'
  THROW(err.NE.0)

  IncZeeman = 0_jpim
  IF (coef%IncZeeman) IncZeeman = 1_jpim
  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)IncZeeman, '! Zeeman flag'
  THROW(err.NE.0)

  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(1x,a,t20,a)', iostat=err)TRIM(gas_name(coef%fmv_gas_id(i))), '! Gas identification'
    THROW(err.NE.0)

    WRITE (file_id, '(1x,3i4,t20,a)', iostat=err)coef%fmv_var(i), coef%fmv_coe(i), coef%fmv_lvl(i), &
                    '! Variables/predictors  levels (pressure/absorber)'
    THROW(err.NE.0)
  ENDDO


  IF (coef%readme_srf(1) .NE. 'xxxx') THEN
    section = 'README_SPECTRAL_RESPONSE_FUNCTION'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    DO i = 1, SIZE(coef%readme_srf)
      IF (coef%readme_srf(i) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)TRIM(coef%readme_srf(i))
      THROW(err.NE.0)
    ENDDO
  ENDIF


  section = 'FILTER_FUNCTIONS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel number (from instrument original description)'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel status'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Central wavenumber'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Band correction coefficients (offset, slope)'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gamma correction factor'
  THROW(err.NE.0)

  DO i = 1, coef%fmv_chn
    WRITE (file_id, '(1x,i5,1x,i4,4(1x,e18.10))', iostat=err) &
           coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)
    THROW(err.NE.0)
  ENDDO


  section = 'FUNDAMENTAL_CONSTANTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Units of constants for spectral radiance'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! First radiation constant (mW/(m2.sr.cm-4))'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Second radiation constant (cm.K)'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,f14.1,t30,a)', iostat=err)speedl, '! Speed of light (cm/s)'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,1p,e16.9,0p,f11.7,t30,a)', iostat=err)coef%fc_planck_c1, coef%fc_planck_c2, '! Planck constants'
  THROW(err.NE.0)

  WRITE (file_id, '(1x,f10.1,t30,a)', iostat=err)coef%fc_sat_height, '! Nominal satellite height (km)'
  THROW(err.NE.0)


  IF (ASSOCIATED(coef%ss_solar_spectrum)) THEN
    IF (ANY(coef%ss_val_chn > 0)) THEN
      section = 'SOLAR_SPECTRUM'
      WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)TRIM(section)
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)' ! '
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)' ! Channel number'
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)' ! Type of channel (0 => thermal; 1 => thermal+solar; 2 => solar)'
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'
      THROW(err.NE.0)

      WRITE (file_id, '(a)', iostat=err)' ! Solar spectrum'
      THROW(err.NE.0)

      DO i = 1, coef%fmv_chn
        WRITE (file_id, '(2i5,2e19.10)', iostat=err)i, coef%ss_val_chn(i), coef%ff_cwn(i), coef%ss_solar_spectrum(i)
        THROW(err.NE.0)
      ENDDO
    ENDIF
  ENDIF


  IF (ANY(coef%pw_val_chn > 0)) THEN
    section = 'PLANCK_WEIGHTED'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Planck-weighted flag (1 => yes; 0 => no)'
    THROW(err.NE.0)

    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5)', iostat=err)i, coef%pw_val_chn(i)
      THROW(err.NE.0)
    ENDDO
  ENDIF


  section = 'GAZ_UNITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gas concentrations can be expressed in '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! volume mixing ratio (ppmv)'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! specific concentration (kg/kg)'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
    THROW(err.NE.0)

    WRITE (file_id, '(1x,i4,t20,"! ",a)', iostat=err)2, 'ppmv'  ! v11 value for ppmv dry
    THROW(err.NE.0)
  ENDDO


  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
    section = 'FASTEM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! S. English fast generic millimetre wave ocean emissivity model'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Polarisation of each channel', &
                ' !       MPOL=0: Average of vertical and horizontal polarisation ie 0.5(H+V)', &
                ' !       MPOL=1: Nominal vertical at nadir rotating with view angle QV', &
                ' !       MPOL=2: Nominal horizontal at nadir rotating with view angle QH', &
                ' !       MPOL=3: Vertical V', &
                ' !       MPOL=4: Horizontal H', &
                ' !       MPOL=5: +45 minus -45 (3rd stokes vector) S3', &
                ' !       MPOL=6: Left circular - right circular (4th stokes vector) S4'
    THROW(err.NE.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)5_jpim, '   ! Version number'
    THROW(err.NE.0)

    WRITE (file_id, '(20i3)', iostat=err)(coef%fastem_polar(i), i = 1, coef%fmv_chn)
    THROW(err.NE.0)
  ENDIF


  IF (coef%ssirem_ver >= 1) THEN
    section = 'SSIREM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! 5 coefficients for emissivity model SSIREM'
    THROW(err.NE.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)coef%ssirem_ver, '   ! Version number'
    THROW(err.NE.0)

    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(1x,i5,3f12.7,2f4.1)', iostat=err)i, coef%ssirem_a0(i), coef%ssirem_a1(i), &
             coef%ssirem_a2(i), coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)
      THROW(err.NE.0)
    ENDDO
  ENDIF


  IF (coef%fmv_model_ver >= 9) THEN
    section = 'TRANSMITTANCE_TRESHOLD'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Apply transmittance threshold (1 => yes; 0 => no)'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Transmittance threshold'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Transmittance value'
    THROW(err.NE.0)

    ALLOCATE(tt_val_chn(coef%fmv_chn), tt_a0(coef%fmv_chn), tt_a1(coef%fmv_chn), STAT=err)
    THROW(err.NE.0)
    tt_val_chn = 0
    tt_a0 = 0._jprb
    tt_a1 = 0._jprb
    IF (ASSOCIATED(coef%tt_val_chn)) tt_val_chn = coef%tt_val_chn
    IF (ASSOCIATED(coef%tt_a0)) tt_a0 = coef%tt_a0
    IF (ASSOCIATED(coef%tt_a1)) tt_a1 = coef%tt_a1

    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5,3e19.10)', iostat=err)     &
          i, tt_val_chn(i), coef%ff_cwn(i), tt_a0(i), tt_a1(i)
      THROW(err.NE.0)
    ENDDO
    DEALLOCATE(tt_val_chn, tt_a0, tt_a1)
  ENDIF


  IF (ASSOCIATED(coef%woc_waopc_ow)) THEN
    section = 'WATER_OPTICAL_CONSTANT'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Ocean water (real and imaginary part)'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Fresh water (real and imaginary part)'
    THROW(err.NE.0)

    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(i5,e19.10,2(" (",e17.10,",",e17.10,")"))', iostat=err) &
             i, coef%ff_cwn(i), coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)
      THROW(err.NE.0)
    ENDDO
  ENDIF


  section = 'REFERENCE_PROFILE'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Reference pressure (hPa), reference temperature (K) and'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! reference volume mixing ratio (ppmv) for each gas'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Note that mixing ratio is "missing" for mixed gases'
  THROW(err.NE.0)

  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
    THROW(err.NE.0)

    DO l = 1, coef%fmv_lvl(i)
      WRITE (file_id, '(1x,f9.4,2x,f7.3,1x,e13.6)')coef%ref_prfl_p(l), coef%ref_prfl_t(l, i), coef%ref_prfl_mr(l, i)
      THROW(err.NE.0)
    ENDDO
  ENDDO


  section = 'PROFILE_LIMITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Reference pressure (hPa), temperature max and min (K) and'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! volume mixing ratio max and min (ppmv) for each gas'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' !      Temperature'
  THROW(err.NE.0)

  DO l = 1, coef%fmv_lvl(1)
    WRITE (file_id, '(1x,f9.4,2(1x,f7.2))', iostat=err)&
      coef%lim_prfl_p(l), &
      coef%env_prfl_tmax(l) * (1._jprb + tscale_v11), &
      coef%env_prfl_tmin(l) * (1._jprb - tscale_v11)
    THROW(err.NE.0)
  ENDDO

  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
    THROW(err.NE.0)

    DO l = 1, coef%fmv_lvl(i)
      IF (coef%fmv_gas_id(i) == gas_id_mixed) THEN
        WRITE (file_id, '(1x,f9.4,2x,e12.4,e12.4)', iostat=err)     &
          coef%lim_prfl_p(l), &
          coef%env_prfl_gmax(l, i), &
          coef%env_prfl_gmin(l, i)
      ELSE
        WRITE (file_id, '(1x,f9.4,2x,e12.4,e12.4)', iostat=err)     &
          coef%lim_prfl_p(l), &
          coef%env_prfl_gmax(l, i) * (1._jprb + gscale_v11), &
          coef%env_prfl_gmin(l, i) * (1._jprb - gscale_v11)
      ENDIF
      THROW(err.NE.0)
    ENDDO
  ENDDO


  section = 'FAST_COEFFICIENTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! '
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Transmission coefficients'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)' ! Order of the gases:'
  THROW(err.NE.0)

  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
    THROW(err.NE.0)
  ENDDO

  DO l = 1, coef%fmv_gas
    IF (coef%fmv_var(l) > 0_jpim) THEN
      gas_id = coef%fmv_gas_id(l)
      gas_pos = coef%fmv_gas_pos(gas_id)
      ncoef = coef%fmv_coe(gas_pos)

      WRITE (file_id, '(a)', iostat=err) gas_name(gas_id)
      THROW(err.NE.0)

      ALLOCATE(coef_temp(coef%nlayers, coef%fmv_chn, ncoef), stat = err)
      THROW(err.NE.0)

      DO i = 1, coef%fmv_chn
        IF (ASSOCIATED(coef%thermal(i)%gasarray(gas_pos)%coef)) THEN
          coef_temp(:,i,:) = TRANSPOSE(coef%thermal(i)%gasarray(gas_pos)%coef)
        ELSE
          coef_temp(:,i,:) = 0._jprb
        ENDIF
      ENDDO

      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          (((coef_temp(i,j,k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, ncoef)
      THROW(err.NE.0)

      DEALLOCATE(coef_temp)
    ENDIF
  ENDDO


  IF (coef%solarcoef) THEN
    section = 'SOLAR_FAST_COEFFICIENTS'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Transmission coefficients'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Order of the gases:'
    THROW(err.NE.0)

    DO i = 1, coef%fmv_gas
      WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
      THROW(err.NE.0)
    ENDDO

    DO l = 1, coef%fmv_gas
      IF (coef%fmv_var(l) > 0_jpim) THEN
        gas_id = coef%fmv_gas_id(l)
        gas_pos = coef%fmv_gas_pos(gas_id)
        ncoef = coef%fmv_coe(gas_pos)

        WRITE (file_id, '(a)', iostat=err) gas_name(gas_id)
        THROW(err.NE.0)

        ALLOCATE(coef_temp(coef%nlayers, coef%fmv_chn, ncoef), stat = err)
        THROW(err.NE.0)

        DO i = 1, coef%fmv_chn
          IF (ASSOCIATED(coef%solar(i)%gasarray(gas_pos)%coef)) THEN
            coef_temp(:,i,:) = TRANSPOSE(coef%solar(i)%gasarray(gas_pos)%coef)
          ELSE
            coef_temp(:,i,:) = 0._jprb
          ENDIF
        ENDDO

        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
            (((coef_temp(i,j,k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, ncoef)
        THROW(err.NE.0)

        DEALLOCATE(coef_temp, stat = err)
        THROW(err.NE.0)
      ENDIF
    ENDDO
  ENDIF


  IF (ASSOCIATED(coef%ws_npoint)) THEN
    section = 'WAVE_SPECTRUM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! '
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Number of points'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Point number'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err)' ! Wave spectrum'
    THROW(err.NE.0)

    WRITE (file_id,  * , iostat=err)coef%ws_nomega
    THROW(err.NE.0)

    DO i = 1, coef%ws_nomega
      WRITE (file_id, '(f10.3,f12.5)', iostat=err)coef%ws_npoint(i), coef%ws_k_omega(i)
      THROW(err.NE.0)
    ENDDO
  ENDIF


!   IF (coef%nltecoef) THEN
!     section = 'NLTE_RADIANCE_COEFS'
!     WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'  
!     WRITE (file_id, '(a)', iostat=err) TRIM(section)  
!     WRITE (file_id, '(a)', iostat=err) ' ! '
!     WRITE (file_id, '(a)', iostat=err) ' ! Radiance coefficients'
!     WRITE (file_id, '(a)', iostat=err) ' ! ncoef nsol nsat nchan start_chan'
!     WRITE (file_id, '(a)', iostat=err) ' ! sol_zen_angle'
!     WRITE (file_id, '(a)', iostat=err) ' ! sat_zen_angle'
!     WRITE (file_id, '(a)', iostat=err) &
!       ' ! nsol * nsat * nchan * [co c1 c2] coefficients'
!     THROW(err .NE. 0)
! 
!     WRITE (file_id, '(5(i6))', iostat=err) coef%nlte_coef%ncoef, &
!       coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
!       coef%nlte_coef%start_chan
!     THROW(err .NE. 0)
! 
!     WRITE(FMT_sol,'(a,i1,a)') '(',coef%nlte_coef%nsol,'(f6.0))'
!     IF (coef%nlte_coef%nsat < 10) THEN
!       WRITE(FMT_sat,'(a,i1,a)') '(',coef%nlte_coef%nsat,'(f6.3))'
!     ELSE
!       WRITE(FMT_sat,'(a,i2,a)') '(',coef%nlte_coef%nsat,'(f6.3))'
!     ENDIF
!     WRITE(FMT_coef,'(a,i1,a)') '(',coef%nlte_coef%ncoef,'(e15.7))'
! 
!     WRITE (file_id, TRIM(FMT_sol), iostat=err) coef%nlte_coef%sol_zen_angle
!     THROW(err .NE. 0)
!     WRITE (file_id, TRIM(FMT_sat), iostat=err) coef%nlte_coef%sec_sat
!     THROW(err .NE. 0)
! 
!     DO ichan = 1, coef%nlte_coef%nchan
!       DO isol = 1, coef%nlte_coef%nsol
!         DO isat = 1, coef%nlte_coef%nsat
!           WRITE(file_id, TRIM(FMT_coef), iostat=ERR) &
!             coef%nlte_coef%coef(:, isat, isol, ichan)
!           THROW(err .NE. 0)
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDIF


  IF (coef%pmc_shift) THEN
    section = 'PRESSURE_MODULATED_CELL'
    WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'  
    WRITE (file_id, '(a)', iostat=err) TRIM(section)  
    WRITE (file_id, '(a)', iostat=err) ' ! '

    WRITE (file_id, '(a)', iostat=err) ' ! Cell length (cm)' 
    WRITE (file_id, '(a)', iostat=err) ' ! Nominal cell pressures (hPa) - as used for the above fast coefficients' 
    WRITE (file_id, '(a)', iostat=err) ' ! Temperature of cell (K) - fixed'
    WRITE (file_id, '(a)', iostat=err) ' ! Ratio gamma(co2)/gamma(air) - band-averaged co2 halfwidth ratio'
    WRITE (file_id, '(a)', iostat=err) ' ! Number of layers used  - may be a subset of coef%nlayers'
    WRITE (file_id, '(a)', iostat=err) ' ! Number of variables used'
    WRITE (file_id, '(a)', iostat=err) ' ! Coefficients - lev, ((coef(lev,chan,var), var=1,nvar), chan=1,nchan)'
    THROW(err .NE. 0)

    WRITE (file_id, *, iostat=err)  coef%pmc_lengthcell
    THROW(err .NE. 0)

    IF(coef%fmv_chn < 10) THEN  
      WRITE(FMT_pnom,'(a,i1,a)') '(',coef%fmv_chn,'(f6.2,1x))'
    ELSE
      WRITE(FMT_pnom,'(a,i2,a)') '(',coef%fmv_chn,'(f6.2,1x))'
    ENDIF
    WRITE (file_id, TRIM(FMT_pnom), iostat=err) (coef%pmc_pnominal(i),i=1,coef%fmv_chn)
    THROW(err .NE. 0)
    WRITE (file_id, *, iostat=err)  coef%pmc_tempcell
    THROW(err .NE. 0)
    WRITE (file_id, *, iostat=err)  coef%pmc_betaplus1
    THROW(err .NE. 0)
    WRITE (file_id, *, iostat=err)  coef%pmc_nlay
    THROW(err .NE. 0)
    WRITE (file_id, *, iostat=err)  coef%pmc_nvar
    THROW(err .NE. 0)
    IF(coef%fmv_chn < 10) THEN  
      WRITE(FMT_coef,'(a,i1,a)') '(',coef%pmc_nvar*coef%fmv_chn,'(e17.8))'
    ELSE
      WRITE(FMT_coef,'(a,i2,a)') '(',coef%pmc_nvar*coef%fmv_chn,'(e17.8))'
    ENDIF
    DO i = 1, coef%pmc_nlay
      WRITE (file_id, TRIM(FMT_coef), iostat=ERR)     &
          ((coef%pmc_coef(i, j, k), k = 1, coef%pmc_nvar), j = 1, coef%fmv_chn)
      THROW(err .NE. 0)
    ENDDO
  ENDIF


  section = 'END'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err)TRIM(section)
  THROW(err.NE.0)

  INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
