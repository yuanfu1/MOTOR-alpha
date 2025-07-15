! Description:
!> @file
!!   Read a v10/v11 ASCII coefficient file into the v12 rttov_coef structure.
!
!> @brief
!!   Read a v10/v11 ASCII coefficient file into the v12 rttov_coef structure.
!!
!! @details
!!   There is no optional channel selection: do this before or after conversion.
!!
!! @param[out]    err             status on exit
!! @param[in]     file_id         logical unit for input rtcoef file
!! @param[in,out] coef            RTTOV optical depth coefficient structure
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
SUBROUTINE rttov11_read_ascii_coef(err, file_id, coef)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb, jplm
  USE rttov_const, ONLY :  &
        sensor_id_hi,           &
        sensor_id_mw,           &
        sensor_id_ir,           &
        sensor_id_po,           &
        gas_mass,               &
        mair,                   &
        gas_unit_specconc,      &
        gas_id_mixed,           &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_wvcont,          &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        tscale_v11,             &
        gscale_v11,             &
        sensor_name,            &
        ngases_max,             &
        gas_name,               &
        lensection
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)   :: err
  INTEGER(KIND=jpim), INTENT(IN)    :: file_id
  TYPE(rttov_coef),   INTENT(INOUT) :: coef
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_deletecomment.interface"
#include "rttov_cmpuc.interface"
#include "rttov_findnextsection.interface"
#include "rttov_nullify_coef.interface"

  INTEGER(KIND=jpim) :: io_status
  REAL   (KIND=jprb) :: pres
  INTEGER(KIND=jpim) :: i, j, k, l, n, inczeeman
!   INTEGER(KIND=jpim) :: isat, isol, ichan
  INTEGER(KIND=jpim) :: itmp
  REAL   (KIND=jprb) :: rtmp
  REAL   (KIND=jprb) :: fc_speedl
  INTEGER(KIND=jpim) :: fastem_ver

  REAL   (KIND=jprb), POINTER :: coeffsarray(:, :, :)
  INTEGER(KIND=jpim), POINTER :: gaz_units(:) => NULL()

  CHARACTER(LEN=36)         :: input_string
  CHARACTER(LEN=32)         :: gas_type
  CHARACTER(LEN=lensection) :: section
  LOGICAL(KIND=jplm)        :: found

!- End of header --------------------------------------------------------
  TRY

  CALL rttov_nullify_coef(coef)

  coef%iremis_version = 0

  readfile : DO
    CALL rttov_findnextsection(file_id, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (TRIM(section))

    CASE ('IDENTIFICATION')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef%id_platform, coef%id_sat, coef%id_inst
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id, '(a)', iostat=err)coef%id_common_name
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)input_string
      THROWM(err.NE.0, 'io status while reading section '//section)

      SELECT CASE (input_string)
      CASE (sensor_name(sensor_id_ir))
        coef%id_sensor = sensor_id_ir
      CASE (sensor_name(sensor_id_mw))
        coef%id_sensor = sensor_id_mw
      CASE (sensor_name(sensor_id_hi))
        coef%id_sensor = sensor_id_hi
      CASE (sensor_name(sensor_id_po))
        coef%id_sensor = sensor_id_po
      CASE DEFAULT
        coef%id_sensor = sensor_id_ir
      END SELECT

      READ (file_id,  * , iostat=err)coef%id_comp_lvl
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (coef%id_comp_lvl < 10 .OR. coef%id_comp_lvl > 11) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "Version of coefficient file is incompatible with this subroutine (must be RTTOV v10/v11)")
      ENDIF
      coef%id_comp_lvl = 12

      READ (file_id, '(a)', iostat=err)coef%id_creation
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef%id_creation_date
      THROWM(err.NE.0, 'io status while reading section '//section)


    CASE ('README_SPECTRAL_RESPONSE_FUNCTION')
      CALL rttov_skipcommentline(file_id, err)
      LoopRSRF : DO i = 1, SIZE(coef%readme_srf)
        READ (file_id, '(a)', iostat=err)coef%readme_srf(i)
        DO j = 1, LEN(coef%readme_srf(i))

          SELECT CASE (coef%readme_srf(i)(j:j))
          CASE ('!')
            EXIT LoopRSRF
          CASE (' ')
            CYCLE
          CASE DEFAULT
            EXIT
          END SELECT

        ENDDO
      ENDDO LoopRSRF
      IF (i .LE. SIZE(coef%readme_srf)) coef%readme_srf(i) = 'xxxx'


    CASE ('LINE-BY-LINE')
      CALL rttov_skipcommentline(file_id, err)

      LoopLBL : DO i = 1, SIZE(coef%line_by_line)
        READ (file_id, '(a)', iostat=err)coef%line_by_line(i)

        DO j = 1, LEN(coef%line_by_line(i))

          SELECT CASE (coef%line_by_line(i)(j:j))
          CASE ('!')
            EXIT LoopLBL
          CASE (' ')
            CYCLE
          CASE DEFAULT
            EXIT
          END SELECT

        ENDDO
      ENDDO LoopLBL
      IF (i .LE. SIZE(coef%line_by_line)) coef%line_by_line(i) = 'xxxx'


    CASE ('FAST_MODEL_VARIABLES')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef%fmv_model_def
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef%fmv_model_ver
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%fmv_ori_nchn
      THROWM(err.NE.0,"io status while reading section "//section)

      coef%fmv_chn = coef%fmv_ori_nchn

      READ (file_id,  * , iostat=err)coef%fmv_gas
      THROWM(err.NE.0,"io status while reading section "//section)

      coef%id_comp_pc = 0_jpim
      READ (file_id,  * , iostat=err)coef%id_comp_pc
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)inczeeman
      THROWM(err.NE.0,"io status while reading section "//section)
      coef%inczeeman = (inczeeman == 1)

      ALLOCATE (coef%fmv_gas_id(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_gas_id")

      ALLOCATE (coef%fmv_gas_pos(ngases_max), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_gas_pos")

      ALLOCATE (coef%fmv_var(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_var")

      ALLOCATE (coef%fmv_lvl(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_lvl")

      ALLOCATE (coef%fmv_coe(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_coe")

      coef%fmv_gas_id(:)  = 0_jpim
      coef%fmv_gas_pos(:) = 0_jpim
      coef%fmv_var(:)     = 0_jpim
      coef%fmv_lvl(:)     = 0_jpim
      coef%fmv_coe(:)     = 0_jpim

      DO n = 1, coef%fmv_gas
! gas id. number i gas_id list (fmv_gas)
        READ (file_id, '(a)', iostat=err)gas_Type
        THROWM(err.NE.0,"io status while reading section "//section)

        CALL rttov_deletecomment(gas_Type)
        found = .FALSE.

        DO i = 1, ngases_max
          IF (rttov_cmpuc(gas_Type, gas_name(i))) THEN
            coef%fmv_gas_id(n) = i
            found = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF (.NOT. found) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0,Trim(gas_Type))
          RETURN
        ENDIF

! store also the index of this gas in the identification list
! so fmv_gas_pos(1) will give position of MxG in the file
        coef%fmv_gas_pos(coef%fmv_gas_id(n)) = n
! number of variables/predictors by gaz
! number of levels(pres/absorber) by gaz

        READ (file_id,  * , iostat=err)coef%fmv_var(n), coef%fmv_coe(n), coef%fmv_lvl(n)
        THROWM(err.NE.0,"io status while reading section "//section)

! Transfer information to some "classical" variables with more common names
! Note that the number of levels is taken from the Mixed Gases line

        SELECT CASE (coef%fmv_gas_id(n))
        CASE (gas_id_mixed)
          coef%nmixed  = coef%fmv_var(n)
          coef%ncmixed = coef%fmv_coe(n)
          coef%nlevels = coef%fmv_lvl(n)
        CASE (gas_id_watervapour)
          coef%nwater  = coef%fmv_var(n)
          coef%ncwater = coef%fmv_coe(n)
        CASE (gas_id_ozone)
          coef%nozone  = coef%fmv_var(n)
          coef%ncozone = coef%fmv_coe(n)
        CASE (gas_id_wvcont)
          coef%nwvcont  = coef%fmv_var(n)
          coef%ncwvcont = coef%fmv_coe(n)
        CASE (gas_id_co2)
          coef%nco2  = coef%fmv_var(n)
          coef%ncco2 = coef%fmv_coe(n)
        CASE (gas_id_n2o)
          coef%nn2o  = coef%fmv_var(n)
          coef%ncn2o = coef%fmv_coe(n)
        CASE (gas_id_co)
          coef%nco  = coef%fmv_var(n)
          coef%ncco = coef%fmv_coe(n)
        CASE (gas_id_ch4)
          coef%nch4  = coef%fmv_var(n)
          coef%ncch4 = coef%fmv_coe(n)
        END SELECT
      ENDDO

      coef%nlayers = coef%nlevels - 1

      ALLOCATE (gaz_units(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of gaz_units")
      gaz_units(:) = 1  ! v11 value for kg/kg


    CASE ('GAZ_UNITS')
! the array has already been allocated and initialised
! to specific concentration (kg/kg)
!
! This section needs one input line per gaz
! in the same order as the gaz list defined inside
!
! This is defining the units used for the sections
! REFERENCE_PROFILE and PROFILE_LIMITS
!
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

        READ (file_id,  * , iostat=err)gaz_units(n)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO

!-------------------------------------------------------


    CASE ('FILTER_FUNCTIONS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ff_ori_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_ori_chn")

      ALLOCATE (coef%ff_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_val_chn")

      ALLOCATE (coef%ff_cwn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_cwn")

      ALLOCATE (coef%ff_bco(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_bco")

      ALLOCATE (coef%ff_bcs(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_bcs")

      ALLOCATE (coef%ff_gam(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ff_gam")

      DO i = 1, coef%fmv_chn
      coef%ff_cwn(i) = -999.0_JPRB
        READ (file_id,  * , iostat=err) &
              coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('TRANSMITTANCE_TRESHOLD')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%tt_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_val_chn")

      ALLOCATE (coef%tt_a0(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_a0")

      ALLOCATE (coef%tt_a1(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_a1")

      DO i = 1, coef%fmv_chn
! chan number
! validity of channel
! central wave number
! transmittance treshold
! transmittance value
        READ (file_id,  * , iostat=err) itmp, coef%tt_val_chn(i), rtmp, coef%tt_a0(i), coef%tt_a1(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('PLANCK_WEIGHTED')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%pw_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of pw_val_chn")

      DO i = 1, coef%fmv_chn
        READ (file_id,  * , iostat=err)itmp, coef%pw_val_chn(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('SOLAR_SPECTRUM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ss_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ss_val_chn")

      ALLOCATE (coef%ss_solar_spectrum(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ss_solar_spectrum")

      DO i = 1, coef%fmv_chn
! chan number
! validity of channel
! central wave number
! solar spectrum
        READ (file_id,  * , iostat=err)itmp, coef%ss_val_chn(i), rtmp, coef%ss_solar_spectrum(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('WATER_OPTICAL_CONSTANT')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of woc_waopc_ow")

      ALLOCATE (coef%woc_waopc_fw(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of woc_waopc_fw")

      DO i = 1, coef%fmv_chn
! chan number
! central wave number
! ocean water optical constants(real and imaginary part)
! fresh water optical constants(real and imaginary part)
        READ (file_id,  * , iostat=err)itmp, rtmp, coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('WAVE_SPECTRUM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%ws_nomega
      ALLOCATE (coef%ws_k_omega(coef%ws_nomega), STAT = err)
      THROWM(err.NE.0, "allocation of ws_k_omega")

      ALLOCATE (coef%ws_npoint(coef%ws_nomega), STAT = err)
      THROWM(err.NE.0, "allocation of ws_npoint")

      DO i = 1, coef%ws_nomega
        READ (file_id,  * , iostat=err)coef%ws_npoint(i), coef%ws_k_omega(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO


    CASE ('FUNDAMENTAL_CONSTANTS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)fc_speedl
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%fc_planck_c1, coef%fc_planck_c2
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%fc_sat_height
      THROWM(err.NE.0,"io status while reading section "//section)


    CASE ('FASTEM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err) fastem_ver
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of fastem_polar")

      READ (file_id,  * , iostat=err)(coef%fastem_polar(i), i = 1, coef%fmv_chn)
      THROWM(err.NE.0,"io status while reading section "//section)


!-------------------------------------------------------
    CASE ('SSIREM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%ssirem_ver
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ssirem_a0(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ssirem_a0")

      ALLOCATE (coef%ssirem_a1(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ssirem_a1")

      ALLOCATE (coef%ssirem_a2(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ssirem_a2")

      ALLOCATE (coef%ssirem_xzn1(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ssirem_xzn1")

      ALLOCATE (coef%ssirem_xzn2(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ssirem_xzn2")

      DO i = 1, coef%fmv_chn
! original chan number
! constant coef
! first order coef
! second order coef
! 1st exponent on zenith angle
! 2nd exponent on zenith angle
        READ (file_id,  * , iostat=err)itmp, coef%ssirem_a0(i), coef%ssirem_a1(i), coef%ssirem_a2(i), &
              coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO

!-------------------------------------------------------
    CASE ('REFERENCE_PROFILE')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ref_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_p")

      ALLOCATE (coef%ref_prfl_t(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_t")

      ALLOCATE (coef%ref_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_mr")

      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

        DO i = 1, coef%nlevels
          READ (file_id,  * , iostat=err)pres, coef%ref_prfl_t(i, n), coef%ref_prfl_mr(i, n)
          THROWM(err.NE.0,"io status while reading section "//section)
          IF (coef%fmv_gas_id(n) == gas_id_mixed) coef%ref_prfl_p(i) = pres
        ENDDO
      ENDDO

!-------------------------------------------------------
    CASE ('PROFILE_LIMITS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%lim_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = err)
      THROWM(err.NE.0, "allocation of lim_prfl_p")

      ALLOCATE (coef%env_prfl_tmax(coef%fmv_lvl(gas_id_mixed)), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_tmax")

      ALLOCATE (coef%env_prfl_tmin(coef%fmv_lvl(gas_id_mixed)), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_tmin")

      ALLOCATE (coef%env_prfl_gmin(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_gmin")

      ALLOCATE (coef%env_prfl_gmax(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_gmax")

      DO l = 1, coef%nlevels
! pressure  (hPa)       (levels)
! max temperature (K)   (levels)
! min temperature (K)   (levels)
        READ (file_id,  * , iostat=err)coef%lim_prfl_p(l), coef%env_prfl_tmax(l), coef%env_prfl_tmin(l)
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDDO

      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

        DO i = 1, coef%nlevels
! max specific concentration (kg/kg) (levels, gases)
! min specific concentration (kg/kg) (levels, gases)
! or
! max volume mixing r (ppmv) (levels, gases)
! min volume mixing r (ppmv) (levels, gases)
! according to
! units specified in GAZ_UNITS section (default is specific concentration (kg/kg))
          READ (file_id,  * , iostat=err)pres, coef%env_prfl_gmax(i, n), coef%env_prfl_gmin(i, n)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO
      ENDDO


!-------------------------------------------------------
    CASE ('FAST_COEFFICIENTS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%thermal(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of thermal fast coefs")
      DO i = 1, coef%fmv_chn
        ALLOCATE (coef%thermal(i)%gasarray(coef%fmv_gas), STAT = err)
        THROWM(err.NE.0, "allocation of gasarray")
        NULLIFY (coef%thermal(i)%mixedgas)
        NULLIFY (coef%thermal(i)%watervapour)
        NULLIFY (coef%thermal(i)%ozone)
        NULLIFY (coef%thermal(i)%wvcont)
        NULLIFY (coef%thermal(i)%co2)
        NULLIFY (coef%thermal(i)%n2o)
        NULLIFY (coef%thermal(i)%co)
        NULLIFY (coef%thermal(i)%ch4)
      ENDDO

! loop on gases
      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

! read dummy string of gas name or filename of the sub_coefficient file
        READ (file_id,  * , iostat=err)input_string
        THROWM(err.NE.0,"io status while reading section "//section)

        ! test existence of gas
        IF (coef%fmv_coe(n) > 0_jpim) THEN
          CALL read_ascii_fast_coef(err, coef%thermal, n, coef%fmv_gas_id(n), coef%fmv_coe(n))
          THROWM(err.NE.0,"error occurred inside read_ascii_fast_coef while reading gas "//gas_name(n))
        ENDIF

      ENDDO


!-------------------------------------------------------
    CASE ('SOLAR_FAST_COEFFICIENTS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

! For instruments with solar-affected Planck-weighted channels, a second coefficient section
! is present. For pure-thermal channels these coefs will be zero. For all solar-affected channels
! these coefs will be non-PW. If this section is not present, the solar coefs structure will
! point to the thermal coefs.

      ALLOCATE (coef%solar(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of solar fast coefs")
      DO i = 1, coef%fmv_chn
        ALLOCATE (coef%solar(i)%gasarray(coef%fmv_gas), STAT = err)
        THROWM(err.NE.0, "allocation of gasarray")
        NULLIFY (coef%solar(i)%mixedgas)
        NULLIFY (coef%solar(i)%watervapour)
        NULLIFY (coef%solar(i)%ozone)
        NULLIFY (coef%solar(i)%wvcont)
        NULLIFY (coef%solar(i)%co2)
        NULLIFY (coef%solar(i)%n2o)
        NULLIFY (coef%solar(i)%co)
        NULLIFY (coef%solar(i)%ch4)
      ENDDO

      coef%solarcoef = .TRUE.   ! Solar coefs present

! loop on gases
      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

! read dummy string of gas name
        READ (file_id,  * , iostat=err)input_string
        THROWM(err.NE.0,"io status while reading section "//section)

        ! test existence of gas
        IF (coef%fmv_coe(n) > 0_jpim) THEN
          CALL read_ascii_fast_coef(err, coef%solar, n, coef%fmv_gas_id(n), coef%fmv_coe(n))
          THROWM(err.NE.0,"error occurred inside read_ascii_fast_coef while reading gas "//gas_name(n))
        ENDIF
      ENDDO


    CASE ('NLTE_RADIANCE_COEFS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err .NE. 0, "io status while reading section "//section)

      CYCLE readfile

!       coef%nltecoef = .TRUE. ! nlte coef present
!       ALLOCATE (coef%nlte_coef, STAT = err)
! 
!       READ(file_id, *, iostat=err) coef%nlte_coef%ncoef, coef%nlte_coef%nsol, &
!                                    coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
!                                    coef%nlte_coef%start_chan
!       THROWM(err.NE.0, "io status while reading NLTE dimensions")
! 
!       NULLIFY (coef%nlte_coef%coef)
!       NULLIFY (coef%nlte_coef%sol_zen_angle, coef%nlte_coef%cos_sol)
!       NULLIFY (coef%nlte_coef%sat_zen_angle, coef%nlte_coef%sec_sat)
! 
!       IF (coef%nltecoef) THEN
! 
!         ALLOCATE(coef%nlte_coef%sol_zen_angle(coef%nlte_coef%nsol), STAT = err)
!         THROWM(err.NE.0, "allocation of NLTE solar zenith angle array")
!         ALLOCATE(coef%nlte_coef%sec_sat(coef%nlte_coef%nsat), STAT = err)
!         THROWM(err.NE.0, "allocation of NLTE satellite zenith angle array")
!         ALLOCATE(coef%nlte_coef%coef(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
!                  coef%nlte_coef%nsol, coef%nlte_coef%nchan), STAT = err)
!         THROWM(err.NE.0, "allocation of NLTE coef array")
! 
!         READ(file_id,*,iostat=err) (coef%nlte_coef%sol_zen_angle(i), &
!           i = 1, coef%nlte_coef%nsol)
!         THROWM(err.NE.0, "io status while reading NLTE solar zenith angles")
! 
!         READ(file_id,*,iostat=err)(coef%nlte_coef%sec_sat(i), &
!           i = 1, coef%nlte_coef%nsat)
!         THROWM(err.NE.0, "io status while reading NLTE satellite zenith angles")
! 
!         DO ichan = 1, coef%nlte_coef%nchan
!           DO isol = 1, coef%nlte_coef%nsol
!             DO isat = 1, coef%nlte_coef%nsat
!               READ(file_id,*,iostat=err) &
!                 (coef%nlte_coef%coef(i, isat, isol, ichan), &
!                 i = 1, coef%nlte_coef%ncoef)
!               THROWM(err.NE.0, "io status while reading NLTE coefs")
!             ENDDO
!           ENDDO
!         ENDDO
! 
!       ENDIF


    CASE ('PRESSURE_MODULATED_CELL')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      coef%pmc_shift = .TRUE. ! pmc coef present

      READ (file_id,  * , iostat=err)coef%pmc_lengthcell
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%pmc_pnominal(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of coef%pmc_pnominal array")
      READ (file_id,  * , iostat=err) (coef%pmc_pnominal(i),i=1,coef%fmv_chn)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%pmc_tempcell
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%pmc_betaplus1
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%pmc_nlay
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%pmc_nvar
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%pmc_coef(coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar), STAT = err)
      THROWM(err.NE.0, "allocation of coef%pmc_coef array")
      coeffsarray => coef%pmc_coef

      READ (file_id,  * , iostat=err) &
           (((coeffsarray(i, j, k), k = 1, coef%pmc_nvar), j = 1, coef%fmv_ori_nchn), i = 1, coef%pmc_nlay)
      THROWM(err.NE.0,"io status while reading section "//section)


    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile


  ! Convert limits read from file to profile envelope and take care of gas units (inc. background profile)
  ! Conversion of gas mixing ratio units (kg/kg wet -> ppmv dry) if needed
  DO n = 1, coef%fmv_gas
    IF (n /= gas_id_mixed) THEN
      IF (gaz_units(n) == gas_unit_specconc) THEN
        coef%ref_prfl_mr(:,n)   = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%ref_prfl_mr(:,n) / (1._jprb - coef%ref_prfl_mr(:,n))
        coef%env_prfl_gmin(:,n) = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%env_prfl_gmin(:,n) / (1._jprb - coef%env_prfl_gmin(:,n))
        coef%env_prfl_gmax(:,n) = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%env_prfl_gmax(:,n) / (1._jprb - coef%env_prfl_gmax(:,n))
      ENDIF
      coef%env_prfl_gmax(:,n) = coef%env_prfl_gmax(:,n) / (1._jprb + gscale_v11)
      coef%env_prfl_gmin(:,n) = coef%env_prfl_gmin(:,n) / (1._jprb - gscale_v11)
    ENDIF
  ENDDO

  coef%env_prfl_tmax = coef%env_prfl_tmax / (1._jprb + tscale_v11)
  coef%env_prfl_tmin = coef%env_prfl_tmin / (1._jprb - tscale_v11)

  ALLOCATE (coef%bkg_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = err)
  THROWM(err.NE.0, "allocation of bkg_prfl_mr")
  coef%bkg_prfl_mr = coef%ref_prfl_mr

  IF (ASSOCIATED(gaz_units)) DEALLOCATE(gaz_units)

  CATCH

CONTAINS
  SUBROUTINE read_ascii_fast_coef(err, fast_coef, gas_pos, gas_id, ncoef)

    USE rttov_types, ONLY : rttov_fast_coef
    USE rttov_fast_coef_utils_mod, ONLY : set_pointers, nullify_pointers

    INTEGER(jpim),         INTENT(OUT)          :: err
    TYPE(rttov_fast_coef), INTENT(INOUT)        :: fast_coef(:)
    INTEGER(jpim),         INTENT(IN)           :: ncoef
    INTEGER(jpim),         INTENT(IN)           :: gas_pos
    INTEGER(jpim),         INTENT(IN)           :: gas_id

    REAL(jprb), ALLOCATABLE :: coef_temp(:,:,:)
    INTEGER(jpim) :: i, j, k, chn

!This is transposed from pre-RTTOV 11.3 order
!Coefs will be transposed for more efficent access in memory
!but remain same on disk (for pre-RTTOV 11.3 coefs)

    ALLOCATE(coef_temp(ncoef, coef%nlayers, coef%fmv_ori_nchn), stat = err)
    IF (err.NE.0) RETURN

    READ (file_id,  * , iostat=err)     &
        (((coef_temp(k, i, j),     i = 1, coef%nlayers), &
      j = 1, coef%fmv_ori_nchn), &
      k = 1, ncoef)
    IF (err.NE.0) RETURN

    DO i = 1, coef%fmv_chn
      chn = i

      ! Allocate space only for non-zero coefs
      IF (ANY(coef_temp(:, :, chn) /= 0._jprb)) THEN
        ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef, coef%nlayers), STAT = err)
        IF (err.NE.0) RETURN

        fast_coef(i)%gasarray(gas_pos)%coef = coef_temp(:, :, chn)
        CALL set_pointers(fast_coef(i), gas_pos, gas_id)
      ELSE
        NULLIFY (fast_coef(i)%gasarray(gas_pos)%coef)
        CALL nullify_pointers(fast_coef(i), gas_id)
      ENDIF

    ENDDO

    DEALLOCATE(coef_temp, stat = err)

  END SUBROUTINE read_ascii_fast_coef

END SUBROUTINE rttov11_read_ascii_coef
