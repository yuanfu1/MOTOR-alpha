! Description:
!> @file
!!   Read an ASCII coefficient file, optionally extracting a subset of channels.
!
!> @brief
!!   Read an ASCII coefficient file, optionally extracting a subset of channels.
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
SUBROUTINE rttov_read_ascii_coef(err, coef, file_id, channels, lbl)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
        version_compatible_min, &
        version_compatible_max, &
        sensor_id_hi,           &
        sensor_id_mw,           &
        sensor_id_ir,           &
        sensor_id_po,           &
        gas_id_mixed,           &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_wvcont,          &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        gas_id_so2,             &
        sensor_name,            &
        ngases_max,             &
        gas_name,               &
        lensection
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)            :: err
  TYPE(rttov_coef),   INTENT(INOUT)          :: coef
  INTEGER(KIND=jpim), INTENT(IN)             :: file_id
  INTEGER(KIND=jpim), INTENT(IN),   OPTIONAL :: channels(:)
  LOGICAL(KIND=jplm), INTENT(IN),   OPTIONAL :: lbl
!INTF_END
#include "rttov_opencoeff.interface"
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_deletecomment.interface"
#include "rttov_cmpuc.interface"
#include "rttov_findnextsection.interface"
#include "rttov_nullify_coef.interface"

  LOGICAL(KIND=jplm) :: for_output
  INTEGER(KIND=jpim) :: file_id_coef
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: io_status
  REAL   (KIND=jprb) :: pres
  INTEGER(KIND=jpim) :: i, j, k, l, n, inczeeman, itmp
  INTEGER(KIND=jpim) :: isat, isol, ichan
  INTEGER(KIND=jpim) :: nlte_count, nlte_start, nlte_file_nchan
  INTEGER(KIND=jpim), ALLOCATABLE :: nlte_chans(:)
  COMPLEX(KIND=jprb), POINTER :: values_c_0   (:)
  COMPLEX(KIND=jprb), POINTER :: values_c_1   (:)
  REAL   (KIND=jprb), POINTER :: values0      (:)
  REAL   (KIND=jprb), POINTER :: values1      (:)
  REAL   (KIND=jprb), POINTER :: values2      (:)
  REAL   (KIND=jprb), POINTER :: values3      (:)
  REAL   (KIND=jprb), POINTER :: values4      (:)
  REAL   (KIND=jprb), POINTER :: iremis_values(:,:)
  REAL   (KIND=jprb), POINTER :: nlte_values  (:,:,:,:)
  INTEGER(KIND=jpim), POINTER :: ivalues0     (:)
  INTEGER(KIND=jpim), POINTER :: ivalues1     (:)
  REAL   (KIND=jprb), POINTER :: coeffsarray  (:, :, :)

  CHARACTER(LEN=36)         :: input_string
  CHARACTER(LEN=32)         :: gas_type
  CHARACTER(LEN=lensection) :: section
  LOGICAL(KIND=jplm)        :: found
  LOGICAL(KIND=jplm)        :: lbl1
!- End of header --------------------------------------------------------

  TRY

  lbl1 = .FALSE.
  IF(PRESENT(lbl)) lbl1 = lbl

  all_channels = .NOT. PRESENT(channels)

  CALL rttov_nullify_coef(coef)


  readfile : DO
    CALL rttov_findnextsection(file_id, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (TRIM(section))

    CASE ('IDENTIFICATION')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef%id_platform, coef%id_sat, coef%id_inst
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id, '(1x,a)', iostat=err)coef%id_Common_name
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)input_string
      THROWM(err.NE.0, 'io status while reading section '//section)


      SELECT CASE (input_string)
      CASE (sensor_name(sensor_id_ir))
        coef%id_sensor = sensor_id_ir! Infrared
      CASE (sensor_name(sensor_id_mw))
        coef%id_sensor = sensor_id_mw! Micro Wave
      CASE (sensor_name(sensor_id_hi))
        coef%id_sensor = sensor_id_hi! High resolution
      CASE (sensor_name(sensor_id_po))
        coef%id_sensor = sensor_id_po! Polarimetric
      CASE DEFAULT
        coef%id_sensor = sensor_id_ir
      END SELECT

      READ (file_id,  * , iostat=err)coef%id_comp_lvl
      THROWM(err.NE.0, 'io status while reading section '//section)

! Error if the compatibility version of the coefficient file
! is not in the range defined by the constant module

      IF (coef%id_comp_lvl < version_compatible_min .OR. &
          coef%id_comp_lvl > version_compatible_max) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "Version of coefficient file is incompatible with RTTOV library")
      ENDIF

      READ (file_id, '(1x,a)', iostat=err)coef%id_creation
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

! fast model variables definition
      READ (file_id,  * , iostat=err)coef%fmv_model_def
      THROWM(err.NE.0, 'io status while reading section '//section)

! fast model variables version
      READ (file_id,  * , iostat=err)coef%fmv_model_ver
      THROWM(err.NE.0,"io status while reading section "//section)

! number of channels stored
      READ (file_id,  * , iostat=err)coef%fmv_ori_nchn
      THROWM(err.NE.0,"io status while reading section "//section)

! number of levels
      IF (coef%fmv_model_ver > 9) THEN
        READ (file_id,  * , iostat=err)coef%nlevels
        THROWM(err.NE.0,"io status while reading section "//section)
      ENDIF

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

! number of gases in file
      READ (file_id,  * , iostat=err)coef%fmv_gas
      THROWM(err.NE.0,"io status while reading section "//section)

! Flag to check compatibility with principal component regression file
      coef%id_comp_pc = 0_jpim
      READ (file_id,  * , iostat=err)coef%id_comp_pc
      THROWM(err.NE.0,"io status while reading section "//section)

! Flag for Zeeman coefficients
      READ (file_id,  * , iostat=err)inczeeman
      THROWM(err.NE.0,"io status while reading section "//section)
      coef%inczeeman = (inczeeman == 1)

! allocate arrays of FAST_MODEL_VARIABLES section
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

      ALLOCATE (coef%fmv_ncorr(coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of fmv_ncorr")

      coef%fmv_gas_id(:)  = 0_jpim
      coef%fmv_gas_pos(:) = 0_jpim
      coef%fmv_var(:)     = 0_jpim
      coef%fmv_lvl(:)     = 0_jpim
      coef%fmv_coe(:)     = 0_jpim
      coef%fmv_ncorr(:)   = 0_jpim

      DO n = 1, coef%fmv_gas
! gas id. number i gas_id list (fmv_gas)
        READ (file_id, '(a)', iostat=err)gas_type
        THROWM(err.NE.0,"io status while reading section "//section)

        CALL Rttov_deletecomment(gas_type)
        found = .FALSE.

        DO i = 1, ngases_max
          IF (rttov_cmpuc(gas_type, gas_name(i))) THEN
            coef%fmv_gas_id(n) = i
            found              = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF (.NOT. found) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0,TRIM(gas_type))
          RETURN
        ENDIF

! store also the index of this gas in the identification list
! so fmv_gas_pos(1) will give position of MxG in the file
        coef%fmv_gas_pos(coef%fmv_gas_id(n)) = n
! number of variables/predictors by gas
! number of levels(pres/absorber) by gas
        IF (coef%fmv_model_ver <= 9) THEN
          READ (file_id,  * , iostat=err)coef%fmv_var(n), coef%fmv_coe(n), coef%nlevels
          THROWM(err.NE.0,"io status while reading section "//section)
        ELSE
          READ (file_id,  * , iostat=err)coef%fmv_var(n), coef%fmv_coe(n), coef%fmv_ncorr(n)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDIF

! Transfer information to some "classical" variables with more common names
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

      coef%fmv_lvl(:) = coef%nlevels
      coef%nlayers = coef%nlevels - 1


    CASE ('FILTER_FUNCTIONS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

!allocate FILTER_FUNCTIONS section  array size is fmv_chn
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


      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
          coef%ff_cwn(i) = -999.0_jprb
          READ (file_id,  * , iostat=err)     &
              coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iv0")

        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iv1")

        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v1")

        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v2")

        ALLOCATE (values3(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v3")

        DO i = 1, coef%fmv_ori_nchn
          READ (file_id,  * , iostat=err)ivalues0(i), ivalues1(i), values0(i), values1(i), values2(i), values3(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

        coef%ff_ori_chn(:) = ivalues0(channels(:))
        coef%ff_val_chn(:) = ivalues1(channels(:))
        coef%ff_cwn(:)     = values0(channels(:))
        coef%ff_bco(:)     = values1(channels(:))
        coef%ff_bcs(:)     = values2(channels(:))
        coef%ff_gam(:)     = values3(channels(:))
        DEALLOCATE (ivalues0, STAT = err)
        THROWM(err.NE.0, "deallocation of iv0")

        DEALLOCATE (ivalues1, STAT = err)
        THROWM(err.NE.0, "deallocation of iv1")

        DEALLOCATE (values0, STAT = err)
        THROWM(err.NE.0, "deallocation of v0")

        DEALLOCATE (values1, STAT = err)
        THROWM(err.NE.0, "deallocation of v1")

        DEALLOCATE (values2, STAT = err)
        THROWM(err.NE.0, "deallocation of v2")

        DEALLOCATE (values3, STAT = err)
        THROWM(err.NE.0, "deallocation of v3")
      ENDIF


    CASE ('TRANSMITTANCE_TRESHOLD')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%tt_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_val_chn")

      ALLOCATE (coef%tt_a0(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_a0")

      ALLOCATE (coef%tt_a1(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of tt_a1")

      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
! chan number
! validity of channel
! transmittance treshold
! transmittance value
          READ (file_id,  * , iostat=err)     &
              itmp, coef%tt_val_chn(i), coef%tt_a0(i), coef%tt_a1(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

      ELSE
        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iv1")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v1")

        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v2")

        DO i = 1, coef%fmv_ori_nchn
          READ (file_id,  * , iostat=err)itmp, ivalues1(i), values1(i), values2(i)

          THROWM(err.NE.0,"io status while reading section "//section)

        ENDDO

        coef%tt_val_chn(:) = ivalues1(channels(:))
        coef%tt_a0(:)      = values1(channels(:))
        coef%tt_a1(:)      = values2(channels(:))

        DEALLOCATE (ivalues1, STAT = err)
        THROWM(err.NE.0, "deallocation of iv1")

        DEALLOCATE (values1, STAT = err)
        THROWM(err.NE.0, "deallocation of v1")

        DEALLOCATE (values2, STAT = err)
        THROWM(err.NE.0, "deallocation of v2")
      ENDIF


    CASE ('PLANCK_WEIGHTED')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%pw_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of pw_val_chn")

      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
          READ (file_id,  * , iostat=err)itmp, coef%pw_val_chn(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

      ELSE

        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iv1")

        DO i = 1, coef%fmv_ori_nchn
          READ (file_id,  * , iostat=err)itmp, ivalues1(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

        coef%pw_val_chn(:)        = ivalues1(channels(:))

        DEALLOCATE (ivalues1, STAT = err)
        THROWM(err.NE.0, "deallocation of iv1")
      ENDIF


    CASE ('SOLAR_SPECTRUM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ss_val_chn(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ss_val_chn")

      ALLOCATE (coef%ss_solar_spectrum(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of ss_solar_spectrum")

      IF (coef%fmv_model_ver > 9) THEN
        ALLOCATE (coef%ss_rayleigh_ext(coef%fmv_chn), STAT = err)
        THROWM(err.NE.0, "allocation of ss_rayleigh_ext")
      ENDIF

      IF (all_channels) THEN

        IF (coef%fmv_model_ver > 9) THEN
          DO i = 1, coef%fmv_chn
            READ (file_id,  * , iostat=err)itmp, coef%ss_val_chn(i), coef%ss_solar_spectrum(i), coef%ss_rayleigh_ext(i)
            THROWM(err.NE.0,"io status while reading section "//section)
          ENDDO
        ELSE
          DO i = 1, coef%fmv_chn
            READ (file_id,  * , iostat=err)itmp, coef%ss_val_chn(i), coef%ss_solar_spectrum(i)
            THROWM(err.NE.0,"io status while reading section "//section)
          ENDDO
        ENDIF

      ELSE
        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iv1")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of v1")

        IF (coef%fmv_model_ver > 9) THEN
          ALLOCATE (values2(coef%fmv_ori_nchn), STAT = err)
          THROWM(err.NE.0, "allocation of v2")

          DO i = 1, coef%fmv_ori_nchn
            READ (file_id,  * , iostat=err)itmp, ivalues1(i), values1(i), values2(i)
            THROWM(err.NE.0,"io status while reading section "//section)
          ENDDO

          coef%ss_rayleigh_ext(:) = values2(channels(:))

          DEALLOCATE (values2, STAT = err)
          THROWM(err.NE.0, "deallocation of v2")
        ELSE
          DO i = 1, coef%fmv_ori_nchn
            READ (file_id,  * , iostat=err)itmp, ivalues1(i), values1(i)
            THROWM(err.NE.0,"io status while reading section "//section)
          ENDDO
        ENDIF

        coef%ss_val_chn(:)        = ivalues1(channels(:))
        coef%ss_solar_spectrum(:) = values1(channels(:))

        DEALLOCATE (ivalues1, STAT = err)
        THROWM(err.NE.0, "deallocation of iv1")

        DEALLOCATE (values1, STAT = err)
        THROWM(err.NE.0, "deallocation of v1")
      ENDIF


    CASE ('WATER_OPTICAL_CONSTANT')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of woc_waopc_ow")

      ALLOCATE (coef%woc_waopc_fw(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of woc_waopc_fw")

      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
! chan number
! ocean water optical constants(real and imaginary part)
! fresh water optical constants(real and imaginary part)
          READ (file_id,  * , iostat=err)itmp, coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

      ELSE
        ALLOCATE (values_c_0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values_c_0")

        ALLOCATE (values_c_1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values_c_1")

        DO i = 1, coef%fmv_ori_nchn
          READ (file_id,  * , iostat=err)itmp, values_c_0(i), values_c_1(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

        coef%woc_waopc_ow(:) = values_c_0(channels(:))
        coef%woc_waopc_fw(:) = values_c_1(channels(:))

        DEALLOCATE (values_c_0, STAT = err)
        THROWM(err.NE.0, "deallocation of values_c_0")

        DEALLOCATE (values_c_1, STAT = err)
        THROWM(err.NE.0, "deallocation of values_c_1")
      ENDIF


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

! Planck constants
      READ (file_id,  * , iostat=err)coef%fc_planck_c1, coef%fc_planck_c2
      THROWM(err.NE.0,"io status while reading section "//section)

! satellite nominal altitude (km)
      READ (file_id,  * , iostat=err)coef%fc_sat_height
      THROWM(err.NE.0,"io status while reading section "//section)

    CASE ('FASTEM')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of fastem_polar")

! polarisation of each channel
      IF (all_channels) THEN
        READ (file_id,  * , iostat=err)(coef%fastem_polar(i), i = 1, coef%fmv_chn)
        THROWM(err.NE.0,"io status while reading section "//section)

        IF (ANY(coef%fastem_polar == 7)) THEN
          CALL rttov_skipcommentline(file_id, err)
          THROWM(err.NE.0,"io status while reading section "//section)

          ALLOCATE (coef%pol_phi(coef%fmv_chn), STAT = err)
          THROWM(err.NE.0, "allocation of pol_phi")

          READ (file_id,  * , iostat=err)(coef%pol_phi(i), i = 1, coef%fmv_chn)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDIF
      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of ivalues0")

        READ (file_id,  * , iostat=err)(ivalues0(i), i = 1, coef%fmv_ori_nchn)
        THROWM(err.NE.0,"io status while reading section "//section)
        coef%fastem_polar(:) = ivalues0(channels(:))
        DEALLOCATE (ivalues0, STAT = err)
        THROWM(err.NE.0, "deallocation of ivalues0")

        IF (ANY(coef%fastem_polar == 7)) THEN
          CALL rttov_skipcommentline(file_id, err)
          THROWM(err.NE.0,"io status while reading section "//section)

          ALLOCATE (coef%pol_phi(coef%fmv_chn), STAT = err)
          THROWM(err.NE.0, "allocation of pol_phi")
          ALLOCATE (values0(coef%fmv_ori_nchn), STAT = err)
          THROWM(err.NE.0, "allocation of values0")

          READ (file_id,  * , iostat=err)(values0(i), i = 1, coef%fmv_ori_nchn)
          THROWM(err.NE.0,"io status while reading section "//section)
          coef%pol_phi(:) = values0(channels(:))
          DEALLOCATE (values0, STAT = err)
          THROWM(err.NE.0, "deallocation of values0")
        ENDIF
      ENDIF


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

      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
! original chan number
! constant coef
! first order coef
! second order coef
! 1st exponent on zenith angle
! 2nd exponent on zenith angle
          READ (file_id,  * , iostat=err)itmp, coef%ssirem_a0(i), coef%ssirem_a1(i), coef%ssirem_a2(i),      &
              coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

      ELSE
        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values1")

        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values2")

        ALLOCATE (values3(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values3")

        ALLOCATE (values4(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values4")

        DO i = 1, coef%fmv_ori_nchn
          READ (file_id,  * , iostat=err)itmp, values0(i), values1(i), values2(i), values3(i), values4(i)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO

        coef%ssirem_a0(:)   = values0(channels(:))
        coef%ssirem_a1(:)   = values1(channels(:))
        coef%ssirem_a2(:)   = values2(channels(:))
        coef%ssirem_xzn1(:) = values3(channels(:))
        coef%ssirem_xzn2(:) = values4(channels(:))

        DEALLOCATE (values0, STAT = err)
        THROWM(err.NE.0, "deallocation of values0")

        DEALLOCATE (values1, STAT = err)
        THROWM(err.NE.0, "deallocation of values1")

        DEALLOCATE (values2, STAT = err)
        THROWM(err.NE.0, "deallocation of values2")

        DEALLOCATE (values3, STAT = err)
        THROWM(err.NE.0, "deallocation of values3")

        DEALLOCATE (values4, STAT = err)
        THROWM(err.NE.0, "deallocation of values4")
      ENDIF


!-------------------------------------------------------
    CASE ('IR_SEA_EMIS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%iremis_version
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%iremis_ncoef
      THROWM(err.NE.0,"io status while reading section "//section)

      READ (file_id,  * , iostat=err)coef%iremis_angle0, coef%iremis_tskin0
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%iremis_coef(coef%iremis_ncoef,coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of iremis_coef")

      IF (all_channels) THEN
        READ (file_id,  * , iostat=err)coef%iremis_coef(:,:)
        THROWM(err.NE.0,"io status while reading section "//section)
      ELSE
        ALLOCATE (iremis_values(coef%iremis_ncoef,coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of iremis_values")

        READ (file_id,  * , iostat=err)iremis_values
        THROWM(err.NE.0,"io status while reading section "//section)
        coef%iremis_coef(:,:) = iremis_values(:,channels(:))

        DEALLOCATE (iremis_values, STAT = err)
        THROWM(err.NE.0, "deallocation of iremis_values")
      ENDIF


!-------------------------------------------------------
    CASE ('REFERENCE_PROFILE')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%ref_prfl_p(coef%nlevels), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_p")

      ALLOCATE (coef%ref_prfl_t(coef%nlevels, coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_t")

      ALLOCATE (coef%ref_prfl_mr(coef%nlevels, coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_mr")

      ALLOCATE (coef%bkg_prfl_mr(coef%nlevels, coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of ref_prfl_mr")

      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

        DO i = 1, coef%nlevels
          READ (file_id,  * , iostat=err)pres, coef%ref_prfl_t(i, n), coef%ref_prfl_mr(i, n), coef%bkg_prfl_mr(i, n)
          THROWM(err.NE.0,"io status while reading section "//section)
          IF (coef%fmv_gas_id(n) == gas_id_mixed) coef%ref_prfl_p(i) = pres
        ENDDO
      ENDDO

!-------------------------------------------------------
    CASE ('PROFILE_ENVELOPE')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      ALLOCATE (coef%lim_prfl_p(coef%nlevels), STAT = err)
      THROWM(err.NE.0, "allocation of lim_prfl_p")

      ALLOCATE (coef%env_prfl_tmax(coef%nlevels), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_tmax")

      ALLOCATE (coef%env_prfl_tmin(coef%nlevels), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_tmin")

      ALLOCATE (coef%env_prfl_gmin(coef%nlevels, coef%fmv_gas), STAT = err)
      THROWM(err.NE.0, "allocation of env_prfl_gmin")

      ALLOCATE (coef%env_prfl_gmax(coef%nlevels, coef%fmv_gas), STAT = err)
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
! max volume mixing r (ppmv) (levels, gases)
! min volume mixing r (ppmv) (levels, gases)
          READ (file_id,  * , iostat=err)pres, coef%env_prfl_gmax(i, n), coef%env_prfl_gmin(i, n)
          THROWM(err.NE.0,"io status while reading section "//section)
        ENDDO
      ENDDO


!-------------------------------------------------------
    CASE ('FAST_COEFFICIENTS', 'COEF_SUB_FILES')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

! If section is COEF_SUB_FILES then coefficients for each gaz is stored
! in separate files.
! This possibility could be used to store very large coefficient files
! (large number of channels or gases)
! Section contains 1 line per gaz in the same order as the
! FAST_MODEL_VARIABLES section
! line indicates the filename of the coefficient for that gas
!
! No verification is done on the file.
! header lines starting with "!" sign are skipped

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

! loop on gases
      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0,"io status while reading section "//section)

! read dummy string of gas name or filename of the sub_coefficient file
        READ (file_id,  * , iostat=err)input_string
        THROWM(err.NE.0,"io status while reading section "//section)

! Case of Sub coefficient files
! Open the file and skip the header

        IF (TRIM(section) == 'COEF_SUB_FILES') THEN
          file_id_coef = 0
          for_output   = .FALSE.
          CALL rttov_opencoeff( &
                  err,          &
                  input_string, &
                  file_id_coef, &
                  for_output)
          THROWM(err.NE.0,"Error opening sub_coef file")

          CALL rttov_skipcommentline(file_id_coef, err)
          THROWM(err.NE.0,"io status while reading section "//section)
        ELSE
          file_id_coef = file_id
        ENDIF

        ! test existence of gas
        IF (coef%fmv_coe(n) > 0_jpim) THEN
          CALL read_ascii_fast_coef(err, coef%thermal, n, coef%fmv_gas_id(n), &
                                    coef%fmv_coe(n), coef%fmv_model_ver, lbl1)
          THROWM(err.NE.0,"error occurred in read_ascii_fast_coef while reading gas "//gas_name(n))

          IF (coef%fmv_ncorr(n) > 0_jpim) THEN
            CALL read_ascii_fast_coef(err, coef%thermal_corr, n, coef%fmv_gas_id(n), &
                                      coef%fmv_ncorr(n), coef%fmv_model_ver, lbl1)
            THROWM(err.NE.0,"error occurred in read_ascii_fast_coef while reading gas "//gas_name(n))
          ENDIF
        ENDIF

! For COEF_SUB_FILES close the intermediate coef file
        IF (TRIM(section) == 'COEF_SUB_FILES') THEN
          CLOSE (unit=file_id_coef)
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
          CALL read_ascii_fast_coef(err, coef%solar, n, coef%fmv_gas_id(n), &
                                    coef%fmv_coe(n), coef%fmv_model_ver, lbl1)
          THROWM(err.NE.0,"error occurred inside read_ascii_fast_coef while reading gas "//gas_name(n))

          IF (coef%fmv_ncorr(n) > 0_jpim) THEN
            CALL read_ascii_fast_coef(err, coef%solar_corr, n, coef%fmv_gas_id(n), &
                                      coef%fmv_ncorr(n), coef%fmv_model_ver, lbl1)
            THROWM(err.NE.0,"error occurred inside read_ascii_fast_coef while reading gas "//gas_name(n))
          ENDIF
        ENDIF
      ENDDO


    CASE ('NLTE_RADIANCE_COEFS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err .NE. 0, "io status while reading section "//section)

      coef%nltecoef = .TRUE. ! nlte coef present
      ALLOCATE (coef%nlte_coef, STAT = err)

      READ(file_id, *, iostat=err) coef%nlte_coef%ncoef, coef%nlte_coef%nsol, &
                                   coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
                                   coef%nlte_coef%start_chan
      THROWM(err.NE.0, "io status while reading NLTE dimensions")

      NULLIFY (coef%nlte_coef%coef)
      NULLIFY (coef%nlte_coef%sol_zen_angle, coef%nlte_coef%cos_sol)
      NULLIFY (coef%nlte_coef%sat_zen_angle, coef%nlte_coef%sec_sat)

      ! For any monotonic channel selection we must find those selected channels
      ! which lie within the range of NLTE channels in the coef file. This
      ! constitutes another contiguous block of channels in the coef structure.
      IF (.NOT. all_channels) THEN
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
          IF (channels(i) >= coef%nlte_coef%start_chan .AND. &
              channels(i) < coef%nlte_coef%start_chan + coef%nlte_coef%nchan) THEN
            nlte_count = nlte_count + 1
            nlte_chans(nlte_count) = channels(i) - coef%nlte_coef%start_chan + 1
            IF (nlte_count == 1) nlte_start = i
          ENDIF
        ENDDO

        IF (nlte_count > 0) THEN
          ! Reset NLTE channel variables according to input channel list
          nlte_file_nchan           = coef%nlte_coef%nchan
          coef%nlte_coef%start_chan = nlte_start
          coef%nlte_coef%nchan      = nlte_count
        ELSE
          ! No NLTE channels selected so no need for NLTE coef section
          coef%nlte_coef%start_chan = 0_jpim
          coef%nlte_coef%nchan = 0_jpim
          DEALLOCATE(coef%nlte_coef, nlte_chans)
          coef%nltecoef = .FALSE.
        ENDIF
      ENDIF

      IF (coef%nltecoef) THEN

        ALLOCATE(coef%nlte_coef%sol_zen_angle(coef%nlte_coef%nsol), STAT = err)
        THROWM(err.NE.0, "allocation of NLTE solar zenith angle array")
        ALLOCATE(coef%nlte_coef%sec_sat(coef%nlte_coef%nsat), STAT = err)
        THROWM(err.NE.0, "allocation of NLTE satellite zenith angle array")
        ALLOCATE(coef%nlte_coef%coef(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
                 coef%nlte_coef%nsol, coef%nlte_coef%nchan), STAT = err)
        THROWM(err.NE.0, "allocation of NLTE coef array")

        READ(file_id,*,iostat=err) (coef%nlte_coef%sol_zen_angle(i), &
          i = 1, coef%nlte_coef%nsol)
        THROWM(err.NE.0, "io status while reading NLTE solar zenith angles")

        READ(file_id,*,iostat=err)(coef%nlte_coef%sec_sat(i), &
          i = 1, coef%nlte_coef%nsat)
        THROWM(err.NE.0, "io status while reading NLTE satellite zenith angles")

        IF (all_channels) THEN
          DO ichan = 1, coef%nlte_coef%nchan
            DO isol = 1, coef%nlte_coef%nsol
              DO isat = 1, coef%nlte_coef%nsat
                READ(file_id,*,iostat=err) &
                  (coef%nlte_coef%coef(i, isat, isol, ichan), &
                  i = 1, coef%nlte_coef%ncoef)
                THROWM(err.NE.0, "io status while reading NLTE coefs")
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ALLOCATE(nlte_values(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
                   coef%nlte_coef%nsol, nlte_file_nchan), STAT = err)
          THROWM(err.NE.0, "allocation of nlte_values array")
          DO ichan = 1, nlte_file_nchan
            DO isol = 1, coef%nlte_coef%nsol
              DO isat = 1, coef%nlte_coef%nsat
                READ(file_id,*,iostat=err) &
                  (nlte_values(i, isat, isol, ichan), &
                  i = 1, coef%nlte_coef%ncoef)
                THROWM(err.NE.0, "io status while reading NLTE coefs")
              ENDDO
            ENDDO
          ENDDO
          coef%nlte_coef%coef(:,:,:,:) = nlte_values(:,:,:,nlte_chans(1:nlte_count))
          DEALLOCATE(nlte_chans, nlte_values)
        ENDIF

      ENDIF


    CASE ('PRESSURE_MODULATED_CELL')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0,"io status while reading section "//section)

      coef%pmc_shift = .TRUE. ! pmc coef present

! cell length (cm)
      READ (file_id,  * , iostat=err)coef%pmc_lengthcell
      THROWM(err.NE.0,"io status while reading section "//section)

! nominal channel cell pressures (hPa)
      ALLOCATE (coef%pmc_pnominal(coef%fmv_chn), STAT = err)
      THROWM(err.NE.0, "allocation of coef%pmc_pnominal array")
      IF (all_channels) THEN
        READ (file_id,  * , iostat=err) (coef%pmc_pnominal(i), i = 1, coef%fmv_chn)
        THROWM(err.NE.0,"io status while reading section "//section)
      ELSE
        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = err)
        THROWM(err.NE.0, "allocation of values1")
        READ (file_id,  * , iostat=err) (values1(i), i = 1, coef%fmv_ori_nchn)
        coef%pmc_pnominal(:) = values1(channels(:))
        DEALLOCATE (values1, STAT = err)
        THROWM(err.NE.0, "deallocation of values1")
      ENDIF

! cell temperature (K)
      READ (file_id,  * , iostat=err)coef%pmc_tempcell
      THROWM(err.NE.0,"io status while reading section "//section)

! gamma(co2)/gamma(air)
! ratio of band-averaged self- and air-broadened halfwidths of CO2 lines
      READ (file_id,  * , iostat=err)coef%pmc_betaplus1
      THROWM(err.NE.0,"io status while reading section "//section)

! Number of layers used
      READ (file_id,  * , iostat=err)coef%pmc_nlay
      THROWM(err.NE.0,"io status while reading section "//section)

! Number of variables used
      READ (file_id,  * , iostat=err)coef%pmc_nvar
      THROWM(err.NE.0,"io status while reading section "//section)

! Cell pressure coefficients
      ALLOCATE (coef%pmc_coef(coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar), STAT = err)
      THROWM(err.NE.0, "allocation of coef%pmc_coef array")
      IF (all_channels) THEN
        coeffsarray => coef%pmc_coef
      ELSE
        ALLOCATE (coeffsarray(coef%pmc_nlay, coef%fmv_ori_nchn, coef%pmc_nvar), STAT = err)
        THROWM(err.NE.0, "allocation of coeffsarray")
      ENDIF

      READ (file_id_coef,  * , iostat=err)     &
            (((coeffsarray(i, j, k), k = 1, coef%pmc_nvar), j = 1, coef%fmv_ori_nchn), i = 1, coef%pmc_nlay)
      THROWM(err.NE.0,"io status while reading section "//section)

      IF (.NOT. all_channels) THEN
        coef%pmc_coef(:,:,:) = coeffsarray(:, channels(:), :)
        DEALLOCATE (coeffsarray, STAT = err)
        THROWM(err.NE.0, "deallocation of coeffsarray")
      ENDIF


    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH

CONTAINS
  SUBROUTINE read_ascii_fast_coef(err, fast_coef, gas_pos, gas_id, ncoef, version, lbl_mode)

    USE rttov_types, ONLY : rttov_fast_coef
    USE rttov_fast_coef_utils_mod, ONLY : set_pointers, nullify_pointers

    INTEGER(jpim),         INTENT(OUT)          :: err
    TYPE(rttov_fast_coef), INTENT(INOUT)        :: fast_coef(:)
    INTEGER(jpim),         INTENT(IN)           :: gas_pos
    INTEGER(jpim),         INTENT(IN)           :: gas_id
    INTEGER(jpim),         INTENT(IN)           :: ncoef
    INTEGER(jpim),         INTENT(IN)           :: version
    LOGICAL(jplm),         INTENT(IN)           :: lbl_mode

    REAL(jprb), ALLOCATABLE :: coef_temp(:,:,:)
    INTEGER(jpim) :: i, j, k, chn

    TRY
    ALLOCATE(coef_temp(ncoef,coef%nlayers,coef%fmv_ori_nchn), stat = err)
    THROW(err.NE.0)

    IF (version > 9) THEN
      ! New coefs are stored in the order we want them in
      READ (file_id_coef,  * , iostat=err) coef_temp
    ELSE
      ! This is transposed from pre-RTTOV 11.3 order
      ! Coefs will be transposed for more efficent access in memory
      ! but remain same on disk (for pre-RTTOV 11.3 coefs)
      READ (file_id_coef, *, iostat=err) &
        (((coef_temp(k,i,j), i = 1, coef%nlayers), j = 1, coef%fmv_ori_nchn), k = 1, ncoef)
    ENDIF
    THROW(err.NE.0)

    DO i = 1, coef%fmv_chn
      chn = i
      IF (.NOT. all_channels) chn = channels(i)

      ! Allocate space only for non-zero coefs
      IF (lbl_mode .OR. ANY(coef_temp(:,:,chn) /= 0._jprb)) THEN
        ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef,coef%nlayers), STAT = err)
        THROW(err.NE.0)

        fast_coef(i)%gasarray(gas_pos)%coef = coef_temp(:,:,chn)
        CALL set_pointers(fast_coef(i), gas_pos, gas_id)
      ELSE
        NULLIFY (fast_coef(i)%gasarray(gas_pos)%coef)
        CALL nullify_pointers(fast_coef(i), gas_id)
      ENDIF
    ENDDO

    DEALLOCATE(coef_temp, stat = err)
    CATCH
  END SUBROUTINE read_ascii_fast_coef

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

END SUBROUTINE rttov_read_ascii_coef
