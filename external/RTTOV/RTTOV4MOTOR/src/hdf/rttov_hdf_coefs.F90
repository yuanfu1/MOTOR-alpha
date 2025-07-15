! Description:
!> @file
!!   Subroutines for HDF5 I/O of coefs structure
!
!> @brief
!!   Subroutines for HDF5 I/O of coefs structure
!!
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
     MODULE RTTOV_HDF_COEFS
#include "throw.h"
      USE RTTOV_TYPES
      USE RTTOV_HDF_MOD
      USE HDF5, only : HID_T
      USE H5LT

      USE RTTOV_CONST, ONLY :  &
        & gas_id_mixed,           &
        & gas_id_watervapour,     &
        & gas_id_ozone,           &
        & gas_id_wvcont,          &
        & gas_id_co2,             &
        & gas_id_n2o,             &
        & gas_id_co,              &
        & gas_id_ch4,             &
        & gas_id_so2,             &
        & mfasis_maxzenangle

      USE RTTOV_FAST_COEF_UTILS_MOD, ONLY : set_pointers
!
      IMPLICIT NONE
#include "rttov_errorreport.interface"
#include "rttov_alloc_phfn_int.interface"

      PRIVATE
      PUBLIC :: RTTOV_HDF_COEF_WH,       RTTOV_HDF_COEF_RH
      PUBLIC :: RTTOV_HDF_SCCLDCOEF_WH,  RTTOV_HDF_SCCLDCOEF_RH
      PUBLIC :: RTTOV_HDF_SCAERCOEF_WH,  RTTOV_HDF_SCAERCOEF_RH
      PUBLIC :: RTTOV_HDF_PCCOEF_WH,     RTTOV_HDF_PCCOEF_RH
      PUBLIC :: RTTOV_HDF_MFASISCOEF_WH, RTTOV_HDF_MFASISCOEF_RH

      CONTAINS

!> Write an optical depth coefficient structure to HDF5 file
!! param[in]  x             optical depth coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
      SUBROUTINE RTTOV_HDF_COEF_WH(X,LUN,ERR,COMPRESS,FORCE_DOUBLE)
USE RTTOV_HDF_RTTOV_COEF_IO
USE RTTOV_HDF_RTTOV_FAST_COEF_IO
USE RTTOV_HDF_RTTOV_NLTE_COEF_IO

      TYPE(RTTOV_COEF),INTENT(IN)    ::X
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS
      LOGICAL,INTENT(IN),OPTIONAL    ::FORCE_DOUBLE

!
      INTEGER(HID_T) :: G_ID_SUB

      TYPE(rttov_fast_coef_hdf_io) :: FAST_COEF_temp
!
TRY

        CALL RTTOV_HDF_RTTOV_COEF_WH( x, LUN, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE COEF")

        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV coefficient structure" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

        CALL MKPAR( LUN, "THERMAL", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP THERMAL")

! Move fast_coef structure in to temp for writing

        IF (x%ncmixed > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%mixedgas, gas_id_mixed)
        IF (x%ncwater > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%watervapour, gas_id_watervapour)
        IF (x%ncozone > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%ozone, gas_id_ozone)
        IF (x%ncwvcont > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%wvcont, gas_id_wvcont)
        IF (x%ncco2 > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%co2, gas_id_co2)
        IF (x%ncn2o > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%n2o, gas_id_n2o)
        IF (x%ncco > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%co, gas_id_co)
        IF (x%ncch4 > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%ch4, gas_id_ch4)
        IF (x%ncso2 > 0_jpim ) &
          CALL move_fast_coef_hdf(x, x%thermal(:), FAST_COEF_temp%so2, gas_id_so2)

        CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
        THROWM(ERR.NE.0,"CANNOT WRITE THERMAL COEFFICIENTS")

        IF (x%ncmixed > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
        IF (x%ncwater > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
        IF (x%ncozone > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
        IF (x%ncwvcont > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
        IF (x%ncco2 > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
        IF (x%ncn2o > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
        IF (x%ncco > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
        IF (x%ncch4 > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
        IF (x%ncso2 > 0_jpim ) &
          DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")

        IF (x%fmv_model_ver > 9) THEN
          CALL MKPAR( LUN, "THERMAL_CORR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP THERMAL_CORR")

          IF (x%nccmixed > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%mixedgas, gas_id_mixed, corr = .TRUE._jplm)
          IF (x%nccwater > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%watervapour, gas_id_watervapour, corr = .TRUE._jplm)
          IF (x%nccozone > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%ozone, gas_id_ozone, corr = .TRUE._jplm)
          IF (x%nccwvcont > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%wvcont, gas_id_wvcont, corr = .TRUE._jplm)
          IF (x%nccco2 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%co2, gas_id_co2, corr = .TRUE._jplm)
          IF (x%nccn2o > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%n2o, gas_id_n2o, corr = .TRUE._jplm)
          IF (x%nccco > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%co, gas_id_co, corr = .TRUE._jplm)
          IF (x%nccch4 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%ch4, gas_id_ch4, corr = .TRUE._jplm)
          IF (x%nccso2 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%thermal_corr(:), FAST_COEF_temp%so2, gas_id_so2, corr = .TRUE._jplm)

          CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
          THROWM(ERR.NE.0,"CANNOT WRITE THERMAL_CORR COEFFICIENTS")

          IF (x%nccmixed > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
          IF (x%nccwater > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
          IF (x%nccozone > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
          IF (x%nccwvcont > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
          IF (x%nccco2 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
          IF (x%nccn2o > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
          IF (x%nccco > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
          IF (x%nccch4 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
          IF (x%nccso2 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL_CORR")
        ENDIF

        IF( ASSOCIATED( x%SOLAR ) .AND. x%SOLARCOEF ) THEN
          CALL MKPAR( LUN, "SOLAR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP SOLAR")

          IF (x%ncmixed > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%mixedgas, gas_id_mixed)
          IF (x%ncwater > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%watervapour, gas_id_watervapour)
          IF (x%ncozone > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%ozone, gas_id_ozone)
          IF (x%ncwvcont > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%wvcont, gas_id_wvcont)
          IF (x%ncco2 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%co2, gas_id_co2)
          IF (x%ncn2o > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%n2o, gas_id_n2o)
          IF (x%ncco > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%co, gas_id_co)
          IF (x%ncch4 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%ch4, gas_id_ch4)
          IF (x%ncso2 > 0_jpim ) &
            CALL move_fast_coef_hdf(x, x%solar(:), FAST_COEF_temp%so2, gas_id_so2)

          CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
          THROWM(ERR.NE.0,"CANNOT WRITE SOLAR COEFFICIENTS")

          IF (x%ncmixed > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
          IF (x%ncwater > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
          IF (x%ncozone > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
          IF (x%ncwvcont > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
          IF (x%ncco2 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
          IF (x%ncn2o > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
          IF (x%ncco > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
          IF (x%ncch4 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
          IF (x%ncso2 > 0_jpim ) &
            DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")

          IF (x%fmv_model_ver > 9) THEN
            CALL MKPAR( LUN, "SOLAR_CORR", G_ID_SUB, ERR )
            THROWM(ERR.NE.0,"CANNOT CREATE GROUP SOLAR_CORR")

            IF (x%nccmixed > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%mixedgas, gas_id_mixed, corr = .TRUE._jplm)
            IF (x%nccwater > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%watervapour, gas_id_watervapour, corr = .TRUE._jplm)
            IF (x%nccozone > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%ozone, gas_id_ozone, corr = .TRUE._jplm)
            IF (x%nccwvcont > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%wvcont, gas_id_wvcont, corr = .TRUE._jplm)
            IF (x%nccco2 > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%co2, gas_id_co2, corr = .TRUE._jplm)
            IF (x%nccn2o > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%n2o, gas_id_n2o, corr = .TRUE._jplm)
            IF (x%nccco > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%co, gas_id_co, corr = .TRUE._jplm)
            IF (x%nccch4 > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%ch4, gas_id_ch4, corr = .TRUE._jplm)
            IF (x%nccso2 > 0_jpim ) &
              CALL move_fast_coef_hdf(x, x%solar_corr(:), FAST_COEF_temp%so2, gas_id_so2, corr = .TRUE._jplm)

            CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
            THROWM(ERR.NE.0,"CANNOT WRITE SOLAR_CORR COEFFICIENTS")

            IF (x%nccmixed > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
            IF (x%nccwater > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
            IF (x%nccozone > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
            IF (x%nccwvcont > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
            IF (x%nccco2 > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
            IF (x%nccn2o > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
            IF (x%nccco > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
            IF (x%nccch4 > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
            IF (x%nccso2 > 0_jpim ) &
              DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)

            CALL H5GCLOSE_F( G_ID_SUB, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR_CORR")
          ENDIF
        ENDIF

        IF( x%NLTECOEF ) THEN
          CALL MKPAR( LUN, "NLTE_COEF", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP NLTE_COEF")

          CALL RTTOV_HDF_RTTOV_NLTE_COEF_WH( x%NLTE_COEF, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
          THROWM(ERR.NE.0,"CANNOT WRITE COEF%NLTE_COEF")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")
        ENDIF

CATCH
      END SUBROUTINE

!> Write a cloud coefficient structure to HDF5 file
!! param[in]  x             RTTOV optical property coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
    SUBROUTINE RTTOV_HDF_SCCLDCOEF_WH(X,LUN,ERR,COMPRESS,FORCE_DOUBLE)
      TYPE(RTTOV_COEF_SCATT), INTENT(IN)           :: X
      INTEGER(HID_T),         INTENT(IN)           :: LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)          :: ERR
      LOGICAL,                INTENT(IN), OPTIONAL :: COMPRESS
      LOGICAL,                INTENT(IN), OPTIONAL :: FORCE_DOUBLE

      INTEGER(HID_T) :: G_ID_SUB
      TRY

      CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
      "This is a RTTOV cloud coefficient structure" // &
      CHAR(0), ERR )
      THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

      IF ( X%OPTP_WCL_OPAC%NCHAN > 0 ) THEN
        CALL MKPAR( LUN, "WATER_CLOUD_OPAC", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP WATER_CLOUD_OPAC")

        CALL RTTOV_HDF_OPTP_WH( X%OPTP_WCL_OPAC, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
        THROWM(ERR.NE.0,"CANNOT WRITE WATER_CLOUD_OPAC COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATER_CLOUD_OPAC")
      ENDIF

      IF ( X%OPTP_WCL_DEFF%NCHAN > 0 ) THEN
        CALL MKPAR( LUN, "WATER_CLOUD_DEFF", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP WATER_CLOUD_DEFF")

        CALL RTTOV_HDF_OPTP_WH( X%OPTP_WCL_DEFF, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
        THROWM(ERR.NE.0,"CANNOT WRITE WATER_CLOUD_DEFF COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATER_CLOUD_DEFF")
      ENDIF

      IF ( X%OPTP_ICL_BAUM%NCHAN > 0 ) THEN
        CALL MKPAR( LUN, "ICE_CLOUD_BAUM", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP ICE_CLOUD_BAUM")

        CALL RTTOV_HDF_OPTP_WH( X%OPTP_ICL_BAUM, G_ID_SUB, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE )
        THROWM(ERR.NE.0,"CANNOT WRITE ICE_CLOUD_BAUM COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP ICE_CLOUD_BAUM")
      ENDIF

      CATCH
    END SUBROUTINE

!> Write an aerosol coefficient structure to HDF5 file
!! param[in]  x             RTTOV optical property coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
    SUBROUTINE RTTOV_HDF_SCAERCOEF_WH(X,LUN,ERR,COMPRESS,FORCE_DOUBLE)
      TYPE(RTTOV_COEF_SCATT), INTENT(IN)           :: X
      INTEGER(HID_T),         INTENT(IN)           :: LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)          :: ERR
      LOGICAL,                INTENT(IN), OPTIONAL :: COMPRESS
      LOGICAL,                INTENT(IN), OPTIONAL :: FORCE_DOUBLE
      TRY

      CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        & "This is a RTTOV aerosol coefficient structure" // &
        & CHAR(0), ERR )
      THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

      CALL RTTOV_HDF_OPTP_WH(X%OPTP_AER,LUN,ERR,COMPRESS,FORCE_DOUBLE)
      THROWM(ERR.NE.0,"CANNOT WRITE AEROSOL OPTICAL PROPERTIES")

      CATCH
    END SUBROUTINE



!> Write a cloud/aerosol optical property structure to HDF5 file
!! param[in]  x             RTTOV optical property coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
    SUBROUTINE RTTOV_HDF_OPTP_WH(X,LUN,ERR,COMPRESS,FORCE_DOUBLE)
      USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE

      TYPE(RTTOV_OPTP),       INTENT(IN)           :: X
      INTEGER(HID_T),         INTENT(IN)           :: LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)          :: ERR
      LOGICAL,                INTENT(IN), OPTIONAL :: COMPRESS
      LOGICAL,                INTENT(IN), OPTIONAL :: FORCE_DOUBLE

      CHARACTER(LEN=LENSH)  :: SNAME
      CHARACTER(LEN=LENSH)  :: GNAME
      INTEGER(KIND=JPIM)    :: I
      CHARACTER(LEN=4), POINTER :: C1(:)
      INTEGER(HID_T) :: G_ID_SUB
      TRY

      sname='VERSION'
      call write_array_hdf(lun,sname,&
        & 'File format version',&
        & err,i0=x%version )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      sname='ID'
      call write_array_hdf(lun,sname,&
        & 'Contents ID',&
        & err,i0=x%id )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      sname='NCHAN'
      call write_array_hdf(lun,sname,&
        & 'Number of channels for which optical parameters are stored',&
        & err,i0=x%nchan )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      sname='NCHAN_PHA'
      call write_array_hdf(lun,sname,&
        & 'Number of channels for which phase function values are stored',&
        & err,i0=x%nchan_pha)
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      sname='NTYPES'
      call write_array_hdf(lun,sname,&
        & 'Number of aerosols components',&
        & err,i0=x%ntypes )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      sname='MAXNMOM'
      call write_array_hdf(lun,sname,&
        & 'Maximum number of Legendre coefficients for phase functions',&
        & err,i0=x%maxnmom )
      THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

      if(x%nchan_pha > 0)then
        sname='NPHANGLE'
        call write_array_hdf(lun,sname,&
          & 'Number of phase angles for phase functions',&
          & err,i0=x%nphangle )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif

      if(associated(x%chan_pha))then
        sname='CHAN_PHA'
        call write_array_hdf(lun,sname,&
          & 'The solar channel indexes',&
          & err,i1=x%chan_pha )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif

      if(associated(x%phangle))then
        sname='PHANGLE'
        call write_array_hdf(lun,sname,&
          & 'The phase function angle grid',&
          & err,r1=x%phangle )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
      endif

      IF (x%ntypes > 1) THEN
        sname='NAMES'
        allocate(c1(x%ntypes))
        c1 = x%data(:)%name
        call write_array_hdf(lun,sname,&
          & 'Particle group names',&
          & err,c1=c1)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        deallocate(c1)
      ENDIF

      DO I = 1, x%ntypes

        IF (x%ntypes > 1) THEN
          CALL RTTOV_UPPER_CASE(GNAME, x%data(i)%name)
        ELSE
          GNAME = 'DATA'
        ENDIF

        CALL MKPAR( LUN, TRIM(GNAME), G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//TRIM(GNAME))

        sname='NRELHUM'
        call write_array_hdf(g_id_sub,sname,&
          & 'Number of relative humidity values',&
          & err,i0=x%data(i)%nrelhum )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='RELHUM'
        call write_array_hdf(g_id_sub,sname,&
          & 'Relative humidity values',&
          & err,r1=x%data(i)%relhum , units = 'percent')
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='CONFAC'
        call write_array_hdf(g_id_sub,sname,&
          & 'Conversion factor to particle density',&
          & err,r0=x%data(i)%confac , units = '[particle.cm^-3]/[g.m^-3]')
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='NDEFF'
        call write_array_hdf(g_id_sub,sname,&
          & 'Number of effective diameter values',&
          & err,i0=x%data(i)%ndeff )
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='DEFF'
        call write_array_hdf(g_id_sub,sname,&
          & 'Effective diameter values',&
          & err,r1=x%data(i)%deff , units = 'microns')
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='ABS'
        call write_array_hdf(g_id_sub,sname,&
          & 'Absorption (humidity, deff, channels)',&
          & err,r3=x%data(i)%abs(:,:,:) , units = 'm-1',compress=compress,force_double=force_double)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='SCA'
        call write_array_hdf(g_id_sub,sname,&
          & 'Scattering (humidity, deff, channels)',&
          & err,r3=x%data(i)%sca(:,:,:) , units = 'm-1',compress=compress,force_double=force_double)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        sname='BPR'
        if(associated(x%data(i)%bpr))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Back scattering factor (humidity, deff, channels)',&
            & err,r3=x%data(i)%bpr(:,:,:) ,compress=compress,force_double=force_double)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif

        sname='NMOM'
        if(associated(x%data(i)%nmom))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Number of Legendre coefficients (humidity, channels)',&
            & err,i2=x%data(i)%nmom ,compress=compress)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif

        sname='LEGCOEF'
        if(associated(x%data(i)%legcoef))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Phase function Legendre coefficients (1:maxnmom+1, humidity, deff, channels)',&
            & err,r4=x%data(i)%legcoef ,compress=compress,force_double=force_double)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif

        sname='PHA'
        if(associated(x%data(i)%pha))then
          call write_array_hdf(g_id_sub,sname,&
            & 'Phase functions for solar channels (nphangle, humidity, deff, channels)',&
            & err,r4=x%data(i)%pha ,compress=compress,force_double=force_double)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        endif

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))

      ENDDO

      CATCH
    END SUBROUTINE


!> Write PC-RTTOV coefficient structure to HDF5 file
!! param[in]  x             PC coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
      SUBROUTINE RTTOV_HDF_PCCOEF_WH(X,LUN,ERR,COMPRESS,FORCE_DOUBLE)
USE RTTOV_HDF_RTTOV_COEF_PCC_IO
USE RTTOV_HDF_RTTOV_COEF_PCC1_IO
USE RTTOV_HDF_RTTOV_COEF_PCC2_IO

      TYPE(RTTOV_COEF_PCCOMP),INTENT(IN)    ::X
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS
      LOGICAL,INTENT(IN),OPTIONAL    ::FORCE_DOUBLE

      CHARACTER(LEN=LENSH)  :: GNAME, GNAME2
      INTEGER(KIND=JPIM)    :: I, J
!
      INTEGER(HID_T) :: G_ID_SUB, G_ID_SUB2, G_ID_SUB3
!
TRY

        CALL RTTOV_HDF_RTTOV_COEF_PCC_WH(X,LUN,ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
        THROWM(ERR.NE.0,"CANNOT WRITE PC COEF")

        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV coefficient structure PCCOMP" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

        ! PCREG structure
        CALL MKPAR( LUN, "PCREG", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG")

        DO J = 1, X%FMV_PC_BANDS
          WRITE(GNAME,'(I2.2)') J
          CALL MKPAR( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG/"//TRIM(GNAME))

          DO I = 1, X%FMV_PC_SETS(J)

            WRITE(GNAME2,'(I2.2)') I
            CALL MKPAR( G_ID_SUB2, TRIM(GNAME2), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CREATE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_WH( X%PCREG(j,i), G_ID_SUB3, ERR , COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
            THROWM(ERR.NE.0,"CANNOT WRITE PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL H5GCLOSE_F( G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

          END DO

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG")

        ! EIGEN structure
        CALL MKPAR( LUN, "EIGEN", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP EIGEN")

        DO I = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') I
          CALL MKPAR( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP EIGEN/"//TRIM(GNAME))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_WH( X%EIGEN(i), G_ID_SUB2, ERR, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
          THROWM(ERR.NE.0,"CANNOT WRITE EIGEN/"//TRIM(GNAME))

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN")


CATCH
      END SUBROUTINE

!> Write MFASIS LUT structure to HDF5 file
!! param[in]  x             MFASIS LUT structure
!! param[in]  file_type     CLD or AER for cloud or aerosol structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  compress      if true will apply internal HDF5 compression, optional
!! param[in]  force_double  if true all real values are stored as H5T_NATIVE_DOUBLE, optional
      SUBROUTINE RTTOV_HDF_MFASISCOEF_WH(X, FILE_TYPE, LUN,ERR,COMPRESS,FORCE_DOUBLE)
        USE RTTOV_CONST, ONLY: ERRORSTATUS_FATAL

        TYPE(RTTOV_COEF_MFASIS), INTENT(IN)    :: X
        CHARACTER(LEN=3),        INTENT(IN)    :: FILE_TYPE
        INTEGER(HID_T),          INTENT(IN)    :: LUN
        INTEGER(KIND=JPIM),      INTENT(OUT)   :: ERR
        LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS
        LOGICAL,INTENT(IN),OPTIONAL    ::FORCE_DOUBLE

        CHARACTER(LEN=LENSH) :: GNAME, SNAME
        INTEGER(KIND=JPIM)   :: I
        INTEGER(HID_T)       :: G_ID_SUB

        TRY

        SNAME = 'NCHANNELS_COEF'
        CALL write_array_hdf(LUN,SNAME,&
        'Number of channels in corresponding rtcoef file', &
        ERR,I0=X%NCHANNELS_COEF)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'VERSION'
        CALL write_array_hdf(LUN,SNAME,&
        'Version number',&
        ERR,I0=X%VERSION)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'FILE_TYPE'
        CALL write_array_hdf(LUN,SNAME,&
        'Type of MFASIS file: 1 => cloud, 2 => aerosol',&
        ERR,I0=X%FILE_TYPE)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'README_LUT'
        CALL write_array_hdf(LUN,SNAME,&
        'Readme for LUTs',&
        ERR,C1=X%README_LUT)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'NCHANNELS'
        CALL write_array_hdf(LUN,SNAME,&
        'Number of channels supported by MFASIS',&
        ERR,I0=X%NCHANNELS)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'CHANNEL_LIST'
        CALL write_array_hdf(LUN,SNAME,&
        'List of channels for which LUTs are stored',&
        ERR,I1=X%CHANNEL_LIST)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'CHANNEL_LUT_INDEX'
        CALL write_array_hdf(LUN,SNAME,&
        'Index into channel_list for each channel in rtcoef file',&
        ERR,I1=X%CHANNEL_LUT_INDEX)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'NDIMS'
        CALL write_array_hdf(LUN,SNAME,&
        'Number of dimensions in LUTs',&
        ERR,I0=X%NDIMS)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        SNAME = 'NPARTICLES'
        CALL write_array_hdf(LUN,SNAME,&
        'Number of particle types included in LUTs',&
        ERR,I0=X%NPARTICLES)
        THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

        IF (FILE_TYPE == 'AER') THEN
          SNAME = 'AER_TYPES'
          CALL write_array_hdf(LUN,SNAME,&
          'Aerosol types included in LUTs',&
          ERR,I1=X%AER_TYPES)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        ELSE
          SNAME = 'CLW_SCHEME'
          CALL write_array_hdf(LUN,SNAME,&
          'CLW scheme used for training LUT',&
          ERR,I0=X%CLW_SCHEME)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = 'ICE_SCHEME'
          CALL write_array_hdf(LUN,SNAME,&
          'ICE scheme used for training LUT',&
          ERR,I0=X%ICE_SCHEME)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        ENDIF

        IF (X%VERSION > 0) THEN
          SNAME = 'MAXZENANGLE'
          CALL write_array_hdf(LUN,SNAME,&
          'Mazimum zenith angle used in LUT training',&
          ERR,R0=X%MAXZENANGLE)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
        ENDIF

        DO I = 1, X%NDIMS
          WRITE(GNAME, '(A,I3.3)') 'DIM', I

          CALL MKPAR( LUN, GNAME, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//GNAME)

          SNAME = TRIM(GNAME)//'/NAME'
          CALL write_array_hdf(LUN,SNAME,&
          'Dimension name',&
          ERR,C0=X%LUT_AXES(I)%NAME)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = TRIM(GNAME)//'/DIM_TYPE'
          CALL write_array_hdf(LUN,SNAME,&
          'Type of dimension',&
          ERR,I0=X%LUT_AXES(I)%DIM_TYPE)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = TRIM(GNAME)//'/NVALUES'
          CALL write_array_hdf(LUN,SNAME,&
          'Size of dimension',&
          ERR,I0=X%LUT_AXES(I)%NVALUES)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = TRIM(GNAME)//'/VALUES'
          CALL write_array_hdf(LUN,SNAME,&
          'Dimension values',&
          ERR,R1=X%LUT_AXES(I)%VALUES)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))
        ENDDO

        DO I = 1, X%NCHANNELS
          WRITE(GNAME, '(A,I6.6)') 'CH', I

          CALL MKPAR( LUN, GNAME, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CREATE GROUP "//GNAME)

          SNAME = TRIM(GNAME)//'/NLUTS'
          CALL write_array_hdf(LUN,SNAME,&
          'Number of LUTs for channel '//TRIM(GNAME),&
          ERR,I0=X%LUT(I)%NLUTS)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = TRIM(GNAME)//'/QINT'
          CALL write_array_hdf(LUN,SNAME,&
          'Water vapour values for each LUT',&
          ERR,R2=X%LUT(I)%QINT)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          SNAME = TRIM(GNAME)//'/DATA'
          CALL write_array_hdf(LUN,SNAME,&
          'MFASIS LUT(s) for channel '//TRIM(GNAME),&
          ERR,R2=X%LUT(I)%DATA, COMPRESS=COMPRESS, FORCE_DOUBLE=FORCE_DOUBLE)
          THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))
        ENDDO

        CATCH
      END SUBROUTINE RTTOV_HDF_MFASISCOEF_WH


!> Read an optical depth coefficient structure from HDF5 file
!! param[out] x             optical depth coefficient structure
!! param[in]  lun           file ID of HDF5 file
!! param[out] err           return status
!! param[in]  lbl           set this to true only if reading optical depth coefs
!!                            from LBL code (default false), optional
      SUBROUTINE RTTOV_HDF_COEF_RH(X,LUN,ERR,LBL)

USE RTTOV_HDF_RTTOV_COEF_IO
USE RTTOV_HDF_RTTOV_FAST_COEF_IO
USE RTTOV_HDF_RTTOV_NLTE_COEF_IO
USE RTTOV_CONST, ONLY : &
  version_compatible_min, &
  version_compatible_max

      TYPE(RTTOV_COEF),INTENT(OUT)    ::X
      INTEGER(HID_T),INTENT(IN)       ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT)  ::ERR
      LOGICAL(KIND=JPLM),INTENT(IN)   ::LBL

      LOGICAL :: LEXT
      INTEGER(HID_T) :: G_ID_SUB
      INTEGER(KIND=JPIM) :: n, i

      TYPE(rttov_fast_coef_hdf_io) :: FAST_COEF_temp
TRY

        ! First check the file format version
        CALL read_array_hdf(LUN,'ID_COMP_LVL',ERR,i0=x%ID_COMP_LVL)
        THROWM(ERR.NE.0,"CANNOT READ "//trim('ID_COMP_LVL'))
        IF (x%id_comp_lvl < version_compatible_min .OR. &
            x%id_comp_lvl > version_compatible_max) THEN
          err = errorstatus_fatal
          THROWM(ERR.NE.0, "Version of coefficient file is incompatible with RTTOV library")
        ENDIF

        CALL RTTOV_HDF_RTTOV_COEF_RH( x, LUN, ERR )
        THROWM(ERR.NE.0,"CANNOT READ COEF")

        ! No channel selection with HDF5 so number of channels
        ! in file is same as number of channels extracted
        x%fmv_ori_nchn = x%fmv_chn

        IF (x%fmv_model_ver <= 9) THEN
          IF (.NOT. ASSOCIATED(x%fmv_ncorr)) ALLOCATE(x%fmv_ncorr(x%fmv_gas))
          x%fmv_ncorr = 0
          x%nlevels = x%fmv_lvl(1)
        ELSE
          IF (.NOT. ASSOCIATED(x%fmv_lvl)) ALLOCATE(x%fmv_lvl(x%fmv_gas))
          x%fmv_lvl = x%nlevels
        ENDIF

        DO n = 1, x%fmv_gas
          SELECT CASE (x%fmv_gas_id(n))
          CASE (gas_id_mixed)
            x%nmixed  = x%fmv_var(n)
            x%ncmixed = x%fmv_coe(n)
            x%nccmixed = x%fmv_ncorr(n)
          CASE (gas_id_watervapour)
            x%nwater  = x%fmv_var(n)
            x%ncwater = x%fmv_coe(n)
            x%nccwater = x%fmv_ncorr(n)
          CASE (gas_id_ozone)
            x%nozone  = x%fmv_var(n)
            x%ncozone = x%fmv_coe(n)
            x%nccozone = x%fmv_ncorr(n)
          CASE (gas_id_wvcont)
            x%nwvcont  = x%fmv_var(n)
            x%ncwvcont = x%fmv_coe(n)
            x%nccwvcont = x%fmv_ncorr(n)
          CASE (gas_id_co2)
            x%nco2  = x%fmv_var(n)
            x%ncco2 = x%fmv_coe(n)
            x%nccco2 = x%fmv_ncorr(n)
          CASE (gas_id_n2o)
            x%nn2o  = x%fmv_var(n)
            x%ncn2o = x%fmv_coe(n)
            x%nccn2o = x%fmv_ncorr(n)
          CASE (gas_id_co)
            x%nco  = x%fmv_var(n)
            x%ncco = x%fmv_coe(n)
            x%nccco = x%fmv_ncorr(n)
          CASE (gas_id_ch4)
            x%nch4  = x%fmv_var(n)
            x%ncch4 = x%fmv_coe(n)
            x%nccch4 = x%fmv_ncorr(n)
          CASE (gas_id_so2)
            x%nso2  = x%fmv_var(n)
            x%ncso2 = x%fmv_coe(n)
            x%nccso2 = x%fmv_ncorr(n)
          END SELECT
        END DO
        x%NLAYERS = x%NLEVELS - 1

        NULLIFY(FAST_COEF_temp%mixedgas)
        NULLIFY(FAST_COEF_temp%watervapour)
        NULLIFY(FAST_COEF_temp%ozone)
        NULLIFY(FAST_COEF_temp%wvcont)
        NULLIFY(FAST_COEF_temp%co2)
        NULLIFY(FAST_COEF_temp%n2o)
        NULLIFY(FAST_COEF_temp%co)
        NULLIFY(FAST_COEF_temp%ch4)
        NULLIFY(FAST_COEF_temp%so2)

        CALL H5GOPEN_F( LUN, "THERMAL", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP THERMAL")

!DAR read ALL gases in to temp structure that looks like old fast_coef
        CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ COEF%THERMAL")

!DAR add code to move read in data to coef structure
        ALLOCATE (x%thermal(x%fmv_chn), STAT = ERR)
        THROWM(ERR.NE.0, "allocation of thermal fast coefs")
        DO i = 1, x%fmv_chn
          ALLOCATE (x%thermal(i)%gasarray(x%fmv_gas), STAT = ERR)
          THROWM(ERR.NE.0, "allocation of gasarray")
          CALL nullify_gas_coef_pointers(x%thermal(i))
        ENDDO

        IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
          CALL move_hdf_fast_coef(x, x%thermal, FAST_COEF_temp%mixedgas, gas_id_mixed, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
          CALL move_hdf_fast_coef(x, x%thermal, FAST_COEF_temp%watervapour, gas_id_watervapour, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
          CALL move_hdf_fast_coef(x, x%thermal, FAST_COEF_temp%ozone, gas_id_ozone, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
          CALL move_hdf_fast_coef(x, x%thermal, FAST_COEF_temp%wvcont, gas_id_wvcont, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
          CALL move_hdf_fast_coef(x, x%thermal,FAST_COEF_temp%co2, gas_id_co2, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
          CALL move_hdf_fast_coef(x, x%thermal,FAST_COEF_temp%n2o, gas_id_n2o, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
          CALL move_hdf_fast_coef(x, x%thermal,FAST_COEF_temp%co, gas_id_co, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
          CALL move_hdf_fast_coef(x, x%thermal,FAST_COEF_temp%ch4, gas_id_ch4, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
        ENDIF
        IF (ASSOCIATED(FAST_COEF_temp%so2)) THEN
          CALL move_hdf_fast_coef(x, x%thermal,FAST_COEF_temp%so2, gas_id_so2, lbl_mode = LBL)
          DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)
        ENDIF

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")

        IF (x%fmv_model_ver > 9) THEN
          CALL H5GOPEN_F( LUN, "THERMAL_CORR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP THERMAL_CORR")

  !DAR read ALL gases in to temp structure that looks like old fast_coef
          CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ COEF%THERMAL_CORR")

  !DAR add code to move read in data to coef structure
          ALLOCATE (x%thermal_corr(x%fmv_chn), STAT = ERR)
          THROWM(ERR.NE.0, "allocation of thermal fast coefs")
          DO i = 1, x%fmv_chn
            ALLOCATE (x%thermal_corr(i)%gasarray(x%fmv_gas), STAT = ERR)
            THROWM(ERR.NE.0, "allocation of gasarray")
            CALL nullify_gas_coef_pointers(x%thermal_corr(i))
          ENDDO

          IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr, FAST_COEF_temp%mixedgas, gas_id_mixed, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr, FAST_COEF_temp%watervapour, gas_id_watervapour, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr, FAST_COEF_temp%ozone, gas_id_ozone, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr, FAST_COEF_temp%wvcont, gas_id_wvcont, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr,FAST_COEF_temp%co2, gas_id_co2, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr,FAST_COEF_temp%n2o, gas_id_n2o, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr,FAST_COEF_temp%co, gas_id_co, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr,FAST_COEF_temp%ch4, gas_id_ch4, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%so2)) THEN
            CALL move_hdf_fast_coef(x, x%thermal_corr,FAST_COEF_temp%so2, gas_id_so2, &
                                    lbl_mode = LBL, corr = .TRUE._jplm)
            DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)
          ENDIF

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL_CORR")
        ENDIF

        x%PMC_SHIFT = ASSOCIATED(x%PMC_COEF)

        x%SOLARCOEF = .FALSE.
        NULLIFY( x%SOLAR)

        CALL H5LEXISTS_F( LUN, "SOLAR", LEXT, ERR )
        THROWM(ERR.NE.0,"CANNOT TEST SOLAR ")
        IF( LEXT ) THEN

          CALL H5GOPEN_F( LUN, "SOLAR", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP SOLAR")

!DAR read ALL gases in to temp structure that looks like old fast_coef
          CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ COEF%SOLAR")

!DAR add code to move read in data to coef structure
          ALLOCATE (x%solar(x%fmv_chn), STAT = ERR)
          THROWM( ERR .NE. 0, "ALLOCATION OF SOLAR FAST COEFS")
          DO i = 1, x%fmv_chn
            ALLOCATE (x%solar(i)%gasarray(x%fmv_gas), STAT = ERR)
            THROWM(ERR.NE.0, "allocation of gasarray")
            CALL nullify_gas_coef_pointers(x%solar(i))
          ENDDO

          IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%mixedgas, gas_id_mixed, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%watervapour, gas_id_watervapour, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%ozone, gas_id_ozone, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%wvcont, gas_id_wvcont, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%co2, gas_id_co2, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%n2o, gas_id_n2o, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%co, gas_id_co, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%ch4, gas_id_ch4, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
          ENDIF
          IF (ASSOCIATED(FAST_COEF_temp%so2)) THEN
            CALL move_hdf_fast_coef(x, x%solar, FAST_COEF_temp%so2, gas_id_so2, lbl_mode = LBL)
            DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)
          ENDIF

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")

          x%SOLARCOEF = .TRUE.

          IF (x%fmv_model_ver > 9) THEN
            CALL H5GOPEN_F( LUN, "SOLAR_CORR", G_ID_SUB, ERR )
            THROWM(ERR.NE.0,"CANNOT OPEN GROUP SOLAR_CORR")

    !DAR read ALL gases in to temp structure that looks like old fast_coef
            CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
            THROWM(ERR.NE.0,"CANNOT READ COEF%SOLAR_CORR")

    !DAR add code to move read in data to coef structure
            ALLOCATE (x%solar_corr(x%fmv_chn), STAT = ERR)
            THROWM(ERR.NE.0, "allocation of solar fast coefs")
            DO i = 1, x%fmv_chn
              ALLOCATE (x%solar_corr(i)%gasarray(x%fmv_gas), STAT = ERR)
              THROWM(ERR.NE.0, "allocation of gasarray")
              CALL nullify_gas_coef_pointers(x%solar_corr(i))
            ENDDO

            IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr, FAST_COEF_temp%mixedgas, gas_id_mixed, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr, FAST_COEF_temp%watervapour, gas_id_watervapour, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr, FAST_COEF_temp%ozone, gas_id_ozone, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr, FAST_COEF_temp%wvcont, gas_id_wvcont, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr,FAST_COEF_temp%co2, gas_id_co2, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr,FAST_COEF_temp%n2o, gas_id_n2o, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr,FAST_COEF_temp%co, gas_id_co, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr,FAST_COEF_temp%ch4, gas_id_ch4, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
            ENDIF
            IF (ASSOCIATED(FAST_COEF_temp%so2)) THEN
              CALL move_hdf_fast_coef(x, x%solar_corr,FAST_COEF_temp%so2, gas_id_so2, &
                                      lbl_mode = LBL, corr = .TRUE._jplm)
              DEALLOCATE(FAST_COEF_temp%so2)         ; NULLIFY(FAST_COEF_temp%so2)
            ENDIF

            CALL H5GCLOSE_F( G_ID_SUB, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR_CORR")
          ENDIF

        ENDIF


        x%NLTECOEF = .FALSE.
        NULLIFY( x%NLTE_COEF)

        CALL H5LEXISTS_F( LUN, "NLTE_COEF", LEXT, ERR )
        THROWM(ERR.NE.0,"CANNOT TEST NLTE_COEF ")
        IF( LEXT ) THEN

          CALL H5GOPEN_F( LUN, "NLTE_COEF", G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP NLTE_COEF")

          ALLOCATE (x%NLTE_COEF, STAT = ERR)
          THROWM( ERR .NE. 0, "ALLOCATION OF NLTE_COEF FAST COEFS")

          CALL RTTOV_HDF_RTTOV_NLTE_COEF_RH( x%NLTE_COEF, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ COEF%NLTE_COEF")

          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")

          x%NLTECOEF = .TRUE.

        ENDIF

CATCH
      END SUBROUTINE

!> Read a cloud coefficient structure from HDF5 file
!! param[inout]  x             RTTOV optical property coefficient structure
!! param[in]     z             optical depth coefficient structure
!! param[in]     lun           file ID of HDF5 file
!! param[out]    err           return status
    SUBROUTINE RTTOV_HDF_SCCLDCOEF_RH(X, Z, LUN, ERR)
      TYPE(RTTOV_COEF_SCATT), INTENT(INOUT) :: X
      TYPE(RTTOV_COEF),       INTENT(IN)    :: Z
      INTEGER(HID_T),         INTENT(IN)    :: LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)   :: ERR
      LOGICAL  :: LEXT
      INTEGER(HID_T) :: G_ID_SUB
      TRY

      CALL H5LEXISTS_F( LUN, "WATER_CLOUD_OPAC", LEXT, ERR )
      THROWM(ERR.NE.0,"CANNOT TEST WATER_CLOUD_OPAC")
      IF( LEXT ) THEN
        CALL H5GOPEN_F( LUN, "WATER_CLOUD_OPAC", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP WATER_CLOUD_OPAC")

        CALL RTTOV_HDF_OPTP_RH( X%OPTP_WCL_OPAC, Z, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ WATER_CLOUD_OPAC COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATER_CLOUD_OPAC")
      ENDIF

      CALL H5LEXISTS_F( LUN, "WATER_CLOUD_DEFF", LEXT, ERR )
      THROWM(ERR.NE.0,"CANNOT TEST WATER_CLOUD_DEFF")
      IF( LEXT ) THEN
        CALL H5GOPEN_F( LUN, "WATER_CLOUD_DEFF", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP WATER_CLOUD_DEFF")

        CALL RTTOV_HDF_OPTP_RH( X%OPTP_WCL_DEFF, Z, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ WATER_CLOUD_DEFF COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP WATER_CLOUD_DEFF")
      ENDIF

      CALL H5LEXISTS_F( LUN, "ICE_CLOUD_BAUM", LEXT, ERR )
      THROWM(ERR.NE.0,"CANNOT TEST ICE_CLOUD_BAUM")
      IF( LEXT ) THEN
        CALL H5GOPEN_F( LUN, "ICE_CLOUD_BAUM", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP ICE_CLOUD_BAUM")

        CALL RTTOV_HDF_OPTP_RH( X%OPTP_ICL_BAUM, Z, G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ ICE_CLOUD_BAUM COEFS")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP ICE_CLOUD_BAUM")
      ENDIF

      CATCH
    END SUBROUTINE

!> Read an aerosol coefficient structure from HDF5 file
!! param[inout]  x             RTTOV optical property coefficient structure
!! param[in]     z             optical depth coefficient structure
!! param[in]     lun           file ID of HDF5 file
!! param[out]    err           return status
    SUBROUTINE RTTOV_HDF_SCAERCOEF_RH(X, Z, LUN, ERR)
      TYPE(RTTOV_COEF_SCATT), INTENT(INOUT) :: X
      TYPE(RTTOV_COEF),       INTENT(IN)    :: Z
      INTEGER(HID_T),         INTENT(IN)    :: LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)   :: ERR
      TRY

      CALL RTTOV_HDF_OPTP_RH( X%OPTP_AER, Z, LUN, ERR )
      THROWM(ERR.NE.0,"CANNOT READ AEROSOL COEFS")

      CATCH
    END SUBROUTINE

!> Read a cloud/aerosol optical property structure from HDF5 file
!! param[inout]  x             RTTOV optical property coefficient structure
!! param[in]     z             optical depth coefficient structure
!! param[in]     lun           file ID of HDF5 file
!! param[out]    err           return status
    SUBROUTINE RTTOV_HDF_OPTP_RH(X, Z, LUN,ERR)
      USE RTTOV_UNIX_ENV, ONLY : RTTOV_UPPER_CASE
      USE RTTOV_CONST, ONLY : SCAERCLD_VERSION_COMPATIBLE_MIN, &
                              SCAERCLD_VERSION_COMPATIBLE_MAX
      TYPE(RTTOV_OPTP),       INTENT(INOUT) :: X
      TYPE(RTTOV_COEF),       INTENT(IN)    :: Z
      INTEGER(HID_T),         INTENT(IN)    ::LUN
      INTEGER(KIND=JPIM),     INTENT(OUT)   ::ERR

      CHARACTER(LEN=LENSH)  ::SNAME
      CHARACTER(LEN=LENSH)  ::GNAME
      INTEGER(KIND=JPIM)    :: I
      LOGICAL :: LEXT
      CHARACTER(LEN=4),   POINTER :: C1(:)
      INTEGER(HID_T) :: G_ID_SUB
      TRY

      sname='VERSION'
      call read_array_hdf(lun,sname,err,i0=x%version)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))

      IF (x%version < scaercld_version_compatible_min .OR. &
          x%version > scaercld_version_compatible_max) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "Version of sccld/scaer file is incompatible with RTTOV library")
      ENDIF

      sname='ID'
      call read_array_hdf(lun,sname,err,i0=x%id)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))

      sname='NCHAN'
      call read_array_hdf(lun,sname,err,i0=x%nchan)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))

      IF (z%fmv_ori_nchn /= x%nchan) THEN
          err = errorstatus_fatal
          THROWM(err.ne.0,"Incompatible channels between rtcoef and sccoef files")
      ENDIF

      sname='NCHAN_PHA'
      call read_array_hdf(lun,sname,err,i0=x%nchan_pha)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))

      sname='NTYPES'
      call read_array_hdf(lun,sname,err,i0=x%ntypes)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))

      sname='MAXNMOM'
      call read_array_hdf(lun,sname,err,i0=x%maxnmom)
      THROWM(err.ne.0,"CANNOT READ "//trim(sname))


      ! Sort out the solar channels/phase functions
      IF (x%nchan_pha > 0) THEN

        sname='NPHANGLE'
        call read_array_hdf(lun,sname,err,i0=x%nphangle)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        ! Always get solar channel indexes from HDF file
        sname='CHAN_PHA'
        call read_array_hdf(lun,sname,err,pi1=x%chan_pha)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        ! Create map from extracted channel list into pha array
        ALLOCATE (x%chan_pha_index(z%fmv_chn))
        x%chan_pha_index(:) = 0
        x%chan_pha_index(x%chan_pha(1:x%nchan_pha)) = &
              (/ (i, i = 1, x%nchan_pha) /)

        sname='PHANGLE'
        call read_array_hdf(lun,sname,err,pr1=x%phangle)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        CALL rttov_alloc_phfn_int(err, x%phangle, x%phfn_int, 1_jpim)
        THROWM(err.NE.0, "initialisation of optp%phfn_int")

      ENDIF

      ALLOCATE (x%data(x%ntypes), STAT = ERR)
      THROWM( ERR .NE. 0, "allocation of optp%data")

      DO i = 1, x%ntypes
        CALL rttov_hdf_nullify_optp_data(x%data(i))
      ENDDO

      IF (x%ntypes > 1) THEN
        sname='NAMES'
        call read_array_hdf(lun,sname,err,pc1=c1)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))
        x%data(:)%name = c1
        deallocate(c1)
      ENDIF

      DO I = 1, x%ntypes

        IF (x%ntypes > 1) THEN
          CALL RTTOV_UPPER_CASE(GNAME, x%data(i)%name)
        ELSE
          GNAME = 'DATA'
        ENDIF

        CALL H5GOPEN_F( LUN, TRIM(GNAME), G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))

        sname='NRELHUM'
        call read_array_hdf(G_ID_SUB,sname,err,i0=x%data(i)%nrelhum)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='RELHUM'
        call read_array_hdf(G_ID_SUB,sname,err,pr1=x%data(i)%relhum)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='CONFAC'
        call read_array_hdf(G_ID_SUB,sname,err,r0=x%data(i)%confac)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='NDEFF'
        call read_array_hdf(G_ID_SUB,sname,err,i0=x%data(i)%ndeff)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='DEFF'
        call read_array_hdf(G_ID_SUB,sname,err,pr1=x%data(i)%deff)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='ABS'
        call read_array_hdf(G_ID_SUB,sname,err,pr3=x%data(i)%abs)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='SCA'
        call read_array_hdf(G_ID_SUB,sname,err,pr3=x%data(i)%sca)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        sname='BPR'
        call read_array_hdf(G_ID_SUB,sname,err,pr3=x%data(i)%bpr)
        THROWM(err.ne.0,"CANNOT READ "//trim(sname))

        IF (x%maxnmom > 0) THEN
          sname='NMOM'
          call read_array_hdf(G_ID_SUB,sname,err,pi2=x%data(i)%nmom)
          THROWM(err.ne.0,"CANNOT READ "//trim(sname))

          sname='LEGCOEF'
          call read_array_hdf(G_ID_SUB,sname,err,pr4=x%data(i)%legcoef)
          THROWM(err.ne.0,"CANNOT READ "//trim(sname))
        endif

        sname='PHA'
        call h5lexists_f( G_ID_SUB, sname, lext, err )
        if( lext ) then
          call read_array_hdf(G_ID_SUB,sname,err,pr4=x%data(i)%pha)
          THROWM(err.ne.0,"CANNOT READ "//trim(sname))
        endif

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(GNAME))
      ENDDO

      CATCH
    END SUBROUTINE

!> Read PC-RTTOV coefficient structure from HDF5 file
!! param[inout]  x             PC coefficient structure
!! param[in]     y             optical depth coefficient structure
!! param[in]     lun           file ID of HDF5 file
!! param[out]    err           return status
      SUBROUTINE RTTOV_HDF_PCCOEF_RH(X, Y,LUN,ERR)
USE RTTOV_CONST, ONLY: ERRORSTATUS_FATAL, SENSOR_ID_HI, &
                       PCCOEF_VERSION_COMPATIBLE_MIN, &
                       PCCOEF_VERSION_COMPATIBLE_MAX

USE RTTOV_HDF_RTTOV_COEF_PCC_IO
USE RTTOV_HDF_RTTOV_COEF_PCC1_IO
USE RTTOV_HDF_RTTOV_COEF_PCC2_IO

      TYPE(RTTOV_COEF_PCCOMP),INTENT(INOUT)  ::X
      TYPE(RTTOV_COEF),INTENT(IN)            ::Y
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR

      CHARACTER(LEN=LENSH)  :: GNAME, GNAME2
      INTEGER(KIND=JPIM)    :: I, J
!
      INTEGER(HID_T) :: G_ID_SUB, G_ID_SUB2, G_ID_SUB3
!
TRY

        ! First check the file format version
        CALL read_array_hdf(LUN,'FMV_PC_COMP_PC',ERR,i0=x%FMV_PC_COMP_PC)
        THROWM(ERR.NE.0,"CANNOT READ "//trim('FMV_PC_COMP_PC'))
        IF (x%fmv_pc_comp_pc < pccoef_version_compatible_min .OR. &
            x%fmv_pc_comp_pc > pccoef_version_compatible_max) THEN
          err = errorstatus_fatal
          THROWM(ERR.NE.0, "Version of PC coef file is incompatible with RTTOV library")
        ENDIF

        IF (y%id_sensor == sensor_id_hi) THEN
          IF (y%id_comp_pc /= x%fmv_pc_comp_pc) THEN
            err = errorstatus_fatal
            THROWM( ERR .NE. 0, "Version of PC coef file is incompatible with RTTOV regression file")
          ENDIF
        ELSE        
          err = errorstatus_fatal
          THROWM(err.NE.0, "PC-RTTOV only compatible with hi-res sounders")
        ENDIF

        CALL RTTOV_HDF_RTTOV_COEF_PCC_RH(X,LUN,ERR)
        THROWM(ERR.NE.0,"CANNOT READ COEF")

        !PCREG structure
        CALL H5GOPEN_F( LUN, "PCREG", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP PCREG")

        ALLOCATE (X%PCREG(X%FMV_PC_BANDS,X%FMV_PC_MSETS), STAT = ERR)
        THROWM( ERR .NE. 0, "allocation of coef_pccomp %pcreg")

        DO J = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') J
          CALL H5GOPEN_F( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP PCREG/"//TRIM(GNAME))

          DO I = 1, X%FMV_PC_SETS(J)

            WRITE(GNAME2,'(I2.2)') I
            CALL H5GOPEN_F( G_ID_SUB2, TRIM(GNAME2), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT OPEN GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_INIT(X%PCREG(j,i))

            CALL RTTOV_HDF_RTTOV_COEF_PCC1_RH( X%PCREG(j,i), G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT READ PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

            CALL H5GCLOSE_F( G_ID_SUB3, ERR )
            THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME)//'/'//TRIM(GNAME2))

          END DO

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP PCREG")

        !EIGEN structure
        CALL H5GOPEN_F( LUN, "EIGEN", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(GNAME))

        ALLOCATE (X%EIGEN(X%FMV_PC_BANDS), STAT = ERR)
        THROWM( ERR .NE. 0, "allocation of coef_pccomp %eigen")


        DO I = 1, X%FMV_PC_BANDS

          WRITE(GNAME,'(I2.2)') I
          CALL H5GOPEN_F( G_ID_SUB, TRIM(GNAME), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP EIGEN/"//TRIM(GNAME))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_INIT(X%EIGEN(i))

          CALL RTTOV_HDF_RTTOV_COEF_PCC2_RH( X%EIGEN(i), G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT READ EIGEN/"//TRIM(GNAME))

          CALL H5GCLOSE_F( G_ID_SUB2, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN/"//TRIM(GNAME))

        END DO

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP EIGEN")

CATCH
      END SUBROUTINE

!> Read MFASIS LUT structure from HDF5 file
!! param[inout]  x             MFASIS LUT structure
!! param[in]     file_type     CLD or AER for cloud or aerosol structure
!! param[in]     lun           file ID of HDF5 file
!! param[out]    err           return status
      SUBROUTINE RTTOV_HDF_MFASISCOEF_RH(X, FILE_TYPE, LUN,ERR)
        USE RTTOV_CONST, ONLY: ERRORSTATUS_FATAL
        USE RTTOV_CONST, ONLY : MFASISLUT_VERSION_COMPATIBLE_MIN, &
                                MFASISLUT_VERSION_COMPATIBLE_MAX

        TYPE(RTTOV_COEF_MFASIS), INTENT(INOUT) :: X
        CHARACTER(LEN=3),        INTENT(IN)    :: FILE_TYPE
        INTEGER(HID_T),          INTENT(IN)    :: LUN
        INTEGER(KIND=JPIM),      INTENT(OUT)   :: ERR

        CHARACTER(LEN=LENSH) :: GNAME
        INTEGER(KIND=JPIM)   :: I

        TRY

        CALL read_array_hdf(LUN,'VERSION',ERR,I0=X%VERSION)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        IF (x%version < mfasislut_version_compatible_min .OR. &
            x%version > mfasislut_version_compatible_max) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "Version of MFASIS LUT file is incompatible with RTTOV library")
        ENDIF

        CALL read_array_hdf(LUN,'NCHANNELS_COEF',ERR,I0=X%NCHANNELS_COEF)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'FILE_TYPE',ERR,I0=X%FILE_TYPE)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'README_LUT',ERR,C1=X%README_LUT)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'NCHANNELS',ERR,I0=X%NCHANNELS)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'CHANNEL_LIST',ERR,PI1=X%CHANNEL_LIST)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'CHANNEL_LUT_INDEX',ERR,PI1=X%CHANNEL_LUT_INDEX)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'NDIMS',ERR,I0=X%NDIMS)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        CALL read_array_hdf(LUN,'NPARTICLES',ERR,I0=X%NPARTICLES)
        THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

        IF (FILE_TYPE == 'AER') THEN
          X%CLW_SCHEME = 0
          X%ICE_SCHEME = 0
          CALL read_array_hdf(LUN,'AER_TYPES',ERR,PI1=X%AER_TYPES)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
        ELSE
          NULLIFY(X%AER_TYPES)
          CALL read_array_hdf(LUN,'CLW_SCHEME',ERR,I0=X%CLW_SCHEME)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
          CALL read_array_hdf(LUN,'ICE_SCHEME',ERR,I0=X%ICE_SCHEME)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
        ENDIF

        IF (x%version > 0) THEN
          CALL read_array_hdf(LUN,'MAXZENANGLE',ERR,R0=X%MAXZENANGLE)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
        ELSE ! coef_mfasis%version = 0
          X%MAXZENANGLE = mfasis_maxzenangle
        ENDIF

        ALLOCATE(X%LUT_AXES(X%NDIMS), STAT=ERR)
        THROWM(ERR.NE.0,"Failure allocating MFASIS LUT")

        DO I = 1, X%NDIMS
          WRITE(GNAME, '(A,I3.3)') 'DIM', I

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/NAME',ERR,C0=X%LUT_AXES(I)%NAME)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/DIM_TYPE',ERR,I0=X%LUT_AXES(I)%DIM_TYPE)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/NVALUES',ERR,I0=X%LUT_AXES(I)%NVALUES)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/VALUES',ERR,PR1=X%LUT_AXES(I)%VALUES)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
        ENDDO

        ALLOCATE(X%LUT(X%NCHANNELS), STAT=ERR)
        THROWM(ERR.NE.0,"Failure allocating MFASIS LUT")

        DO I = 1, X%NCHANNELS
          WRITE(GNAME, '(A,I6.6)') 'CH', I

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/NLUTS',ERR,I0=X%LUT(I)%NLUTS)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/QINT',ERR,PR2=X%LUT(I)%QINT)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")

          CALL read_array_hdf(LUN,TRIM(GNAME)//'/DATA',ERR,PR2=X%LUT(I)%DATA)
          THROWM(ERR.NE.0,"Failure reading MFASIS LUT")
        ENDDO

        CATCH
      END SUBROUTINE RTTOV_HDF_MFASISCOEF_RH

!> Nullify/initialise an RTTOV aerosol/cloud optical property data structure
!! param[inout]  data             aerosol/cloud optical property data structure
      SUBROUTINE rttov_hdf_nullify_optp_data(data)
        USE parkind1, ONLY : jprb
        USE rttov_types, ONLY : rttov_optp_data
        IMPLICIT NONE
        TYPE(rttov_optp_data), INTENT(INOUT) :: data

        data%name = ''
        data%confac = 0._jprb
        data%nrelhum = 0
        data%ndeff = 0

        NULLIFY(data%relhum, data%deff, &
                data%abs, data%sca, data%bpr, &
                data%nmom, data%legcoef, data%pha)

      END SUBROUTINE rttov_hdf_nullify_optp_data

      !> Transfer fast coefficients from format in HDF5 files to coefficient structure
      !! @param[inout]    coef          optical depth coefficient structure
      !! @param[inout]    fast_coef     fast coefficient structure in coef
      !! @param[in]       gas_in        fast coefficients read from HDF5 file
      !! @param[in]       gas_id        gas ID of coefficients to transfer
      !! @param[in]       lbl_mode      set this to true only if reading optical depth coefs from LBL code
      !!                                  (forces reading of all coefs even where they are all zero)
      !! @param[in]       corr          flag to indicate correction coefs or gas opdep coefs
      SUBROUTINE move_hdf_fast_coef(coef, fast_coef, gas_in, gas_id, lbl_mode, corr)

        TYPE(rttov_coef), INTENT(INOUT) :: coef
        TYPE(rttov_fast_coef), INTENT(INOUT) :: fast_coef(:)
        REAL(jprb), INTENT(IN) :: gas_in(:,:,:)
        INTEGER(jpim), INTENT(IN) :: gas_id
        LOGICAL(jplm), INTENT(IN) :: lbl_mode
        LOGICAL(jplm), INTENT(IN), OPTIONAL :: corr

        INTEGER(jpim) :: gas_pos, ncoef
        INTEGER(jpim) :: i

        gas_pos = coef%fmv_gas_pos(gas_id)
        ncoef = coef%fmv_coe(gas_pos)
        IF (PRESENT(corr)) THEN
          IF (corr) ncoef = coef%fmv_ncorr(gas_pos)
        ENDIF

        IF (coef%fmv_model_ver <= 9) THEN
          DO i = 1, coef%fmv_chn
            IF (lbl_mode .OR. ANY(gas_in(:,i,:) /= 0._jprb)) THEN
              ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef, coef%nlayers))
              fast_coef(i)%gasarray(gas_pos)%coef = TRANSPOSE(gas_in(:,i,:))
              CALL set_pointers(fast_coef(i), gas_pos, gas_id)
            ELSE
              NULLIFY(fast_coef(i)%gasarray(gas_pos)%coef)
            ENDIF
          ENDDO
        ELSE
          DO i = 1, coef%fmv_chn
            IF (lbl_mode .OR. ANY(gas_in(:,:,i) /= 0._jprb)) THEN
              ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef, coef%nlayers))
              fast_coef(i)%gasarray(gas_pos)%coef = gas_in(:,:,i)
              CALL set_pointers(fast_coef(i), gas_pos, gas_id)
            ELSE
              NULLIFY(fast_coef(i)%gasarray(gas_pos)%coef)
            ENDIF
          ENDDO
        ENDIF

      END SUBROUTINE move_hdf_fast_coef

      !> Transfer fast coefficients from coefficient structure to format for HDF5 files
      !! @param[in]       coef          optical depth coefficient structure
      !! @param[in]       fast_coef     fast coefficient structure in coef
      !! @param[inout]    gas_in        fast coefficients to write to HDF5 file
      !! @param[in]       gas_id        gas ID of coefficients to transfer
      !! @param[in]       corr          flag to indicate correction coefs or gas opdep coefs
      SUBROUTINE move_fast_coef_hdf(coef, fast_coef, gas_in, gas_id, corr)
        TYPE(rttov_coef), INTENT(IN) :: coef
        TYPE(rttov_fast_coef), INTENT(IN) :: fast_coef(:)
        REAL(jprb), INTENT(INOUT), POINTER :: gas_in(:,:,:)
        INTEGER(jpim), INTENT(IN) :: gas_id
        LOGICAL(jplm), INTENT(IN), OPTIONAL :: corr

        INTEGER(jpim) :: i, gas_pos, ncoef

        gas_pos = coef%fmv_gas_pos(gas_id)
        ncoef = coef%fmv_coe(gas_pos)
        IF (PRESENT(corr)) THEN
          IF (corr) ncoef = coef%fmv_ncorr(gas_pos)
        ENDIF

        IF (coef%fmv_model_ver <= 9) THEN
          ALLOCATE(gas_in(coef%nlayers,coef%fmv_chn,ncoef))
          DO i = 1, coef%fmv_chn
            IF (ASSOCIATED(fast_coef(i)%gasarray(gas_pos)%coef)) THEN
              gas_in(:,i,:) = TRANSPOSE(fast_coef(i)%gasarray(gas_pos)%coef)
            ELSE
              gas_in(:,i,:) = 0._jprb
            ENDIF
          ENDDO
        ELSE
          ALLOCATE(gas_in(ncoef,coef%nlayers,coef%fmv_chn))
          DO i = 1, coef%fmv_chn
            IF (ASSOCIATED(fast_coef(i)%gasarray(gas_pos)%coef)) THEN
              gas_in(:,:,i) = fast_coef(i)%gasarray(gas_pos)%coef
            ELSE
              gas_in(:,:,i) = 0._jprb
            ENDIF
          ENDDO
        ENDIF

      END SUBROUTINE move_fast_coef_hdf

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

    END MODULE RTTOV_HDF_COEFS

