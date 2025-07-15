! Description:
!> @file
!!   Write coefficient file(s) from RTTOV coefficients structure.
!
!> @brief
!!   Write coefficient file(s) from RTTOV coefficients structure.
!!
!! @details
!!   This is called by the coefficient generation software and the
!!   coefficient conversion and other similar subroutines.
!!
!!   This subroutine always writes an optical depth "rtcoef" coefficient
!!   file. It can also optionally write additional RTTOV input files
!!   containing cloud or aerosol optical properties, MFASIS look-up
!!   tables (LUTs) or PC-RTTOV coefficients.
!!
!!   The output format(s) must be specified for each type of file
!!   being written: they do not have to be in the same format.
!!
!!   Coefficients can be written by specifying the filename(s)
!!   explicitly, by passing in pre-opened logical unit(s) or by
!!   specifying the instrument triplet in which case RTTOV constructs
!!   the filename(s) automatically from the platform and instrument
!!   names in rttov_const.F90. In this latter case the files are written
!!   to the current directory.
!!
!!   Cloud coefficient (sccldcoef) files will only be written if the
!!   addclouds option is true and the user_cld_opt_param option is false.
!!
!!   The MFASIS cloud LUT file will only be written if the addclouds option
!!   is true, the user_cld_opt_param option is false, and MFASIS is specified
!!   in the vis_scatt_model option.
!!
!!   Aerosol coefficient (scaercoef) files will only be written if the
!!   addaerosl option is true and the user_aer_opt_param option is false.
!!
!!   PC coefficient (pccoef) files will only be written if the addpc option
!!   is true.
!!
!! @param[out]    err                status on exit
!! @param[in]     coefs              RTTOV coefs structure
!! @param[in]     opts               options to configure the simulations
!! @param[in]     form_coef          format of rtcoef file, optional
!! @param[in]     form_scaer         format of IR aerosol scaer file, optional
!! @param[in]     form_sccld         format of IR cloud scaer file, optional
!! @param[in]     form_mfasis_cld    format of MFASIS cloud LUT file, optional
!! @param[in]     form_pccoef        format of PC-RTTOV pccoef file, optional
!! @param[in]     file_coef          file name of rtcoef file, optional
!! @param[in]     file_scaer         file name of IR aerosol scaer file, optional
!! @param[in]     file_sccld         file name of IR cloud scaer file, optional
!! @param[in]     file_mfasis_cld    file name of MFASIS cloud LUT file, optional
!! @param[in]     file_pccoef        file name of PC-RTTOV pccoef file, optional
!! @param[in]     file_id_coef       logical unit for rtcoef file, optional
!! @param[in]     file_id_scaer      logical unit for IR aerosol scaer file, optional
!! @param[in]     file_id_sccld      logical unit for IR cloud scaer file, optional
!! @param[in]     file_id_mfasis_cld logiucal unit of MFASIS cloud LUT file, optional
!! @param[in]     file_id_pccoef     logical unit for PC-RTTOV pccoef file, optional
!! @param[in]     instrument         (platform,satellite,instrument) ID triplet, optional
!! @param[in]     compress           HDF5 only: if true larger data sets will be compressed, optional
!! @param[in]     force_double       HDF5 only: if true large real data sets to be written as 64-bit,
!!                                   if false they will be written as 32-bit,
!!                                   if omitted the real kind specified when the HDF5 lib was opened is used, optional
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
SUBROUTINE rttov_write_coefs(err, coefs, opts,   &
                             form_coef,          &
                             form_scaer,         &
                             form_sccld,         &
                             form_mfasis_cld,    &
                             form_pccoef,        &
                             file_coef,          &
                             file_scaer,         &
                             file_sccld,         &
                             file_mfasis_cld,    &
                             file_pccoef,        &
                             file_id_coef,       &
                             file_id_scaer,      &
                             file_id_sccld,      &
                             file_id_mfasis_cld, &
                             file_id_pccoef,     &
                             instrument,         &
                             compress,           &
                             force_double)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb, jplm
  USE rttov_const, ONLY : vis_scatt_mfasis
  USE rttov_coef_io_mod, ONLY : getlun, closelun
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err
  TYPE(rttov_coefs),  INTENT(IN)  :: coefs
  TYPE(rttov_options),INTENT(IN)  :: opts
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_coef
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_scaer
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_sccld
!   CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_mfasis_aer
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_mfasis_cld
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form_pccoef
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_coef
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_scaer
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_sccld
!   CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_mfasis_aer
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_mfasis_cld
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: file_pccoef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_coef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_scaer
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_sccld
!   INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_mfasis_aer
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_mfasis_cld
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_pccoef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: instrument(3)
  LOGICAL,            INTENT(IN), OPTIONAL :: compress
  LOGICAL,            INTENT(IN), OPTIONAL :: force_double
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"

#include "rttov_write_ascii_coef.interface"
#include "rttov_write_ascii_pccoef.interface"
#include "rttov_write_ascii_scaercoef.interface"
#include "rttov_write_ascii_sccldcoef.interface"
#include "rttov_write_ascii_mfasis_file.interface"
#include "rttov_write_binary_coef.interface"
#include "rttov_write_binary_pccoef.interface"
#include "rttov_write_binary_scaercoef.interface"
#include "rttov_write_binary_sccldcoef.interface"
#include "rttov_write_binary_mfasis_file.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

  INTEGER(KIND=jpim) :: file_id_coef1
  INTEGER(KIND=jpim) :: file_id_scaer1
  INTEGER(KIND=jpim) :: file_id_sccld1
!   INTEGER(KIND=jpim) :: file_id_mfasis_aer1
  INTEGER(KIND=jpim) :: file_id_mfasis_cld1
  INTEGER(KIND=jpim) :: file_id_pccoef1

  CHARACTER(LEN=256) :: file_coef1
  CHARACTER(LEN=256) :: file_scaer1
  CHARACTER(LEN=256) :: file_sccld1
!   CHARACTER(LEN=256) :: file_mfasis_aer1
  CHARACTER(LEN=256) :: file_mfasis_cld1
  CHARACTER(LEN=256) :: file_pccoef1

  character(LEN=32) :: form_coef1
  character(LEN=32) :: form_scaer1
  character(LEN=32) :: form_sccld1
!   CHARACTER(LEN=32) :: form_mfasis_aer1
  CHARACTER(LEN=32) :: form_mfasis_cld1
  character(LEN=32) :: form_pccoef1

  LOGICAL :: create
  REAL(KIND=jprb)   :: ZHOOK_HANDLE

  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 0_jpim, ZHOOK_HANDLE)

  CALL getlun(err, file_id_coef1, file_coef1, form_coef1, .TRUE._jplm, "rtcoef", &
                   file_id_coef, file_coef, form_coef, instrument )
  THROW(err.NE.0)

  IF (rttov_cmpuc(form_coef1,"unformatted")) THEN
    CALL rttov_write_binary_coef(err, coefs%coef, file_id_coef1, opts%config%verbose)
    THROW(err.NE.0)
    CALL closelun(err, file_id_coef1, file_id_coef)
    THROW(err.NE.0)
  ELSE IF (rttov_cmpuc(form_coef1,"formatted")) THEN
    CALL rttov_write_ascii_coef(err, coefs%coef, file_id_coef1, opts%config%verbose)
    THROW(err.NE.0)
    CALL closelun(err, file_id_coef1, file_id_coef)
    THROW(err.NE.0)
  ELSE IF (rttov_cmpuc(form_coef1,"hdf5")) THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
    CALL rttov_hdf_save(err, file_coef1, "/COEF", create=.TRUE., coef=coefs%coef, &
                        compress=compress, force_double=force_double)
    THROW(err.NE.0)
#endif
  ELSE
    err = errorstatus_fatal
    THROWM(err.NE.0,"Unknown format "//TRIM(form_coef1))
  ENDIF

  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN

    IF (PRESENT(instrument) .OR. PRESENT(file_scaer) .OR. PRESENT(file_id_scaer)) THEN

      CALL getlun(err, file_id_scaer1, file_scaer1, form_scaer1, .TRUE._jplm, "scaercoef", &
                       file_id_scaer, file_scaer, form_scaer, instrument )
      THROW(err.NE.0)
      IF (rttov_cmpuc(form_scaer1,"unformatted")) THEN
        CALL rttov_write_binary_scaercoef(err, coefs%coef, coefs%coef_scatt, file_id_scaer1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_scaer1, file_id_scaer)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_scaer1,"formatted")) THEN
        CALL rttov_write_ascii_scaercoef(err, coefs%coef, coefs%coef_scatt, file_id_scaer1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_scaer1, file_id_scaer)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_scaer1,"hdf5")) THEN
#ifndef _RTTOV_HDF
        err = errorstatus_fatal
        THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
        create = .NOT. (TRIM(file_coef) .EQ. TRIM(file_scaer) .AND. TRIM(file_coef) .NE. "")
        CALL rttov_hdf_save(err, file_scaer, "/SCAER", create=create, &
                            scaercoef=coefs%coef_scatt, &
                            compress=compress, force_double=force_double)
        THROW(err.NE.0)
#endif
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,"Unknown format "//TRIM(form_scaer1))
      ENDIF
    ENDIF

!     IF (opts%rt_ir%vis_scatt_model == vis_scatt_mfasis .AND. &
!         (PRESENT(instrument) .OR. PRESENT(file_mfasis_aer) .OR. PRESENT(file_id_mfasis_aer))) THEN
!       CALL getlun(err, file_id_mfasis_aer1, file_mfasis_aer1, form_mfasis_aer1, .TRUE._jplm, "rttov_mfasis_aer", &
!                        file_id_mfasis_aer, file_mfasis_aer, form_mfasis_aer, instrument)
!       THROW(err.NE.0)
!       IF (rttov_cmpuc(form_mfasis_aer1,"unformatted")) THEN
!         CALL rttov_write_binary_mfasis_file(err, coefs%coef_mfasis_aer, &
!                                             file_id_mfasis_aer1, opts%config%verbose)
!         THROW(err.NE.0)
!         CALL closelun(err, file_id_mfasis_aer1, file_id_mfasis_aer)
!         THROW(err.NE.0)
!       ELSE IF (rttov_cmpuc(form_mfasis_aer1,"formatted")) THEN
!         CALL rttov_write_ascii_mfasis_file(err, coefs%coef, coefs%coef_mfasis_aer, &
!                                            file_id_mfasis_aer1, opts%config%verbose)
!         THROW(err.NE.0)
!         CALL closelun(err, file_id_mfasis_aer1, file_id_mfasis_aer)
!         THROW(err.NE.0)
!       ELSE IF (rttov_cmpuc(form_mfasis_aer1,"hdf5")) THEN
! #ifndef _RTTOV_HDF
!         err = errorstatus_fatal
!     THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
! #else
!         create = .NOT. (TRIM(file_coef) .EQ. TRIM(file_mfasis_aer) .AND. TRIM(file_coef) .NE. "")
!         CALL rttov_hdf_save(err, file_mfasis_aer, "/MFASIS_AER", create=create, &
!                     mfasislutaer=coefs%coef_mfasis_aer, &
!                     compress=compress, force_double=force_double)
!         THROW(err.NE.0)
! #endif
!       ELSE
!         err = errorstatus_fatal
!         THROWM(err.NE.0,"Unknown format "//TRIM(form_mfasis_aer1))
!       ENDIF
!     ENDIF
  ENDIF

  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN

    IF (PRESENT(instrument) .OR. PRESENT(file_sccld) .OR. PRESENT(file_id_sccld)) THEN

      CALL getlun(err, file_id_sccld1, file_sccld1, form_sccld1, .TRUE._jplm, "sccldcoef", &
                       file_id_sccld, file_sccld, form_sccld, instrument)
      THROW(err.NE.0)
      IF (rttov_cmpuc(form_sccld1,"unformatted")) THEN
        CALL rttov_write_binary_sccldcoef(err, coefs%coef, coefs%coef_scatt, file_id_sccld1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_sccld1, file_id_sccld)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_sccld1,"formatted")) THEN
        CALL rttov_write_ascii_sccldcoef(err, coefs%coef, coefs%coef_scatt, file_id_sccld1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_sccld1, file_id_sccld)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_sccld1,"hdf5")) THEN
#ifndef _RTTOV_HDF
        err = errorstatus_fatal
        THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
        create = .NOT. (TRIM(file_coef) .EQ. TRIM(file_sccld) .AND. TRIM(file_coef) .NE. "")
        CALL rttov_hdf_save(err, file_sccld, "/SCCLD", create=create, &
                            sccldcoef=coefs%coef_scatt, &
                            compress=compress, force_double=force_double)
        THROW(err.NE.0)
#endif
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,"Unknown format "//TRIM(form_sccld1))
      ENDIF
    ENDIF

    IF (opts%rt_ir%vis_scatt_model == vis_scatt_mfasis .AND. &
        (PRESENT(instrument) .OR. PRESENT(file_mfasis_cld) .OR. PRESENT(file_id_mfasis_cld))) THEN

      CALL getlun(err, file_id_mfasis_cld1, file_mfasis_cld1, form_mfasis_cld1, .TRUE._jplm, "rttov_mfasis_cld", &
                       file_id_mfasis_cld, file_mfasis_cld, form_mfasis_cld, instrument)
      THROW(err.NE.0)
      IF (rttov_cmpuc(form_mfasis_cld1,"unformatted")) THEN
        CALL rttov_write_binary_mfasis_file(err, coefs%coef_mfasis_cld, &
                                            file_id_mfasis_cld1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_mfasis_cld1, file_id_mfasis_cld)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_mfasis_cld1,"formatted")) THEN
        CALL rttov_write_ascii_mfasis_file(err, coefs%coef, coefs%coef_mfasis_cld, &
                                           file_id_mfasis_cld1, opts%config%verbose)
        THROW(err.NE.0)
        CALL closelun(err, file_id_mfasis_cld1, file_id_mfasis_cld)
        THROW(err.NE.0)
      ELSE IF (rttov_cmpuc(form_mfasis_cld1,"hdf5")) THEN
#ifndef _RTTOV_HDF
        err = errorstatus_fatal
        THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
        create = .NOT. (TRIM(file_coef) .EQ. TRIM(file_mfasis_cld) .AND. TRIM(file_coef) .NE. "")
        CALL rttov_hdf_save(err, file_mfasis_cld, "/MFASIS_CLD", create=create, &
                    mfasislutcld=coefs%coef_mfasis_cld, &
                    compress=compress, force_double=force_double)
        THROW(err.NE.0)
#endif
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,"Unknown format "//TRIM(form_mfasis_cld1))
      ENDIF
    ENDIF
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .TRUE._jplm, "pccoef", &
                     file_id_pccoef, file_pccoef, form_pccoef, instrument)
    THROW(err.NE.0)
    IF (rttov_cmpuc(form_pccoef1,"unformatted")) THEN
      CALL rttov_write_binary_pccoef(err, coefs%coef_pccomp, file_id_pccoef1, opts%config%verbose)
      THROW(err.NE.0)
      CALL closelun(err, file_id_pccoef1, file_id_pccoef)
      THROW(err.NE.0)
    ELSE IF (rttov_cmpuc(form_pccoef1,"formatted")) THEN
      CALL rttov_write_ascii_pccoef(err,  coefs%coef_pccomp, file_id_pccoef1, opts%config%verbose)
      THROW(err.NE.0)
      CALL closelun(err, file_id_pccoef1, file_id_pccoef)
      THROW(err.NE.0)
    ELSE IF (rttov_cmpuc(form_pccoef1,"hdf5")) THEN
#ifndef _RTTOV_HDF
      err = errorstatus_fatal
      THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
      create = .NOT. (TRIM(file_coef) .EQ. TRIM(file_pccoef) .AND. TRIM(file_coef) .NE. "")
      CALL rttov_hdf_save(err, file_pccoef, "/PC", create=create, &
                  pccoef=coefs%coef_pccomp, compress=compress, force_double=force_double)
      THROW(err.NE.0)
#endif
    ELSE
      err = errorstatus_fatal
      THROWM(err.NE.0,"Unknown format "//TRIM(form_pccoef1))
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)

END SUBROUTINE
