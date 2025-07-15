! Description:
!> @file
!!   Read a binary cloud coefficient file, optionally extracting a subset of channels.
!
!> @brief
!!   Read a binary cloud coefficient file, optionally extracting a subset of channels.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!!   Note that after reading a subset of channels RTTOV will identify them by
!!   indexes 1...SIZE(channels), not by the original channel numbers.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in,out] coef_scatt      RTTOV cloud/aerosol coefficient structure
!! @param[in]     file_id         logical unit for input sccldcoef file
!! @param[in]     channels        list of channels to read, optional
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
SUBROUTINE rttov_read_binary_sccldcoef(err, coef, coef_scatt, file_id, channels)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_scatt
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : rttov_magic_string, rttov_magic_number, errorstatus_fatal
  USE parkind1, ONLY : jprb, jplm
  USE rttov_cldaer_io_mod, ONLY : read_binary_optp
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),     INTENT(OUT)          :: err
  TYPE(rttov_coef),       INTENT(IN)           :: coef
  TYPE(rttov_coef_scatt), INTENT(INOUT)        :: coef_scatt
  INTEGER(KIND=jpim),     INTENT(IN)           :: file_id
  INTEGER(KIND=jpim),     INTENT(IN), OPTIONAL :: channels(:)
!INTF_END
#include "rttov_errorreport.interface"

  LOGICAL(KIND=jplm) :: section_present
  CHARACTER(LEN=16)  :: bin_check_string
  REAL(KIND=jprb)    :: bin_check_number
  REAL(KIND=jprb)    :: bin_check_value
  CHARACTER(LEN=10)  :: filetype_check_string
!- End of header --------------------------------------------------------
  TRY

  READ (file_id, iostat=err) bin_check_string, bin_check_number
  THROWM(err.NE.0,'io status while reading header')

  ! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal
  THROWM(err.NE.0,'Wrong header string in file')

  ! Verification of single/double precision using a 5 digit number
  ! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._jprb - ABS(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_jprb .OR. bin_check_value < 0.99_jprb) err = errorstatus_fatal
  THROWM(err.NE.0,'File created with a different real precision (R4<->R8)')

  ! Verify this is a cloud file (to avoid mixing up with aerosol files)
  READ (file_id, iostat=err) filetype_check_string
  IF (TRIM(filetype_check_string) /= 'sccld_coef') err = errorstatus_fatal
  THROWM(err.NE.0,'This is not a cloud coefficient file')

  ! OPAC water clouds
  READ (file_id, iostat=err) section_present
  IF (err == 0 .AND. section_present) THEN
    CALL read_binary_optp(err, file_id, coef, coef_scatt%optp_wcl_opac, channels)
    THROW(err.NE.0)
  ENDIF

  ! Deff water clouds
  READ (file_id, iostat=err) section_present
  IF (err == 0 .AND. section_present) THEN
    CALL read_binary_optp(err, file_id, coef, coef_scatt%optp_wcl_deff, channels)
    THROW(err.NE.0)
  ENDIF

  ! Baum ice clouds
  READ (file_id, iostat=err) section_present
  IF (err == 0 .AND. section_present) THEN
    CALL read_binary_optp(err, file_id, coef, coef_scatt%optp_icl_baum, channels)
    THROW(err.NE.0)
  ENDIF

  CATCH
END SUBROUTINE rttov_read_binary_sccldcoef
