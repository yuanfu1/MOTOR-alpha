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
!! @param[in]     coef_scatt      RTTOV cloud/aerosol coefficient structure
!! @param[in]     file_id         logical unit for output sccldcoef file
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
SUBROUTINE rttov_write_binary_sccldcoef(err, coef, coef_scatt, file_id, verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_scatt
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : rttov_magic_string, rttov_magic_number
  USE rttov_cldaer_io_mod, ONLY : write_binary_optp
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),     INTENT(OUT)          :: err
  TYPE(rttov_coef),       INTENT(IN)           :: coef
  TYPE(rttov_coef_scatt), INTENT(IN)           :: coef_scatt
  INTEGER(KIND=jpim),     INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm),     INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"

  LOGICAL(KIND=jplm) :: section_present
  LOGICAL(KIND=jplm) :: lverbose
  CHARACTER(LEN=80)  :: errMessage
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")') file_id
    INFO(errMessage)
  END IF

  ! Write a string that could be displayed
  ! Write a real number to be able to check single/double precision
  WRITE (file_id, iostat=err) rttov_magic_string, rttov_magic_number
  THROW(err.NE.0)

  ! Write an identifying string for cloud files (to avoid mixing up with aerosol files)
  WRITE (file_id, iostat=err) 'sccld_coef'
  THROW(err.NE.0)

  ! OPAC water clouds
  section_present = (coef_scatt%optp_wcl_opac%ntypes > 0)
  WRITE (file_id, iostat=err) section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    CALL write_binary_optp(err, file_id, coef, coef_scatt%optp_wcl_opac)
    THROW(err.NE.0)
  ENDIF

  ! Deff water clouds
  section_present = (coef_scatt%optp_wcl_deff%ntypes > 0)
  WRITE (file_id, iostat=err) section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    CALL write_binary_optp(err, file_id, coef, coef_scatt%optp_wcl_deff)
    THROW(err.NE.0)
  ENDIF

  ! Baum ice clouds
  section_present = (coef_scatt%optp_icl_baum%ntypes > 0)
  WRITE (file_id, iostat=err) section_present
  THROW(err.NE.0)

  IF (section_present) THEN
    CALL write_binary_optp(err, file_id, coef, coef_scatt%optp_icl_baum)
    THROW(err.NE.0)
  ENDIF

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE
