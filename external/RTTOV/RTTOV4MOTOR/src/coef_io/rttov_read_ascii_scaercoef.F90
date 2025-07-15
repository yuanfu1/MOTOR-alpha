! Description:
!> @file
!!   Read an ASCII aerosol coefficient file, optionally extracting a subset of channels.
!
!> @brief
!!   Read an ASCII aerosol coefficient file, optionally extracting a subset of channels.
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
!! @param[in]     file_id         logical unit for input scaercoef file
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
SUBROUTINE rttov_read_ascii_scaercoef(err, coef, coef_scatt, file_id, channels)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_scatt
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : lensection, errorstatus_fatal
  USE rttov_cldaer_io_mod, ONLY : read_ascii_optp
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),          INTENT(OUT)          :: err
  TYPE(rttov_coef),       INTENT(IN)           :: coef
  TYPE(rttov_coef_scatt), INTENT(INOUT)        :: coef_scatt
  INTEGER(jpim),          INTENT(IN)           :: file_id
  INTEGER(jpim),          INTENT(IN), OPTIONAL :: channels(:)
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_findnextsection.interface"

  INTEGER(jpim) :: io_status
  CHARACTER(LEN=lensection) :: section
!- End of header --------------------------------------------------------
  TRY

  readfile : DO
    CALL rttov_findnextsection(file_id, io_status, section)
    IF (io_status < 0) EXIT !end-of-file

    SELECT CASE (TRIM(section))

    CASE ('WATER_CLOUD_OPAC', 'WATER_CLOUD_DEFF', 'ICE_CLOUD_BAUM')
      ! Aerosol/cloud files have very similar formats which could cause confusion
      ! (reported by v12 beta tester). This check helps prevent that.
      err = errorstatus_fatal
      THROWM(err.NE.0, "Trying to read cloud coefficient file as an aerosol file")

    CASE ('AEROSOLS')

      CALL read_ascii_optp(err, file_id, section, coef, coef_scatt%optp_aer, channels)
      THROW(err.NE.0)

    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_scaercoef
