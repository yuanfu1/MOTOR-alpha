! Description:
!> @file
!!   Write an ASCII aerosol coefficient file.
!
!> @brief
!!   Write an ASCII aerosol coefficient file.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in]     coef_scatt      RTTOV cloud/aerosol coefficient structure
!! @param[in]     file_id         logical unit for output scaercoef file
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
SUBROUTINE rttov_write_ascii_scaercoef(err, coef, coef_scatt, file_id, verbose)
!INTF_OFF
#include "throw.h"
  USE rttov_cldaer_io_mod, ONLY : write_ascii_optp
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_scatt
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim),     INTENT(OUT)          :: err
  TYPE(rttov_coef),       INTENT(IN)           :: coef
  TYPE(rttov_coef_scatt), INTENT(IN)           :: coef_scatt
  INTEGER(KIND=jpim),     INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm),     INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"

  LOGICAL(KIND=jplm) :: lverbose
  CHARACTER(LEN=32)  :: section
  CHARACTER(LEN=80)  :: errMessage
  CHARACTER(LEN=*), PARAMETER :: routinename = 'rttov_write_ascii_scaercoef'
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")') file_id
    INFO(errMessage)
  END IF

  WRITE (file_id, '(a)', iostat=err) ' ! RTTOV coefficient file '//TRIM(coef%id_common_name)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' ! Automatic creation by subroutine '//routinename
  THROW(err.NE.0)


  section = 'AEROSOLS'
  WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' !'
  THROW(err.NE.0)

  CALL write_ascii_optp(err, file_id, coef, coef_scatt%optp_aer, 'aerosol', routinename)
  THROW(err.NE.0)

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
