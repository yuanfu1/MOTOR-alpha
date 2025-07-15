! Description:
!> @file
!!   Construct an RTTOV coefficient filename based on the platform,
!!   satellite and instrument ID triplet.
!
!> @brief
!!   Construct an RTTOV coefficient filename based on the platform,
!!   satellite and instrument ID triplet.
!!
!! @details
!!   This subroutine constructs default RTTOV filenames from the
!!   satellite/instrument IDs.
!!
!!   The instrument argument is a triplet consisting of the platform
!!   ID, satellite ID and instrument ID.
!!
!!   For example, for MSG-3 SEVIRI would the triplet would be (/12, 3, 21/)
!!   where the platform and instrument IDs are taken from rttov_const.
!!
!!   The filetype argument consists of the filename prefix: rtcoef, scaercoef,
!!   sccldcoef or pccoef.
!!
!!   The resulting filename is of the form "rtcoef_platform_satellite_inst",
!!   for example.
!!
!!   Note that no file extension is added. Also note that some RTTOV
!!   coefficient filenames include additional information (e.g. tags to denote
!!   Zeeman-compatible coefficients, PC-compatible optical depth coefficients,
!!   coefficients with shifted spectral responses/passbands, and so on). These
!!   are not accounted for.
!!
!! @param[out]     err              status on exit
!! @param[in]      instrument       instrument triplet (platform, satellite, instrument)
!! @param[in]      filetype         string containing file prefix
!! @param[out]     coeffname        output filename excluding file extension
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
SUBROUTINE rttov_coeffname(err, instrument, filetype, coeffname)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : &
      nplatforms,         &
      ninst,              &
      inst_name,          &
      platform_name
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err            ! return code
  INTEGER(KIND=jpim), INTENT(IN)  :: instrument(3)  ! (platform, sat_id, inst) numbers
  CHARACTER(LEN=*),   INTENT(IN)  :: filetype       ! file type e.g. "rtcoef", "pccoef", etc
  CHARACTER(LEN=*),   INTENT(OUT) :: coeffname      ! filename of the coefficient file (excluding extension)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: platform
  INTEGER(KIND=jpim) :: sat_id
  INTEGER(KIND=jpim) :: inst
  CHARACTER(LEN=2)   :: ch_sat_id
!- End of header --------------------------------------------------------

  TRY

  coeffname = 'no_name'
  err = errorstatus_success

  ! Expand instrument triplet
  platform = instrument(1)
  sat_id   = instrument(2)
  inst     = instrument(3)

  ! Test sat_id and convert to string
  IF (sat_id < 10 .AND. sat_id >= 0) THEN
    ! one digit
    WRITE(ch_sat_id,'(i1)') sat_id
  ELSE IF (sat_id >= 10 .AND. sat_id < 99) THEN
    ! two digits
    WRITE(ch_sat_id,'(i2)') sat_id
  ELSE
    ! ERROR and exit
    err = sat_id
    THROWM(err .NE. 0,"invalid sat_id")
  ENDIF

  ! Test platform number
  IF (platform <= 0 .OR. platform > nplatforms) err = errorstatus_fatal
  THROWM(err .NE. 0,"invalid platform number")

  ! Test instrument number  (0 is HIRS)
  IF (inst < 0 .OR. inst > ninst) err = errorstatus_fatal
  THROWM(err .NE. 0,"invalid instrument number")

  ! Create the file name e.g. "rtcoef_platform_satellite_inst"
  coeffname = TRIM(filetype) // '_'      // &
         & TRIM(platform_name(platform)) // &
         & '_'                           // &
         & TRIM(ch_sat_id)               // &
         & '_'                           // &
         & TRIM(inst_name(inst))

  CATCH
END SUBROUTINE rttov_coeffname
