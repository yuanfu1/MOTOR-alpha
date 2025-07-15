! Description:
!> @file
!!   Write out fatal error and informational messages to the logical unit
!!   specified by rttov_errorhandling.
!
!> @brief
!!   Write out fatal error and informational messages to the logical unit
!!   specified by rttov_errorhandling.
!!
!! @details
!!   If rttov_errorhandling was not previously called to set the output unit
!!   the RTTOV default unit is used.
!!
!! @param[in] error_status      errorstatus_success or errorstatus_fatal (defined in rttov_const)
!! @param[in] error_message     message to output
!! @param[in] name_of_routine   name of subroutine calling this one, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_errorreport(error_status, error_message, name_of_routine)
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : errorstatus_success
  USE rttov_global, ONLY : error_unit
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),           INTENT(IN) :: error_status
  CHARACTER(LEN=*),             INTENT(IN) :: error_message
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN) :: name_of_routine
!INTF_END

  CHARACTER(LEN=8)   :: date
  CHARACTER(LEN=10)  :: time
  CHARACTER(LEN=21)  :: datetime
  !- End of header --------------------------------------------------------

  CALL DATE_AND_TIME(date, time)

  WRITE(datetime,"(1X,a4,2('/',a2),2x,2(a2,':'),a2)")   &
              & date(1:4), date(5:6), date(7:8),        &
              & time(1:2), time(3:4), time(5:6)

  IF (error_status == errorstatus_success) THEN

    ! Standard warning or informational message
    IF (PRESENT(name_of_routine)) THEN
      WRITE(error_unit,"(3a)") TRIM(datetime), "  ", TRIM(name_of_routine)
      WRITE(error_unit,"(5X,A)") TRIM(error_message)
    ELSE
      WRITE(error_unit,"(3a)") TRIM(datetime), "  ", TRIM(error_message)
    END IF

  ELSE

    ! fatal error message
    IF (PRESENT(name_of_routine)) THEN
      WRITE(error_unit,"(3a)") TRIM(datetime), "  fatal error in module ", TRIM(name_of_routine)
    ELSE
      WRITE(error_unit,"(2a)") TRIM(datetime), "  fatal error"
    END IF
    WRITE(error_unit,"(5X,A)") TRIM(error_message)

  END IF

END SUBROUTINE rttov_errorreport
