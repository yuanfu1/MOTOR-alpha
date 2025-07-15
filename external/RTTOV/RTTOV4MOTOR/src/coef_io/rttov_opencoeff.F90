! Description:
!> @file
!!   Open a file with the given filename
!
!> @brief
!!   Open a file with the given filename
!!
!! @details
!!   Opens sequential formatted (default) or unformatted file with given
!!   filename for input (default) or output. You can specify a logical unit
!!   to use or set file_id to a value of zero or less and the subroutine
!!   will use the first free logical unit.
!!
!! @param[in]        err          status on exit
!! @param[in]        coeffname    full path of file to open
!! @param[in,out]    file_id      logical unit to use, if <=0 opens file on first free unit
!! @param[in]        for_output   if true file is opened for writing, optional (default: open for reading)
!! @param[in]        lbinary      if true file is assumed sequential unformatted, optional (default:
!!                                assume sequential formatted)
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
SUBROUTINE rttov_opencoeff (&
       & err,        &
       & coeffname,  &
       & file_id,    &
       & for_output, &
       & lbinary     )
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  CHARACTER(LEN=*),   INTENT(IN)           :: coeffname
  INTEGER(KIND=jpim), INTENT(INOUT)        :: file_id
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: for_output
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: lbinary
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=8)   :: file_status
  CHARACTER(LEN=8)   :: file_action
  CHARACTER(LEN=16)  :: file_form
  INTEGER(KIND=jpim) :: file_output  ! 1 for output; 0 for input
  LOGICAL            :: file_open
  LOGICAL            :: existence
  INTEGER(KIND=jpim) :: file_unit
  !- End of header --------------------------------------------------------
  TRY

  file_unit = file_id

  ! Consider file_id argument to determine unit for open
  ! Be careful of the following loop for searching
  ! the first free logical unit. It has been observed that
  ! with some high level compiler options it can have some
  ! side effect, like returning file_id with 0 value.
  IF (file_id <= 0) THEN

    ! Find first free logical unit
    file_unit = 9
    file_open = .TRUE.
    DO
      file_unit = file_unit + 1
      INQUIRE(file_unit, OPENED=file_open)
      IF (.NOT. file_open) EXIT
    ENDDO
  ENDIF

  IF (file_id <= 0 .AND. file_unit >= 9) THEN
    file_id = file_unit
  ENDIF

  ! Consider lbinary option to create the option
  file_form = 'formatted'
  IF (PRESENT(lbinary)) THEN
    IF (lbinary) file_form = 'unformatted'
  ENDIF

  ! Access mode
  file_output = 0
  IF (PRESENT(for_output)) THEN
    IF (for_output) file_output = 1
  ENDIF


  ! Check data file existence

  INQUIRE(FILE=coeffname, EXIST=existence)
  IF (file_output == 0) THEN
    ! If data file does not exist, return an error
    IF (.NOT. existence) err = errorstatus_fatal
    THROWM(err .NE. 0, "Coefficient file"//TRIM(coeffname)//" not found")

    ! Set OPEN keywords for reading
    file_status = 'OLD   '
    file_action = 'READ '
  ELSE
    ! If data file does exist, output a warning message
    IF (existence) THEN
      INFO("Coefficient file"//TRIM(coeffname)//" will be overwritten")
    ENDIF

    ! Set OPEN keywords for writing
    file_status = 'REPLACE'
    file_action = 'WRITE'
  ENDIF


  ! Open the data file

  OPEN (file_id, FILE = coeffname,  &
        STATUS = TRIM(file_status), &
        ACTION = TRIM(file_action), &
        ACCESS = 'SEQUENTIAL',      &
        FORM   = TRIM(file_form),   &
        IOSTAT = err)
  THROWM(err .NE. 0, "Error opening "//TRIM(coeffname))

CATCH
END SUBROUTINE rttov_opencoeff
