! Description:
!> @file
!!   Skip comment lines (starting '!') in an ASCII coefficient file.
!
!> @brief
!!   Skip comment lines (starting '!') in an ASCII coefficient file.
!!
!! @details
!!   File pointer is left at the beginning of the line following the comments.
!!
!! @param[in]        fileunit     open fileunit for coefficient file
!! @param[out]       readstatus   returns iostat from most recent READ
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
SUBROUTINE rttov_skipcommentline(fileunit, readstatus)

  USE parkind1, ONLY : jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)  :: fileunit
  INTEGER(KIND=jpim), INTENT(OUT) :: readstatus
!INTF_END

  CHARACTER(LEN=80) :: line
!- End of header --------------------------------------------------------

  readfile: DO

     READ(UNIT=fileunit, fmt='(a)', IOSTAT=readstatus) line
     IF (readstatus /= 0) EXIT

     line = ADJUSTL(line)
     IF (line(1:1) == '!' .OR. line == '') THEN
        CYCLE !skip blank/comment lines
     ELSE
        !reposition file at the start of the line and exit
        BACKSPACE(fileunit)
        EXIT readfile
     ENDIF

  ENDDO readfile

END SUBROUTINE rttov_skipcommentline
