! Description:
!> @file
!!   Read ASCII coefficient file until the next section is reached.
!
!> @brief
!!   Read ASCII coefficient file until the next section is reached.
!!
!! @details
!!   File pointer is left at the beginning of the section.
!!
!! @param[in]        fileunit     open fileunit for coefficient file
!! @param[out]       readstatus   returns iostat from most recent READ
!! @param[out]       section      returns the name of the section that was found (if any)
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
SUBROUTINE rttov_findnextsection(fileunit, readstatus, section)

  USE rttov_const, ONLY : lensection
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : nsections, section_types
  USE parkind1, ONLY : jplm
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),        INTENT(IN)  :: fileunit
  INTEGER(KIND=jpim),        INTENT(OUT) :: readstatus
  CHARACTER(LEN=lensection), INTENT(OUT) :: section
!INTF_END

  INTEGER(KIND=jpim) :: i
  LOGICAL(KIND=jplm) :: sectionfound
  CHARACTER(LEN=80)  :: line
!- End of header --------------------------------------------------------

  section = ''
  sectionfound = .FALSE.

  readfile: DO

     READ(UNIT=fileunit, FMT='(a)', IOSTAT=readstatus) line
     IF (readstatus /= 0) EXIT

     line = ADJUSTL(line)
     IF (line(1:1) == '!' .OR. line == '') THEN
        CYCLE !skip blank/comment lines
     ELSE IF (.NOT. sectionfound) THEN
        !check for a section name
        DO i = 1, nsections
           IF (section_types(i) == line) THEN
              sectionfound = .TRUE.
              section = section_types(i)
              EXIT
           ENDIF
        ENDDO
     ELSE
        !reposition file at the start of the line and exit
        BACKSPACE(fileunit)
        EXIT readfile
     ENDIF
  ENDDO readfile

END SUBROUTINE rttov_findnextsection
