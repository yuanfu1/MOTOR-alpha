! Description:
!> @file
!!   Compare two strings after removing all spaces and converting them to upper
!!   case.
!
!> @brief
!!   Compare two strings after removing all spaces and converting them to upper
!!   case.
!!
!! @details
!!   Returns TRUE if strings are identical.
!!
!! @param[in]    string1          first string
!! @param[in]    string2          second string
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
FUNCTION rttov_cmpuc(string1, string2)

  USE parkind1, ONLY: jplm
!INTF_OFF
  USE parkind1, ONLY: jpim
!INTF_ON
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: string1
  CHARACTER(LEN=*), INTENT(IN) :: string2
  LOGICAL(KIND=jplm) :: rttov_cmpuc
!INTF_END

  CHARACTER(LEN=LEN(string1)) :: wstr1  ! working string 1
  CHARACTER(LEN=LEN(string2)) :: wstr2  ! working string 2
  INTEGER(KIND=jpim) :: pos             ! position of space character
  INTEGER(KIND=jpim) :: cur_char        ! ASCII indice for current character
  INTEGER(KIND=jpim) :: amin            ! ASCII indice for 'a'
  INTEGER(KIND=jpim) :: amaj            ! ASCII indice for 'A'
  INTEGER(KIND=jpim) :: zmin            ! ASCII indice for 'z'
  INTEGER(KIND=jpim) :: i               ! loop indice
!- End of header --------------------------------------------------------

  amin = ICHAR('a')
  zmin = ICHAR('z')
  amaj = ICHAR('A')

  ! reduce string 1
  wstr1 = string1
  ! remove all spaces
  DO
    pos = INDEX(wstr1(1:LEN_TRIM(wstr1)), ' ')
    IF (pos == 0) EXIT
    wstr1(pos:) = wstr1(pos+1:)
  ENDDO

  ! reduce string 2
  wstr2 = string2
  ! remove all spaces
  DO
    pos = INDEX(wstr2(1:LEN_TRIM(wstr2)), ' ')
    IF (pos == 0) EXIT
    wstr2(pos:) = wstr2(pos+1:)
  ENDDO

  ! upcase string 1
  DO i = 1, LEN(wstr1)
    cur_char = ICHAR(wstr1(i:i))
    IF (cur_char >= amin .AND. cur_char <= zmin) THEN
      wstr1(i:i) = CHAR(cur_char + amaj - amin)
    ENDIF
  ENDDO

  ! upcase string 2
  DO i = 1, LEN(wstr2)
    cur_char = ICHAR(wstr2(i:i))
    IF (cur_char >= amin .AND. cur_char <= zmin) THEN
      wstr2(i:i) = CHAR(cur_char + amaj - amin)
    ENDIF
  ENDDO

  ! compare the 2 working strings
  rttov_cmpuc = wstr2 .EQ. wstr1

END FUNCTION rttov_cmpuc
