! Description:
!> @file
!!   Delete comments (starting '!') from a character string.
!
!> @brief
!!   Delete comments (starting '!') from a character string.
!!
!! @details
!!   Comments begin '!'. Every character from '!' (if present) to the end of
!!   the string is replaced by a space. Finally string is adjusted left.
!!
!! @param[in,out]    string       string to modify
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
SUBROUTINE rttov_deletecomment(string)
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON

  IMPLICIT NONE
  CHARACTER(LEN=*) , INTENT(INOUT) :: string
!INTF_END

  CHARACTER(LEN=1)   :: comment = '!' ! character for starting comment
  INTEGER(KIND=jpim) :: pos_mark      ! position of character '!' in current string
  INTEGER(KIND=jpim) :: length        ! length of string
  INTEGER(KIND=jpim) :: i             ! loop indice
!- End of header

  ! find position of comment character in string
  pos_mark = SCAN(string, comment)
  length = LEN(string)

  ! if comment is present, replace comment by spaces
  IF (pos_mark > 0) THEN
    DO i = pos_mark, length
      string(i:i) = ' '
    ENDDO
  ENDIF

  ! Adjust left string
  string = ADJUSTL(string)

END SUBROUTINE rttov_deletecomment
