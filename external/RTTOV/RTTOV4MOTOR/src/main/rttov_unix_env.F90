! Description:
!> @file
!!   Wraps useful functions.
!
!> @brief
!!   Wraps useful functions.
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
MODULE rttov_unix_env

  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rttov_getenv, rttov_iargc, rttov_getarg, rttov_exit, &
            rttov_dirname, rttov_basename, rttov_upper_case, rttov_lower_case, &
            rttov_countlines, rttov_countwords, rttov_isalpha, rttov_isdigit, &
            rttov_date_and_time, rttov_cpu_time, rttov_banner, rttov_unlink, &
            rttov_get_file

CONTAINS

  !> Get value of specified environment variable
  !! @param[in]       key     name of environment variable
  !! @param[out]      val     value of environment variable
  SUBROUTINE rttov_getenv(key, val)
    CHARACTER(LEN=*), INTENT(IN) :: key
    CHARACTER(LEN=*), INTENT(OUT) :: val
    CALL GET_ENVIRONMENT_VARIABLE(key, val)
  END SUBROUTINE

  !> Return the number of commandline arguments
  INTEGER(KIND=jpim) FUNCTION rttov_iargc()
    rttov_iargc = COMMAND_ARGUMENT_COUNT()
  END FUNCTION

  !> Return argument with given index
  !! @param[in]       key     index of argument
  !! @param[out]      val     value of argument
  SUBROUTINE rttov_getarg(key, val)
    INTEGER(KIND=jpim), INTENT(IN)  :: key
    CHARACTER(LEN=*),   INTENT(OUT) :: val
    CALL GET_COMMAND_ARGUMENT(INT(key, SELECTED_INT_KIND(9)), val)
  END SUBROUTINE

  !> Exit program with 0/1 return status (using STOP)
  !! @param[in]       status     return status
  SUBROUTINE rttov_exit(status)
    INTEGER(KIND=jpim), INTENT(IN) :: status
    IF (status == 0) THEN
      STOP 0
    ELSE
      STOP 1
    ENDIF
  END SUBROUTINE

  !> Return name of directory containing the given path
  !! @param[in]       path     path to file (or directory)
  CHARACTER(LEN=256) FUNCTION rttov_dirname(path)
    CHARACTER(LEN=*), INTENT(IN) :: path

    INTEGER(KIND=jpim) :: i
    rttov_dirname = ""
    i = LEN(TRIM(path)) - 1
    DO
      IF (i .LE. 0) RETURN
      IF (path(i:i) .EQ. '/') EXIT
      i = i - 1
    ENDDO
    rttov_dirname = path(1:i)
  END FUNCTION

  !> Return filename excluding path
  !! @param[in]       path     full path to file
  CHARACTER(LEN=256) FUNCTION rttov_basename(path)
    CHARACTER(LEN=*), INTENT(IN) :: path

    INTEGER(KIND=jpim) :: i
    rttov_basename = ""
    i = LEN(TRIM(path)) - 1
    DO
      IF (i .LE. 0) THEN
        i = 0
        EXIT
      ENDIF
      IF (path(i:i) .EQ. '/') EXIT
      i = i - 1
    ENDDO
    rttov_basename = path(i+1:)
  END FUNCTION

  !> Convert string to lower case
  !! @param[out]      ous     output string in lower case
  !! @param[in]       ins     input string
  ELEMENTAL SUBROUTINE rttov_lower_case(ous, ins)
    ! convert a word to lower case
    CHARACTER(LEN=*), INTENT(OUT) :: ous
    CHARACTER(LEN=*), INTENT(IN)  :: ins
    INTEGER :: i, ic, nlen
    nlen = LEN(ins)
    ous = ''
    DO i = 1, nlen
      ic = ICHAR(ins(i:i))
      IF (ic >= 65 .AND. ic < 90) THEN
        ous(i:i) = CHAR(ic+32)
      ELSE
        ous(i:i) = ins(i:i)
      ENDIF
    ENDDO
  END SUBROUTINE rttov_lower_case

  !> Convert string to upper case
  !! @param[out]      ous     output string in upper case
  !! @param[in]       ins     input string
  ELEMENTAL SUBROUTINE rttov_upper_case(ous, ins)
    ! convert a word to upper case
    CHARACTER(LEN=*), INTENT(OUT) :: ous
    CHARACTER(LEN=*), INTENT(IN)  :: ins
    INTEGER :: i, ic, nlen
    nlen = LEN(ins)
    ous = ''
    DO i = 1, nlen
      ic = ICHAR(ins(i:i))
      IF (ic >= 97 .AND. ic < 122) THEN
        ous(i:i) = CHAR(ic-32)
      ELSE
        ous(i:i) = ins(i:i)
      ENDIF
    ENDDO
  END SUBROUTINE rttov_upper_case

  !> Returns true if single character c is in [a-z,A-Z]
  !! @param[in]       c     character to test
  LOGICAL(KIND=jplm) FUNCTION rttov_isalpha(c)
    CHARACTER, INTENT(IN) :: c

    rttov_isalpha = ((c.GE.'a') .AND. (c.LE.'z')) &
               .OR. ((c.GE.'A') .AND. (c.LE.'Z'))
  END FUNCTION

  !> Returns true if single character c is in [0-9]
  !! @param[in]       c     character to test
  LOGICAL(KIND=jplm) FUNCTION rttov_isdigit(c)
    CHARACTER, INTENT(IN) :: c

    rttov_isdigit = (c.GE.'0') .AND. (c.LE.'9')
  END FUNCTION

  !> Wraps date_and_time
  !! @param[out]       vl     returned date and time
  SUBROUTINE rttov_date_and_time(vl)
    INTEGER(KIND=jpim), INTENT(OUT) :: vl(8)

    INTEGER :: vlx(8)

    CALL DATE_AND_TIME(values = vlx)
    vl = vlx
  END SUBROUTINE

  !> Wraps cpu_time
  !! @param[out]       t     returned CPU time
  SUBROUTINE rttov_cpu_time(t)
    REAL, INTENT(OUT) :: t
    CALL CPU_TIME(t)
  END SUBROUTINE

  !> Counts the number of lines in a file
  !! @param[out]       nlines     number of lines
  !! @param[in]        f          path to file
  !! @param[out]       err        return status
  SUBROUTINE rttov_countlines(nlines, f, err)
    INTEGER(KIND=jpim), INTENT(OUT) :: nlines
    CHARACTER(LEN=*),   INTENT(IN)  :: f
    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=32) :: str

    nlines = 0
    OPEN(77, FILE=f, ERR=888)
    DO
      READ(77, *, ERR=888, END=777) str
      nlines = nlines + 1
    ENDDO
    777 CONTINUE
    CLOSE(77)

    RETURN
    888 CONTINUE
    err = 1_jpim
  END SUBROUTINE

  !> Counts the number of words in a character string
  !! @param[in]        s          input string
  INTEGER(KIND=jpim) FUNCTION rttov_countwords(s)
    CHARACTER(LEN=*), INTENT(IN) :: s
    INTEGER(KIND=jpim) :: n, i, l
    LOGICAL(KIND=jplm) :: in
    n = 0_jpim
    in = .FALSE.
    l = LEN(TRIM(s))
    DO i = 1, l
      IF (s(i:i) .EQ. ' ') THEN
        in = .FALSE.
      ELSEIF (.NOT. in) THEN
        n = n + 1
        in = .TRUE.
      ENDIF
    ENDDO
    rttov_countwords = n
  END FUNCTION

  !> Tries to find file fic under directory path
  !! Iterates on parent directories until found
  !! @param[in]   path   input directory
  !! @param[in]   fic    file name searched
  !! @param[out]  full   full pathname of file found
  !! @param[out]  err    return status
  RECURSIVE SUBROUTINE rttov_get_file(path, fic, full, err)
    CHARACTER(LEN=*),   INTENT(IN)  :: path
    CHARACTER(LEN=*),   INTENT(IN)  :: fic
    CHARACTER(LEN=*),   INTENT(OUT) :: full
    INTEGER(KIND=jpim), INTENT(OUT) :: err
    LOGICAL            :: existence
    CHARACTER(LEN=256) :: p
    INTEGER(KIND=jpim) :: i

    err = 0_jpim

    i = LEN(TRIM(path)) - 1
    IF (i .LE. 0) THEN
      full = ""
      err = 1_jpim
      RETURN
    ENDIF

    full = TRIM(path)//"/"//TRIM(fic)
    INQUIRE(FILE=full, EXIST=existence)
    IF (existence) THEN
      RETURN
    ELSE
      p = rttov_dirname(path)
      CALL rttov_get_file(p, fic, full, err)
    ENDIF
  END SUBROUTINE

  !> Prints out an RTTOV banner in ASCII to the specified open lun. Text to
  !! print can be specified
  !! @param[in]        lun        logical unit open for writing
  !! @param[in]        text       text to write
  !! @param[out]       err        status on exit
  SUBROUTINE rttov_banner(lun, text, err)
    INTEGER(KIND=jpim), INTENT(IN)  :: lun
    CHARACTER(LEN=*),   INTENT(IN)  :: text(:)
    INTEGER(KIND=jpim), INTENT(OUT) :: err

    INTEGER(KIND=jpim), PARAMETER :: xcol = 8_jpim
    CHARACTER(LEN=*),   PARAMETER :: cstar = '**********'
    INTEGER(KIND=jpim) :: ln, il, nl, ic, nc1, nc2, kc

    err = 0
    ln = LEN(text(1))
    nl = SIZE(text)

    DO ic = 1, xcol
      WRITE(lun, '(a)', ADVANCE='no', err=888) cstar
    ENDDO
    WRITE(lun, *, err=888)

    kc = xcol * 10 - 2 - ln

    IF (kc .LT. 0) THEN
      ln = xcol * 10 - 2
      kc = 0
    ENDIF

    nc1 = kc / 2
    nc2 = kc - nc1

    DO il = 1, nl + 2
      WRITE(lun, '(a)', ADVANCE='no', err=888) '*'
      IF ((il .EQ. 1) .OR. (il .EQ. nl + 2)) THEN
        DO ic = 1, xcol * 10 - 2
          WRITE(lun, '(" ")', ADVANCE='no', err=888)
        ENDDO
      ELSE
        DO ic = 1, nc1
          WRITE(lun, '(" ")', ADVANCE='no', err=888)
        ENDDO
        WRITE(lun, '(a)', ADVANCE='no', err=888) text(il-1)(1:ln)
        DO ic = 1, nc2
          WRITE(lun, '(" ")', ADVANCE='no', err=888)
        ENDDO
      ENDIF
      WRITE(lun, '(a)', ADVANCE='no', err=888) '*'
      WRITE(lun, *, err=888)
    ENDDO

    DO ic = 1, xcol
      WRITE(lun, '(a)', ADVANCE='no', err=888) cstar
    ENDDO
    WRITE(lun, *, err=888)

    RETURN
    888  CONTINUE
    err = 1_jpim
  END SUBROUTINE

  !> Delete (unlink) a file
  !! @param[in]        f          full path to file
  !! @param[out]       err        status on exit
  SUBROUTINE rttov_unlink(f, err)
    CHARACTER(LEN=*),   INTENT(IN)  :: f
    INTEGER(KIND=jpim), INTENT(OUT) :: err
    err = 0
    OPEN(77, FILE=f, IOSTAT=err)
    IF (err .NE. 0) RETURN
    CLOSE(77, STATUS='delete', IOSTAT=err)
  END SUBROUTINE

END MODULE
