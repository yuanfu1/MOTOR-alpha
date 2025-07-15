! Description:
!> @file
!!   Defines functions to handle program arguments
!
!> @brief
!!   Defines functions to handle program arguments
!!
!! @details
!!   These functions and subroutines are used internally to
!!   process arguments passed to RTTOV executables on the
!!   commandline.
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
MODULE rttov_getoptions

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_iargc, rttov_getarg, &
    rttov_basename, rttov_countwords, rttov_getenv,     &
    rttov_isalpha, rttov_isdigit, rttov_exit

  IMPLICIT NONE

  INTERFACE GetOption
    MODULE PROCEDURE GetOptionS, GetOptionSL, &
                     GetOptionI, GetOptionIL, &
                     GetOptionR, GetOptionRL, &
                     GetOptionB

  END INTERFACE

  !! @todo : list with fixed size

  PRIVATE
  PUBLIC :: GetOption, InitOptions, InitOptionsStr, CheckOptions, AddGroup

  INTEGER, PARAMETER :: argsizemax = 256

  CHARACTER(LEN=argsizemax), POINTER :: myargs(:) => NULL()
  LOGICAL(KIND=jplm), POINTER :: check_args(:) => NULL()
  LOGICAL(KIND=jplm) :: lhelp  = .FALSE., lshell = .FALSE.

  CHARACTER(LEN=1056) :: message_opt = ""
  INTEGER(KIND=jpim) :: sepind(0:1024)

  TYPE rttov_opt
    CHARACTER(LEN=32)   :: key, type
    CHARACTER(LEN=1024) :: use
    LOGICAL(KIND=jplm)  :: group = .FALSE.
  END TYPE

  INTEGER(KIND=jpim) :: nopt_seen
  TYPE(rttov_opt), POINTER :: opt_seen(:) => NULL()

CONTAINS

SUBROUTINE AddGroup(use)
  CHARACTER(LEN=*), INTENT(IN) :: use

  CALL init_opt_seen()
  nopt_seen = nopt_seen + 1
  CALL grow_opt_seen()

  opt_seen(nopt_seen)%group = .TRUE.
  opt_seen(nopt_seen)%use = use

END SUBROUTINE

CHARACTER(LEN=argsizemax) FUNCTION get_env_opt(key)
  CHARACTER(LEN=*), INTENT(IN) :: key

  CHARACTER(LEN=argsizemax) :: key_env, val_env
  INTEGER(KIND=jpim) :: i, n
  CHARACTER :: c

  key_env = key(3:)

  n = LEN(TRIM(key_env))
  DO i = 1, n
    c = key_env(i:i)
    IF ((.NOT.rttov_isalpha(c)) .AND. &
        (.NOT.rttov_isdigit(c)) .AND. &
        (c .NE. '_')) THEN
      key_env(i:i) = '_'
    ENDIF
  ENDDO

  val_env = ""
  CALL rttov_getenv('rttov_opt_'//TRIM(key_env), val_env)
  !PRINT *, " key = ", TRIM(key_env), " val = ", TRIM(val_env)
  get_env_opt = val_env

END FUNCTION

SUBROUTINE mygetarg(i, s)
  INTEGER(KIND=jpim), INTENT(IN)  :: i
  CHARACTER(LEN=*),   INTENT(OUT) :: s

  IF (ASSOCIATED(myargs)) THEN
    IF (i .LE. UBOUND(myargs, 1)) THEN
      s = myargs(i)
    ELSE
      s = ""
    ENDIF
  ELSE
    CALL rttov_getarg(i, s)
  ENDIF
END SUBROUTINE

INTEGER FUNCTION myiargc()
  INTEGER :: n
  IF (ASSOCIATED(myargs)) THEN
    n = UBOUND(myargs, 1)
  ELSE
    n = rttov_iargc()
  ENDIF
  myiargc = n
END FUNCTION

SUBROUTINE addopt_shell(key, type, mnd, use)
  CHARACTER(LEN=*),   INTENT(IN) :: key, type, use
  LOGICAL(KIND=jplm), INTENT(IN) :: mnd
  OPTIONAL :: use, mnd

  CHARACTER(LEN=argsizemax) :: str
  INTEGER :: nn, n, n1, i1, i2, k
  CHARACTER(LEN=argsizemax), POINTER :: myargs1(:)

  myargs1 => NULL()

  IF (PRESENT(use)) WRITE(*, '("> ",A)') TRIM(use)
  IF (PRESENT(mnd)) THEN
    IF (mnd) WRITE(*, *) "[Mandatory]"
  ENDIF
  WRITE(*, *) "* Option: [", type, "]", " ", TRIM(key)
  READ(*, '(A)') str

! PRINT *, "str = ",TRIM(str)
  IF (TRIM(str) .NE. "") THEN
    IF (type .EQ. 'flag') THEN
      nn = 0
    ELSE
      nn = rttov_countwords(str)
    ENDIF
    n  = UBOUND(myargs, 1)
    n1 = n + nn + 1

!
! realloc myargs
!
    ALLOCATE(myargs1(0:n1))
    myargs1(0:n) = myargs(0:n)
    DEALLOCATE(myargs)
    myargs => myargs1
    myargs(n+1) = key

!
! parse argument list
!
    IF (type .NE. 'flag') THEN
      k = 1
      i1 = 1
      loop_i1 : DO
        DO
          IF (i1 .GT. LEN(str)) EXIT loop_i1
          IF (str(i1:i1) .NE. ' ') EXIT
          i1 = i1+1
        ENDDO
        i2 = i1+1
        DO
          IF (i2 .GT. LEN(str)) EXIT
          IF (str(i2:i2) .EQ. ' ') EXIT
          i2 = i2+1
        ENDDO
!PRINT *, i1, i2
        myargs(n+1+k) = str(i1:i2-1)
!PRINT *, k, TRIM(myargs(n+1+k))
        k = k+1
        i1 = i2+1
      ENDDO loop_i1
    ENDIF
  ENDIF

END SUBROUTINE

SUBROUTINE init_opt_seen()

  IF (.NOT. ASSOCIATED(opt_seen)) THEN
    nopt_seen = 0
    ALLOCATE(opt_seen(32))
  ENDIF

END SUBROUTINE

SUBROUTINE grow_opt_seen()
  INTEGER(KIND=jpim) :: n
  TYPE(rttov_opt), POINTER :: opt_seen1(:)

  n = SIZE(opt_seen)
  IF (nopt_seen .GE. n) THEN ! realloc data
    opt_seen1 => opt_seen
    ALLOCATE(opt_seen(2 * n))
    opt_seen(1:nopt_seen) = opt_seen1(1:nopt_seen)
    DEALLOCATE(opt_seen1)
  ENDIF

END SUBROUTINE

SUBROUTINE addopt(key, type, use)
  CHARACTER(LEN=*), INTENT(IN) :: key, type, use
  OPTIONAL :: use

  CALL init_opt_seen()

  nopt_seen = nopt_seen + 1

  CALL grow_opt_seen()

  opt_seen(nopt_seen)%key  = key
  opt_seen(nopt_seen)%type = type

  IF (PRESENT(use)) THEN
    opt_seen(nopt_seen)%use  = use
  ELSE
    opt_seen(nopt_seen)%use  = ''
  ENDIF

END SUBROUTINE

SUBROUTINE InitOptions(Message)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: Message
  INTEGER(KIND=jpim) :: n, i
  CHARACTER(LEN=32) :: str
  n = rttov_iargc()

  ALLOCATE(myargs(0:n))
  DO i = 0, n
    CALL rttov_getarg(i, myargs(i))
  ENDDO

  IF (PRESENT(Message)) THEN
  message_opt = Message
  ELSE
  message_opt = ""
  ENDIF

  IF (n .EQ. 0) THEN
    lhelp = .TRUE.
    RETURN
  ENDIF

  IF (n .EQ. 1) THEN
    CALL mygetarg(1_jpim, str)
    IF (TRIM(str) .EQ. '--help') THEN
      lhelp = .TRUE.
      RETURN
    ELSE IF (TRIM(str) .EQ. '--shell') THEN
      lshell = .TRUE.
      RETURN
    ENDIF
  ENDIF

  lhelp = .FALSE.
  ALLOCATE(check_args(n))
  check_args = .FALSE.

END SUBROUTINE

SUBROUTINE InitOptionsStr(fromstr, Message)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: Message
  CHARACTER(LEN=*), INTENT(IN) :: fromstr

  INTEGER(KIND=jpim) :: n, i
  CHARACTER(LEN=32) :: str
  CHARACTER(LEN=LEN(fromstr)) :: nstr

  sepind(0)=0

  nstr = ADJUSTL(fromstr)
  CALL nseparators(' ', TRIM(nstr), n, sepind(1:1000))

  ALLOCATE(myargs(0:n))
  myargs(0) = "fromstring"
  DO i = 1, n
    myargs(i)= ADJUSTL(nstr(sepind(i-1_jpim)+1:sepind(i)))
  ENDDO

  IF (PRESENT(Message)) THEN
   message_opt = Message
  ELSE
   message_opt = ""
  ENDIF

  IF (n .EQ. 1) THEN
    CALL mygetarg(1_jpim, str)
    IF (TRIM(str) .EQ. '--help') THEN
      lhelp = .TRUE.
      RETURN
    ELSE IF (TRIM(str) .EQ. '--shell') THEN
      lshell = .TRUE.
      RETURN
    ENDIF
  ENDIF

  lhelp = .FALSE.
  ALLOCATE(check_args(n))
  check_args = .FALSE.
END SUBROUTINE

SUBROUTINE CheckOptions()
  INTEGER(KIND=jpim) :: i, n, is, ns, ks
  CHARACTER(LEN=argsizemax) :: opt, prog
  LOGICAL(KIND=jplm) :: pb
  CHARACTER(LEN=10) :: fmt

  INTEGER(KIND=jpim), PARAMETER :: lgline = 130
  INTEGER(KIND=jpim) :: halflgline = lgline/2
  CHARACTER(LEN=lgline+14) :: buf
  CHARACTER(LEN=3) :: buffmt

  CHARACTER(LEN=lgline*20) line_opt
  INTEGER(KIND=jpim) :: sepind(lgline*20)
  INTEGER(KIND=jpim) :: nbsep
  INTEGER(KIND=jpim) :: l, i1, i2

  CALL mygetarg(0_jpim, prog)

  IF (lhelp) THEN
    PRINT *, "Program: ", TRIM(rttov_basename(prog))
    IF (TRIM(message_opt) .NE. "") THEN
      CALL nseparators(';', message_opt, nbsep, sepind)
      i1 = 0_jpim
      i2 = 0_jpim
      DO l = 1, nbsep
        i1 = i2 + 1_jpim
        i2 = sepind(l)
        IF (message_opt(i2:i2) .EQ. ";") THEN
          line_opt=message_opt(i1:i2-1)
        ELSE
          line_opt=message_opt(i1:i2)
        ENDIF
        ns = LEN(line_opt)

        DO is = 1, ns / lgline
          ks = LEN(TRIM(line_opt(1+(is-1)*lgline:is*lgline)))
          IF (ks .GT. 0) THEN
            IF (is .EQ. 1) THEN
              WRITE(*, '("    ")', advance = 'no')
            ELSE
              WRITE(*, '("  > ")', advance = 'no')
            ENDIF
            WRITE(fmt, '("(A",I3.3,")")') ks
            WRITE(*, fmt) TRIM(line_opt(1+(is-1)*lgline:is*lgline))
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    WRITE(buffmt, '(i3)') LEN(buf)

    DO i = 1, nopt_seen

      IF (opt_seen(i)%group) THEN
        WRITE(*, *)
        IF (TRIM(opt_seen(i)%use) .NE. "") &
          WRITE(*, *) '* '//TRIM(opt_seen(i)%use)
        cycle
      ENDIF

      buf = ""

      WRITE(buf, '(A32,"   ",A15)') &
          TRIM(opt_seen(i)%key), &
          TRIM(opt_seen(i)%type)

      IF (TRIM(opt_seen(i)%use) .NE. '') THEN
        CALL nseparators(';',opt_seen(i)%use, nbsep, sepind)
        i1 = 0_jpim
        i2 = 0_jpim
        DO l = 1, nbsep
          i1 = i2 + 1_jpim
          i2 = sepind(l)

          IF (opt_seen(i)%use(i2:i2) .EQ. ";") THEN
            line_opt=opt_seen(i)%use(i1:i2-1)
          ELSE
            line_opt=opt_seen(i)%use(i1:i2)
          ENDIF
          ns = LEN(line_opt)

          DO is = 1, ns / halflgline
            ks = LEN(TRIM(line_opt(1+(is-1)*halflgline:is*halflgline)))
            IF (ks .GT. 0) THEN
              IF (is .EQ. 1 .AND. l .EQ. 1) THEN
                buf = TRIM(buf)//" :   "//TRIM(line_opt(1+(is-1)*halflgline:is*halflgline))
              ELSEIF (is.EQ.1) THEN
!                     000000000011111111112222222222333333333344444444445555555555
!                     012345678901234567890123456789012345678901234567890123456789
                buf = "                                                       "&
                      //TRIM(line_opt(1+(is-1)*halflgline:is*halflgline))
              ELSE
                buf = "                                                     > "&
                      //TRIM(line_opt(1+(is-1)*halflgline:is*halflgline))
              ENDIF
              WRITE(*, '(A'//TRIM(buffmt)//')') buf
            ENDIF
          ENDDO
        ENDDO
      ELSE
        WRITE(*, '(A'//TRIM(buffmt)//')') buf
        WRITE(*, *)
      ENDIF

    ENDDO
    stop
  ELSE IF (ASSOCIATED(check_args)) THEN
    n = SIZE(check_args)
    pb = .FALSE.
    DO i = 1, n
      IF (.NOT. check_args(i)) THEN
        CALL mygetarg(i, opt)
        IF (opt(1:2) .EQ. '--') THEN
          PRINT *, 'Invalid option: ', TRIM(opt)
          pb = .TRUE.
          check_args(i) = .TRUE.
        ENDIF
      ENDIF
    ENDDO

    DO i = 1, n
      IF (.NOT. check_args(i)) THEN
        CALL mygetarg(i, opt)
        PRINT *, 'Garbage in options:`', TRIM(opt), "'"
        pb = .TRUE.
        EXIT
      ENDIF
    ENDDO

    IF (pb) CALL rttov_exit(1_jpim)

    DEALLOCATE(check_args)
  ELSE IF (lshell) THEN
    OPEN(77, file = TRIM(prog)//'.sh', form = 'formatted')
    WRITE(77, '("#!/bin/sh")')
    WRITE(77, *)
    WRITE(77, '(A)', advance = 'no') TRIM(prog)
    n = UBOUND(myargs, 1)
    DO i = 1, n
      IF (myargs(i) .EQ. '--shell') cycle
      IF (myargs(i)(1:2) .EQ. '--') THEN
        WRITE(77, '(" \")')
        WRITE(77, '("    ")', advance = 'no')
      ENDIF
      WRITE(77, '(" ",A)', advance = 'no') TRIM(myargs(i))
    ENDDO
    WRITE(77, *)
    CLOSE(77)
  ENDIF

  IF (ASSOCIATED(opt_seen)) DEALLOCATE(opt_seen)
  IF (ASSOCIATED(myargs)) DEALLOCATE(myargs)
END SUBROUTINE

SUBROUTINE Check_mnd(key, mnd, use)
  CHARACTER(LEN=*),             INTENT(IN) :: key
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN) :: use
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(IN) :: mnd

  CHARACTER(LEN=argsizemax) :: prog

  IF (PRESENT(mnd)) THEN
    IF (mnd) THEN
      CALL mygetarg(0_jpim, prog)
      WRITE(*, '("PROGRAM: ",(a))') TRIM(prog)
      WRITE(*, '("ERROR:   Option `",(a),"'' is mandatory")') TRIM(key)
      IF (PRESENT(use)) WRITE(*, '("         ",(a)," : ",(a))') TRIM(key), TRIM(use)
      CALL rttov_exit(1_jpim)
    ENDIF
  ENDIF
END SUBROUTINE

SUBROUTINE FindArgIndex(key, i, n)
  CHARACTER(LEN=*),   INTENT(IN)  :: key
  INTEGER(KIND=jpim), INTENT(OUT) :: i, n
  CHARACTER(LEN=argsizemax) :: arg

  n = myiargc()
  DO i = 1, n
    CALL mygetarg(i, arg)
    IF (TRIM(ADJUSTL(arg)) .EQ. TRIM(ADJUSTL(key))) RETURN
  ENDDO
  i = -1_jpim
END SUBROUTINE

SUBROUTINE FindNextArgIndex(i, j)
  INTEGER(KIND=jpim), INTENT(IN)  :: i
  INTEGER(KIND=jpim), INTENT(OUT) :: j

  CHARACTER(LEN=argsizemax) :: arg
  INTEGER(KIND=jpim) :: n

  n = myiargc()
  DO j = i+1, n
    CALL mygetarg(j, arg)
    IF (arg(1:2) .EQ. '--') EXIT
  ENDDO
END SUBROUTINE

SUBROUTINE GetOptionS(key, val, mnd, use, exists)
  CHARACTER(LEN=*),   INTENT(IN)              :: key
  CHARACTER(LEN=*),   INTENT(INOUT)           :: val
  LOGICAL(KIND=jplm), INTENT(IN),    OPTIONAL :: mnd
  CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL :: use
  LOGICAL(KIND=jplm), INTENT(INOUT), OPTIONAL :: exists

  INTEGER(KIND=jpim) :: i, n
  CHARACTER(LEN=argsizemax) :: arg
  LOGICAL(KIND=jplm) :: lshell1
  LOGICAL(KIND=jplm) :: found

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'string', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'string', mnd, use)
  ENDIF

  CALL FindArgIndex(key, i, n)

  found = (0 .LT. i) .AND. (i .LT. n)

  IF (found) THEN
    IF (ASSOCIATED(check_args)) THEN
      check_args(i)   = .TRUE.
      check_args(i+1) = .TRUE.
    ENDIF
    CALL mygetarg(i+1_jpim, val)
  ELSE
    arg = get_env_opt(key)
    found = arg .NE. ""
    IF (found) val = arg
  ENDIF

  IF (.NOT. found) &
    CALL check_mnd(key, mnd, use)

  IF (PRESENT(exists)) exists = found

  lshell = lshell1

END SUBROUTINE

SUBROUTINE GetOptionI(key, val, mnd, use, exists)
  CHARACTER(LEN=*),   INTENT(IN)              :: key
  INTEGER(KIND=jpim), INTENT(INOUT)           :: val
  LOGICAL(KIND=jplm), INTENT(IN),    OPTIONAL :: mnd
  CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL :: use
  LOGICAL(KIND=jplm), INTENT(INOUT), OPTIONAL :: exists

  CHARACTER(LEN=argsizemax) :: sval
  INTEGER :: err
  LOGICAL(KIND=jplm) :: lshell1

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'INTEGER', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'INTEGER', mnd, use)
  ENDIF

  sval = ""
  CALL GetOptionS(key, sval, mnd, use, exists)
  IF (TRIM(sval) .NE. "") THEN
    READ(sval, *, iostat = err) val
    IF (err .NE. 0) THEN
      PRINT *, "Error while parsing option "//TRIM(key)
      CALL rttov_exit(1_jpim)
    ENDIF
  ENDIF

  lshell = lshell1

END SUBROUTINE

SUBROUTINE GetOptionR(key, val, mnd, use, exists)
  CHARACTER(LEN=*),   INTENT(IN)              :: key
  REAL(KIND=jprb),    INTENT(INOUT)           :: val
  LOGICAL(KIND=jplm), INTENT(IN),    OPTIONAL :: mnd
  CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL :: use
  LOGICAL(KIND=jplm), INTENT(INOUT), OPTIONAL :: exists

  CHARACTER(LEN=argsizemax) :: sval
  INTEGER :: err
  LOGICAL(KIND=jplm) :: lshell1

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'real', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'real', mnd, use)
  ENDIF

  sval = ""
  CALL GetOptionS(key, sval, mnd, use, exists)
  IF (TRIM(sval) .NE. "") THEN
    READ(sval, *, iostat = err) val
    IF (err .NE. 0) THEN
      PRINT *, "Error while parsing option "//TRIM(key)
      CALL rttov_exit(1_jpim)
    ENDIF
  ENDIF

  lshell = lshell1

END SUBROUTINE

SUBROUTINE ReadASLFromString(val, sval)
  CHARACTER(LEN=*), INTENT(OUT) :: val(:)
  CHARACTER(LEN=*), INTENT(IN) :: sval

  INTEGER(KIND=jpim) :: i, j, k, n

  n = LEN(sval)
  i = 1
  k = 1
  do1 : DO
    DO
      IF (i .GT. n) EXIT do1
      IF (sval(i:i) .NE. ' ') EXIT
      i = i + 1
    ENDDO
    j = i
    DO
      IF (j .GT. n) EXIT
      IF (sval(j:j) .EQ. ' ') EXIT
      j = j + 1
    ENDDO
    val(k) = sval(i:j-1)
    i = j
    k = k + 1
  ENDDO do1
END SUBROUTINE

SUBROUTINE ReadSLFromString(val, sval)
  CHARACTER(LEN=*), POINTER    :: val(:)
  CHARACTER(LEN=*), INTENT(IN) :: sval

  INTEGER(KIND=jpim) :: n

  n = rttov_countwords(sval)
  ALLOCATE(val(n))

  CALL ReadASLFromString(val, sval)
END SUBROUTINE

SUBROUTINE ReadSLFromFile(val, sval)
  CHARACTER(LEN=*), POINTER    :: val(:)
  CHARACTER(LEN=*), INTENT(IN) :: sval

  INTEGER(KIND=jpim) :: k, n
  INTEGER(KIND=jpim) :: ioerr
  CHARACTER(LEN=4096) :: buffer

  OPEN(77, file = TRIM(sval), form = 'formatted', status = 'old', iostat = ioerr)
  IF (ioerr .NE. 0) THEN
    PRINT '( "Could not open ",A, " for reading")', TRIM(sval)
    CALL rttov_exit(1_jpim)
  ENDIF
  n = 0_jpim
  DO
    READ(77, '(A)', end = 500) buffer
    n = n + Rttov_CountWords(buffer)
  ENDDO

  500 CONTINUE

  REWIND(77)

  ALLOCATE(val(n))
  k = 1
  DO
    READ(77, '(A)', end = 600) buffer
    n = rttov_countwords(buffer)
    CALL ReadASLFromString(val(k:k+n-1), buffer)
    k = k + n
  ENDDO

  600 CONTINUE

  CLOSE(77)
END SUBROUTINE

SUBROUTINE GetOptionSL(key, val, mnd, use, exists)

  CHARACTER(LEN=*),             INTENT(IN)    :: key
  CHARACTER(LEN=*),             POINTER       :: val(:)
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(IN)    :: mnd
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)    :: use
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(INOUT) :: exists

  INTEGER(KIND=jpim) :: i, j, k, n
  CHARACTER(LEN=argsizemax) :: arg
  CHARACTER(LEN=argsizemax) :: sval
  LOGICAL(KIND=jplm) :: lshell1
  LOGICAL(KIND=jplm) :: found

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'string-list', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'string-list', mnd, use)
  ENDIF

  CALL FindArgIndex(key, i, n)

  found = i >= 0

  IF (found) THEN

    CALL FindNextArgIndex(i, j)

    ALLOCATE(val(j - i - 1))

    IF (ASSOCIATED(check_args)) &
      check_args(i) = .TRUE.

    DO k = i+1, j-1
      IF (ASSOCIATED(check_args)) &
        check_args(k)   = .TRUE.
      CALL mygetarg(k, arg)
      val(k-i) = arg
    ENDDO

    IF ((SIZE(val) .EQ. 1).AND.(val(1)(1:7).EQ.'file://')) THEN
      arg = val(1)(8:)
      DEALLOCATE(val)
      CALL ReadSLFromFile(val, arg)
    ENDIF

  ENDIF

  IF (.NOT. found) THEN
    sval = get_env_opt(key)
    found = sval .NE. ""
    IF (found) CALL ReadSLFromString(val, sval)
  ENDIF

  IF (.NOT. found) CALL check_mnd(key, mnd, use)

  IF (PRESENT(exists)) exists = found

  lshell = lshell1

END SUBROUTINE

SUBROUTINE GetOptionIL(key, val, mnd, use, exists)
  CHARACTER(LEN=*),             INTENT(IN)    :: key
  INTEGER(KIND=jpim),           POINTER       :: val(:)
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(IN)    :: mnd
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)    :: use
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(INOUT) :: exists

  CHARACTER(LEN=argsizemax), POINTER :: sval(:) => NULL()
  INTEGER(KIND=jpim) :: i, n
  INTEGER :: err
  LOGICAL(KIND=jplm) :: lshell1

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'integer-list', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'integer-list', mnd, use)
  ENDIF

  CALL GetOptionSL(key, sval, mnd, use, exists)

  IF (.NOT. ASSOCIATED(sval)) goto 999

  n = SIZE(sval)
  ALLOCATE(val(n))
  DO i = 1, n
    READ(sval(i), *, iostat = err) val(i)
    IF (err .NE. 0) THEN
      PRINT *, "Error while parsing option "//TRIM(key)
      CALL rttov_exit(1_jpim)
    ENDIF
  ENDDO

  DEALLOCATE(sval)

  999 CONTINUE
  lshell = lshell1
END SUBROUTINE

SUBROUTINE GetOptionRL(key, val, mnd, use, exists)
  CHARACTER(LEN=*),             INTENT(IN)    :: key
  REAL(KIND=jprb),              POINTER       :: val(:)
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(IN)    :: mnd
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)    :: use
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(INOUT) :: exists

  CHARACTER(LEN=argsizemax), POINTER :: sval(:) => NULL()
  INTEGER(KIND=jpim) :: i, n
  INTEGER :: err
  LOGICAL(KIND=jplm) :: lshell1

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, 'real-list', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, 'real-list', mnd, use)
  ENDIF

  CALL GetOptionSL(key, sval, mnd, use, exists)

  IF (.NOT. ASSOCIATED(sval)) goto 999

  n = SIZE(sval)
  ALLOCATE(val(n))
  DO i = 1, n
    READ(sval(i), *, iostat = err) val(i)
    IF (err .NE. 0) THEN
      PRINT *, "Error while parsing option "//TRIM(key)
      CALL rttov_exit(1_jpim)
    ENDIF
  ENDDO

  DEALLOCATE(sval)

  999 CONTINUE
  lshell = lshell1
END SUBROUTINE

SUBROUTINE GetOptionB(key, val, use, exists)
  CHARACTER(LEN=*),             INTENT(IN)    :: key
  LOGICAL(KIND=jplm),           INTENT(INOUT) :: val
  CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)    :: use
  LOGICAL(KIND=jplm), OPTIONAL, INTENT(INOUT) :: exists

  LOGICAL(KIND=jplm) :: lshell1
  LOGICAL(KIND=jplm) :: found
  CHARACTER(LEN=argsizemax) :: sval
  INTEGER(KIND=jpim) :: i, n

  lshell1 = lshell

  IF (lhelp) THEN
    CALL addopt(key, '(flag)', use)
    RETURN
  ELSE IF (lshell) THEN
    lshell = .FALSE.
    CALL addopt_shell(key, '(flag)', .FALSE._jplm, use)
  ENDIF

  CALL FindArgIndex(key, i, n)
  found = i > 0
  IF (found .AND. ASSOCIATED(check_args)) THEN
    check_args(i)   = .TRUE.
    val = .TRUE.
  ELSE
    sval = get_env_opt(key)
    IF (sval .NE. "") READ(sval, *) val
  ENDIF

  IF (PRESENT(exists)) exists = found

  lshell = lshell1
END SUBROUTINE

SUBROUTINE nseparators(sep, line, nbsep, sepind)
  CHARACTER(LEN=1),   INTENT(IN)  :: sep
  CHARACTER(LEN=*),   INTENT(IN)  :: line
  INTEGER(KIND=jpim), INTENT(OUT) :: nbsep
  INTEGER(KIND=jpim), INTENT(OUT) :: sepind(:)

  INTEGER(KIND=jpim) :: l, last, ns

  nbsep     = 0_jpim
  sepind(:) = 0_jpim

  IF (TRIM(line) .NE. "") THEN
    ns = LEN(line)

    l         = 1_jpim
    last      = 0_jpim

    DO WHILE (l .GT. 0_jpim)
      l = INDEX(line(last+1_jpim:ns), sep)
      IF (l .GT. 0) THEN
        nbsep = nbsep + 1_jpim
        sepind(nbsep) = l + last
        last = sepind(nbsep)
        IF (nbsep .GT. 1) THEN
          IF (last .EQ. sepind(nbsep-1)+1) THEN
            sepind(nbsep) = 0
            nbsep = nbsep -1
            !last = sepind(nbsep)
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    IF (last .NE. ns) THEN
      nbsep = nbsep + 1_jpim
      sepind(nbsep) = ns
    ENDIF
  ENDIF
END SUBROUTINE

END MODULE
