!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis

MODULE time_utils

  IMPLICIT NONE

  INTEGER, PRIVATE  :: days_in_month(12)
  DATA days_in_month/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE geth_idts(ndate, odate, idts)

    IMPLICIT NONE

    !  From 2 input mdates ('YYYY-MM-DD HH:MM:SS.ffff'),
    !  compute the time difference.

    !  on entry     -  ndate  -  the new hdate.
    !                  odate  -  the old hdate.

    !  on exit      -  idts    -  the change in time in seconds.

    CHARACTER(LEN=*), INTENT(INOUT) :: ndate, odate
    INTEGER, INTENT(OUT)   :: idts

    !  Local Variables

    !  yrnew    -  indicates the year associated with "ndate"
    !  yrold    -  indicates the year associated with "odate"
    !  monew    -  indicates the month associated with "ndate"
    !  moold    -  indicates the month associated with "odate"
    !  dynew    -  indicates the day associated with "ndate"
    !  dyold    -  indicates the day associated with "odate"
    !  hrnew    -  indicates the hour associated with "ndate"
    !  hrold    -  indicates the hour associated with "odate"
    !  minew    -  indicates the minute associated with "ndate"
    !  miold    -  indicates the minute associated with "odate"
    !  scnew    -  indicates the second associated with "ndate"
    !  scold    -  indicates the second associated with "odate"
    !  i        -  loop counter
    !  mday     -  a list assigning the number of days in each month

    CHARACTER(LEN=24) :: tdate
    INTEGER :: olen, nlen
    INTEGER :: yrnew, monew, dynew, hrnew, minew, scnew
    INTEGER :: yrold, moold, dyold, hrold, miold, scold
    INTEGER :: mday(12), i, newdys, olddys
    LOGICAL :: npass, opass
    INTEGER :: isign

    IF (odate .GT. ndate) THEN
      isign = -1
      tdate = ndate
      ndate = odate
      odate = tdate
    ELSE
      isign = 1
    END IF

    !  Assign the number of days in a months

    mday(1) = 31
    mday(2) = 28
    mday(3) = 31
    mday(4) = 30
    mday(5) = 31
    mday(6) = 30
    mday(7) = 31
    mday(8) = 31
    mday(9) = 30
    mday(10) = 31
    mday(11) = 30
    mday(12) = 31

    !  Break down old hdate into parts

    hrold = 0
    miold = 0
    scold = 0
    olen = LEN(odate)

    READ (odate(1:4), '(I4)') yrold
    READ (odate(6:7), '(I2)') moold
    READ (odate(9:10), '(I2)') dyold
    IF (olen .GE. 13) THEN
      READ (odate(12:13), '(I2)') hrold
      IF (olen .GE. 16) THEN
        READ (odate(15:16), '(I2)') miold
        IF (olen .GE. 19) THEN
          READ (odate(18:19), '(I2)') scold
        END IF
      END IF
    END IF

    !  Break down new hdate into parts

    hrnew = 0
    minew = 0
    scnew = 0
    nlen = LEN(ndate)

    READ (ndate(1:4), '(I4)') yrnew
    READ (ndate(6:7), '(I2)') monew
    READ (ndate(9:10), '(I2)') dynew
    IF (nlen .GE. 13) THEN
      READ (ndate(12:13), '(I2)') hrnew
      IF (nlen .GE. 16) THEN
        READ (ndate(15:16), '(I2)') minew
        IF (nlen .GE. 19) THEN
          READ (ndate(18:19), '(I2)') scnew
        END IF
      END IF
    END IF

    !  Check that the dates make sense.

    npass = .TRUE.
    opass = .TRUE.

    !  Check that the month of NDATE makes sense.

    IF ((monew .GT. 12) .OR. (monew .LT. 1)) THEN
      PRINT *, 'GETH_IDTS:  Month of NDATE = ', monew
      npass = .FALSE.
    END IF

    !  Check that the month of ODATE makes sense.

    IF ((moold .GT. 12) .OR. (moold .LT. 1)) THEN
      PRINT *, 'GETH_IDTS:  Month of ODATE = ', moold
      opass = .FALSE.
    END IF

    !  Check that the day of NDATE makes sense.

    IF (monew .NE. 2) THEN
      ! ...... For all months but February
      IF ((dynew .GT. mday(monew)) .OR. (dynew .LT. 1)) THEN
        PRINT *, 'GETH_IDTS:  Day of NDATE = ', dynew
        npass = .FALSE.
      END IF
    ELSE IF (monew .EQ. 2) THEN
      ! ...... For February
      IF ((dynew .GT. nfeb(yrnew)) .OR. (dynew .LT. 1)) THEN
        PRINT *, 'GETH_IDTS:  Day of NDATE = ', dynew
        npass = .FALSE.
      END IF
    END IF

    !  Check that the day of ODATE makes sense.

    IF (moold .NE. 2) THEN
      ! ...... For all months but February
      IF ((dyold .GT. mday(moold)) .OR. (dyold .LT. 1)) THEN
        PRINT *, 'GETH_IDTS:  Day of ODATE = ', dyold
        opass = .FALSE.
      END IF
    ELSE IF (moold .EQ. 2) THEN
      ! ....... For February
      IF ((dyold .GT. nfeb(yrold)) .OR. (dyold .LT. 1)) THEN
        PRINT *, 'GETH_IDTS:  Day of ODATE = ', dyold
        opass = .FALSE.
      END IF
    END IF

    !  Check that the hour of NDATE makes sense.

    IF ((hrnew .GT. 23) .OR. (hrnew .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  Hour of NDATE = ', hrnew
      npass = .FALSE.
    END IF

    !  Check that the hour of ODATE makes sense.

    IF ((hrold .GT. 23) .OR. (hrold .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  Hour of ODATE = ', hrold
      opass = .FALSE.
    END IF

    !  Check that the minute of NDATE makes sense.

    IF ((minew .GT. 59) .OR. (minew .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  Minute of NDATE = ', minew
      npass = .FALSE.
    END IF

    !  Check that the minute of ODATE makes sense.

    IF ((miold .GT. 59) .OR. (miold .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  Minute of ODATE = ', miold
      opass = .FALSE.
    END IF

    !  Check that the second of NDATE makes sense.

    IF ((scnew .GT. 59) .OR. (scnew .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  SECOND of NDATE = ', scnew
      npass = .FALSE.
    END IF

    !  Check that the second of ODATE makes sense.

    IF ((scold .GT. 59) .OR. (scold .LT. 0)) THEN
      PRINT *, 'GETH_IDTS:  Second of ODATE = ', scold
      opass = .FALSE.
    END IF

    IF (.NOT. npass) THEN
      PRINT *, 'Screwy NDATE: ', ndate(1:nlen)
      STOP 'ndate_2'
    END IF

    IF (.NOT. opass) THEN
      PRINT *, 'Screwy ODATE: ', odate(1:olen)
      STOP 'odate_1'
    END IF

    !  Date Checks are completed.  Continue.

    !  Compute number of days from 1 January ODATE, 00:00:00 until ndate
    !  Compute number of hours from 1 January ODATE, 00:00:00 until ndate
    !  Compute number of minutes from 1 January ODATE, 00:00:00 until ndate

    newdys = 0
    DO i = yrold, yrnew - 1
      newdys = newdys + (365 + (nfeb(i) - 28))
    END DO

    IF (monew .GT. 1) THEN
      mday(2) = nfeb(yrnew)
      DO i = 1, monew - 1
        newdys = newdys + mday(i)
      END DO
      mday(2) = 28
    END IF

    newdys = newdys + dynew - 1

    !  Compute number of hours from 1 January ODATE, 00:00:00 until odate
    !  Compute number of minutes from 1 January ODATE, 00:00:00 until odate

    olddys = 0

    IF (moold .GT. 1) THEN
      mday(2) = nfeb(yrold)
      DO i = 1, moold - 1
        olddys = olddys + mday(i)
      END DO
      mday(2) = 28
    END IF

    olddys = olddys + dyold - 1

    !  Determine the time difference in seconds

    idts = (newdys - olddys) * 86400
    idts = idts + (hrnew - hrold) * 3600
    idts = idts + (minew - miold) * 60
    idts = idts + (scnew - scold)

    IF (isign .EQ. -1) THEN
      tdate = ndate
      ndate = odate
      odate = tdate
      idts = idts * isign
    END IF

  END SUBROUTINE geth_idts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE geth_newdate(ndate, odate, idt)

    IMPLICIT NONE

    !  From old date ('YYYY-MM-DD HH:MM:SS.ffff') and
    !  delta-time, compute the new date.

    !  on entry     -  odate  -  the old hdate.
    !                  idt    -  the change in time

    !  on exit      -  ndate  -  the new hdate.

    INTEGER, INTENT(IN)           :: idt
    CHARACTER(LEN=*), INTENT(OUT) :: ndate
    CHARACTER(LEN=*), INTENT(IN)  :: odate

    !  Local Variables

    !  yrold    -  indicates the year associated with "odate"
    !  moold    -  indicates the month associated with "odate"
    !  dyold    -  indicates the day associated with "odate"
    !  hrold    -  indicates the hour associated with "odate"
    !  miold    -  indicates the minute associated with "odate"
    !  scold    -  indicates the second associated with "odate"

    !  yrnew    -  indicates the year associated with "ndate"
    !  monew    -  indicates the month associated with "ndate"
    !  dynew    -  indicates the day associated with "ndate"
    !  hrnew    -  indicates the hour associated with "ndate"
    !  minew    -  indicates the minute associated with "ndate"
    !  scnew    -  indicates the second associated with "ndate"

    !  mday     -  a list assigning the number of days in each month

    !  i        -  loop counter
    !  nday     -  the integer number of days represented by "idt"
    !  nhour    -  the integer number of hours in "idt" after taking out
    !              all the whole days
    !  nmin     -  the integer number of minutes in "idt" after taking out
    !              all the whole days and whole hours.
    !  nsec     -  the integer number of minutes in "idt" after taking out
    !              all the whole days, whole hours, and whole minutes.

    INTEGER :: nlen, olen
    INTEGER :: yrnew, monew, dynew, hrnew, minew, scnew, frnew
    INTEGER :: yrold, moold, dyold, hrold, miold, scold, frold
    INTEGER :: mday(12), nday, nhour, nmin, nsec, nfrac, i, ifrc
    LOGICAL :: opass
    CHARACTER(LEN=10) :: hfrc
    CHARACTER(LEN=1) :: sp
    ! INTEGER, EXTERNAL :: nfeb  ! in the same module now

    !  Assign the number of days in a months

    mday(1) = 31
    mday(2) = 28
    mday(3) = 31
    mday(4) = 30
    mday(5) = 31
    mday(6) = 30
    mday(7) = 31
    mday(8) = 31
    mday(9) = 30
    mday(10) = 31
    mday(11) = 30
    mday(12) = 31

    !  Break down old hdate into parts

    hrold = 0
    miold = 0
    scold = 0
    frold = 0
    olen = LEN(odate)
    IF (olen .GE. 11) THEN
      sp = odate(11:11)
    ELSE
      sp = ' '
    END IF

    !  Use internal READ statements to convert the CHARACTER string
    !  date into INTEGER components.

    READ (odate(1:4), '(I4)') yrold
    READ (odate(6:7), '(I2)') moold
    READ (odate(9:10), '(I2)') dyold
    IF (olen .GE. 13) THEN
      READ (odate(12:13), '(I2)') hrold
      IF (olen .GE. 16) THEN
        READ (odate(15:16), '(I2)') miold
        IF (olen .GE. 19) THEN
          READ (odate(18:19), '(I2)') scold
          IF (olen .GT. 20) THEN
            READ (odate(21:olen), '(I2)') frold
          END IF
        END IF
      END IF
    END IF

    !  Set the number of days in February for that year.

    mday(2) = nfeb(yrold)

    !  Check that ODATE makes sense.

    opass = .TRUE.

    !  Check that the month of ODATE makes sense.

    IF ((moold .GT. 12) .OR. (moold .LT. 1)) THEN
      WRITE (*, *) 'GETH_NEWDATE:  Month of ODATE = ', moold
      opass = .FALSE.
    END IF

    !  Check that the day of ODATE makes sense.

    IF ((dyold .GT. mday(moold)) .OR. (dyold .LT. 1)) THEN
      WRITE (*, *) 'GETH_NEWDATE:  Day of ODATE = ', dyold
      opass = .FALSE.
    END IF

    !  Check that the hour of ODATE makes sense.

    IF ((hrold .GT. 23) .OR. (hrold .LT. 0)) THEN
      WRITE (*, *) 'GETH_NEWDATE:  Hour of ODATE = ', hrold
      opass = .FALSE.
    END IF

    !  Check that the minute of ODATE makes sense.

    IF ((miold .GT. 59) .OR. (miold .LT. 0)) THEN
      WRITE (*, *) 'GETH_NEWDATE:  Minute of ODATE = ', miold
      opass = .FALSE.
    END IF

    !  Check that the second of ODATE makes sense.

    IF ((scold .GT. 59) .OR. (scold .LT. 0)) THEN
      WRITE (*, *) 'GETH_NEWDATE:  Second of ODATE = ', scold
      opass = .FALSE.
    END IF

    !  Check that the fractional part  of ODATE makes sense.

    IF (.NOT. opass) THEN
      WRITE (*, *) 'GETH_NEWDATE: Crazy ODATE: ', odate(1:olen), olen
      STOP 'odate_3'
    END IF

    !  Date Checks are completed.  Continue.

    !  Compute the number of days, hours, minutes, and seconds in idt

    IF (olen .GT. 20) THEN !idt should be in fractions of seconds
      ifrc = olen - 20
      ifrc = 10**ifrc
      nday = ABS(idt) / (86400 * ifrc)
      nhour = MOD(ABS(idt), 86400 * ifrc) / (3600 * ifrc)
      nmin = MOD(ABS(idt), 3600 * ifrc) / (60 * ifrc)
      nsec = MOD(ABS(idt), 60 * ifrc) / (ifrc)
      nfrac = MOD(ABS(idt), ifrc)
    ELSE IF (olen .EQ. 19) THEN  !idt should be in seconds
      ifrc = 1
      nday = ABS(idt) / 86400 ! Integer number of days in delta-time
      nhour = MOD(ABS(idt), 86400) / 3600
      nmin = MOD(ABS(idt), 3600) / 60
      nsec = MOD(ABS(idt), 60)
      nfrac = 0
    ELSE IF (olen .EQ. 16) THEN !idt should be in minutes
      ifrc = 1
      nday = ABS(idt) / 1440 ! Integer number of days in delta-time
      nhour = MOD(ABS(idt), 1440) / 60
      nmin = MOD(ABS(idt), 60)
      nsec = 0
      nfrac = 0
    ELSE IF (olen .EQ. 13) THEN !idt should be in hours
      ifrc = 1
      nday = ABS(idt) / 24 ! Integer number of days in delta-time
      nhour = MOD(ABS(idt), 24)
      nmin = 0
      nsec = 0
      nfrac = 0
    ELSE IF (olen .EQ. 10) THEN !idt should be in days
      ifrc = 1
      nday = ABS(idt) / 24 ! Integer number of days in delta-time
      nhour = 0
      nmin = 0
      nsec = 0
      nfrac = 0
    ELSE
      WRITE (*, '(''GETH_NEWDATE: Strange length for ODATE: '', i3)') &
        olen
      WRITE (*, *) odate(1:olen)
      STOP 'odate_4'
    END IF

    IF (idt .GE. 0) THEN

      frnew = frold + nfrac
      IF (frnew .GE. ifrc) THEN
        frnew = frnew - ifrc
        nsec = nsec + 1
      END IF

      scnew = scold + nsec
      IF (scnew .GE. 60) THEN
        scnew = scnew - 60
        nmin = nmin + 1
      END IF

      minew = miold + nmin
      IF (minew .GE. 60) THEN
        minew = minew - 60
        nhour = nhour + 1
      END IF

      hrnew = hrold + nhour
      IF (hrnew .GE. 24) THEN
        hrnew = hrnew - 24
        nday = nday + 1
      END IF

      dynew = dyold
      monew = moold
      yrnew = yrold
      DO i = 1, nday
        dynew = dynew + 1
        IF (dynew .GT. mday(monew)) THEN
          dynew = dynew - mday(monew)
          monew = monew + 1
          IF (monew .GT. 12) THEN
            monew = 1
            yrnew = yrnew + 1
            ! If the year changes, recompute the number of days in February
            mday(2) = nfeb(yrnew)
          END IF
        END IF
      END DO

    ELSE IF (idt .LT. 0) THEN

      frnew = frold - nfrac
      IF (frnew .LT. 0) THEN
        frnew = frnew + ifrc
        nsec = nsec - 1
      END IF

      scnew = scold - nsec
      IF (scnew .LT. 00) THEN
        scnew = scnew + 60
        nmin = nmin + 1
      END IF

      minew = miold - nmin
      IF (minew .LT. 00) THEN
        minew = minew + 60
        nhour = nhour + 1
      END IF

      hrnew = hrold - nhour
      IF (hrnew .LT. 00) THEN
        hrnew = hrnew + 24
        nday = nday + 1
      END IF

      dynew = dyold
      monew = moold
      yrnew = yrold
      DO i = 1, nday
        dynew = dynew - 1
        IF (dynew .EQ. 0) THEN
          monew = monew - 1
          IF (monew .EQ. 0) THEN
            monew = 12
            yrnew = yrnew - 1
            ! If the year changes, recompute the number of days in February
            mday(2) = nfeb(yrnew)
          END IF
          dynew = mday(monew)
        END IF
      END DO
    END IF

    !  Now construct the new mdate

    nlen = LEN(ndate)

    IF (nlen .GT. 20) THEN
      WRITE (ndate(1:19), 19) yrnew, monew, dynew, hrnew, minew, scnew
      WRITE (hfrc, '(I10)') frnew + 1000000000
      ndate = ndate(1:19)//'.'//hfrc(31 - nlen:10)

    ELSE IF (nlen .EQ. 19 .OR. nlen .EQ. 20) THEN
      WRITE (ndate(1:19), 19) yrnew, monew, dynew, hrnew, minew, scnew
19    FORMAT(I4, '-', I2.2, '-', I2.2, '_', I2.2, ':', I2.2, ':', I2.2)
      IF (nlen .EQ. 20) ndate = ndate(1:19)//'.'

    ELSE IF (nlen .EQ. 16) THEN
      WRITE (ndate, 16) yrnew, monew, dynew, hrnew, minew
16    FORMAT(I4, '-', I2.2, '-', I2.2, '_', I2.2, ':', I2.2)

    ELSE IF (nlen .EQ. 13) THEN
      WRITE (ndate, 13) yrnew, monew, dynew, hrnew
13    FORMAT(I4, '-', I2.2, '-', I2.2, '_', I2.2)

    ELSE IF (nlen .EQ. 10) THEN
      WRITE (ndate, 10) yrnew, monew, dynew
10    FORMAT(I4, '-', I2.2, '-', I2.2)

    END IF

    IF (olen .GE. 11) ndate(11:11) = sp

  END SUBROUTINE geth_newdate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION nfeb(year) RESULT(num_days)

    ! Compute the number of days in February for the given year

    IMPLICIT NONE

    INTEGER :: year
    INTEGER :: num_days

    num_days = 28 ! By default, February has 28 days ...
    IF (MOD(year, 4) .EQ. 0) THEN
      num_days = 29  ! But every four years, it has 29 days ...
      IF (MOD(year, 100) .EQ. 0) THEN
        num_days = 28  ! Except every 100 years, when it has 28 days ...
        IF (MOD(year, 400) .EQ. 0) THEN
          num_days = 29  ! Except every 400 years, when it has 29 days.
        END IF
      END IF
    END IF

  END FUNCTION nfeb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE split_date_char(date, century_year, month, day, hour, minute, second)

    IMPLICIT NONE

    !  Input data.

    CHARACTER(LEN=19), INTENT(IN) :: date

    !  Output data.

    INTEGER, INTENT(OUT) :: century_year, month, day, hour, minute, second

    READ (date, FMT='(    I4.4)') century_year
    READ (date, FMT='( 5X,I2.2)') month
    READ (date, FMT='( 8X,I2.2)') day
    READ (date, FMT='(11X,I2.2)') hour
    READ (date, FMT='(14X,I2.2)') minute
    READ (date, FMT='(17X,I2.2)') second

  END SUBROUTINE split_date_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION compute_day_of_year(century_year, month, day) RESULT(day_of_year)

    IMPLICIT NONE
    INTEGER                     :: century_year, month, day
    INTEGER                     :: day_of_year
    INTEGER                     :: m

    days_in_month(2) = nfeb(century_year)
    IF (month .EQ. 1) THEN
      day_of_year = day
    ELSE
      day_of_year = 0
      DO m = 1, month - 1
        day_of_year = day_of_year + days_in_month(m)
      END DO
      day_of_year = day_of_year + day
    END IF

  END FUNCTION compute_day_of_year
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION hms_to_wrf_time(hour, minute, second) RESULT(seconds_utc)
    IMPLICIT NONE
    INTEGER                     :: hour, minute, second
    REAL                        :: seconds_utc
    seconds_utc = FLOAT(hour * 3600 + minute * 60 + second)
  END FUNCTION hms_to_wrf_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrf_date_to_ymd(wrf_date, century_year, month, day)

    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: wrf_date
    INTEGER, INTENT(OUT)               :: century_year
    INTEGER, INTENT(OUT)               :: month
    INTEGER, INTENT(OUT)               :: day
    INTEGER :: day_of_year, m, total_days

    century_year = wrf_date / 1000
    day_of_year = MOD(wrf_date, 1000)
    days_in_month(2) = nfeb(century_year)
    total_days = 0
    month_loop: DO m = 1, 12
      total_days = total_days + days_in_month(m)
      IF (total_days .LT. day_of_year) THEN
        CYCLE month_loop
      ELSE
        month = m
        day = days_in_month(m) - (total_days - day_of_year)
        EXIT month_loop
      END IF
    END DO month_loop
  END SUBROUTINE wrf_date_to_ymd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrf_time_to_hms(wrf_time, hour, minute, second)
    IMPLICIT NONE
    REAL, INTENT(IN)                    :: wrf_time
    INTEGER, INTENT(OUT)                :: hour, minute, second
    INTEGER                             :: leftover_seconds
    hour = FLOOR(wrf_time) / 3600
    leftover_seconds = MOD(FLOOR(wrf_time), 3600)
    minute = leftover_seconds / 60
    second = MOD(leftover_seconds, 60)
  END SUBROUTINE wrf_time_to_hms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrf_to_mm5_date(wrf_date, wrf_time, mm5_date)

    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: wrf_date
    REAL, INTENT(IN)                :: wrf_time
    CHARACTER(LEN=19), INTENT(OUT)     :: mm5_date
    INTEGER :: century_year, month, day, hour, minute, second

    CALL wrf_date_to_ymd(wrf_date, century_year, month, day)
    CALL wrf_time_to_hms(wrf_time, hour, minute, second)
    WRITE (mm5_date, 91) century_year, month, day, hour, minute, second
91  FORMAT(I4.4, '-', I2.2, '-', I2.2, '_', I2.2, ':', I2.2, ':', I2.2)

  END SUBROUTINE wrf_to_mm5_date

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mm5_to_wrf_date(mm5_date, wrf_date, wrf_time)

    IMPLICIT NONE
    CHARACTER(LEN=19), INTENT(IN)      :: mm5_date
    INTEGER, INTENT(OUT)                :: wrf_date
    REAL, INTENT(OUT)                   :: wrf_time
    INTEGER :: century_year, month, day, hour, minute, second
    INTEGER :: day_of_year

    CALL split_date_char(mm5_date, century_year, month, day, hour, &
                         minute, second)
    wrf_time = hms_to_wrf_time(hour, minute, second)
    day_of_year = compute_day_of_year(century_year, month, day)
    wrf_date = century_year * 1000 + day_of_year

  END SUBROUTINE mm5_to_wrf_date

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mm5_to_lapstime(mm5_date, lapstime)

    IMPLICIT NONE
    CHARACTER(LEN=24), INTENT(IN)   :: mm5_date
    INTEGER, INTENT(OUT)            :: lapstime
    INTEGER                         :: year, month, day, hour, minute, second
    INTEGER                         :: jday
    INTEGER                         :: total_years, leap_days, total_days
    INTEGER                         :: total_hours

    CALL split_date_char(mm5_date(1:19), year, month, day, hour, minute, second)
    jday = compute_day_of_year(year, month, day)

    total_years = year - 1960
    leap_days = NINT((FLOAT(total_years) + 1.) / 4.)
    total_days = total_years * 365 + leap_days + jday - 1
    total_hours = total_days * 24 + hour
    lapstime = total_hours * 3600 + minute * 60
    RETURN
  END SUBROUTINE mm5_to_lapstime
END MODULE time_utils

