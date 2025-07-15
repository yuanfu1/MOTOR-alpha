!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Project datetime-fortran on https://github.com/wavebitscientific/datetime-fortran
! AUTOHR(S)         : ZAlexander Gavrikov
!                     Bjoern Hendrik Fock
!                     Izaak Beekman
!                     Marco Galli
!                     Mark Carter
!                     Michael Hirsch
!                     Milan Curcic
!                     Oleksandr Sasha Huziy
!                     Paul Leopardi
!                     Stefano Zaghi
!                     Tom Canich
!                     Wadud Miah
! VERSION           : V 0.0
! LICENSE           : MIT license
! HISTORY           :
!!--------------------------------------------------------------------------------------------------

MODULE datetime_m

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64, REAL32, REAL64, &
                                                                              stderr => ERROR_UNIT
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: datetime, timedelta, clock
  PUBLIC :: date2num
  PUBLIC :: datetimeRange
  PUBLIC :: daysInMonth
  PUBLIC :: daysInYear
  PUBLIC :: isLeapYear
  PUBLIC :: num2date
  PUBLIC :: strptime
  PUBLIC :: tm2date
  PUBLIC :: tm_struct
  PUBLIC :: c_strftime
  PUBLIC :: c_strptime

  REAL(REAL64), PARAMETER :: zero = 0._REAL64, one = 1._REAL64

  ! Constant multipliers to transform a number of some time unit to another
  REAL(REAL64), PARAMETER :: d2h = 24._REAL64 ! day -> hour
  REAL(REAL64), PARAMETER :: h2d = one / d2h ! hour -> day
  REAL(REAL64), PARAMETER :: d2m = d2h * 60._REAL64 ! day -> minute
  REAL(REAL64), PARAMETER :: m2d = one / d2m ! minute -> day
  REAL(REAL64), PARAMETER :: m2h = one / 60 ! minute -> hour
  REAL(REAL64), PARAMETER :: s2d = m2d / 60 ! second -> day
  REAL(REAL64), PARAMETER :: d2s = 86400._REAL64 ! day -> second
  REAL(REAL64), PARAMETER :: h2s = 3600._REAL64 ! hour -> second
  REAL(REAL64), PARAMETER :: s2h = one / h2s ! second -> hour
  REAL(REAL64), PARAMETER :: m2s = 60._REAL64 ! minute -> second
  REAL(REAL64), PARAMETER :: s2m = one / m2s ! second -> minute

  INTEGER, PARAMETER :: MAXSTRLEN = 99 ! maximum string length for strftime

  TYPE :: datetime

    PRIVATE

    INTEGER :: year = 1 ! year [1-HUGE(year)]
    INTEGER :: month = 1 ! month in year [1-12]
    INTEGER :: day = 1 ! day in month [1-31]
    INTEGER :: hour = 0 ! hour in day [0-23]
    INTEGER :: minute = 0 ! minute in hour [0-59]
    INTEGER :: second = 0 ! second in minute [0-59]
    INTEGER :: millisecond = 0 ! milliseconds in second [0-999]
    REAL(REAL64) :: tz = 0 ! timezone offset from UTC [hours]

  CONTAINS

    ! getter functions
    PROCEDURE, PASS(self), PUBLIC :: getYear
    PROCEDURE, PASS(self), PUBLIC :: getMonth
    PROCEDURE, PASS(self), PUBLIC :: getDay
    PROCEDURE, PASS(self), PUBLIC :: getHour
    PROCEDURE, PASS(self), PUBLIC :: getMinute
    PROCEDURE, PASS(self), PUBLIC :: getSecond
    PROCEDURE, PASS(self), PUBLIC :: getMillisecond
    PROCEDURE, PASS(self), PUBLIC :: getTz

    ! public methods
    PROCEDURE, PASS(self), PUBLIC :: isocalendar
    PROCEDURE, PASS(self), PUBLIC :: isoformat
    PROCEDURE, PASS(self), PUBLIC :: isValid
    PROCEDURE, NOPASS, PUBLIC :: now
    PROCEDURE, PASS(self), PUBLIC :: secondsSinceEpoch
    PROCEDURE, PASS(self), PUBLIC :: strftime
    PROCEDURE, PASS(self), PUBLIC :: tm
    PROCEDURE, PASS(self), PUBLIC :: tzOffset
    PROCEDURE, PASS(self), PUBLIC :: isoweekday
    PROCEDURE, PASS(self), PUBLIC :: isoweekdayLong
    PROCEDURE, PASS(self), PUBLIC :: isoweekdayShort
    PROCEDURE, PASS(self), PUBLIC :: utc
    PROCEDURE, PASS(self), PUBLIC :: weekday
    PROCEDURE, PASS(self), PUBLIC :: weekdayLong
    PROCEDURE, PASS(self), PUBLIC :: weekdayShort
    PROCEDURE, PASS(self), PUBLIC :: yearday

    ! private methods
    PROCEDURE, PASS(self), PRIVATE :: addMilliseconds
    PROCEDURE, PASS(self), PRIVATE :: addSeconds
    PROCEDURE, PASS(self), PRIVATE :: addMinutes
    PROCEDURE, PASS(self), PRIVATE :: addHours
    PROCEDURE, PASS(self), PRIVATE :: addDays

    ! operator overloading procedures
    PROCEDURE, PASS(d0), PRIVATE :: datetime_plus_timedelta
    PROCEDURE, PASS(d0), PRIVATE :: timedelta_plus_datetime
    PROCEDURE, PASS(d0), PRIVATE :: datetime_minus_datetime
    PROCEDURE, PASS(d0), PRIVATE :: datetime_minus_timedelta
    PROCEDURE, PASS(d0), PRIVATE :: datetime_eq
    PROCEDURE, PASS(d0), PRIVATE :: datetime_neq
    PROCEDURE, PASS(d0), PRIVATE :: datetime_gt
    PROCEDURE, PASS(d0), PRIVATE :: datetime_ge
    PROCEDURE, PASS(d0), PRIVATE :: datetime_lt
    PROCEDURE, PASS(d0), PRIVATE :: datetime_le

    GENERIC :: OPERATOR(+) => datetime_plus_timedelta, &
      timedelta_plus_datetime
    GENERIC :: OPERATOR(-) => datetime_minus_datetime, &
      datetime_minus_timedelta
    GENERIC :: OPERATOR(==) => datetime_eq
    GENERIC :: OPERATOR(/=) => datetime_neq
    GENERIC :: OPERATOR(>) => datetime_gt
    GENERIC :: OPERATOR(>=) => datetime_ge
    GENERIC :: OPERATOR(<) => datetime_lt
    GENERIC :: OPERATOR(<=) => datetime_le

  END TYPE datetime

  INTERFACE datetime
    MODULE PROCEDURE :: datetime_constructor
  END INTERFACE datetime

  TYPE :: timedelta
    PRIVATE

    INTEGER :: days = 0
    INTEGER :: hours = 0
    INTEGER :: minutes = 0
    INTEGER :: seconds = 0
    INTEGER :: milliseconds = 0

  CONTAINS

    PROCEDURE, PASS(self), PUBLIC :: getDays
    PROCEDURE, PASS(self), PUBLIC :: getHours
    PROCEDURE, PASS(self), PUBLIC :: getMinutes
    PROCEDURE, PASS(self), PUBLIC :: getSeconds
    PROCEDURE, PASS(self), PUBLIC :: getMilliseconds

    PROCEDURE, PUBLIC :: total_seconds

    PROCEDURE, PRIVATE :: timedelta_plus_timedelta
    PROCEDURE, PRIVATE :: timedelta_minus_timedelta
    PROCEDURE, PRIVATE :: unary_minus_timedelta
    PROCEDURE, PRIVATE :: timedelta_eq
    PROCEDURE, PRIVATE :: timedelta_neq
    PROCEDURE, PRIVATE :: timedelta_gt
    PROCEDURE, PRIVATE :: timedelta_ge
    PROCEDURE, PRIVATE :: timedelta_lt
    PROCEDURE, PRIVATE :: timedelta_le

    GENERIC :: OPERATOR(+) => timedelta_plus_timedelta
    GENERIC :: OPERATOR(-) => timedelta_minus_timedelta, unary_minus_timedelta
    GENERIC :: OPERATOR(==) => timedelta_eq
    GENERIC :: OPERATOR(/=) => timedelta_neq
    GENERIC :: OPERATOR(>) => timedelta_gt
    GENERIC :: OPERATOR(>=) => timedelta_ge
    GENERIC :: OPERATOR(<) => timedelta_lt
    GENERIC :: OPERATOR(<=) => timedelta_le

  END TYPE timedelta

  INTERFACE timedelta
    MODULE PROCEDURE :: timedelta_constructor
  END INTERFACE timedelta

  TYPE, BIND(c) :: tm_struct
    ! Derived type for compatibility with C and C++ struct tm.
    ! Enables calling strftime and strptime using iso_c_binding.
    ! See http://www.cplusplus.com/reference/ctime/tm for reference.
    INTEGER(C_INT) :: tm_sec = 0 ! Seconds [0-60] (1 leap second)
    INTEGER(C_INT) :: tm_min = 0 ! Minutes [0-59]
    INTEGER(C_INT) :: tm_hour = 0 ! Hours [0-23]
    INTEGER(C_INT) :: tm_mday = 0 ! Day [1-31]
    INTEGER(C_INT) :: tm_mon = 0 ! Month [0-11]
    INTEGER(C_INT) :: tm_year = 0 ! Year - 1900
    INTEGER(C_INT) :: tm_wday = 0 ! Day of week [0-6]
    INTEGER(C_INT) :: tm_yday = 0 ! Days in year [0-365]
    INTEGER(C_INT) :: tm_isdst = 0 ! DST [-1/0/1]
  END TYPE tm_struct

  INTERFACE

    FUNCTION c_strftime(str, slen, FORMAT, tm) BIND(c, name='strftime') RESULT(rc)
      ! Returns a formatted time string, given input time struct and format.
      ! See https://www.cplusplus.com/reference/ctime/strftime for reference.
      IMPORT :: C_CHAR, C_INT, tm_struct
      CHARACTER(kind=C_CHAR), INTENT(out) :: str(*) ! result string
      INTEGER(C_INT), VALUE, INTENT(in) :: slen ! string length
      CHARACTER(kind=C_CHAR), INTENT(in) :: FORMAT(*) ! time format
      TYPE(tm_struct), INTENT(in) :: tm ! tm_struct instance
      INTEGER(C_INT) :: rc ! return code
    END FUNCTION c_strftime

    FUNCTION c_strptime(str, FORMAT, tm) BIND(c, name='strptime') RESULT(rc)
      ! Interface to POSIX strptime.
      ! Returns a time struct object based on the input time string str and format.
      ! See http://man7.org/linux/man-pages/man3/strptime.3.html for reference.
      IMPORT :: C_CHAR, C_INT, tm_struct
      CHARACTER(kind=C_CHAR), INTENT(in) :: str(*) ! input string
      CHARACTER(kind=C_CHAR), INTENT(in) :: FORMAT(*) ! time format
      TYPE(tm_struct), INTENT(out) :: tm ! result tm_struct
      INTEGER(C_INT) :: rc ! return code
    END FUNCTION c_strptime

  END INTERFACE

  TYPE :: clock
    TYPE(datetime) :: startTime
    TYPE(datetime) :: stopTime
    TYPE(datetime) :: currentTime
    TYPE(timedelta) :: tickInterval
    LOGICAL :: alarm = .FALSE.
    LOGICAL :: started = .FALSE.
    LOGICAL :: stopped = .FALSE.
  CONTAINS
    PROCEDURE :: reset
    PROCEDURE :: tick
  END TYPE clock

CONTAINS

  PURE ELEMENTAL SUBROUTINE reset(self)
    ! Resets the clock to its start time.
    CLASS(clock), INTENT(in out) :: self
    self%currentTime = self%startTime
    self%started = .FALSE.
    self%stopped = .FALSE.
  END SUBROUTINE reset

  PURE ELEMENTAL SUBROUTINE tick(self)
    ! Increments the currentTime of the clock instance by one tickInterval.
    CLASS(clock), INTENT(in out) :: self
    IF (self%stopped) RETURN
    IF (.NOT. self%started) THEN
      self%started = .TRUE.
      self%currentTime = self%startTime
    END IF
    self%currentTime = self%currentTime + self%tickInterval
    IF (self%currentTime >= self%stopTime) self%stopped = .TRUE.
  END SUBROUTINE tick

  PURE ELEMENTAL TYPE(datetime) FUNCTION datetime_constructor( &
    year, month, day, hour, minute, second, millisecond, tz)
    ! Constructor function for the `datetime` class.
    INTEGER, INTENT(in), OPTIONAL :: year, month, day, hour, minute, second, millisecond
    REAL(REAL64), INTENT(in), OPTIONAL :: tz ! timezone offset in hours

    datetime_constructor%year = 1
    IF (PRESENT(year)) datetime_constructor%year = year

    datetime_constructor%month = 1
    IF (PRESENT(month)) datetime_constructor%month = month

    datetime_constructor%day = 1
    IF (PRESENT(day)) datetime_constructor%day = day

    datetime_constructor%hour = 0
    IF (PRESENT(hour)) datetime_constructor%hour = hour

    datetime_constructor%minute = 0
    IF (PRESENT(minute)) datetime_constructor%minute = minute

    datetime_constructor%second = 0
    IF (PRESENT(second)) datetime_constructor%second = second

    datetime_constructor%millisecond = 0
    IF (PRESENT(millisecond)) datetime_constructor%millisecond = millisecond

    datetime_constructor%tz = 0
    IF (PRESENT(tz)) datetime_constructor%tz = tz

  END FUNCTION datetime_constructor

  PURE ELEMENTAL INTEGER FUNCTION getYear(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getYear = self%year
  END FUNCTION getYear

  PURE ELEMENTAL INTEGER FUNCTION getMonth(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getMonth = self%month
  END FUNCTION getMonth

  PURE ELEMENTAL INTEGER FUNCTION getDay(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getDay = self%day
  END FUNCTION getDay

  PURE ELEMENTAL INTEGER FUNCTION getHour(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getHour = self%hour
  END FUNCTION getHour

  PURE ELEMENTAL INTEGER FUNCTION getMinute(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getMinute = self%minute
  END FUNCTION getMinute

  PURE ELEMENTAL INTEGER FUNCTION getSecond(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getSecond = self%second
  END FUNCTION getSecond

  PURE ELEMENTAL INTEGER FUNCTION getMillisecond(self)
    ! Returns the year component
    CLASS(datetime), INTENT(in) :: self
    getMillisecond = self%millisecond
  END FUNCTION getMillisecond

  PURE ELEMENTAL REAL(REAL64) FUNCTION getTz(self)
    ! Returns the timezone offset component
    CLASS(datetime), INTENT(in) :: self
    getTz = self%tz
  END FUNCTION getTz

  PURE ELEMENTAL SUBROUTINE addMilliseconds(self, ms)
    ! Adds an integer number of milliseconds to self. Called by `datetime`
    ! addition (`+`) and subtraction (`-`) operators.
    CLASS(datetime), INTENT(in out) :: self
    INTEGER, INTENT(in) :: ms
    self%millisecond = self%millisecond + ms
    DO
      IF (self%millisecond >= 1000) THEN
        CALL self%addSeconds(self%millisecond / 1000)
        self%millisecond = MOD(self%millisecond, 1000)
      ELSE IF (self%millisecond < 0) THEN
        CALL self%addSeconds(self%millisecond / 1000 - 1)
        self%millisecond = MOD(self%millisecond, 1000) + 1000
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE addMilliseconds

  PURE ELEMENTAL SUBROUTINE addSeconds(self, s)
    ! Adds an integer number of seconds to self. Called by `datetime`
    ! addition (`+`) and subtraction (`-`) operators.
    CLASS(datetime), INTENT(in out) :: self
    INTEGER, INTENT(in) :: s
    self%second = self%second + s
    DO
      IF (self%second >= 60) THEN
        CALL self%addMinutes(self%second / 60)
        self%second = MOD(self%second, 60)
      ELSE IF (self%second < 0) THEN
        CALL self%addMinutes(self%second / 60 - 1)
        self%second = MOD(self%second, 60) + 60
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE addSeconds

  PURE ELEMENTAL SUBROUTINE addMinutes(self, m)
    ! Adds an integer number of minutes to self. Called by `datetime`
    ! addition (`+`) and subtraction (`-`) operators.
    CLASS(datetime), INTENT(in out) :: self
    INTEGER, INTENT(in) :: m
    self%minute = self%minute + m
    DO
      IF (self%minute >= 60) THEN
        CALL self%addHours(self%minute / 60)
        self%minute = MOD(self%minute, 60)
      ELSE IF (self%minute < 0) THEN
        CALL self%addHours(self%minute / 60 - 1)
        self%minute = MOD(self%minute, 60) + 60
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE addMinutes

  PURE ELEMENTAL SUBROUTINE addHours(self, h)
    ! Adds an integer number of hours to self. Called by `datetime`
    ! addition (`+`) and subtraction (`-`) operators.
    CLASS(datetime), INTENT(in out) :: self
    INTEGER, INTENT(in) :: h
    self%hour = self%hour + h
    DO
      IF (self%hour >= 24) THEN
        CALL self%addDays(self%hour / 24)
        self%hour = MOD(self%hour, 24)
      ELSE IF (self%hour < 0) THEN
        CALL self%addDays(self%hour / 24 - 1)
        self%hour = MOD(self%hour, 24) + 24
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE addHours

  PURE ELEMENTAL SUBROUTINE addDays(self, d)
    ! Adds an integer number of dayss to self. Called by `datetime`
    ! addition (`+`) and subtraction (`-`) operators.
    CLASS(datetime), INTENT(in out) :: self
    INTEGER, INTENT(in) :: d
    INTEGER :: daysInCurrentMonth
    self%day = self%day + d
    DO
      daysInCurrentMonth = daysInMonth(self%month, self%year)
      IF (self%day > daysInCurrentMonth) THEN
        self%day = self%day - daysInCurrentMonth
        self%month = self%month + 1
        IF (self%month > 12) THEN
          self%year = self%year + self%month / 12
          self%month = MOD(self%month, 12)
        END IF
      ELSE IF (self%day < 1) THEN
        self%month = self%month - 1
        IF (self%month < 1) THEN
          self%year = self%year + self%month / 12 - 1
          self%month = 12 + MOD(self%month, 12)
        END IF
        self%day = self%day + daysInMonth(self%month, self%year)
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE addDays

  PURE ELEMENTAL CHARACTER(23) FUNCTION isoformat(self, sep)
    ! Returns character string with time in ISO 8601 format.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER, INTENT(in), OPTIONAL :: sep
    CHARACTER :: separator

    separator = 'T'
    IF (PRESENT(sep)) separator = sep

    ! TODO below is a bit cumbersome and was implemented
    ! at a time before the interface to strftime. Now we
    ! could do something like:
    !
    ! isoformat = self % strftime('%Y-%m-%d'//separator//'%H:%M:%S')
    !
    isoformat = int2str(self%year, 4)//'-'// &
                int2str(self%month, 2)//'-'// &
                int2str(self%day, 2)//separator// &
                int2str(self%hour, 2)//':'// &
                int2str(self%minute, 2)//':'// &
                int2str(self%second, 2)//'.'// &
                int2str(self%millisecond, 3)

  END FUNCTION isoformat

  PURE ELEMENTAL LOGICAL FUNCTION isValid(self)
    ! Checks whether the `datetime` instance has valid component values.
    ! Returns `.true.` if the `datetime` instance is valid, and `.false.`
    ! otherwise.
    CLASS(datetime), INTENT(in) :: self

    ! assume valid
    isValid = .TRUE.

    IF (self%year < 1) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%month < 1 .OR. self%month > 12) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%day < 1 .OR. &
        self%day > daysInMonth(self%month, self%year)) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%hour < 0 .OR. self%hour > 23) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%minute < 0 .OR. self%minute > 59) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%second < 0 .OR. self%second > 59) THEN
      isValid = .FALSE.
      RETURN
    END IF

    IF (self%millisecond < 0 .OR. self%millisecond > 999) THEN
      isValid = .FALSE.
      RETURN
    END IF

  END FUNCTION isValid

  TYPE(datetime) FUNCTION now()
    ! Returns a `datetime` instance with current time.
    CHARACTER(5) :: zone
    INTEGER :: values(8)
    INTEGER :: hour, minute

    ! Obtain local machine time zone information
    CALL DATE_AND_TIME(zone=zone, values=values)

    READ (zone(1:3), '(i3)') hour
    READ (zone(4:5), '(i2)') minute

    now = datetime(year=values(1), month=values(2), day=values(3), &
                   hour=values(5), minute=values(6), second=values(7), &
                   millisecond=values(8))

    now%tz = hour + minute * m2h

  END FUNCTION now

  PURE ELEMENTAL INTEGER FUNCTION weekday(self)
    ! Returns the day of the week calculated using Zeller's congruence.
    ! Returned value is an integer scalar in the range [0-6], such that:
    !
    ! 0: Sunday
    ! 1: Monday
    ! 2: Tuesday
    ! 3: Wednesday
    ! 4: Thursday
    ! 5: Friday
    ! 6: Saturday
    CLASS(datetime), INTENT(in) :: self
    INTEGER :: year, month, j, k

    year = self%year
    month = self%month

    IF (month <= 2) THEN
      month = month + 12
      year = year - 1
    END IF

    j = year / 100
    k = MOD(year, 100)

    weekday = MOD(self%day + ((month + 1) * 26) / 10 + k + k / 4 + j / 4 + 5 * j, 7) - 1

    IF (weekday < 0) weekday = 6

  END FUNCTION weekday

  PURE ELEMENTAL INTEGER FUNCTION isoweekday(self)
    ! Returns the day of the week per ISO 8601 returned from weekday().
    ! Returned value is an integer scalar in the range [1-7].
    CLASS(datetime), INTENT(in) :: self
    isoweekday = self%weekday()
    IF (isoweekday == 0) isoweekday = 7
  END FUNCTION isoweekday

  PURE ELEMENTAL CHARACTER(9) FUNCTION weekdayLong(self)
    ! Returns the full name of the day of the week.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER(9), PARAMETER :: &
      days(*) = ['Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', &
                 'Thursday ', 'Friday   ', 'Saturday ']
    weekdayLong = days(self%weekday() + 1)
  END FUNCTION weekdayLong

  PURE ELEMENTAL CHARACTER(9) FUNCTION isoweekdayLong(self)
    ! Returns the full name of the day of the week for ISO 8601
    ! ordered weekdays.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER(9), PARAMETER :: &
      days(7) = ['Monday   ', 'Tuesday  ', 'Wednesday', 'Thursday ', &
                 'Friday   ', 'Saturday ', 'Sunday   ']
    isoweekdayLong = days(self%isoweekday())
  END FUNCTION isoweekdayLong

  PURE ELEMENTAL CHARACTER(3) FUNCTION weekdayShort(self)
    ! Returns the short (3-letter) name of the day of the week.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER(3), PARAMETER :: days(7) = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
    weekdayShort = days(self%weekday() + 1)
  END FUNCTION weekdayShort

  PURE ELEMENTAL CHARACTER(3) FUNCTION isoweekdayShort(self)
    ! Returns the short (3-letter) name of the day of the week
    ! based on ISO 8601 ordering.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER(3), PARAMETER :: days(7) = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
    isoweekdayShort = days(self%isoweekday())
  END FUNCTION isoweekdayShort

  FUNCTION isocalendar(self)
    ! Returns an array of 3 integers, year, week number, and week day,
    ! as defined by ISO 8601 week date. Essentially a wrapper around C
    ! `strftime` function.
    CLASS(datetime), INTENT(in) :: self
    INTEGER :: isocalendar(3)
    INTEGER :: year, week, wday, rc
    CHARACTER(20) :: string

    rc = c_strftime(string, LEN(string), '%G %V %u'//C_NULL_CHAR, self%tm())

    READ (string(1:4), '(i4)') year
    READ (string(6:7), '(i2)') week
    READ (string(9:9), '(i1)') wday

    isocalendar = [year, week, wday]

  END FUNCTION isocalendar

  INTEGER(INT64) FUNCTION secondsSinceEpoch(self)
    ! Returns an integer number of seconds since the UNIX Epoch (1 Jan 1970).
    ! Since Windows does not have strftime('%s'), we implement this using
    ! datetime itself.
    CLASS(datetime), INTENT(in) :: self
    TYPE(timedelta) :: delta
    TYPE(datetime) :: this_time, unix_time

    this_time = datetime(self%year, self%month, self%day, &
                         self%hour, self%minute, self%second)
    unix_time = datetime(1970, 1, 1, 0, 0, 0)
    delta = this_time - unix_time
    secondsSinceEpoch = delta%total_seconds()

  END FUNCTION secondsSinceEpoch

  FUNCTION strftime(self, FORMAT)
    ! Wrapper around C and C++ `strftime` function.
    CLASS(datetime), INTENT(in) :: self
    CHARACTER(*), INTENT(in)  :: FORMAT
    CHARACTER(:), ALLOCATABLE :: strftime
    INTEGER :: rc
    CHARACTER(MAXSTRLEN) :: resultString
    resultString = ""
    rc = c_strftime(resultString, LEN(resultString), TRIM(FORMAT)//C_NULL_CHAR, &
                    self%tm())
    IF (rc == 0) WRITE (stderr, *) "datetime:strftime failure, format: ", TRIM(FORMAT)
    strftime = resultString(1:LEN_TRIM(resultString) - 1)  !< strip null
  END FUNCTION strftime

  PURE ELEMENTAL TYPE(tm_struct) FUNCTION tm(self)
    ! Returns a `tm_struct` instance of the current `datetime`.
    CLASS(datetime), INTENT(in) :: self
    tm%tm_sec = self%second
    tm%tm_min = self%minute
    tm%tm_hour = self%hour
    tm%tm_mday = self%day
    tm%tm_mon = self%month - 1
    tm%tm_year = self%year - 1900
    tm%tm_wday = self%weekday()
    tm%tm_yday = self%yearday() - 1
    tm%tm_isdst = -1
  END FUNCTION tm

  PURE ELEMENTAL CHARACTER(5) FUNCTION tzOffset(self)
    ! Returns a character string with timezone offset in hours from UTC,
    ! in format +/-[hh][mm].
    CLASS(datetime), INTENT(in) :: self
    INTEGER :: hours, minutes

    IF (self%tz < 0) THEN
      tzOffset(1:1) = '-'
    ELSE
      tzOffset(1:1) = '+'
    END IF

    hours = INT(ABS(self%tz))
    minutes = NINT((ABS(self%tz) - hours) * 60)

    IF (minutes == 60) THEN
      minutes = 0
      hours = hours + 1
    END IF

    WRITE (tzOffset(2:5), '(2i2.2)') hours, minutes

  END FUNCTION tzOffset

  PURE ELEMENTAL TYPE(datetime) FUNCTION utc(self)
    ! Returns the `datetime` instance at Coordinated Universal Time (UTC).
    CLASS(datetime), INTENT(in) :: self
    INTEGER :: hours, minutes, sgn
    hours = INT(ABS(self%tz))
    minutes = NINT((ABS(self%tz) - hours) * 60)
    sgn = INT(SIGN(one, self%tz))
    utc = self - timedelta(hours=sgn * hours, minutes=sgn * minutes)
    utc%tz = 0
  END FUNCTION utc

  PURE ELEMENTAL INTEGER FUNCTION yearday(self)
    ! Returns the integer day of the year (ordinal date).
    CLASS(datetime), INTENT(in) :: self
    INTEGER :: month
    yearday = 0
    DO month = 1, self%month - 1
      yearday = yearday + daysInMonth(month, self%year)
    END DO
    yearday = yearday + self%day
  END FUNCTION yearday

  PURE ELEMENTAL FUNCTION datetime_plus_timedelta(d0, t) RESULT(d)
    ! Adds a `timedelta` instance to a `datetime` instance, and returns a
    ! new `datetime` instance. Overloads the operator `+`.
    CLASS(datetime), INTENT(in) :: d0
    CLASS(timedelta), INTENT(in) :: t
    TYPE(datetime) :: d

    INTEGER :: milliseconds, seconds, minutes, hours, days

    d = datetime(year=d0%getYear(), &
                 month=d0%getMonth(), &
                 day=d0%getDay(), &
                 hour=d0%getHour(), &
                 minute=d0%getMinute(), &
                 second=d0%getSecond(), &
                 millisecond=d0%getMillisecond(), &
                 tz=d0%getTz())

    milliseconds = t%getMilliseconds()
    seconds = t%getSeconds()
    minutes = t%getMinutes()
    hours = t%getHours()
    days = t%getDays()

    IF (milliseconds /= 0) CALL d%addMilliseconds(milliseconds)
    IF (seconds /= 0) CALL d%addSeconds(seconds)
    IF (minutes /= 0) CALL d%addMinutes(minutes)
    IF (hours /= 0) CALL d%addHours(hours)
    IF (days /= 0) CALL d%addDays(days)

  END FUNCTION datetime_plus_timedelta

  PURE ELEMENTAL FUNCTION timedelta_plus_datetime(t, d0) RESULT(d)
    ! Adds a `timedelta` instance to a `datetime` instance, and returns a
    ! new `datetime` instance. Overloads the operator `+`.
    CLASS(timedelta), INTENT(in) :: t
    CLASS(datetime), INTENT(in) :: d0
    TYPE(datetime) :: d
    d = d0 + t
  END FUNCTION timedelta_plus_datetime

  PURE ELEMENTAL FUNCTION datetime_minus_timedelta(d0, t) RESULT(d)
    ! Subtracts a `timedelta` instance from a `datetime` instance and
    ! returns a new `datetime` instance. Overloads the operator `-`.
    CLASS(datetime), INTENT(in) :: d0
    CLASS(timedelta), INTENT(in) :: t
    TYPE(datetime) :: d
    d = d0 + (-t)
  END FUNCTION datetime_minus_timedelta

  PURE ELEMENTAL FUNCTION datetime_minus_datetime(d0, d1) RESULT(t)
    ! Subtracts a `datetime` instance from another `datetime` instance,
    ! and returns a `timedelta` instance. Overloads the operator `-`.
    CLASS(datetime), INTENT(in) :: d0, d1
    TYPE(timedelta) :: t
    REAL(REAL64) :: daysDiff
    INTEGER :: days, hours, minutes, seconds, milliseconds
    INTEGER :: sign_

    daysDiff = date2num(d0) - date2num(d1)

    IF (daysDiff < 0) THEN
      sign_ = -1
      daysDiff = ABS(daysDiff)
    ELSE
      sign_ = 1
    END IF

    days = INT(daysDiff)
    hours = INT((daysDiff - days) * d2h)
    minutes = INT((daysDiff - days - hours * h2d) * d2m)
    seconds = INT((daysDiff - days - hours * h2d - minutes * m2d) * d2s)
    milliseconds = NINT((daysDiff - days - hours * h2d - minutes * m2d &
                         - seconds * s2d) * d2s * 1E3_REAL64)

    t = timedelta(sign_ * days, sign_ * hours, sign_ * minutes, sign_ * seconds, &
                  sign_ * milliseconds)

  END FUNCTION datetime_minus_datetime

  PURE ELEMENTAL LOGICAL FUNCTION datetime_gt(d0, d1) RESULT(res)
    ! `datetime` comparison operator that returns `.true.` if `d0` is
    ! greater than `d1` and `.false.` otherwise. Overloads the
    ! operator `>`.
    CLASS(datetime), INTENT(in) :: d0, d1
    TYPE(datetime) :: d0_utc, d1_utc

    ! Convert to UTC before making comparison
    d0_utc = d0%utc()
    d1_utc = d1%utc()

    ! Compare years
    IF (d0_utc%year > d1_utc%year) THEN
      res = .TRUE.
    ELSE IF (d0_utc%year < d1_utc%year) THEN
      res = .FALSE.
    ELSE

      ! Compare months
      IF (d0_utc%month > d1_utc%month) THEN
        res = .TRUE.
      ELSE IF (d0_utc%month < d1_utc%month) THEN
        res = .FALSE.
      ELSE

        ! Compare days
        IF (d0_utc%day > d1_utc%day) THEN
          res = .TRUE.
        ELSE IF (d0_utc%day < d1_utc%day) THEN
          res = .FALSE.
        ELSE

          ! Compare hours
          IF (d0_utc%hour > d1_utc%hour) THEN
            res = .TRUE.
          ELSE IF (d0_utc%hour < d1_utc%hour) THEN
            res = .FALSE.
          ELSE

            ! Compare minutes
            IF (d0_utc%minute > d1_utc%minute) THEN
              res = .TRUE.
            ELSE IF (d0_utc%minute < d1_utc%minute) THEN
              res = .FALSE.
            ELSE

              ! Compare seconds
              IF (d0_utc%second > d1_utc%second) THEN
                res = .TRUE.
              ELSE IF (d0_utc%second < d1_utc%second) THEN
                res = .FALSE.
              ELSE

                ! Compare milliseconds
                IF (d0_utc%millisecond > d1_utc%millisecond) THEN
                  res = .TRUE.
                ELSE
                  res = .FALSE.
                END IF

              END IF
            END IF
          END IF
        END IF
      END IF
    END IF

  END FUNCTION datetime_gt

  PURE ELEMENTAL LOGICAL FUNCTION datetime_lt(d0, d1) RESULT(res)
    ! `datetime` comparison operator that returns `.true.` if `d0` is
    ! less than `d1` and `.false.` otherwise. Overloads the operator `<`.
    CLASS(datetime), INTENT(in) :: d0, d1
    res = d1 > d0
  END FUNCTION datetime_lt

  PURE ELEMENTAL LOGICAL FUNCTION datetime_eq(d0, d1) RESULT(res)
    ! `datetime` comparison operator that returns `.true.` if `d0` is
    ! equal to `d1` and `.false.` otherwise. Overloads the operator `==`.
    CLASS(datetime), INTENT(in) :: d0, d1
    TYPE(datetime) :: d0_utc, d1_utc

    ! Convert to UTC before making comparison
    d0_utc = d0%utc()
    d1_utc = d1%utc()

    res = d0_utc%year == d1_utc%year .AND. &
          d0_utc%month == d1_utc%month .AND. &
          d0_utc%day == d1_utc%day .AND. &
          d0_utc%hour == d1_utc%hour .AND. &
          d0_utc%minute == d1_utc%minute .AND. &
          d0_utc%second == d1_utc%second .AND. &
          d0_utc%millisecond == d1_utc%millisecond

  END FUNCTION datetime_eq

  PURE ELEMENTAL LOGICAL FUNCTION datetime_neq(d0, d1) RESULT(res)
    ! `datetime` comparison operator that eturns `.true.` if `d0` is
    ! not equal to `d1` and `.false.` otherwise. Overloads the operator `/=`.
    CLASS(datetime), INTENT(in) :: d0, d1
    res = .NOT. d0 == d1
  END FUNCTION datetime_neq

  PURE ELEMENTAL LOGICAL FUNCTION datetime_ge(d0, d1) RESULT(res)
    ! `datetime` comparison operator. Returns `.true.` if `d0` is greater
    ! than or equal to `d1` and `.false.` otherwise. Overloads the
    ! operator `>=`.
    CLASS(datetime), INTENT(in) :: d0, d1
    res = d0 > d1 .OR. d0 == d1
  END FUNCTION datetime_ge

  PURE ELEMENTAL LOGICAL FUNCTION datetime_le(d0, d1) RESULT(res)
    ! `datetime` comparison operator. Returns `.true.` if `d0` is less
    ! than or equal to `d1`, and `.false.` otherwise. Overloads the
    ! operator `<=`.
    CLASS(datetime), INTENT(in) :: d0, d1
    res = d1 > d0 .OR. d0 == d1
  END FUNCTION datetime_le

  PURE ELEMENTAL LOGICAL FUNCTION isLeapYear(year)
    ! Returns `.true.` if year is leap year and `.false.` otherwise.
    INTEGER, INTENT(in) :: year
    isLeapYear = (MOD(year, 4) == 0 .AND. .NOT. MOD(year, 100) == 0) &
                 .OR. (MOD(year, 400) == 0)
  END FUNCTION isLeapYear

  PURE FUNCTION datetimeRange(d0, d1, t)
    ! Given start and end `datetime` instances `d0` and `d1` and time
    ! increment as `timedelta` instance `t`, returns an array of
    ! `datetime` instances. The number of elements is the number of whole
    ! time increments contained between datetimes `d0` and `d1`.
    TYPE(datetime), INTENT(in) :: d0, d1
    TYPE(timedelta), INTENT(in) :: t
    REAL(REAL64) :: datenum0, datenum1, eps, increment
    TYPE(datetime), ALLOCATABLE :: datetimeRange(:)
    INTEGER :: n, nm
    eps = 1E-10_REAL64
    datenum0 = date2num(d0)
    datenum1 = date2num(d1)
    increment = t%total_seconds() * s2d
    nm = FLOOR((datenum1 - datenum0 + eps) / increment) + 1
    ALLOCATE (datetimeRange(nm))
    DO n = 1, nm
      datetimeRange(n) = num2date(datenum0 + (n - 1) * increment)
    END DO
  END FUNCTION datetimeRange

  PURE ELEMENTAL INTEGER FUNCTION daysInMonth(month, year)
    ! Given integer month and year, returns an integer number
    ! of days in that particular month.
    INTEGER, INTENT(in) :: month, year

    INTEGER, PARAMETER :: days(*) = [31, 28, 31, 30, 31, 30, &
                                     31, 31, 30, 31, 30, 31]

    IF (month < 1 .OR. month > 12) THEN
      ! Should raise an error and abort here, however we want to keep
      ! the pure and elemental attributes. Make sure this function is
      ! called with the month argument in range.
      daysInMonth = 0
      RETURN
    END IF

    IF (month == 2 .AND. isLeapYear(year)) THEN
      daysInMonth = 29
    ELSE
      daysInMonth = days(month)
    END IF

  END FUNCTION daysInMonth

  PURE ELEMENTAL INTEGER FUNCTION daysInYear(year)
    ! Returns the number of days in year.
    INTEGER, INTENT(in) :: year
    IF (isLeapYear(year)) THEN
      daysInYear = 366
    ELSE
      daysInYear = 365
    END IF
  END FUNCTION daysInYear

  PURE ELEMENTAL REAL(REAL64) FUNCTION date2num(d)
    ! Given a datetime instance d, returns number of days since
    ! `0001-01-01 00:00:00`, taking into account the timezone offset.
    TYPE(datetime), INTENT(in) :: d
    TYPE(datetime) :: d_utc
    INTEGER :: year

    ! Convert to UTC first
    d_utc = d%utc()

    ! d_utc % year must be positive:
    IF (d_utc%year < 1) THEN
      date2num = 0
      RETURN
    END IF

    date2num = 0
    DO year = 1, d_utc%year - 1
      date2num = date2num + daysInYear(year)
    END DO

    date2num = date2num &
               + d_utc%yearday() &
               + d_utc%hour * h2d &
               + d_utc%minute * m2d &
               + (d_utc%second + 1E-3_REAL64 * d_utc%millisecond) * s2d

  END FUNCTION date2num

  PURE ELEMENTAL TYPE(datetime) FUNCTION num2date(num)
    ! Given number of days since `0001-01-01 00:00:00`, returns a
    ! correspoding `datetime` instance.
    REAL(REAL64), INTENT(in) :: num
    INTEGER :: year, month, day, hour, minute, second, millisecond
    REAL(REAL64) :: days, totseconds

    ! num must be positive
    IF (num < 0) THEN
      num2date = datetime(1)
      RETURN
    END IF

    days = num

    year = 1
    DO
      IF (INT(days) <= daysInYear(year)) EXIT
      days = days - daysInYear(year)
      year = year + 1
    END DO

    month = 1
    DO
      IF (INT(days) <= daysInMonth(month, year)) EXIT
      days = days - daysInMonth(month, year)
      month = month + 1
    END DO

    day = INT(days)
    totseconds = (days - day) * d2s
    hour = INT(totseconds * s2h)
    minute = INT((totseconds - hour * h2s) * s2m)
    second = INT(totseconds - hour * h2s - minute * m2s)
    millisecond = NINT((totseconds - INT(totseconds)) * 1E3_REAL64)

    num2date = datetime(year, month, day, hour, minute, second, millisecond, tz=zero)

    ! Handle a special case caused by floating-point arithmethic:
    IF (num2date%millisecond == 1000) THEN
      num2date%millisecond = 0
      CALL num2date%addSeconds(1)
    END IF

    IF (num2date%second == 60) THEN
      num2date%second = 0
      CALL num2date%addMinutes(1)
    END IF
    IF (num2date%minute == 60) THEN
      num2date%minute = 0
      CALL num2date%addHours(1)
    END IF
    IF (num2date%hour == 24) THEN
      num2date%hour = 0
      CALL num2date%addDays(1)
    END IF

  END FUNCTION num2date

  TYPE(datetime) FUNCTION strptime(str, FORMAT)
    ! A wrapper function around C/C++ strptime function.
    ! Returns a `datetime` instance.
    CHARACTER(*), INTENT(in) :: str, FORMAT
    INTEGER :: rc
    TYPE(tm_struct) :: tm
    rc = c_strptime(TRIM(str)//C_NULL_CHAR, TRIM(FORMAT)//C_NULL_CHAR, tm)
    strptime = tm2date(tm)
  END FUNCTION strptime

  PURE ELEMENTAL TYPE(datetime) FUNCTION tm2date(ctime)
    ! Given a `tm_struct` instance, returns a corresponding `datetime`
    ! instance.
    TYPE(tm_struct), INTENT(in) :: ctime

    tm2date%millisecond = 0
    tm2date%second = ctime%tm_sec
    tm2date%minute = ctime%tm_min
    tm2date%hour = ctime%tm_hour
    tm2date%day = ctime%tm_mday
    tm2date%month = ctime%tm_mon + 1
    tm2date%year = ctime%tm_year + 1900
    tm2date%tz = 0

  END FUNCTION tm2date

  PURE FUNCTION int2str(i, length)
    ! Converts an integer `i` into a character string of requested length,
    ! pre-pending zeros if necessary.
    INTEGER, INTENT(in) :: i, length
    CHARACTER(length) :: int2str
    CHARACTER(2) :: string
    WRITE (string, '(i2)') length
    WRITE (int2str, '(i'//string//'.'//string//')') i
  END FUNCTION int2str

  PURE ELEMENTAL TYPE(timedelta) FUNCTION timedelta_constructor( &
    days, hours, minutes, seconds, milliseconds)
    ! Constructor function for the `timedelta` class.
    INTEGER, INTENT(in), OPTIONAL :: days, hours, minutes, seconds, milliseconds

    timedelta_constructor%days = 0
    IF (PRESENT(days)) timedelta_constructor%days = days

    timedelta_constructor%hours = 0
    IF (PRESENT(hours)) timedelta_constructor%hours = hours

    timedelta_constructor%minutes = 0
    IF (PRESENT(minutes)) timedelta_constructor%minutes = minutes

    timedelta_constructor%seconds = 0
    IF (PRESENT(seconds)) timedelta_constructor%seconds = seconds

    timedelta_constructor%milliseconds = 0
    IF (PRESENT(milliseconds)) timedelta_constructor%milliseconds = milliseconds

  END FUNCTION timedelta_constructor

  PURE ELEMENTAL INTEGER FUNCTION getDays(self)
    ! Returns the number of days.
    CLASS(timedelta), INTENT(in) :: self
    getDays = self%days
  END FUNCTION getDays

  PURE ELEMENTAL INTEGER FUNCTION getHours(self)
    ! Returns the number of hours.
    CLASS(timedelta), INTENT(in) :: self
    getHours = self%hours
  END FUNCTION getHours

  PURE ELEMENTAL INTEGER FUNCTION getMinutes(self)
    ! Returns the number of minutes.
    CLASS(timedelta), INTENT(in) :: self
    getMinutes = self%minutes
  END FUNCTION getMinutes

  PURE ELEMENTAL INTEGER FUNCTION getSeconds(self)
    ! Returns the number of seconds.
    CLASS(timedelta), INTENT(in) :: self
    getSeconds = self%seconds
  END FUNCTION getSeconds

  PURE ELEMENTAL INTEGER FUNCTION getMilliseconds(self)
    ! Returns the number of milliseconds.
    CLASS(timedelta), INTENT(in) :: self
    getMilliseconds = self%milliseconds
  END FUNCTION getMilliseconds

  PURE ELEMENTAL REAL(REAL64) FUNCTION total_seconds(self)
    ! Returns a total number of seconds contained in a `timedelta`
    ! instance.
    CLASS(timedelta), INTENT(in) :: self
    total_seconds = self%days * 86400._REAL64 &
                    + self%hours * 3600._REAL64 &
                    + self%minutes * 60._REAL64 &
                    + self%seconds &
                    + self%milliseconds * 1E-3_REAL64
  END FUNCTION total_seconds

  PURE ELEMENTAL TYPE(timedelta) FUNCTION timedelta_plus_timedelta(t0, t1) RESULT(t)
    ! Adds two `timedelta` instances together and returns a `timedelta`
    ! instance. Overloads the operator `+`.
    CLASS(timedelta), INTENT(in) :: t0, t1
    t = timedelta(days=t0%days + t1%days, &
                  hours=t0%hours + t1%hours, &
                  minutes=t0%minutes + t1%minutes, &
                  seconds=t0%seconds + t1%seconds, &
                  milliseconds=t0%milliseconds + t1%milliseconds)
  END FUNCTION timedelta_plus_timedelta

  PURE ELEMENTAL TYPE(timedelta) FUNCTION timedelta_minus_timedelta(t0, t1) RESULT(t)
    ! Subtracts a `timedelta` instance from another. Returns a
    ! `timedelta` instance. Overloads the operator `-`.
    CLASS(timedelta), INTENT(in) :: t0, t1
    t = t0 + (-t1)
  END FUNCTION timedelta_minus_timedelta

  PURE ELEMENTAL TYPE(timedelta) FUNCTION unary_minus_timedelta(t0) RESULT(t)
    ! Takes a negative of a `timedelta` instance. Overloads the operator `-`.
    CLASS(timedelta), INTENT(in) :: t0
    t%days = -t0%days
    t%hours = -t0%hours
    t%minutes = -t0%minutes
    t%seconds = -t0%seconds
    t%milliseconds = -t0%milliseconds
  END FUNCTION unary_minus_timedelta

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_eq(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if `td0`
    ! is equal to `td1` and `.false.` otherwise. Overloads the operator
    ! `==`.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() == td1%total_seconds()
  END FUNCTION timedelta_eq

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_neq(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if `td0`
    ! is not equal to `td1` and `.false.` otherwise. Overloads the
    ! operator `/=`.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() /= td1%total_seconds()
  END FUNCTION timedelta_neq

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_gt(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if
    ! `td0` is greater than `td1` and `.false.` otherwise. Overloads the
    ! operator `>`.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() > td1%total_seconds()
  END FUNCTION timedelta_gt

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_ge(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if `td0`
    ! is greater than or equal to `td1` and `.false.` otherwise.
    ! Overloads the operator >=.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() >= td1%total_seconds()
  END FUNCTION timedelta_ge

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_lt(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if `td0`
    ! is less than `td1` and `.false.` otherwise. Overloads the operator
    ! `<`.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() < td1%total_seconds()
  END FUNCTION timedelta_lt

  PURE ELEMENTAL LOGICAL FUNCTION timedelta_le(td0, td1) RESULT(res)
    ! `timedelta` object comparison operator. Returns `.true.` if `td0`
    ! is less than or equal to `td1` and `.false.` otherwise. Overloads
    ! the operator `<=`.
    CLASS(timedelta), INTENT(in) :: td0, td1
    res = td0%total_seconds() <= td1%total_seconds()
  END FUNCTION timedelta_le

END MODULE datetime_m
