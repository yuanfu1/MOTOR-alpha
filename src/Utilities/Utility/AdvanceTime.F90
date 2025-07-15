!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : WRFDA
! VERSION           : V 0.0
! HISTORY           :
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2021-10-8, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2021/11/10, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE AdvanceTime_m
  USE kinds_m
  IMPLICIT NONE

CONTAINS

  FUNCTION isLeapYear(year)
    ! check if year is leapyear
    INTEGER, INTENT(IN) :: year
    LOGICAL :: isLeapYear
    IF (MOD(year, 4) .NE. 0) THEN
      isLeapYear = .FALSE.
    ELSE
      isLeapYear = .TRUE.
      IF (MOD(year, 100) == 0 .AND. MOD(year, 400) .NE. 0) isLeapYear = .FALSE.
    END IF
  END FUNCTION isLeapYear

! Code from
! Clive Page,                         e-mail:  cgp
! Dept of Physics & Astronomy,                 (at) le
! University of Leicester.                     (dot) ac (dot) uk
! http://computer-programming-forum.com/49-fortran/27ae6ab276bc1b48.htm
  SUBROUTINE Time_Unix_to_GMT(uTime_sec, gTime)
    IMPLICIT NONE
    INTEGER uTime_sec, gTime(6)
    ! *uTime_sec input Unix system time, seconds since 1970.0
    ! *gTime output Array:1 = year, 2 = month, 3 = date, 4 = hour, 5 = minute, 6 = secs
    ! *-Author Clive Page, Leicester University, UK.1995 - MAY - 2
    ! Note the MJD algorithm only works from years 1901 to 2099.
    INTEGER mjday, nsecs
    REAL day

    mjday = INT(uTime_sec / 86400 + 40587)
    gTime(1) = 1858 + INT((mjday + 321.51) / 365.25)
    day = AINT(MOD(mjday + 262.25, 365.25)) + 0.5
    gTime(2) = 1 + INT(MOD(day / 30.6 + 2.0, 12.0))
    gTime(3) = 1 + INT(MOD(day, 30.6))
    nsecs = MOD(uTime_sec, 86400)
    gTime(6) = MOD(nsecs, 60)
    nsecs = nsecs / 60
    gTime(5) = MOD(nsecs, 60)
    gTime(4) = nsecs / 60
  END SUBROUTINE

  ! SUBROUTINE Unix_Day_Hour_Sec(year,month,day,hour,minute,second,uDay,uHour,uSec)
  SUBROUTINE Time_GMT_to_Unix(gTime, uTime_sec)
    INTEGER(i_kind), INTENT(IN)  :: gTime(6)
    INTEGER(i_kind), INTENT(OUT) :: uTime_sec
    ! INTEGER(i_kind), INTENT(OUT), OPTIONAL :: uDay, uHour

    INTEGER(i_kind) :: day, month, year, hour, minute, second
    INTEGER(i_kind) :: uDay, uHour
    INTEGER(i_kind) :: ndays, m, nleapyr
    INTEGER(i_kind) :: base_year = 1970
    INTEGER(i_kind) :: days_per_mon(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    year = gTime(1)
    month = gTime(2)
    day = gTime(3)
    hour = gTime(4)
    minute = gTime(5)
    second = gTime(6)

    IF (year < base_year) STOP "YEAR can not be before 1970!"

    ! compute number of leap years fully past since base_year
    nleapyr = (year - base_year + 1) / 4 - (year - base_year + 1) / 100 + (year - base_year + 1) / 400

    ! Count up days in this year
    ndays = 0
    DO m = 1, month - 1
      ndays = ndays + days_per_mon(m)
      IF (isLeapYear(year) .AND. m == 2) ndays = ndays + 1
    END DO

    uDay = day - 1 + ndays + 365 * (year - base_year - nleapyr) + 366 * (nleapyr)
    uHour = hour + uDay * 24
    uTime_sec = second + 60 * minute + uHour * 3600

  END SUBROUTINE Time_GMT_to_Unix

  FUNCTION Get_Time_GMT_to_Unix(gTime) RESULT(uTime_sec)
    INTEGER(i_kind), INTENT(IN)  :: gTime(6)
    REAL(r_kind) :: uTime_sec

    INTEGER(i_kind) :: uTime_sec_int
    CALL Time_GMT_to_Unix(gTime, uTime_sec_int)

    uTime_sec = uTime_sec_int
  END FUNCTION Get_Time_GMT_to_Unix

  ! 根据儒略日计算公历年月日时分秒
  ! mjd: Modified Julian Date
  SUBROUTINE mjd2cal(ju,yr,mt,dy,hr,mn,sc)
    IMPLICIT NONE
    REAL(kind=8)::ju,j0
    INTEGER :: yr,mt,dy,mn,yr0,sc
    REAL(kind=4)::bc,dd,n1,n2,n3,sc0
    INTEGER :: hr

    ju=ju+2451545.0
    ! ju = ju + 2400000.5

    if (ju.lt.1721423.5) then
      bc=1
    else
      bc=0
    end if

    if (ju.lt.2299160.5) then
      j0=floor(ju+0.5)
      dd=ju+0.5-j0
    else
      n1=floor((ju-2342031.5)/36524.25/4)+1
      n2=floor((ju-2378555.5)/36524.25/4)+1
      n3=floor((ju-2415079.5)/36524.25/4)+1
      j0=n1+n2+n3+ju+10
      dd=j0+0.5-floor(j0+0.5)
      j0=floor(j0+0.5)
    end if

    j0=j0+32083
    yr0=ceiling(j0/365.25)-1
    yr=yr0-4800
    dy=j0-floor(yr0*365.25)
    mt=floor((dy-0.6)/30.6)+3
    dy=dy-nint((mt-3)*30.6)

    if (mt.gt.12) then
      mt=mt-12
      yr=yr+1
    end if

    yr=yr-bc
    sc0=nint(dd*86400)

    hr=floor(sc0/3600)
    sc0=sc0-hr*3600
    mn=floor(sc0/60)
    sc=int(sc0-mn*60)

  END SUBROUTINE mjd2cal
  
END MODULE AdvanceTime_m
