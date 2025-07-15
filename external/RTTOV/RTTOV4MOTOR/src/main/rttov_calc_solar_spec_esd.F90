! Description:
!> @file
!!   Compute solar irradiance values scaled by Earth-Sun distance
!
!> @brief
!!   Compute solar irradiance values scaled by Earth-Sun distance
!!
!! @details
!!   The coefficient files contain top-of-atmosphere channel-integrated
!!   solar irradiance values. For each profile with a valid date specified
!!   this subroutine scales the solar irradiance according to the Earth-
!!   Sun distance. If no valid date is specified the values from the
!!   coefficient file are used without modification.
!!
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[out]    solar_spec_esd   output top-of-atmosphere solar irradiance values
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
SUBROUTINE rttov_calc_solar_spec_esd(coef, chanprof, profiles, solar_spec_esd)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_profile
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE rttov_const, ONLY : pi
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_coef),     INTENT(IN)    :: coef
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  REAL(KIND=jprb),      INTENT(OUT)   :: solar_spec_esd(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: nchanprof, nprofiles
  INTEGER(KIND=jpim) :: i, chan, prof
  INTEGER(KIND=jpim) :: year, month, day, day_of_year
  INTEGER(KIND=jpim) :: days(12)
  REAL(KIND=jprb)    :: d
  REAL(KIND=jprb)    :: esdsq_r(SIZE(profiles))  ! Reciprocal of Earth-sun distance (in AU) squared
!- End of header ------------------------------------------------------

  nchanprof = SIZE(chanprof)
  nprofiles = SIZE(profiles)

  days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  DO prof = 1, nprofiles

    year  = profiles(prof) % date(1)
    month = profiles(prof) % date(2)
    day   = profiles(prof) % date(3)

    ! If the date is not valid, assume Earth-sun distance is 1 AU.
    esdsq_r(prof) = 1.0_jprb

    IF (year > 1950) THEN  ! Default profile year is 1950
      IF (month > 0 .AND. month < 13) THEN

        IF (MOD(year, 4_jpim) == 0 .AND. (.NOT. MOD(year, 100_jpim) == 0 .OR. MOD(year, 400_jpim) == 0)) THEN
          days(2) = 29
        ELSE
          days(2) = 28
        ENDIF

        IF (day > 0 .AND. day <= days(month)) THEN

          day_of_year = SUM(days(1:month-1)) + day - 1

          ! JAH - this approximation is common in textbooks/literature.
          d             = 2.0_jprb * pi * day_of_year/365.0_jprb
          esdsq_r(prof) = 1.00011_jprb + 0.034221_jprb * COS(d) + 0.00128_jprb * SIN(d) + &
                          0.000719_jprb * COS(d*2.0_jprb) + 0.000077_jprb * SIN(d*2.0_jprb)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO i = 1, nchanprof
    chan = chanprof(i)%chan
    prof = chanprof(i)%prof

    solar_spec_esd(i) = coef%ss_solar_spectrum(chan) * esdsq_r(prof)
  ENDDO

END SUBROUTINE rttov_calc_solar_spec_esd
