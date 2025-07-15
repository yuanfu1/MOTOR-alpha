! Description:
!> @file
!!   Calculate brightness temperatures from radiances
!
!> @brief
!!   Calculate brightness temperatures from radiances
!!
!! @details
!!   The Planck function is calculated using the channel central wavenumbers
!!   as specified in the coefficient file. To adjust for the finite spectral
!!   bandwidth the temperatures are modified using band correction coefficients
!!   also stored in the coefficient file.
!!
!!   Radiances units are mW/cm-1/ster/m2 and temperature units are Kelvin.
!!   Brightness temperatures are not calculated for channels at wavelengths
!!   below 3um (i.e. with insignificant thermal emission).
!!
!!   Polarimetric sensors are a special case: for the channels representing
!!   the 3rd and 4th Stokes components the output BTs require the corresponding
!!   BTs for the 1st and 2nd components. The calculations assume that where
!!   channels representing the 3rd and 4th components are being simulated the
!!   corresponding 1st and 2nd component channels immediately precede them in
!!   the chanprof structure. This is enforced by rttov_check_options.
!!
!!   Reference:
!!   Sharp, J.C, 1983: A comparison of approximate methods for converting
!!   between radiance and equivalent black body temperature for a radiometer
!!   channel, Met Office technical report (see docs/ directory)
!!
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     thermal          flag to indicate channels with thermal emission
!! @param[in,out] rad              radiance structure
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
SUBROUTINE rttov_calcbt(chanprof, coef, thermal, rad)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_radiance
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_po
  USE parkind1, ONLY : jprb, jpim
  USE rttov_math_mod, ONLY : inv_planck
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
  LOGICAL(KIND=jplm  ), INTENT(IN)    :: thermal(SIZE(chanprof))
  TYPE(rttov_radiance), INTENT(INOUT) :: rad
!INTF_END
  REAL   (KIND=jprb) :: bco(SIZE(coef%ff_bco)), bcs(SIZE(coef%ff_bcs))
  REAL   (KIND=jprb) :: tstore1, tstore2, radtotal, radclear
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: chan, i, pol_id
!- End of header ------------------------------------------------------
  nchanprof = SIZE(chanprof)

  IF (coef%ff_val_bc) THEN
    bco = coef%ff_bco
    bcs = coef%ff_bcs
  ELSE
    bco = 0._jprb
    bcs = 1._jprb
  ENDIF

  IF (coef%id_sensor /= sensor_id_po) THEN! All no polarimetric sensors
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        chan            = chanprof(i)%chan
        CALL inv_planck(coef%planck1(chan), coef%planck2(chan), rad%total(i), tstore1)
        CALL inv_planck(coef%planck1(chan), coef%planck2(chan), rad%clear(i), tstore2)
        rad%bt(i)       = (tstore1 - bco(chan)) / bcs(chan)
        rad%bt_clear(i) = (tstore2 - bco(chan)) / bcs(chan)
      ENDIF
    ENDDO
  ELSE! Special case for polarimetric radiometers
! Add average of 1st two elements of Stokes vector to differences in
! 3rd/4th before conversion to BT and subtract after conversion
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        chan   = chanprof(i)%chan
        pol_id = coef%fastem_polar(chan) + 1_jpim
        IF (pol_id == 6_jpim) THEN! 3rd element
          radtotal        = rad%total(i) + 0.5_jprb * (rad%total(i - 2) + rad%total(i - 1))
          radclear        = rad%clear(i) + 0.5_jprb * (rad%clear(i - 2) + rad%clear(i - 1))
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), radtotal, tstore1)
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), radclear, tstore2)
          rad%bt(i)       = (tstore1 - bco(chan)) / bcs(chan)
          rad%bt_clear(i) = (tstore2 - bco(chan)) / bcs(chan)
          rad%bt(i)       = rad%bt(i) - 0.5_jprb * (rad%bt(i - 2) + rad%bt(i - 1))
          rad%bt_clear(i) = rad%bt_clear(i) - 0.5_jprb * (rad%bt_clear(i - 2) + rad%bt_clear(i - 1))
        ELSE IF (pol_id == 7_jpim) THEN! 4th element
          radtotal        = rad%total(i) + 0.5_jprb * (rad%total(i - 3) + rad%total(i - 2))
          radclear        = rad%clear(i) + 0.5_jprb * (rad%clear(i - 3) + rad%clear(i - 2))
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), radtotal, tstore1)
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), radclear, tstore2)
          rad%bt(i)       = (tstore1 - bco(chan)) / bcs(chan)
          rad%bt_clear(i) = (tstore2 - bco(chan)) / bcs(chan)
          rad%bt(i)       = rad%bt(i) - 0.5_jprb * (rad%bt(i - 3) + rad%bt(i - 2))
          rad%bt_clear(i) = rad%bt_clear(i) - 0.5_jprb * (rad%bt_clear(i - 3) + rad%bt_clear(i - 2))
        ELSE! 1st or 2nd element
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), rad%total(i), tstore1)
          CALL inv_planck(coef%planck1(chan), coef%planck2(chan), rad%clear(i), tstore2)
          rad%bt(i)       = (tstore1 - bco(chan)) / bcs(chan)
          rad%bt_clear(i) = (tstore2 - bco(chan)) / bcs(chan)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE rttov_calcbt
