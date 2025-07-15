! Description:
!> @file
!!   TL of brightness temperature calculation
!
!> @brief
!!   TL of brightness temperature calculation
!!
!! @details
!!  The derivative of the inverse Planck function with respect to radiance is
!!
!!                             C1 * C2 * Nu**4
!! B-1'(R,Nu) = --------------------------------------------- dR
!!                     (    C1 * Nu**3) (     C1 * Nu**3 )**2
!!               R**2 *(1 + ----------) ( Ln( ---------- )
!!                     (         R    ) (        R       )
!!
!! which can be reduced to the following, with
!!  C1 = C1 * Nu**3
!!  C2 = C2 * Nu
!!
!!                  C1 * B-1(R,Nu)**2
!! B-1'(R,Nu) = ----------------------- dR
!!                  C2 * R * (R + C1)
!!
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     thermal          flag to indicate channels with thermal emission
!! @param[in]     rad              radiance structure
!! @param[in,out] rad_tl           radiance and BT perturbations
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
SUBROUTINE rttov_calcbt_tl(chanprof, coef, thermal, rad, rad_tl)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_radiance
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_po
  USE parkind1, ONLY : jprb, jpim
  USE rttov_math_mod, ONLY : inv_planck_tl
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
  LOGICAL(KIND=jplm  ), INTENT(IN)    :: thermal(SIZE(chanprof))
  TYPE(rttov_radiance), INTENT(IN)    :: rad
  TYPE(rttov_radiance), INTENT(INOUT) :: rad_tl
!INTF_END
  REAL   (KIND=jprb) :: bco(SIZE(coef%ff_bco)), bcs(SIZE(coef%ff_bcs))
  REAL   (KIND=jprb) :: tstar    , tstar1   , tstar2
  REAL   (KIND=jprb) :: tstar_tl , tstar1_tl, tstar2_tl
  REAL   (KIND=jprb) :: radtotal , radclear
  INTEGER(KIND=jpim) :: chan     , i, pol_id
  INTEGER(KIND=jpim) :: nchanprof
!- End of header --------------------------------------------------------
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
        chan = chanprof(i)%chan
! Clear+cloudy radiance
! T star for direct model
        tstar              = bco(chan) + bcs(chan) * rad%bt(i)
! TL
        CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), rad%total(i), rad_tl%total(i), tstar, tstar_tl)
        rad_tl%bt(i)       = tstar_tl / bcs(chan)
!clear radiance
! T star for direct model
        tstar              = bco(chan) + bcs(chan) * rad%bt_clear(i)
! TL
        CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), rad%clear(i), rad_tl%clear(i), tstar, tstar_tl)
        rad_tl%bt_clear(i) = tstar_tl / bcs(chan)
      ENDIF
    ENDDO
  ELSE! Special case for polarimetric radiometers
! Add average of 1st two elements of Stokes vector to differences in
! 3rd/4th before conversion to BT and subtract after conversion
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        chan   = chanprof(i)%chan
        pol_id = coef%fastem_polar(chan) + 1_jpim
        IF (pol_id == 6_jpim) THEN! 3rd stoke element
          tstar1             = bco(chan) + bcs(chan) * rad%bt(i)
          tstar1             = tstar1 + bcs(chan) * 0.5_jprb * (rad%bt(i - 2) + rad%bt(i - 1))
          tstar2             = bco(chan) + bcs(chan) * rad%bt_clear(i)
          tstar2             = tstar2 + bcs(chan) * 0.5_jprb * (rad%bt_clear(i - 2) + rad%bt_clear(i - 1))
          radtotal           = rad%total(i) + 0.5 * (rad%total(i - 2) + rad%total(i - 1))
          radclear           = rad%clear(i) + 0.5 * (rad%clear(i - 2) + rad%clear(i - 1))
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), radtotal, rad_tl%total(i), tstar1, tstar1_tl)
          rad_tl%bt(i)       = tstar1_tl / bcs(chan)
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), radclear, rad_tl%clear(i), tstar2, tstar2_tl)
          rad_tl%bt_clear(i) = tstar2_tl / bcs(chan)
        ELSE IF (pol_id == 7_jpim) THEN! 4th stoke element
          tstar1             = bco(chan) + bcs(chan) * rad%bt(i)
          tstar1             = tstar1 + bcs(chan) * 0.5_jprb * (rad%bt(i - 3) + rad%bt(i - 2))
          tstar2             = bco(chan) + bcs(chan) * rad%bt_clear(i)
          tstar2             = tstar2 + bcs(chan) * 0.5_jprb * (rad%bt_clear(i - 3) + rad%bt_clear(i - 2))
          radtotal           = rad%total(i) + 0.5 * (rad%total(i - 3) + rad%total(i - 2))
          radclear           = rad%clear(i) + 0.5 * (rad%clear(i - 3) + rad%clear(i - 2))
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), radtotal, rad_tl%total(i), tstar1, tstar1_tl)
          rad_tl%bt(i)       = tstar1_tl / bcs(chan)
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), radclear, rad_tl%clear(i), tstar2, tstar2_tl)
          rad_tl%bt_clear(i) = tstar2_tl / bcs(chan)
        ELSE! 1st/2nd stoke elements
! Clear+cloudy radiance
! T star for direct model
          tstar              = bco(chan) + bcs(chan) * rad%bt(i)
! TL
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), rad%total(i), rad_tl%total(i), tstar, tstar_tl)
          rad_tl%bt(i)       = tstar_tl / bcs(chan)
!clear radiance
! T star for direct model
          tstar              = bco(chan) + bcs(chan) * rad%bt_clear(i)
! TL
          CALL inv_planck_tl(coef%planck1(chan), coef%planck2(chan), rad%clear(i), rad_tl%clear(i), tstar, tstar_tl)
          rad_tl%bt_clear(i) = tstar_tl / bcs(chan)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE rttov_calcbt_tl
