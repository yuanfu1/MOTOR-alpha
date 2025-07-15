! Description:
!> @file
!!   Calculates level-to-space optical profiles for each channel by
!!   applying the optical depth regression for v7 and v8 predictors.
!!
!> @brief
!!   Calculates level-to-space optical profiles for each channel by
!!   applying the optical depth regression for v7 and v8 predictors.
!!
!! @details
!!   This subroutine calculates optical depths for v7 and v8 predictors
!!   including mixed gases, water vapour, water vapour continuum,
!!   O3 and CO2.
!!
!!   Predictors (calculated previously) are multiplied by coefficients
!!   to obtain layer optical depths which are summed to obtain the
!!   profiles.
!!
!!   This subroutine operates on coefficient layers/levels.
!!
!!   NB The output optical depths are NEGATIVE numbers.
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     predictors      pre-calculated predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients
!! @param[in,out] opdp_path       calculated optical depth profiles
!! @param[out]    opdp_ref        reference layer optical depths (may include positive values)
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
SUBROUTINE rttov_opdep_78( &
            nlayers,    &
            chanprof,   &
            predictors, &
            coef,       &
            fast_coef,  &
            opdp_path,  &
            opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,   &
        rttov_coef,       &
        rttov_fast_coef,  &
        rttov_path_pred,  &
        rttov_opdp_path
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)      , INTENT(IN)         :: nlayers
  TYPE(rttov_chanprof)    , INTENT(IN)         :: chanprof(:)
  TYPE(rttov_path_pred)   , INTENT(IN), TARGET :: predictors(:)
  TYPE(rttov_coef)        , INTENT(IN)         :: coef
  TYPE(rttov_fast_coef)   , INTENT(IN), TARGET :: fast_coef(:)
  TYPE(rttov_opdp_path)   , INTENT(INOUT)      :: opdp_path
  REAL(KIND=jprb)         , INTENT(OUT)        :: opdp_ref(nlayers,SIZE(chanprof))
!INTF_END

  REAL(KIND=jprb) :: t(4), ztemp
  REAL(KIND=jprb) :: od(nlayers) ! layer optical depths
  REAL(KIND=jprb), POINTER :: p(:,:), c(:,:)
  INTEGER(KIND=jpim) :: lev, lay, j, chan, prof, nlevels, ii, nchanprof

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_78', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

!-----------------------------------------
! Calculate layer gaseous optical depths
!-----------------------------------------

  DO j = 1, nchanprof
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof

    t = 0._jprb
!-----------------------------
! Water vapour
!-----------------------------
    IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN

      p => predictors(prof)%watervapour
      c => fast_coef(chan)%watervapour

      IF (coef%ncwater == 15) THEN ! v7 WV predictors
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +    &
                 c(2, lay) * p(2, lay) +    &
                 c(3, lay) * p(3, lay) +    &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +    &
                 c(6, lay) * p(6, lay) +    &
                 c(7, lay) * p(7, lay) +    &
                 c(8, lay) * p(8, lay)

          t(3) = c(9, lay) * p(9, lay) +    &
                 c(10, lay) * p(10, lay) +  &
                 c(11, lay) * p(11, lay) +  &
                 c(12, lay) * p(12, lay)

          t(4) = c(13, lay) * p(13, lay) +  &
                 c(14, lay) * p(14, lay) +  &
                 c(15, lay) * p(15, lay)

          od(lay) = SUM(t(1:4))

        ENDDO
      ELSE ! coef%ncwater == 12      v8 WV predictors
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +   &
                 c(2, lay) * p(2, lay) +   &
                 c(3, lay) * p(3, lay) +   &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +   &
                 c(6, lay) * p(6, lay) +   &
                 c(7, lay) * p(7, lay) +   &
                 c(8, lay) * p(8, lay)

          t(3) = c(9, lay)  * p(9, lay)  + &
                 c(10, lay) * p(10, lay) + &
                 c(11, lay) * p(11, lay) + &
                 c(12, lay) * p(12, lay)

          od(lay) = SUM(t(1:3))
        ENDDO
      ENDIF
    ELSE
      od = 0._jprb
    ENDIF

!--------------------------
! Mixed gases
!--------------------------
    IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
      p => predictors(prof)%mixedgas
      c => fast_coef(chan)%mixedgas

      IF (coef%nmixed == 10) THEN
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +     &
                 c(2, lay) * p(2, lay) +     &
                 c(3, lay) * p(3, lay) +     &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +     &
                 c(6, lay) * p(6, lay) +     &
                 c(7, lay) * p(7, lay) +     &
                 c(8, lay) * p(8, lay)

          od(lay) = od(lay) +                &
                    t(1) + t(2) +            &
                    c(9, lay)  * p(9, lay) + &
                    c(10, lay) * p(10, lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          od(lay) = od(lay) + DOT_PRODUCT(c(:,lay), p(:,lay))
        ENDDO
      ENDIF
    ENDIF

!-------------
! Ozone
!-------------
    IF (coef%nozone > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%ozone)) THEN
        p => predictors(prof)%ozone
        c => fast_coef(chan)%ozone

        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +      &
                 c(2, lay) * p(2, lay) +      &
                 c(3, lay) * p(3, lay) +      &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +      &
                 c(6, lay) * p(6, lay) +      &
                 c(7, lay) * p(7, lay) +      &
                 c(8, lay) * p(8, lay)

          od(lay) = od(lay) + t(1) + t(2) +   &
                    c(9, lay) * p(9, lay)   + &
                    c(10, lay) * p(10, lay) + &
                    c(11, lay) * p(11, lay)
        ENDDO
      ENDIF
    ENDIF

!------------------------------
! Water Vapour Continuum (v8 predictors)
!------------------------------
    IF (coef%nwvcont > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN
        p => predictors(prof)%wvcont
        c => fast_coef(chan)%wvcont

        DO lay = 1, nlayers
          od(lay) = od(lay) + &
            DOT_PRODUCT(c(:, lay), p(:, lay))
        ENDDO
      ENDIF
    ENDIF

!-----------
! CO2 (v8 predictors)
!-----------
    IF (coef%nco2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co2)) THEN
        p => predictors(prof)%co2
        c => fast_coef(chan)%co2

        DO lay = 1, nlayers
          od(lay) = od(lay) + &
            DOT_PRODUCT(c(:, lay), p(:, lay))
        ENDDO
      ENDIF
    ENDIF

!------------------------------------------
! Pressure-modulated cell (pmc) sensors
!------------------------------------------
    IF (coef%pmc_shift) THEN
      DO ii = 1, coef%pmc_nvar
        DO lay = 1, coef%pmc_nlay
          od(lay) = od(lay) + &
            coef%pmc_coef(lay, chan, ii) * predictors(prof)%pmc(ii, lay, chan)
        ENDDO
      ENDDO
    ENDIF

!----------------------------------------
! Calculate level-to-space optical depths
!----------------------------------------
    opdp_ref(:,j) = od(:)
    opdp_path%atm_level(1, j) = 0.0_jprb
    ! Introduce ztemp to stop store followed by load in this recursive loop (DJS)
    ztemp = 0.0_jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ! Note that optical depth in the calculations is negative
      IF (od(lay) < 0._jprb) THEN
        ztemp = ztemp + od(lay)
        opdp_path%atm_level(lev, j) = ztemp
      ELSE
        opdp_path%atm_level(lev, j) = opdp_path%atm_level(lev-1, j)
      ENDIF
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_78', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_78
