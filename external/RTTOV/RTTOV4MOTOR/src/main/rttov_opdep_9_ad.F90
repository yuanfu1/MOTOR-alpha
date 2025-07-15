! Description:
!> @file
!!   AD of optical depth calculation (v9 predictors).
!!
!> @brief
!!   AD of optical depth calculation (v9 predictors).
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in]     predictors      pre-calculated predictors
!! @param[in,out] predictors_ad   increments in predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients (thermal or solar)
!! @param[in,out] opdp_path_ad    optical depth profile increments
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
SUBROUTINE rttov_opdep_9_ad( &
             nlayers,       &
             chanprof,      &
             chanflag,      &
             predictors,    &
             predictors_ad, &
             coef,          &
             fast_coef,     &
             opdp_path_ad,  &
             opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,  &
        rttov_coef,      &
        rttov_fast_coef, &
        rttov_path_pred

  USE parkind1, ONLY : jprb, jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
        rttov9_wv0690_50, &
        rttov9_wv1050_00, &
        rttov9_wv1095_25, &
        rttov9_wv1100_25, &
        rttov9_wv1350_25, &
        rttov9_wv1750_25, &
        rttov9_wv1900_25, &
        rttov9_wv1995_00, &
        rttov9_wv2000_00, &
        rttov9_wv2250_00, &
        rttov9_wv2295_25, &
        rttov9_wv2380_25

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof) , INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors(:)
  TYPE(rttov_path_pred), INTENT(INOUT) :: predictors_ad(:)
  TYPE(rttov_coef)     , INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef(:)
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_path_ad(:,:)
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers,SIZE(chanprof))
!INTF_END

  REAL(KIND=jprb), POINTER :: c(:,:), p_ad(:,:), p(:,:)
  REAL(KIND=jprb)          :: od_ad(nlayers)
  REAL(KIND=jprb)          :: ztemp
  REAL(KIND=jprb)          :: chanx
  INTEGER(KIND=jpim)       :: lev, lay, chan, j, nlevels
  INTEGER(KIND=jpim)       :: prof
  INTEGER(KIND=jpim)       :: nchanprof
  INTEGER(KIND=jpim)       :: map(13)
  REAL   (KIND=jprb)       :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_AD', 0_jpim, ZHOOK_HANDLE)

!----------------------------------------
! Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  nchanprof            = SIZE(chanprof)
  nlevels              = nlayers + 1

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      chanx = coef%ff_cwn(chan)

      ztemp = opdp_path_ad(nlevels, j)
      DO lev = nlevels, 2,  - 1
        lay = lev - 1
        IF (opdp_ref(lay,j) < 0._jprb) THEN
          od_ad(lay) = ztemp
          ztemp = ztemp + opdp_path_ad(lev - 1, j)
        ELSE
          ztemp = ztemp + opdp_path_ad(lev - 1, j)
          od_ad(lay) = 0.0_jprb
        ENDIF
      ENDDO
      opdp_path_ad(1, j) = 0.0_jprb

!-----------
! SO2
!-----------
  IF (coef%nso2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%so2)) THEN !*
          p_ad => predictors_ad(prof)%so2
          c => fast_coef(chan)%so2

          IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
            DO lay = 1, nlayers
              p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
            ENDDO
          ELSE   !use same predictors as WV >= rttov9_wv2380_25
            map(1:6) = (/14,15,18,11,16,19/)
            DO lay = 1, nlayers
              p_ad(1:7, lay) = p_ad(1:7, lay) + c(1:7, lay) * od_ad(lay)
              p_ad(map(1:6), lay) = p_ad(map(1:6), lay) + c(8:13, lay) * od_ad(lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CH4
!-----------
  IF (coef%nch4 > 0) THEN
    IF (ASSOCIATED(fast_coef(chan)%ch4)) THEN
      p_ad => predictors_ad(prof)%ch4
      c => fast_coef(chan)%ch4

      DO lay = 1, nlayers
        p_ad(1:11, lay) = p_ad(1:11, lay) + c(1:11, lay) * od_ad(lay)
      ENDDO

    ENDIF
  ENDIF

!-----------
! CO
!-----------
  IF (coef%nco > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%co)) THEN
          p_ad => predictors_ad(prof)%co
          c => fast_coef(chan)%co

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              p_ad(1:11, lay) = p_ad(1:11, lay) + c(1:11, lay) * od_ad(lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
              p_ad(12:13, lay) = p_ad(12:13, lay) + c(11:12, lay) * od_ad(lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! N2O
!-----------
  IF (coef%nn2o > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%n2o)) THEN
          p_ad => predictors_ad(prof)%n2o
          c => fast_coef(chan)%n2o

          IF (chanx >= rttov9_wv1050_00 .AND. chanx < rttov9_wv1350_25) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO
            IF (coef%nch4 > 0) THEN
              p_ad => predictors_ad(prof)%ch4
              p_ad(1, :) = p_ad(1, :) + c(11, :) * od_ad(:)
              p_ad(10, :) = p_ad(10, :) + c(12, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co; p => predictors(prof)%co
              p_ad(1, :) = p_ad(1, :) + (c(11, :) + c(12, :) * p(14, :)) * od_ad(:)
              p_ad(14, :) = p_ad(14, :) + p(1, :) * c(12, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO
          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              p_ad(1:8, lay) = p_ad(1:8, lay) + c(1:8, lay) * od_ad(lay)
              p_ad(10:13, lay) = p_ad(10:13, lay) + c((/10,9,13,14/), lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co; p => predictors(prof)%co
              p_ad(1, :) = p_ad(1, :) + (c(11, :) + c(12, :) * p(14, :)) * od_ad(:)
              p_ad(14, :) = p_ad(14, :) + p(1, :) * c(12, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co; p => predictors(prof)%co
              p_ad(1, :) = p_ad(1, :) + (c(11, :) + c(12, :) * p(14, :)) * od_ad(:)
              p_ad(14, :) = p_ad(14, :) + p(1, :) * c(12, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              p_ad(1:8, lay) = p_ad(1:8, lay) + c(1:8, lay) * od_ad(lay)
              p_ad(10:13, lay) = p_ad(10:13, lay) + c((/10,9,11,12/), lay) * od_ad(lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CO2
!-----------
  IF (coef%nco2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%co2)) THEN
          p_ad => predictors_ad(prof)%co2
          p => predictors(prof)%co2
          c => fast_coef(chan)%co2

          IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
            map(1:4) = (/1,2,8,9/)
            DO lay = 1, nlayers
              p_ad(map(1:4), lay) = p_ad(map(1:4), lay) + &
                                      c(map(1:4), lay) * od_ad(lay)

              p_ad((/3,7/), lay) = p_ad((/3,7/), lay) + 2._jprb * &
                                 p((/3,7/), lay) * c((/4,7/), lay) * od_ad(lay)

              p_ad(5, lay) = p_ad(5, lay) + (c(5, lay) + p(6, lay) ** 2 * &
                             c(3, lay)) * od_ad(lay)

              p_ad(6, lay) = p_ad(6, lay) + 2._jprb * p(6, lay) * &
                             (c(3, lay) * p(5, lay) + c(6, lay)) * od_ad(lay)

            ENDDO
          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              p_ad(1:9, lay) = p_ad(1:9, lay) + c(1:9, lay) * od_ad(lay)
              p_ad(10, lay) = p_ad(10, lay) + c(11, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co
              p_ad(1, :) = p_ad(1, :) + c(10, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              p_ad(1:9, lay) = p_ad(1:9, lay) + c(1:9, lay) * od_ad(lay)
            ENDDO
          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              p_ad(1:5, lay) = p_ad(1:5, lay) + c(1:5, lay) * od_ad(lay)
              p_ad(7:8,lay) = p_ad(7:8,lay) + c(6:7, lay) * od_ad(lay)
              p_ad(10:15,lay) = p_ad(10:15,lay) + c(9:14, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co
              p_ad(1, :) = p_ad(1, :) + c(8, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
            DO lay = 1, nlayers
              p_ad(1:9, lay) = p_ad(1:9, lay) + c(1:9, lay) * od_ad(lay)
              p_ad(10, lay) = p_ad(10, lay) + c(11, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co
              p_ad(1, :) = p_ad(1, :) + c(10, :) * od_ad(:)
            ENDIF
          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              p_ad(1:5, lay) = p_ad(1:5, lay) + c(1:5, lay) * od_ad(lay)
              p_ad(7:8,lay) = p_ad(7:8,lay) + c(6:7, lay) * od_ad(lay)
              p_ad(10:15,lay) = p_ad(10:15,lay) + c(8:13, lay) * od_ad(lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!------------------------------
! Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN
          p_ad => predictors_ad(prof)%wvcont
          c => fast_coef(chan)%wvcont

          DO lay = 1, nlayers
            p_ad(:, lay) = p_ad(:, lay) + c(:, lay) * od_ad(lay)
          ENDDO
        ENDIF
      ENDIF

!-------------
! Ozone
!-------------
  IF (coef%nozone > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%ozone)) THEN
          p_ad => predictors_ad(prof)%ozone
          c => fast_coef(chan)%ozone

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              p_ad(1:11, lay) = p_ad(1:11, lay) + c(1:11, lay) * od_ad(lay)
            ENDDO
          ELSE
            map = (/1,2,12,4,5,6,7,13,9,10,11,14,15/)
            DO lay = 1, nlayers
              p_ad(map, lay) = p_ad(map, lay) + c(1:13, lay) * od_ad(lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!--------------------
! Water vapour
!--------------------
      IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
          p_ad => predictors_ad(prof)%watervapour
          c => fast_coef(chan)%watervapour

        IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p_ad => predictors_ad(prof)%ch4; p => predictors(prof)%ch4
            p_ad(1, :) = p_ad(1, :) + (c(15, :) * p(3, :) + c(14, :)) * od_ad(:)
            p_ad(3, :) = p_ad(3, :) +  c(15, :) * p(1, :) * od_ad(:)
          ENDIF
        ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p_ad => predictors_ad(prof)%co2
            p_ad(1, :) = p_ad(1, :) + od_ad(:) * c(14, :)
          ENDIF
        ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO

          IF (coef%nco > 0) THEN
            p_ad => predictors_ad(prof)%co
            p_ad(1, :) = p_ad(1, :) + c(15, :) * od_ad(:)
          ENDIF
          IF (coef%nco2 > 0) THEN
            p_ad => predictors_ad(prof)%co2
            p_ad(1, :) = p_ad(1, :) + c(14, :) * od_ad(:)
          ENDIF
        ELSE IF (chanx < rttov9_wv2000_00) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO
        ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
          IF (coef%ss_val_chn(chan) > 0) THEN
            MAP(1:6) = (/14,15,10,11,16,17/)
            DO lay = 1, nlayers
              p_ad(1:7, lay) = p_ad(1:7, lay) + c(1:7, lay) * od_ad(lay)
              p_ad(map(1:6), lay) = p_ad(map(1:6), lay) + &
                                      c(8:13, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co
              p_ad(1, :) = p_ad(1, :) + c(15, :) * od_ad(:)
            ENDIF
            IF (coef%nco2 > 0) THEN
              p_ad => predictors_ad(prof)%co2
              p_ad(1, :) = p_ad(1, :) + c(14, :) * od_ad(:)
            ENDIF
          ELSE
            DO lay = 1, nlayers
              p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_ad => predictors_ad(prof)%co
              p_ad(1, :) = p_ad(1, :) + c(15, :) * od_ad(:)
            ENDIF
            IF (coef%nco2 > 0) THEN
              p_ad => predictors_ad(prof)%co2
              p_ad(1, :) = p_ad(1, :) + c(14, :) * od_ad(:)
            ENDIF
          ENDIF
        ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO

          IF (coef%nco > 0) THEN
            p_ad => predictors_ad(prof)%co
            p_ad(1, :) = p_ad(1, :) + c(15, :) * od_ad(:)
          ENDIF
          IF (coef%nco2 > 0) THEN
            p_ad => predictors_ad(prof)%co2
            p_ad(1, :) = p_ad(1, :) + c(14, :) * od_ad(:)
          ENDIF
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = 1, nlayers
            p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
          ENDDO
        ELSE   !chanx >= rttov9_wv2380_25
          map(1:6) = (/14,15,18,11,16,19/)
          DO lay = 1, nlayers
            p_ad(1:7, lay) = p_ad(1:7, lay) + c(1:7, lay) * od_ad(lay)
            p_ad(map(1:6), lay) = p_ad(map(1:6), lay) + c(8:13, lay) * od_ad(lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p_ad => predictors_ad(prof)%ch4; p => predictors(prof)%ch4

            p_ad(1, :) = p_ad(1, :) + SQRT(p(2, :)) * c(14, :) * od_ad(:)
            p_ad(3, :) = p_ad(3, :) + SQRT(p(2, :)) * c(15, :) * od_ad(:)
            p_ad(2, :) = p_ad(2, :) + 0.5_jprb * &
                         ((c(14, :) * p(1, :) + c(15, :) * p(3, :)) / SQRT(p(2, :))) * &
                         od_ad(:)
          ENDIF
        ENDIF
      ENDIF

!--------------------------
! Mixed gases
!--------------------------
      IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
        p_ad => predictors_ad(prof)%mixedgas
        c => fast_coef(chan)%mixedgas

        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
          (chanx < rttov9_wv2000_00)) THEN
          DO lay = 1, nlayers
            p_ad(1:8, lay) = p_ad(1:8, lay) + c(1:8, lay) * od_ad(lay)
          ENDDO
        ELSE
          DO lay = 1, nlayers
            p_ad(1:10, lay) = p_ad(1:10, lay) + c(1:10, lay) * od_ad(lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_ad
