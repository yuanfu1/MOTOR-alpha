! Description:
!> @file
!!   Calculates level-to-space optical profiles for each channel by
!!   applying the optical depth regression for v9 predictors.
!!
!> @brief
!!   Calculates level-to-space optical profiles for each channel by
!!   applying the optical depth regression for v9 predictors.
!!
!! @details
!!   This subroutine calculates optical depths for v9 predictors
!!   including mixed gases, water vapour, water vapour continuum,
!!   O3, CO2, CO, CH4, N2O and SO2.
!!
!!   Predictors (calculated previously) are multiplied by coefficients
!!   to obtain layer optical depths which are summed to obtain the
!!   profiles. The optical depths for thermal emission calculations
!!   are separate from the optical depths for solar calculations
!!   so this subroutine is called twice where both sources of
!!   radiation are included in the simulations.
!!
!!   This subroutine operates on coefficient layers/levels.
!!
!!   NB The output optical depths are NEGATIVE numbers.
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in]     predictors      pre-calculated predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients (thermal or solar)
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
SUBROUTINE rttov_opdep_9( &
             nlayers,    &
             chanprof,   &
             chanflag,   &
             predictors, &
             coef,       &
             fast_coef,  &
             opdp_path,  &
             opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,   &
        rttov_coef,       &
        rttov_fast_coef,  &
        rttov_path_pred

  USE parkind1, ONLY : jpim, jprb, jplm
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

  INTEGER(KIND=jpim)   , INTENT(IN)         :: nlayers
  TYPE(rttov_chanprof) , INTENT(IN)         :: chanprof(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)         :: chanflag(SIZE(chanprof))
  TYPE(rttov_path_pred), INTENT(IN), TARGET :: predictors(:)
  TYPE(rttov_coef)     , INTENT(IN)         :: coef
  TYPE(rttov_fast_coef), INTENT(IN), TARGET :: fast_coef(:)
  REAL(KIND=jprb)      , INTENT(INOUT)      :: opdp_path(:,:)
  REAL(KIND=jprb)      , INTENT(INOUT)      :: opdp_ref(nlayers,SIZE(chanprof))
!INTF_END

  REAL(KIND=jprb), POINTER :: p(:,:), c(:,:)
  REAL(KIND=jprb)          :: od(nlayers)
  REAL(KIND=jprb)          :: chanx
  REAL(KIND=jprb)          :: t(4), ztemp
  INTEGER(KIND=jpim)       :: lev, lay, chan, j, prof, nlevels, nchanprof
  REAL   (KIND=jprb)       :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

!-----------------------------------------
! Calculate layer gaseous optical depths
!-----------------------------------------

! Summary of which predictors are used (non-solar or solar):
!            WV < 2000.00 Non-solar predictors always used
! 2000.00 <= WV < 2250.00 Solar predictors always used EXCEPT for H2O where it depends on the value of ss_val_chn
! 2250.00 <= WV < 2380.25 Non-solar predictors always used
! 2380.25 <= WV           Solar predictors always used

!--------------------
! Water vapour
!--------------------
  DO j = 1, nchanprof

    t = 0._jprb

    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof

      IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
        p => predictors(prof)%watervapour
        c => fast_coef(chan)%watervapour

        IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) + &
                   c(2, lay) * p(2, lay) + &
                   c(3, lay) * p(3, lay) + &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) + &
                   c(6, lay) * p(6, lay) + &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay)   + &
                   c(10, lay) * p(10, lay) + &
                   c(11, lay) * p(11, lay) + &
                   c(12, lay) * p(12, lay)

            od(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p(13, lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p => predictors(prof)%ch4
            od(:) = od(:) + p(1, :) * (c(14, :) +  c(15, :) * p(3, :))
          ENDIF

        ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay) +       &
                   c(10, lay) * p(10, lay) +     &
                   c(11, lay) * p(11, lay) +     &
                   c(12, lay) * p(12, lay)

          od(lay) = &
            t(1) + t(2) + t(3) + &
            c(13, lay) * p(13, lay)
        ENDDO

        IF (coef%nco2 > 0) THEN
          p => predictors(prof)%co2
          od(:) = od(:) + c(14, :) * p(1, :)
        ENDIF

      ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +      &
                 c(2, lay) * p(2, lay) +       &
                 c(3, lay) * p(3, lay) +       &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +       &
                 c(6, lay) * p(6, lay) +       &
                 c(7, lay) * p(7, lay) + &
                 c(8, lay) * p(8, lay)

          t(3) = c(9, lay) * p(9, lay) +       &
                 c(10, lay) * p(10, lay) +     &
                 c(11, lay) * p(11, lay) +     &
                 c(12, lay) * p(12, lay)

          od(lay) = &
                 t(1) + t(2) + t(3) + &
                 c(13, lay) * p(13, lay)
        ENDDO

        IF (coef%nco2 > 0) THEN
          p => predictors(prof)%co2
          od(:) = od(:) + c(14, :) * p(1, :)
        ENDIF
        IF (coef%nco > 0) THEN
          p => predictors(prof)%co
          od(:) = od(:) + p(1, :) * c(15, :)
        ENDIF
      ELSE IF (chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(1, lay) +      &
                 c(2, lay) * p(2, lay) +       &
                 c(3, lay) * p(3, lay) +       &
                 c(4, lay) * p(4, lay)

          t(2) = c(5, lay) * p(5, lay) +       &
                 c(6, lay) * p(6, lay) +       &
                 c(7, lay) * p(7, lay) + &
                 c(8, lay) * p(8, lay)

          t(3) = c(9, lay) * p(9, lay) +       &
                 c(10, lay) * p(10, lay) +     &
                 c(11, lay) * p(11, lay) +     &
                 c(12, lay) * p(12, lay)

          od(lay) = &
                 t(1) + t(2) + t(3) + &
                 c(13, lay) * p(13, lay)
        ENDDO
      ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
        IF (coef%ss_val_chn(chan) > 0) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(14, lay)

            t(3) = c(9, lay) * p(15, lay) +       &
                   c(10, lay) * p(10, lay) +      &
                   c(11, lay) * p(11, lay) +      &
                   c(12, lay) * p(16, lay)

            od(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p(17, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p => predictors(prof)%co2
            od(:) = od(:) + p(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p => predictors(prof)%co
            od(:) = od(:) + p(1, :) * c(15, :)
          ENDIF
        ELSE
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay) +       &
                   c(10, lay) * p(10, lay) +     &
                   c(11, lay) * p(11, lay) +     &
                   c(12, lay) * p(12, lay)

            od(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p(13, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p => predictors(prof)%co2
            od(:) = od(:) + p(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p => predictors(prof)%co
            od(:) = od(:) + p(1, :) * c(15, :)
          ENDIF
        ENDIF
      ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay) +       &
                   c(10, lay) * p(10, lay) +     &
                   c(11, lay) * p(11, lay) +     &
                   c(12, lay) * p(12, lay)

            od(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p(13, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p => predictors(prof)%co2
            od(:) = od(:) + p(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p => predictors(prof)%co
            od(:) = od(:) + p(1, :) * c(15, :)
          ENDIF
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay) +       &
                   c(10, lay) * p(10, lay) +     &
                   c(11, lay) * p(11, lay) +     &
                   c(12, lay) * p(12, lay)

            od(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p(13, lay)
          ENDDO
        ELSE   !chanx >= rttov9_wv2380_25
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay) +      &
                   c(2, lay) * p(2, lay) +       &
                   c(3, lay) * p(3, lay) +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay) +       &
                   c(6, lay) * p(6, lay) +       &
                   c(7, lay) * p(7, lay) + &
                   c(8, lay) * p(14, lay)

            t(3) = c(9, lay) * p(15, lay) +       &
                   c(10, lay) * p(18, lay) +      & ! not p10
                   c(11, lay) * p(11, lay) +      &
                   c(12, lay) * p(16, lay)

            od(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p(19, lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p => predictors(prof)%ch4
            od(:) = od(:) + &
              SQRT(p(2, :)) * (c(14, :) * p(1, :) + c(15, :) * p(3, :))
          ENDIF
        ENDIF
      ELSE
        od(:) = 0._jprb
      ENDIF

!--------------------------
! Mixed gases
!--------------------------
      IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
        p => predictors(prof)%mixedgas
        c => fast_coef(chan)%mixedgas

        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
          (chanx < rttov9_wv2000_00)) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay)  +      &
                   c(2, lay) * p(2, lay)  +       &
                   c(3, lay) * p(3, lay)  +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay)  +      &
                   c(6, lay) * p(6, lay)  +       &
                   c(7, lay) * p(7, lay)  +       &
                   c(8, lay) * p(8, lay)

            od(lay) = od(lay) + t(1) + t(2)
          ENDDO
        ELSE
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay)  +      &
                   c(2, lay) * p(2, lay)  +       &
                   c(3, lay) * p(3, lay)  +       &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay)  + &
                   c(6, lay) * p(6, lay)  +      &
                   c(7, lay) * p(7, lay)  +       &
                   c(8, lay) * p(8, lay)

            od(lay) = od(lay) + t(1) + t(2) +  &
                      c(9, lay) * p(9, lay)  +       &
                      c(10, lay) * p(10, lay)
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

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  +      &
                     c(2, lay) * p(2, lay)  +       &
                     c(3, lay) * p(3, lay)  +       &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  +      &
                     c(6, lay) * p(6, lay)  +       &
                     c(7, lay) * p(7, lay)  +       &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  +       &
                        c(10, lay) * p(10, lay)  +     &
                        c(11, lay) * p(11, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  +      &
                     c(2, lay) * p(2, lay)  +       &
                     c(3, lay) * p(12, lay)  +       &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  +      &
                     c(6, lay) * p(6, lay)  +       &
                     c(7, lay) * p(7, lay)  +       &
                     c(8, lay) * p(13, lay)

              t(3) = c(9, lay) * p(9, lay)  +       &
                     c(10, lay) * p(10, lay)  +     &
                     c(11, lay) * p(11, lay)  +     &
                     c(12, lay) * p(14, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3) + &
                        c(13, lay) * p(15, lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!------------------------------
! Water Vapour Continuum
!------------------------------
      IF (coef%nwvcont > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN !*
          p => predictors(prof)%wvcont
          c => fast_coef(chan)%wvcont

          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay)  +      &
                   c(2, lay) * p(2, lay)  +      &
                   c(3, lay) * p(3, lay)  +      &
                   c(4, lay) * p(4, lay)
            od(lay) = od(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF

!-----------
! CO2
!-----------
      IF (coef%nco2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%co2)) THEN !*
          p => predictors(prof)%co2
          c => fast_coef(chan)%co2

          IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(5, lay)  * p(6, lay)  ** 2 + &
                     c(4, lay) * p(3, lay)  ** 2

              t(2) = c(5, lay) * p(5, lay)  +       &
                     c(6, lay) * p(6, lay)  ** 2 +  &
                     c(7, lay) * p(7, lay)  ** 2 +  &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) +       &
                        c(9, lay) * p(9, lay)
            ENDDO
          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                c(2, lay) * p(2, lay)  + &
                c(3, lay) * p(3, lay)  + &
                c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(11, lay) * p(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + c(10, :) * p(1, :)
            ENDIF

          ELSE IF (chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)
            ENDDO

          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(7, lay)  + &
                     c(7, lay) * p(8, lay) + &
                     c(9, lay) * p(10, lay)

              t(3) = c(10, lay) * p(11, lay)  +           &
                     c(11, lay) * p(12, lay)  +           &
                     c(12, lay) * p(13, lay)  +           &
                     c(13, lay) * p(14, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3) + &
                        c(14, lay) * p(15, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + c(8, :) * p(1, :)
            ENDIF

          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay) + &
                        c(11, lay) * p(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + c(10, :) * p(1, :)
            ENDIF

          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay) + &
                        c(10, lay) * p(10, lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(7, lay)  + &
                     c(7, lay) * p(8, lay) + &
                     c(8, lay) * p(10, lay)

              t(3) = c(9, lay) * p(11, lay)  +           &
                     c(10, lay) * p(12, lay)  +           &
                     c(11, lay) * p(13, lay)  +           &
                     c(12, lay) * p(14, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3) + &
                        c(13, lay) * p(15, lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! N2O
!-----------
      IF (coef%nn2o > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%n2o)) THEN !*
          p => predictors(prof)%n2o
          c => fast_coef(chan)%n2o

          IF (chanx >= rttov9_wv1050_00 .AND. chanx < rttov9_wv1350_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay)
            ENDDO

            IF (coef%nch4 > 0) THEN
              p => predictors(prof)%ch4
              od(:) = od(:) + &
                c(11, :) * p(1, :) + &
                c(12, :) * p(10, :)
              ENDIF

          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + p(1, :) * ( &
                      c(11, :) +      &
                      c(12, :) * p(14, :))
            ENDIF

          ELSE IF (chanx < rttov9_wv2000_00) THEN

            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay)
            ENDDO

          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              t(3) = c(9, lay) * p(11, lay)  + &
                     c(10, lay) * p(10, lay) + &
                     c(13, lay) * p(12, lay) + &
                     c(14, lay) * p(13, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + p(1, :) * ( &
                      c(11, :) +      &
                      c(12, :) * p(14, :))
            ENDIF

          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN

            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p => predictors(prof)%co
              od(:) = od(:) + p(1, :) * ( &
                      c(11, :) +      &
                      c(12, :) * p(14, :) )
            ENDIF

          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              t(3) = c(9, lay) * p(11, lay)  + &
                     c(10, lay) * p(10, lay) + &
                     c(11, lay) * p(12, lay) + &
                     c(12, lay) * p(13, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CO
!-----------
      IF (coef%nco > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%co)) THEN !*
          p => predictors(prof)%co
          c => fast_coef(chan)%co

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              od(lay) = od(lay) + t(1) + t(2) + &
                        c(9, lay) * p(9, lay)  + &
                        c(10, lay) * p(10, lay) + &
                        c(11, lay) * p(11, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay)  + &
                     c(2, lay) * p(2, lay)  + &
                     c(3, lay) * p(3, lay)  + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay)  + &
                     c(6, lay) * p(6, lay)  + &
                     c(7, lay) * p(7, lay)  + &
                     c(8, lay) * p(8, lay)

              t(3) = c(9, lay) * p(9, lay)  + &
                     c(10, lay) * p(10, lay) + &
                     c(11, lay) * p(12, lay) + &
                     c(12, lay) * p(13, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CH4
!-----------
      IF (coef%nch4 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%ch4)) THEN !*
          p => predictors(prof)%ch4
          c => fast_coef(chan)%ch4

          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(1, lay)  + &
                   c(2, lay) * p(2, lay)  + &
                   c(3, lay) * p(3, lay)  + &
                   c(4, lay) * p(4, lay)

            t(2) = c(5, lay) * p(5, lay)  + &
                   c(6, lay) * p(6, lay)  + &
                   c(7, lay) * p(7, lay)  + &
                   c(8, lay) * p(8, lay)

            t(3) = c(9, lay) * p(9, lay)  + &
                   c(10, lay) * p(10, lay) + &
                   c(11, lay) * p(11, lay)

            od(lay) = od(lay) + t(1) + t(2) + t(3)
          ENDDO
        ENDIF
      ENDIF

!-----------
! SO2
!-----------
      IF (coef%nso2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%so2)) THEN !*
          p => predictors(prof)%so2
          c => fast_coef(chan)%so2

          IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay) + &
                     c(2, lay) * p(2, lay) + &
                     c(3, lay) * p(3, lay) + &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) + &
                     c(6, lay) * p(6, lay) + &
                     c(7, lay) * p(7, lay) + &
                     c(8, lay) * p(8, lay)

              t(3) = c(9, lay) * p(9, lay)   + &
                     c(10, lay) * p(10, lay) + &
                     c(11, lay) * p(11, lay) + &
                     c(12, lay) * p(12, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3) + &
                        c(13, lay) * p(13, lay)
            ENDDO
          ELSE   !use same predictors as WV >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p(1, lay) +      &
                     c(2, lay) * p(2, lay) +       &
                     c(3, lay) * p(3, lay) +       &
                     c(4, lay) * p(4, lay)

              t(2) = c(5, lay) * p(5, lay) +       &
                     c(6, lay) * p(6, lay) +       &
                     c(7, lay) * p(7, lay) + &
                     c(8, lay) * p(14, lay)

              t(3) = c(9, lay) * p(15, lay) +       &
                     c(10, lay) * p(18, lay) +      & ! not p10
                     c(11, lay) * p(11, lay) +      &
                     c(12, lay) * p(16, lay)

              od(lay) = od(lay) + t(1) + t(2) + t(3) + &
                        c(13, lay) * p(19, lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!----------------------------------------
! Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
! single layer
! store (neg) optical depth in reference array on COEF levs ...
! ... for TL, AD and K calculations

! level to space optical depths
! these are path (i.e. cumulated layer) optical depths
      opdp_ref(:,j) = od(:)
      opdp_path(1, j) = 0.0_jprb
!  Introduce ztemp to stop store followed by load in this recursive loop (DJS)
      ztemp = 0.0_jprb
      DO lev = 2, nlevels
        lay = lev - 1
        IF (od(lay) < 0._jprb) THEN
          ztemp = ztemp + od(lay)
          opdp_path(lev, j) = ztemp
        ELSE
          opdp_path(lev, j) = opdp_path(lev - 1, j)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9', 1_jpim, ZHOOK_HANDLE)

END SUBROUTINE rttov_opdep_9
