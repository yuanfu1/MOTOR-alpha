! Description:
!> @file
!!   TL of optical depth calculation (v9 predictors).
!!
!> @brief
!!   TL of optical depth calculation (v9 predictors).
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in]     predictors      pre-calculated predictors
!! @param[in]     predictors_tl   perturbations in predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients (thermal or solar)
!! @param[in,out] opdp_path_tl    calculated optical depth profile perturbations
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
SUBROUTINE rttov_opdep_9_tl( &
             nlayers,       &
             chanprof,      &
             chanflag,      &
             predictors,    &
             predictors_tl, &
             coef,          &
             fast_coef,     &
             opdp_path_tl,  &
             opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,  &
        rttov_coef,      &
        rttov_fast_coef, &
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

  TYPE(rttov_chanprof) , INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors(:)
  TYPE(rttov_path_pred), INTENT(IN)    :: predictors_tl(:)
  TYPE(rttov_coef)     , INTENT(IN)    :: coef
  TYPE(rttov_fast_coef), INTENT(IN)    :: fast_coef(:)
  REAL(KIND=jprb)      , INTENT(INOUT) :: opdp_path_tl(:,:)
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers,SIZE(chanprof))
!INTF_END

  REAL(KIND=jprb), POINTER :: p(:,:), c(:,:), p_tl(:,:)
  REAL(KIND=jprb) :: od_tl(nlayers)
  REAL(KIND=jprb) :: chanx
  INTEGER(KIND=jpim) :: lev, lay, prof, chan, j, nlevels, nchannels
  REAL   (KIND=jprb) :: t(4), ztemp
  REAL   (KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1

!-----------------------------------------
! Calculate layer gaseous optical depths
!-----------------------------------------

!--------------------
! Water vapour
!--------------------
  DO j = 1, nchannels

    t = 0._jprb

    IF (chanflag(j)) THEN
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof

      IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
        p_tl => predictors_tl(prof)%watervapour
        c => fast_coef(chan)%watervapour

        IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) + &
                   c(2, lay) * p_tl(2, lay) + &
                   c(3, lay) * p_tl(3, lay) + &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) + &
                   c(6, lay) * p_tl(6, lay) + &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay)   + &
                   c(10, lay) * p_tl(10, lay) + &
                   c(11, lay) * p_tl(11, lay) + &
                   c(12, lay) * p_tl(12, lay)

            od_tl(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p_tl(13, lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p    => predictors(prof)%ch4
            p_tl => predictors_tl(prof)%ch4
            od_tl(:) = od_tl(:) + &
                       p_tl(1, :) * (c(14, :) +  c(15, :) * p(3, :)) + &
                       p(1, :)    * (c(15, :) * p_tl(3, :))
          ENDIF

        ELSE IF (chanx >= rttov9_wv1900_25 .AND. chanx < rttov9_wv1995_00) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay) +       &
                   c(10, lay) * p_tl(10, lay) +     &
                   c(11, lay) * p_tl(11, lay) +     &
                   c(12, lay) * p_tl(12, lay)

          od_tl(lay) = t(1) + t(2) + t(3) + &
                       c(13, lay) * p_tl(13, lay)
        ENDDO

        IF (coef%nco2 > 0) THEN
          p_tl => predictors_tl(prof)%co2
          od_tl(:) = od_tl(:) + c(14, :) * p_tl(1, :)
        ENDIF

      ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p_tl(1, lay) +      &
                 c(2, lay) * p_tl(2, lay) +       &
                 c(3, lay) * p_tl(3, lay) +       &
                 c(4, lay) * p_tl(4, lay)

          t(2) = c(5, lay) * p_tl(5, lay) +       &
                 c(6, lay) * p_tl(6, lay) +       &
                 c(7, lay) * p_tl(7, lay) + &
                 c(8, lay) * p_tl(8, lay)

          t(3) = c(9, lay) * p_tl(9, lay) +       &
                 c(10, lay) * p_tl(10, lay) +     &
                 c(11, lay) * p_tl(11, lay) +     &
                 c(12, lay) * p_tl(12, lay)

          od_tl(lay) = &
                 t(1) + t(2) + t(3) + &
                 c(13, lay) * p_tl(13, lay)
        ENDDO

        IF (coef%nco2 > 0) THEN
          p_tl => predictors_tl(prof)%co2
          od_tl(:) = od_tl(:) + c(14, :) * p_tl(1, :)
        ENDIF
        IF (coef%nco > 0) THEN
          p_tl => predictors_tl(prof)%co
          od_tl(:) = od_tl(:) + c(15, :) * p_tl(1, :)
        ENDIF
      ELSE IF (chanx < rttov9_wv2000_00) THEN
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p_tl(1, lay) +      &
                 c(2, lay) * p_tl(2, lay) +       &
                 c(3, lay) * p_tl(3, lay) +       &
                 c(4, lay) * p_tl(4, lay)

          t(2) = c(5, lay) * p_tl(5, lay) +       &
                 c(6, lay) * p_tl(6, lay) +       &
                 c(7, lay) * p_tl(7, lay) + &
                 c(8, lay) * p_tl(8, lay)

          t(3) = c(9, lay) * p_tl(9, lay) +       &
                 c(10, lay) * p_tl(10, lay) +     &
                 c(11, lay) * p_tl(11, lay) +     &
                 c(12, lay) * p_tl(12, lay)

          od_tl(lay) = &
                 t(1) + t(2) + t(3) + &
                 c(13, lay) * p_tl(13, lay)
        ENDDO
      ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
        IF (coef%ss_val_chn(chan) > 0) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(14, lay)

            t(3) = c(9, lay) * p_tl(15, lay) +       &
                   c(10, lay) * p_tl(10, lay) +      &
                   c(11, lay) * p_tl(11, lay) +      &
                   c(12, lay) * p_tl(16, lay)

            od_tl(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p_tl(17, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p_tl => predictors_tl(prof)%co2
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p_tl => predictors_tl(prof)%co
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(15, :)
          ENDIF
        ELSE
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay) +       &
                   c(10, lay) * p_tl(10, lay) +     &
                   c(11, lay) * p_tl(11, lay) +     &
                   c(12, lay) * p_tl(12, lay)

            od_tl(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p_tl(13, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p_tl => predictors_tl(prof)%co2
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p_tl => predictors_tl(prof)%co
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(15, :)
          ENDIF
        ENDIF
      ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay) +       &
                   c(10, lay) * p_tl(10, lay) +     &
                   c(11, lay) * p_tl(11, lay) +     &
                   c(12, lay) * p_tl(12, lay)

            od_tl(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p_tl(13, lay)
          ENDDO

          IF (coef%nco2 > 0) THEN
            p_tl => predictors_tl(prof)%co2
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(14, :)
          ENDIF
          IF (coef%nco > 0) THEN
            p_tl => predictors_tl(prof)%co
            od_tl(:) = od_tl(:) + p_tl(1, :) * c(15, :)
          ENDIF
        ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay) +       &
                   c(10, lay) * p_tl(10, lay) +     &
                   c(11, lay) * p_tl(11, lay) +     &
                   c(12, lay) * p_tl(12, lay)

            od_tl(lay) = &
              t(1) + t(2) + t(3) + &
              c(13, lay) * p_tl(13, lay)
          ENDDO
        ELSE   !chanx >= rttov9_wv2380_25
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay) +      &
                   c(2, lay) * p_tl(2, lay) +       &
                   c(3, lay) * p_tl(3, lay) +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay) +       &
                   c(6, lay) * p_tl(6, lay) +       &
                   c(7, lay) * p_tl(7, lay) + &
                   c(8, lay) * p_tl(14, lay)

            t(3) = c(9, lay) * p_tl(15, lay) +       &
                   c(10, lay) * p_tl(18, lay) +      & ! not p10
                   c(11, lay) * p_tl(11, lay) +      &
                   c(12, lay) * p_tl(16, lay)

            od_tl(lay) = t(1) + t(2) + t(3) + &
                      c(13, lay) * p_tl(19, lay)
          ENDDO

          IF (coef%nch4 > 0) THEN
            p    => predictors(prof)%ch4
            p_tl => predictors_tl(prof)%ch4

            od_tl(:) = od_tl(:) + SQRT(p(2, :)) * ( &
              (c(14, :) * p_tl(1, :) + c(15, :) * p_tl(3, :)) + &
              0.5_jprb * p_tl(2, :) * &
                 (c(14, :) * p(1, :) + c(15, :) * p(3, :)) / p(2, :))
          ENDIF
        ENDIF
      ELSE
        od_tl(:) = 0._jprb
      ENDIF

!--------------------------
! Mixed gases
!--------------------------
      IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
        p_tl => predictors_tl(prof)%mixedgas
        c => fast_coef(chan)%mixedgas

        IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
          (chanx < rttov9_wv2000_00)) THEN
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay)  +      &
                   c(2, lay) * p_tl(2, lay)  +       &
                   c(3, lay) * p_tl(3, lay)  +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay)  +      &
                   c(6, lay) * p_tl(6, lay)  +       &
                   c(7, lay) * p_tl(7, lay)  +       &
                   c(8, lay) * p_tl(8, lay)

            od_tl(lay) = od_tl(lay) + t(1) + t(2)
          ENDDO
        ELSE
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay)  +      &
                   c(2, lay) * p_tl(2, lay)  +       &
                   c(3, lay) * p_tl(3, lay)  +       &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay)  + &
                   c(6, lay) * p_tl(6, lay)  +      &
                   c(7, lay) * p_tl(7, lay)  +       &
                   c(8, lay) * p_tl(8, lay)

            od_tl(lay) = od_tl(lay) + t(1) + t(2) +  &
                         c(9, lay) * p_tl(9, lay)  +       &
                         c(10, lay) * p_tl(10, lay)
          ENDDO
        ENDIF
      ENDIF

!-------------
! Ozone
!-------------
      IF (coef%nozone > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%ozone)) THEN
          p_tl => predictors_tl(prof)%ozone
          c => fast_coef(chan)%ozone

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  +      &
                     c(2, lay) * p_tl(2, lay)  +       &
                     c(3, lay) * p_tl(3, lay)  +       &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  +      &
                     c(6, lay) * p_tl(6, lay)  +       &
                     c(7, lay) * p_tl(7, lay)  +       &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  +       &
                           c(10, lay) * p_tl(10, lay)  +     &
                           c(11, lay) * p_tl(11, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  +      &
                     c(2, lay) * p_tl(2, lay)  +       &
                     c(3, lay) * p_tl(12, lay)  +       &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  +      &
                     c(6, lay) * p_tl(6, lay)  +       &
                     c(7, lay) * p_tl(7, lay)  +       &
                     c(8, lay) * p_tl(13, lay)

              t(3) = c(9, lay) * p_tl(9, lay)  +       &
                     c(10, lay) * p_tl(10, lay)  +     &
                     c(11, lay) * p_tl(11, lay)  +     &
                     c(12, lay) * p_tl(14, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(13, lay) * p_tl(15, lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!------------------------------
! Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN !*
          p_tl => predictors_tl(prof)%wvcont
          c => fast_coef(chan)%wvcont

          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay)  +      &
                   c(2, lay) * p_tl(2, lay)  +      &
                   c(3, lay) * p_tl(3, lay)  +      &
                   c(4, lay) * p_tl(4, lay)
            od_tl(lay) = od_tl(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF

!-----------
! CO2
!-----------
  IF (coef%nco2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%co2)) THEN !*
          p    => predictors(prof)%co2
          p_tl => predictors_tl(prof)%co2
          c => fast_coef(chan)%co2

          IF (chanx >= rttov9_wv0690_50 .AND. chanx < rttov9_wv1100_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(5, lay) * p_tl(5, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              t(2) = p(6, lay) * (  &
                       c(3, lay) * (&
                         (2._jprb * p(5, lay) * p_tl(6, lay) + p_tl(5, lay) * p(6, lay))) + &
                       2._jprb * c(6, lay) * p_tl(6, lay))

              t(3) = 2._jprb * ( &
                     c(4, lay) * p(3, lay) * p_tl(3, lay) +  &
                     c(7, lay) * p(7, lay) * p_tl(7, lay))

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(9, lay) * p_tl(9, lay)
            ENDDO
          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(11, lay) * p_tl(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_tl => predictors_tl(prof)%co
              od_tl(:) = od_tl(:) + c(10, :) * p_tl(1, :)
            ENDIF

          ELSE IF (chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)
            ENDDO

          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(7, lay)  + &
                     c(7, lay) * p_tl(8, lay) + &
                     c(9, lay) * p_tl(10, lay)

              t(3) = c(10, lay) * p_tl(11, lay)  +           &
                     c(11, lay) * p_tl(12, lay)  +           &
                     c(12, lay) * p_tl(13, lay)  +           &
                     c(13, lay) * p_tl(14, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(14, lay) * p_tl(15, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_tl => predictors_tl(prof)%co
              od_tl(:) = od_tl(:) + c(8, :) * p_tl(1, :)
            ENDIF

          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay) + &
                           c(11, lay) * p_tl(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p_tl => predictors_tl(prof)%co
              od_tl(:) = od_tl(:) + c(10, :) * p_tl(1, :)
            ENDIF

          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay) + &
                           c(10, lay) * p_tl(10, lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(7, lay) + &
                     c(7, lay) * p_tl(8, lay) + &
                     c(8, lay) * p_tl(10, lay)

              t(3) = c(9, lay) * p_tl(11, lay)  +           &
                     c(10, lay) * p_tl(12, lay)  +           &
                     c(11, lay) * p_tl(13, lay)  +           &
                     c(12, lay) * p_tl(14, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(13, lay) * p_tl(15, lay)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! N2O
!-----------
  IF (coef%nn2o > 0) THEN
    IF (ASSOCIATED(fast_coef(chan)%n2o)) THEN !*
          p_tl => predictors_tl(prof)%n2o
          c => fast_coef(chan)%n2o

          IF (chanx >= rttov9_wv1050_00 .AND. chanx < rttov9_wv1350_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(10, lay) * p_tl(10, lay)
            ENDDO

            IF (coef%nch4 > 0) THEN
              p_tl => predictors_tl(prof)%ch4
              od_tl(:) = od_tl(:) + &
                         c(11, :) * p_tl(1, :) + &
                         c(12, :) * p_tl(10, :)
              ENDIF

          ELSE IF (chanx >= rttov9_wv1995_00 .AND. chanx < rttov9_wv2000_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(10, lay) * p_tl(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p    => predictors(prof)%co
              p_tl => predictors_tl(prof)%co

              od_tl(:) = od_tl(:) + &
                         p(1, :) * (c(12, :) * p_tl(14, :)) + &
                         p_tl(1, :) * (c(11, :) + c(12, :) * p(14, :))
            ENDIF

          ELSE IF (chanx < rttov9_wv2000_00) THEN

            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(10, lay) * p_tl(10, lay)
            ENDDO

          ELSE IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2250_00) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              t(3) = c(9, lay) * p_tl(11, lay)  + &
                     c(10, lay) * p_tl(10, lay) + &
                     c(13, lay) * p_tl(12, lay) + &
                     c(14, lay) * p_tl(13, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3)
            ENDDO

            IF (coef%nco > 0) THEN
              p    => predictors(prof)%co
              p_tl => predictors_tl(prof)%co

              od_tl(:) = od_tl(:) + &
                         p(1, :) * (c(12, :) * p_tl(14, :)) + &
                         p_tl(1, :) * (c(11, :) + c(12, :) * p(14, :))
            ENDIF

          ELSE IF (chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2295_25) THEN

            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                        c(9, lay) * p_tl(9, lay)  + &
                        c(10, lay) * p_tl(10, lay)
            ENDDO

            IF (coef%nco > 0) THEN
              p    => predictors(prof)%co
              p_tl => predictors_tl(prof)%co

              od_tl(:) = od_tl(:) + &
                         p(1, :) * (c(12, :) * p_tl(14, :)) + &
                         p_tl(1, :) * (c(11, :) + c(12, :) * p(14, :))
            ENDIF

          ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2380_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(10, lay) * p_tl(10, lay)
            ENDDO
          ELSE   ! chanx >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              t(3) = c(9, lay) * p_tl(11, lay)  + &
                     c(10, lay) * p_tl(10, lay) + &
                     c(11, lay) * p_tl(12, lay) + &
                     c(12, lay) * p_tl(13, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CO
!-----------
  IF (coef%nco > 0) THEN
    IF (ASSOCIATED(fast_coef(chan)%co)) THEN !*
          p_tl => predictors_tl(prof)%co
          c => fast_coef(chan)%co

          IF ((chanx >= rttov9_wv2250_00 .AND. chanx < rttov9_wv2380_25) .OR. &
            (chanx < rttov9_wv2000_00)) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + &
                           c(9, lay) * p_tl(9, lay)  + &
                           c(10, lay) * p_tl(10, lay) + &
                           c(11, lay) * p_tl(11, lay)
            ENDDO
          ELSE
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay)  + &
                     c(2, lay) * p_tl(2, lay)  + &
                     c(3, lay) * p_tl(3, lay)  + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay)  + &
                     c(6, lay) * p_tl(6, lay)  + &
                     c(7, lay) * p_tl(7, lay)  + &
                     c(8, lay) * p_tl(8, lay)

              t(3) = c(9, lay) * p_tl(9, lay)  + &
                     c(10, lay) * p_tl(10, lay) + &
                     c(11, lay) * p_tl(12, lay) + &
                     c(12, lay) * p_tl(13, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3)
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!-----------
! CH4
!-----------
  IF (coef%nch4 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%ch4)) THEN !*
          p_tl => predictors_tl(prof)%ch4
          c => fast_coef(chan)%ch4

          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(1, lay)  + &
                   c(2, lay) * p_tl(2, lay)  + &
                   c(3, lay) * p_tl(3, lay)  + &
                   c(4, lay) * p_tl(4, lay)

            t(2) = c(5, lay) * p_tl(5, lay)  + &
                   c(6, lay) * p_tl(6, lay)  + &
                   c(7, lay) * p_tl(7, lay)  + &
                   c(8, lay) * p_tl(8, lay)

            t(3) = c(9, lay) * p_tl(9, lay)  + &
                   c(10, lay) * p_tl(10, lay) + &
                   c(11, lay) * p_tl(11, lay)

            od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3)
          ENDDO
        ENDIF
      ENDIF

!-----------
! SO2
!-----------
  IF (coef%nso2 > 0) THEN
        IF (ASSOCIATED(fast_coef(chan)%so2)) THEN !*
          p_tl => predictors_tl(prof)%so2
          c => fast_coef(chan)%so2

          IF (chanx >= rttov9_wv1095_25 .AND. chanx < rttov9_wv1750_25) THEN
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay) + &
                     c(2, lay) * p_tl(2, lay) + &
                     c(3, lay) * p_tl(3, lay) + &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) + &
                     c(6, lay) * p_tl(6, lay) + &
                     c(7, lay) * p_tl(7, lay) + &
                     c(8, lay) * p_tl(8, lay)

              t(3) = c(9, lay) * p_tl(9, lay)   + &
                     c(10, lay) * p_tl(10, lay) + &
                     c(11, lay) * p_tl(11, lay) + &
                     c(12, lay) * p_tl(12, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(13, lay) * p_tl(13, lay)
            ENDDO
          ELSE   !use same predictors as WV >= rttov9_wv2380_25
            DO lay = 1, nlayers
              t(1) = c(1, lay) * p_tl(1, lay) +      &
                     c(2, lay) * p_tl(2, lay) +       &
                     c(3, lay) * p_tl(3, lay) +       &
                     c(4, lay) * p_tl(4, lay)

              t(2) = c(5, lay) * p_tl(5, lay) +       &
                     c(6, lay) * p_tl(6, lay) +       &
                     c(7, lay) * p_tl(7, lay) + &
                     c(8, lay) * p_tl(14, lay)

              t(3) = c(9, lay) * p_tl(15, lay) +       &
                     c(10, lay) * p_tl(18, lay) +      & ! not p10
                     c(11, lay) * p_tl(11, lay) +      &
                     c(12, lay) * p_tl(16, lay)

              od_tl(lay) = od_tl(lay) + t(1) + t(2) + t(3) + &
                           c(13, lay) * p_tl(19, lay)
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
    opdp_path_tl(1, j) = 0.0_jprb
!  Introduce ztemp to stop store followed by load in this recursive loop (DJS)
    ztemp = 0.0_jprb
    DO lev = 2, nlevels
      lay = lev - 1
      IF (opdp_ref(lay,j) < 0._jprb) THEN
        ztemp = ztemp + od_tl(lay)
        opdp_path_tl(lev, j) = ztemp
      ELSE
        opdp_path_tl(lev, j) = opdp_path_tl(lev - 1, j)
      ENDIF
    ENDDO
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_TL', 1_jpim, ZHOOK_HANDLE)

END SUBROUTINE rttov_opdep_9_tl
