! Description:
!> @file
!!   TL of optical depth calculation (v13 predictors).
!!
!> @brief
!!   TL of optical depth calculation (v13 predictors).
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in]     predictors_tl   perturbations in predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients (thermal or solar)
!! @param[in]     fast_coef_corr  correction term coefficients (thermal or solar)
!! @param[in,out] opdp_path_tl    calculated optical depth profile perturbations
!! @param[in]     opdp_ref        reference optical depths (may include negative values)
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
!    Copyright 2019, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_opdep_13_tl( &
             nlayers,        &
             chanprof,       &
             chanflag,       &
             predictors_tl,  &
             coef,           &
             fast_coef,      &
             fast_coef_corr, &
             opdp_path_tl,   &
             opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,   &
        rttov_coef,       &
        rttov_fast_coef,  &
        rttov_path_pred,  &
        rttov_opdp_ref_coef

  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),             INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),      INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),             INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_path_pred),     INTENT(IN)    :: predictors_tl(:)
  TYPE(rttov_coef),          INTENT(IN)    :: coef
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef(:)
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef_corr(:)
  REAL(jprb),                INTENT(INOUT) :: opdp_path_tl(:,:)
  TYPE(rttov_opdp_ref_coef), INTENT(IN)    :: opdp_ref
!INTF_END

  REAL(jprb), POINTER :: p_tl(:,:), c(:,:)
  REAL(jprb)          :: od_tl(nlayers), odtot_tl(nlayers), corr_tl(nlayers)
  REAL(jprb)          :: chanx
  REAL(jprb)          :: t(4), ztemp
  INTEGER(jpim)       :: lev, lay, prof, chan, j, nlevels, nchannels
  REAL(jprb)          :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  nlevels   = nlayers + 1

!-----------------------------------------
! Calculate layer gaseous optical depths
!-----------------------------------------

  DO j = 1, nchannels

    t = 0._jprb

    IF (.NOT. chanflag(j)) CYCLE
    chan  = chanprof(j)%chan
    chanx = coef%ff_cwn(chan)
    prof  = chanprof(j)%prof

!--------------------------
! Mixed gases
!--------------------------
    IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
      p_tl => predictors_tl(prof)%mixedgas
      c => fast_coef(chan)%mixedgas

      DO lay = 1, nlayers
        t(1) = c(1, lay) * p_tl(1, lay)  + &
               c(2, lay) * p_tl(2, lay)  + &
               c(3, lay) * p_tl(3, lay)  + &
               c(4, lay) * p_tl(4, lay)

        t(2) = c(5, lay) * p_tl(5, lay)  + &
               c(6, lay) * p_tl(6, lay)  + &
               c(7, lay) * p_tl(7, lay)  + &
               c(8, lay) * p_tl(8, lay)

        od_tl(lay) = t(1) + t(2) + c(9, lay) * p_tl(9, lay)
      ENDDO

      WHERE (opdp_ref%od_mg_ref(:,j) < 0._jprb) od_tl = 0.
      odtot_tl = od_tl

      IF (ASSOCIATED(fast_coef_corr(chan)%mixedgas)) THEN
        c => fast_coef_corr(chan)%mixedgas

        DO lay = 1, nlayers
          t(1) = c(1, lay) * p_tl(2, lay)  + &
                 c(2, lay) * p_tl(3, lay)  + &
                 c(3, lay) * p_tl(4, lay)  + &
                 c(4, lay) * p_tl(10, lay)

          corr_tl(lay) = t(1)
        ENDDO
      ELSE
        corr_tl = 0._jprb
      ENDIF
    ELSE
      odtot_tl = 0._jprb
      corr_tl = 0._jprb
    ENDIF

!--------------------
! Water vapour
!--------------------
    IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
      p_tl => predictors_tl(prof)%watervapour
      c => fast_coef(chan)%watervapour

      IF (chanx <= 1095._jprb) THEN
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

          od_tl(lay) = t(1) + t(2) + t(3) + c(13, lay) * p_tl(13, lay)
        ENDDO
      ELSE IF (chanx <= 2320._jprb) THEN
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
                       c(13, lay) * p_tl(13, lay) + &
                       c(14, lay) * p_tl(14, lay)
        ENDDO
      ELSE IF (chanx <= 2570._jprb) THEN
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

          od_tl(lay) = t(1) + t(2) + t(3) + c(13, lay) * p_tl(13, lay)
        ENDDO
      ELSE
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
                       c(13, lay) * p_tl(13, lay) + &
                       c(14, lay) * p_tl(14, lay)
        ENDDO
      ENDIF

      WHERE (opdp_ref%od_wv_ref(:,j) < 0._jprb) od_tl = 0.
      odtot_tl = odtot_tl + od_tl

      IF (ASSOCIATED(fast_coef_corr(chan)%watervapour)) THEN
        c => fast_coef_corr(chan)%watervapour

        DO lay = 1, nlayers
          t(1) = c(1, lay) * p_tl(2, lay) + &
                 c(2, lay) * p_tl(4, lay) + &
                 c(3, lay) * p_tl(5, lay) + &
                 c(4, lay) * p_tl(6, lay)
          
          corr_tl(lay) = corr_tl(lay) + t(1) + c(5, lay) * p_tl(15, lay)
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

          od_tl(lay) = t(1) + t(2) + t(3)
        ENDDO

        WHERE (opdp_ref%od_o3_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%ozone)) THEN
          c => fast_coef_corr(chan)%ozone

          DO lay = 1, nlayers
            corr_tl(lay) = corr_tl(lay) + c(1, lay) * p_tl(13, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!------------------------------
! Water Vapour Continuum
!------------------------------
    IF (coef%nwvcont > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN
        p_tl => predictors_tl(prof)%wvcont
        c => fast_coef(chan)%wvcont

        DO lay = 1, nlayers
          od_tl(lay) = c(1, lay) * p_tl(1, lay) + &
                       c(2, lay) * p_tl(2, lay) + &
                       c(3, lay) * p_tl(3, lay) + &
                       c(4, lay) * p_tl(4, lay)
        ENDDO

        WHERE (opdp_ref%od_wvcont_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl
      ENDIF
    ENDIF

!-----------
! CO2
!-----------
    IF (coef%nco2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co2)) THEN
        p_tl => predictors_tl(prof)%co2
        c => fast_coef(chan)%co2

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

          od_tl(lay) = t(1) + t(2) + t(3) + c(13, lay) * p_tl(13, lay)
        ENDDO

        WHERE (opdp_ref%od_co2_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%co2)) THEN
          c => fast_coef_corr(chan)%co2

          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(14, lay) + &
                   c(2, lay) * p_tl(8, lay)  + &
                   c(3, lay) * p_tl(9, lay)

            corr_tl(lay) = corr_tl(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! N2O
!-----------
    IF (coef%nn2o > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%n2o)) THEN
        p_tl => predictors_tl(prof)%n2o
        c => fast_coef(chan)%n2o

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

          od_tl(lay) = t(1) + t(2) + t(3)
        ENDDO

        WHERE (opdp_ref%od_n2o_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%n2o)) THEN
          c => fast_coef_corr(chan)%n2o
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(7, lay)  + &
                   c(2, lay) * p_tl(8, lay)  + &
                   c(3, lay) * p_tl(10, lay) + &
                   c(4, lay) * p_tl(11, lay)

            corr_tl(lay) = corr_tl(lay) + t(1) + c(5, lay) * p_tl(12, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! CO
!-----------
    IF (coef%nco > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co)) THEN
        p_tl => predictors_tl(prof)%co
        c => fast_coef(chan)%co

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

          od_tl(lay) = t(1) + t(2) + t(3) + c(13, lay) * p_tl(13, lay)
        ENDDO

        WHERE (opdp_ref%od_co_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%co)) THEN
          c => fast_coef_corr(chan)%co
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(12, lay) + &
                   c(2, lay) * p_tl(14, lay) + &
                   c(3, lay) * p_tl(15, lay) + &
                   c(4, lay) * p_tl(16, lay)

            corr_tl(lay) = corr_tl(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! CH4
!-----------
    IF (coef%nch4 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%ch4)) THEN
        p_tl => predictors_tl(prof)%ch4
        c => fast_coef(chan)%ch4

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
                 c(11, lay) * p_tl(11, lay)

          od_tl(lay) = t(1) + t(2) + t(3)
        ENDDO

        WHERE (opdp_ref%od_ch4_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%ch4)) THEN
          c => fast_coef_corr(chan)%ch4
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(7, lay)  + &
                   c(2, lay) * p_tl(9, lay)  + &
                   c(3, lay) * p_tl(10, lay) + &
                   c(4, lay) * p_tl(12, lay)

            corr_tl(lay) = corr_tl(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! SO2
!-----------
    IF (coef%nso2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%so2)) THEN
        p_tl => predictors_tl(prof)%so2
        c => fast_coef(chan)%so2

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
                       c(13, lay) * p_tl(13, lay) + &
                       c(14, lay) * p_tl(14, lay)
        ENDDO

        WHERE (opdp_ref%od_so2_ref(:,j) < 0._jprb) od_tl = 0.
        odtot_tl = odtot_tl + od_tl

        IF (ASSOCIATED(fast_coef_corr(chan)%so2)) THEN
          c => fast_coef_corr(chan)%so2
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p_tl(2, lay) + &
                   c(2, lay) * p_tl(4, lay) + &
                   c(3, lay) * p_tl(5, lay) + &
                   c(4, lay) * p_tl(6, lay)
            
            corr_tl(lay) = corr_tl(lay) + t(1) + c(5, lay) * p_tl(15, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!----------------------------------------
! Assemble layer optical depths
!----------------------------------------
    od_tl = odtot_tl + corr_tl
    WHERE (opdp_ref%od_tot_ref(:,j) < 0._jprb) od_tl = 0._jprb

    opdp_path_tl(1,j) = 0._jprb
    ztemp = 0._jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ztemp = ztemp - od_tl(lay)
      opdp_path_tl(lev,j) = ztemp
    ENDDO

  ENDDO ! chanprof
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_13_tl
