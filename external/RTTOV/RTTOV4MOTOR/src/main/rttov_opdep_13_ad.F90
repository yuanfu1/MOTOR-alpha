! Description:
!> @file
!!   AD/K of optical depth calculation (v13 predictors).
!!
!> @brief
!!   AD/K of optical depth calculation (v13 predictors).
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in,out] predictors_ad   increments in predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients (thermal or solar)
!! @param[in]     fast_coef_corr  correction term coefficients (thermal or solar)
!! @param[in,out] opdp_path_ad    optical depth profile increments
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
SUBROUTINE rttov_opdep_13_ad( &
             nlayers,        &
             chanprof,       &
             chanflag,       &
             predictors_ad,  &
             coef,           &
             fast_coef,      &
             fast_coef_corr, &
             opdp_path_ad,   &
             opdp_ref)

  USE rttov_types, ONLY : &
        rttov_chanprof,   &
        rttov_coef,       &
        rttov_fast_coef,  &
        rttov_path_pred,  &
        rttov_opdp_ref_coef

  USE parkind1, ONLY : jprb, jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),             INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),      INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),             INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_path_pred),     INTENT(INOUT) :: predictors_ad(:)
  TYPE(rttov_coef),          INTENT(IN)    :: coef
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef(:)
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef_corr(:)
  REAL(jprb),                INTENT(INOUT) :: opdp_path_ad(:,:)
  TYPE(rttov_opdp_ref_coef), INTENT(IN)    :: opdp_ref
!INTF_END

  REAL(jprb), POINTER :: p_ad(:,:), c(:,:)
  REAL(jprb)          :: od_ad(nlayers), odtot_ad(nlayers), corr_ad(nlayers)
  REAL(jprb)          :: chanx
  REAL(jprb)          :: ztemp
  INTEGER(jpim)       :: lev, lay, prof, chan, j, nlevels, nchanprof
  INTEGER(jpim)       :: map(5), adk
  REAL(jprb)          :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13_AD', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

  IF (SIZE(predictors_ad) == nchanprof) THEN
    adk = 1   ! K
  ELSE
    adk = 0   ! AD
  ENDIF

!----------------------------------------
! Assemble layer optical depths
!----------------------------------------
  DO j = 1, nchanprof
    IF (.NOT. chanflag(j)) CYCLE
    chan = chanprof(j)%chan
    IF (adk == 0) THEN
      prof = chanprof(j)%prof  ! AD
    ELSE
      prof = j                 ! K
    ENDIF
    chanx = coef%ff_cwn(chan)

    ztemp = opdp_path_ad(nlevels,j)
    DO lev = nlevels, 2, -1
      lay = lev - 1
      IF (opdp_ref%od_tot_ref(lay,j) > 0._jprb) THEN
        od_ad(lay) = -ztemp
        ztemp = ztemp + opdp_path_ad(lev-1,j)
      ELSE
        ztemp = ztemp + opdp_path_ad(lev-1,j)
        od_ad(lay) = 0._jprb
      ENDIF
    ENDDO
    opdp_path_ad(1,j) = 0._jprb

    corr_ad  = od_ad
    odtot_ad = od_ad

!-----------
! SO2
!-----------
    IF (coef%nso2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%so2)) THEN
        p_ad => predictors_ad(prof)%so2
        c => fast_coef(chan)%so2

        od_ad = odtot_ad
        WHERE(opdp_ref%od_so2_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:14, lay) = p_ad(1:14, lay) + c(1:14, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%so2)) THEN
          c => fast_coef_corr(chan)%so2

          map(1:5) = (/2,4,5,6,15/)
          DO lay = 1, nlayers
            p_ad(map(1:5), lay) = p_ad(map(1:5), lay) + &
                                  c(1:5, lay) * corr_ad(lay)
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

        od_ad = odtot_ad
        WHERE(opdp_ref%od_ch4_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:11, lay) = p_ad(1:11, lay) + c(1:11, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%ch4)) THEN
          c => fast_coef_corr(chan)%ch4

          map(1:4) = (/7,9,10,12/)
          DO lay = 1, nlayers
            p_ad(map(1:4), lay) = p_ad(map(1:4), lay) + &
                                  c(1:4, lay) * corr_ad(lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! CO
!-----------
    IF (coef%nco > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co)) THEN
        p_ad => predictors_ad(prof)%co
        c => fast_coef(chan)%co

        od_ad = odtot_ad
        WHERE(opdp_ref%od_co_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%co)) THEN
          c => fast_coef_corr(chan)%co

          map(1:4) = (/12,14,15,16/)
          DO lay = 1, nlayers
            p_ad(map(1:4), lay) = p_ad(map(1:4), lay) + &
                                  c(1:4, lay) * corr_ad(lay)
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

        od_ad = odtot_ad
        WHERE(opdp_ref%od_n2o_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:12, lay) = p_ad(1:12, lay) + c(1:12, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%n2o)) THEN
          c => fast_coef_corr(chan)%n2o

          map(1:5) = (/7,8,10,11,12/)
          DO lay = 1, nlayers
            p_ad(map(1:5), lay) = p_ad(map(1:5), lay) + &
                                  c(1:5, lay) * corr_ad(lay)
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
        c => fast_coef(chan)%co2

        od_ad = odtot_ad
        WHERE(opdp_ref%od_co2_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%co2)) THEN
          c => fast_coef_corr(chan)%co2

          map(1:3) = (/14,8,9/)
          DO lay = 1, nlayers
            p_ad(map(1:3), lay) = p_ad(map(1:3), lay) + &
                                  c(1:3, lay) * corr_ad(lay)
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

        od_ad = odtot_ad
        WHERE(opdp_ref%od_wvcont_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:4, lay) = p_ad(1:4, lay) + c(1:4, lay) * od_ad(lay)
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

        od_ad = odtot_ad
        WHERE(opdp_ref%od_o3_ref(:,j) < 0._jprb) od_ad = 0._jprb

        DO lay = 1, nlayers
          p_ad(1:12, lay) = p_ad(1:12, lay) + c(1:12, lay) * od_ad(lay)
        ENDDO

        IF (ASSOCIATED(fast_coef_corr(chan)%ozone)) THEN
          c => fast_coef_corr(chan)%ozone

          DO lay = 1, nlayers
            p_ad(13, lay) = p_ad(13, lay) + c(1, lay) * corr_ad(lay)
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

      od_ad = odtot_ad
      WHERE(opdp_ref%od_wv_ref(:,j) < 0._jprb) od_ad = 0._jprb

      IF (chanx <= 1095._jprb) THEN
        DO lay = 1, nlayers
          p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
        ENDDO
      ELSE IF (chanx <= 2320._jprb) THEN
        DO lay = 1, nlayers
          p_ad(1:14, lay) = p_ad(1:14, lay) + c(1:14, lay) * od_ad(lay)
        ENDDO
      ELSE IF (chanx <= 2570._jprb) THEN
        DO lay = 1, nlayers
          p_ad(1:13, lay) = p_ad(1:13, lay) + c(1:13, lay) * od_ad(lay)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          p_ad(1:14, lay) = p_ad(1:14, lay) + c(1:14, lay) * od_ad(lay)
        ENDDO
      ENDIF

      IF (ASSOCIATED(fast_coef_corr(chan)%watervapour)) THEN
        c => fast_coef_corr(chan)%watervapour

        map(1:5) = (/2,4,5,6,15/)
        DO lay = 1, nlayers
          p_ad(map(1:5), lay) = p_ad(map(1:5), lay) + &
                                c(1:5, lay) * corr_ad(lay)
        ENDDO
      ENDIF
    ENDIF

!--------------------------
! Mixed gases
!--------------------------
    IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
      p_ad => predictors_ad(prof)%mixedgas
      c => fast_coef(chan)%mixedgas

      od_ad = odtot_ad
      WHERE(opdp_ref%od_mg_ref(:,j) < 0._jprb) od_ad = 0._jprb

      DO lay = 1, nlayers
        p_ad(1:9, lay) = p_ad(1:9, lay) + c(1:9, lay) * od_ad(lay)
      ENDDO

      IF (ASSOCIATED(fast_coef_corr(chan)%mixedgas)) THEN
        c => fast_coef_corr(chan)%mixedgas

        map(1:4) = (/2,3,4,10/)
        DO lay = 1, nlayers
          p_ad(map(1:4), lay) = p_ad(map(1:4), lay) + &
                                c(1:4, lay) * corr_ad(lay)
        ENDDO
      ENDIF
    ENDIF
  ENDDO ! chanprof
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_13_ad
