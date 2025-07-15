! Description:
!> @file
!!   AD of optical depth calculation (v7 and v8 predictors).
!!
!> @brief
!!   AD of optical depth calculation (v7 and v8 predictors).
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in,out] predictors_ad   increments in predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       coefficients
!! @param[in,out] opdp_path_ad    optical depth profile increments
!! @param[in]     opdp_ref        reference layer optical depths from direct call
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
SUBROUTINE rttov_opdep_78_ad( &
            nlayers,       &
            chanprof,      &
            predictors_ad, &
            coef,          &
            fast_coef,     &
            opdp_path_ad,  &
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

  INTEGER(KIND=jpim)     , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof(:)
  TYPE(rttov_path_pred)  , INTENT(INOUT) :: predictors_ad(:)
  TYPE(rttov_coef)       , INTENT(IN)    :: coef
  TYPE(rttov_fast_coef)  , INTENT(IN)    :: fast_coef(:)
  TYPE(rttov_opdp_path)  , INTENT(INOUT) :: opdp_path_ad
  REAL(KIND=jprb)        , INTENT(IN)    :: opdp_ref(nlayers,SIZE(chanprof))
!INTF_END

  REAL(KIND=jprb)          :: ztemp
  REAL(KIND=jprb)          :: od_ad(nlayers)
  REAL(KIND=jprb), POINTER :: p_ad(:,:), c(:,:)
  INTEGER(KIND=jpim)       :: lev, lay, chan, prof, j, nlevels, ii, nchanprof

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_78_AD', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

  DO j = 1, nchanprof
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof

!----------------------------------------
! Calculate level-to-space optical depths
!----------------------------------------
    ztemp = opdp_path_ad%atm_level(nlevels, j)
    DO lev = nlevels, 2,  - 1
      lay = lev - 1
      ! Note that optical depth in the calculations is negative
      IF (opdp_ref(lay,j) < 0._jprb) THEN
        od_ad(lay) = ztemp
        ztemp = ztemp + opdp_path_ad%atm_level(lev - 1, j)
      ELSE
        ztemp = ztemp + opdp_path_ad%atm_level(lev - 1, j)
        od_ad(lay) = 0.0_jprb
      ENDIF
    ENDDO
    opdp_path_ad%atm_level(1, j) = 0.0_jprb

!DAR need this again?
    opdp_path_ad%atm_level(2:, j) = 0.0_jprb

!------------------------------------------
! Pressure-modulated cell (pmc) sensors
!------------------------------------------
    IF (coef%pmc_shift) THEN
      DO ii = 1, coef%pmc_nvar
        DO lay = 1, coef%pmc_nlay
          predictors_ad(prof)%pmc(ii, lay, chan) = &
            predictors_ad(prof)%pmc(ii, lay, chan) +              &
            od_ad(lay) * coef%pmc_coef(lay, chan, ii)
        ENDDO
      ENDDO
    ENDIF

!-----------
! CO2 (v8 predictors)
!-----------
    IF (coef%nco2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co2)) THEN
        p_ad => predictors_ad(prof)%co2
        c => fast_coef(chan)%co2

        DO lay = 1, nlayers
          p_ad(:, lay) = p_ad(:, lay) + c(:, lay) * od_ad(lay)
        ENDDO
      ENDIF
    ENDIF

!------------------------------
! Water Vapour Continuum (v8 predictors)
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

        DO lay = 1, nlayers
          p_ad(:, lay) = p_ad(:, lay) + c(:, lay) * od_ad(lay)
        ENDDO
      ENDIF
    ENDIF

!--------------------------
! Mixed gases
!--------------------------
    IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
      p_ad => predictors_ad(prof)%mixedgas
      c => fast_coef(chan)%mixedgas

      DO lay = 1, nlayers
        p_ad(:, lay) = p_ad(:, lay) + c(:, lay) * od_ad(lay)
      ENDDO
    ENDIF

!--------------------
! Water vapour
!--------------------
    IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
      p_ad => predictors_ad(prof)%watervapour
      c => fast_coef(chan)%watervapour

      DO lay = 1, nlayers
        p_ad(:, lay) = p_ad(:, lay) + c(:, lay) * od_ad(lay)
      ENDDO
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_78_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_78_ad
