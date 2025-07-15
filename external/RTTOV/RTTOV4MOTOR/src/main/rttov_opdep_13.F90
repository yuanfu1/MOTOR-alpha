! Description:
!> @file
!!   Calculates level-to-space optical depth profiles for each channel by
!!   applying the optical depth regression for v13 predictors.
!!
!> @brief
!!   Calculates level-to-space optical depth profiles for each channel by
!!   applying the optical depth regression for v13 predictors.
!!
!! @details
!!   This subroutine calculates optical depths for v13 predictors
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
!!   Optical depths for each gas are computed separately, and modified to
!!   ensure they are greater than or equal to zero. In addition, a correction
!!   term is predicted which is added to the summed layer gas optical depths
!!   to account for the error due to the gas optical depths being polychromatic.
!!   The final total layer optical depths are also modified to ensure they are
!!   greater than or equal to zero. Finally, optical depths are accumulated to
!!   obtain level-to-space optical depths.
!!
!!   This subroutine operates on coefficient layers/levels.
!!
!!   NB The output optical depths (opdp_path) are NEGATIVE numbers, but all
!!      other optical depths here are positive numbers.
!!
!! @param[in]     nlayers         number of coefficient layers
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     chanflag        flag to indicate if channel should be processed (thermal/solar)
!! @param[in]     predictors      pre-calculated predictors
!! @param[in]     coef            rttov_coef structure
!! @param[in]     fast_coef       gas optical depth coefficients (thermal or solar)
!! @param[in]     fast_coef_corr  correction term coefficients (thermal or solar)
!! @param[in,out] opdp_path       calculated optical depth profiles
!! @param[in,out] opdp_ref        reference optical depths (may include negative values)
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
SUBROUTINE rttov_opdep_13( &
             nlayers,        &
             chanprof,       &
             chanflag,       &
             predictors,     &
             coef,           &
             fast_coef,      &
             fast_coef_corr, &
             opdp_path,      &
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
  TYPE(rttov_path_pred),     INTENT(IN)    :: predictors(:)
  TYPE(rttov_coef),          INTENT(IN)    :: coef
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef(:)
  TYPE(rttov_fast_coef),     INTENT(IN)    :: fast_coef_corr(:)
  REAL(jprb),                INTENT(INOUT) :: opdp_path(:,:)
  TYPE(rttov_opdp_ref_coef), INTENT(INOUT) :: opdp_ref
!INTF_END

  REAL(jprb), POINTER :: p(:,:), c(:,:)
  REAL(jprb)          :: od(nlayers), odtot(nlayers), corr(nlayers)
  REAL(jprb)          :: chanx
  REAL(jprb)          :: t(4), ztemp
  INTEGER(jpim)       :: lev, lay, j, prof, chan, nlevels, nchanprof
  REAL(jprb)          :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

!-----------------------------------------
! Calculate layer gaseous optical depths
!-----------------------------------------
  DO j = 1, nchanprof

    t = 0._jprb

    IF (.NOT. chanflag(j)) CYCLE
    chan  = chanprof(j)%chan
    chanx = coef%ff_cwn(chan)
    prof  = chanprof(j)%prof

!--------------------------
! Mixed gases
!--------------------------
    IF (ASSOCIATED(fast_coef(chan)%mixedgas)) THEN
      p => predictors(prof)%mixedgas
      c => fast_coef(chan)%mixedgas

      DO lay = 1, nlayers
        t(1) = c(1, lay) * p(1, lay)  + &
               c(2, lay) * p(2, lay)  + &
               c(3, lay) * p(3, lay)  + &
               c(4, lay) * p(4, lay)

        t(2) = c(5, lay) * p(5, lay)  + &
               c(6, lay) * p(6, lay)  + &
               c(7, lay) * p(7, lay)  + &
               c(8, lay) * p(8, lay)

        od(lay) = t(1) + t(2) + c(9, lay) * p(9, lay)
      ENDDO

      ! Store gas ref OD
      opdp_ref%od_mg_ref(:,j) = od

      ! Simple clipping of negative ODs
      WHERE (od < 0._jprb) od = 0._jprb

      ! Accumulate layer optical depths
      odtot = od

      IF (ASSOCIATED(fast_coef_corr(chan)%mixedgas)) THEN
        c => fast_coef_corr(chan)%mixedgas

        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(2, lay)  + &
                 c(2, lay) * p(3, lay)  + &
                 c(3, lay) * p(4, lay)  + &
                 c(4, lay) * p(10, lay)

          corr(lay) = t(1)
        ENDDO
      ELSE
        corr = 0._jprb
      ENDIF
    ELSE
      odtot = 0._jprb
      corr = 0._jprb
    ENDIF

!--------------------
! Water vapour
!--------------------
    IF (ASSOCIATED(fast_coef(chan)%watervapour)) THEN
      p => predictors(prof)%watervapour
      c => fast_coef(chan)%watervapour

      IF (chanx <= 1095._jprb) THEN
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

          od(lay) = t(1) + t(2) + t(3) + c(13, lay) * p(13, lay)
        ENDDO
      ELSE IF (chanx <= 2320._jprb) THEN
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

          od(lay) = t(1) + t(2) + t(3) +      &
                    c(13, lay) * p(13, lay) + &
                    c(14, lay) * p(14, lay)
        ENDDO
      ELSE IF (chanx <= 2570._jprb) THEN
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

          od(lay) = t(1) + t(2) + t(3) + c(13, lay) * p(13, lay)
        ENDDO
      ELSE
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

          od(lay) = t(1) + t(2) + t(3) +      &
                    c(13, lay) * p(13, lay) + &
                    c(14, lay) * p(14, lay)
        ENDDO
      ENDIF

      ! Store gas ref OD
      opdp_ref%od_wv_ref(:,j) = od

      ! Simple clipping of negative ODs
      WHERE (od < 0._jprb) od = 0._jprb

      ! Accumulate layer optical depths
      odtot = odtot + od

      IF (ASSOCIATED(fast_coef_corr(chan)%watervapour)) THEN
        c => fast_coef_corr(chan)%watervapour
        DO lay = 1, nlayers
          t(1) = c(1, lay) * p(2, lay) + &
                 c(2, lay) * p(4, lay) + &
                 c(3, lay) * p(5, lay) + &
                 c(4, lay) * p(6, lay)
          
          corr(lay) = corr(lay) + t(1) + c(5, lay) * p(15, lay)
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

          od(lay) = t(1) + t(2) + t(3)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_o3_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%ozone)) THEN
          c => fast_coef_corr(chan)%ozone
          DO lay = 1, nlayers
            corr(lay) = corr(lay) + c(1, lay) * p(13, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!------------------------------
! Water Vapour Continuum
!------------------------------
    IF (coef%nwvcont > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%wvcont)) THEN
        p => predictors(prof)%wvcont
        c => fast_coef(chan)%wvcont

        DO lay = 1, nlayers
          od(lay) = c(1, lay) * p(1, lay) + &
                    c(2, lay) * p(2, lay) + &
                    c(3, lay) * p(3, lay) + &
                    c(4, lay) * p(4, lay)
          
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_wvcont_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        ! No correction term contribution

      ENDIF
    ENDIF

!-----------
! CO2
!-----------
    IF (coef%nco2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co2)) THEN
        p => predictors(prof)%co2
        c => fast_coef(chan)%co2

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

          od(lay) = t(1) + t(2) + t(3) + c(13, lay) * p(13, lay)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_co2_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%co2)) THEN
          c => fast_coef_corr(chan)%co2
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(14, lay) + &
                   c(2, lay) * p(8, lay)  + &
                   c(3, lay) * p(9, lay)

            corr(lay) = corr(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! N2O
!-----------
    IF (coef%nn2o > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%n2o)) THEN
        p => predictors(prof)%n2o
        c => fast_coef(chan)%n2o

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

          od(lay) = t(1) + t(2) + t(3)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_n2o_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%n2o)) THEN
          c => fast_coef_corr(chan)%n2o
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(7, lay)  + &
                   c(2, lay) * p(8, lay)  + &
                   c(3, lay) * p(10, lay) + &
                   c(4, lay) * p(11, lay)

            corr(lay) = corr(lay) + t(1) + c(5, lay) * p(12, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! CO
!-----------
    IF (coef%nco > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%co)) THEN
        p => predictors(prof)%co
        c => fast_coef(chan)%co

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

          od(lay) = t(1) + t(2) + t(3) + c(13, lay) * p(13, lay)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_co_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%co)) THEN
          c => fast_coef_corr(chan)%co
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(12, lay) + &
                   c(2, lay) * p(14, lay) + &
                   c(3, lay) * p(15, lay) + &
                   c(4, lay) * p(16, lay)

            corr(lay) = corr(lay) + t(1) 
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! CH4
!-----------
    IF (coef%nch4 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%ch4)) THEN
        p => predictors(prof)%ch4
        c => fast_coef(chan)%ch4

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
                 c(11, lay) * p(11, lay)

          od(lay) = t(1) + t(2) + t(3)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_ch4_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%ch4)) THEN
          c => fast_coef_corr(chan)%ch4
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(7, lay)  + &
                   c(2, lay) * p(9, lay)  + &
                   c(3, lay) * p(10, lay) + &
                   c(4, lay) * p(12, lay)

            corr(lay) = corr(lay) + t(1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!-----------
! SO2
!-----------
    IF (coef%nso2 > 0) THEN
      IF (ASSOCIATED(fast_coef(chan)%so2)) THEN
        p => predictors(prof)%so2
        c => fast_coef(chan)%so2

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
                    c(13, lay) * p(13, lay) + &
                    c(14, lay) * p(14, lay)
        ENDDO

        ! Store gas ref OD
        opdp_ref%od_so2_ref(:,j) = od

        ! Simple clipping of negative ODs
        WHERE (od < 0._jprb) od = 0._jprb

        ! Accumulate layer optical depths
        odtot = odtot + od

        IF (ASSOCIATED(fast_coef_corr(chan)%so2)) THEN
          c => fast_coef_corr(chan)%so2
          DO lay = 1, nlayers
            t(1) = c(1, lay) * p(2, lay) + &
                   c(2, lay) * p(4, lay) + &
                   c(3, lay) * p(5, lay) + &
                   c(4, lay) * p(6, lay)
            
            corr(lay) = corr(lay) + t(1) + c(5, lay) * p(15, lay)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

!----------------------------------------
! Assemble layer optical depths
!----------------------------------------
    ! Gas optical depths computed above are positive
    ! - odtot contains total gas layer optical depths and is non-negative
    ! - corr contains the layer optical depth corrections

    ! Calculate total layer optical depths with corrections (these may be negative)
    od = odtot + corr

    ! Store reference total optical depths
    opdp_ref%od_tot_ref(:,j) = od

    ! Simple clipping of negative ODs
    WHERE (od < 0._jprb) od = 0._jprb

    ! od is non-negative: accumulate TOA-level optical depths
    ! Note that opdp_path is expected to be negative by the rest of RTTOV
    opdp_path(1,j) = 0._jprb
    ! Introduce ztemp to stop store followed by load in this recursive loop (DJS)
    ztemp = 0._jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ztemp = ztemp - od(lay)
      opdp_path(lev,j) = ztemp
    ENDDO
  ENDDO ! chanprof
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_13', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_13
