!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.BcPredictors
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/01/17, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/12/25, @GBA-MWF, Shenzhen
!   for adding a safeguard of PRED_CONV to check pressure values.
!   for modifying PRED_CONV for using pressure level from bottom to top as the original is incorrect
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides predictors for bias correction (BC).

MODULE BcPredictors_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE State_m, ONLY: State_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE Satellite_utils_m

CONTAINS

  SUBROUTINE BC_predictors(nobs, configFile, X, pres, t, q, tskin, satzenith, pred)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: nobs
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(IN) ::  pres(:, :), t(:, :), q(:, :), tskin(:), satzenith(:)
    !   TYPE(ObsSet_t), INTENT(IN) :: Y
    CHARACTER(LEN=50), ALLOCATABLE :: bias_predictors(:)
    INTEGER(i_kind) :: npred
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: pred(:, :)
    INTEGER(i_kind) :: ipred, nlevs, nprofs, iobs
    REAL(r_kind) :: layer_base, layer_top
    INTEGER(i_kind) :: istatus

    istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'bias_predictors', bias_predictors)
    npred = SIZE(bias_predictors, 1)
    PRINT *, 'Use bias_predictors: ', bias_predictors

    pred = ZERO
    nlevs = X%sg%vLevel
    nprofs = X%sg%num_icell

    DO ipred = 1, npred
      PRINT *, 'Calculate the predictor: ', TRIM(bias_predictors(ipred))
      SELECT CASE (TRIM(bias_predictors(ipred)))
      CASE ('constant')
        pred(:, ipred) = 1.0
      CASE ('thickness_1000_300')
        layer_base = 1000.0
        layer_top = 300.0
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_thickness(nlevs, pres(iobs, :), t(iobs, :), q(iobs, :), layer_base, layer_top) ! thickness is a function
        END DO
        !  PRINT *, 'thickness_1000_300 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
      CASE ('thickness_200_50')
        layer_base = 200.0
        layer_top = 50.0
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_thickness(nlevs, pres(iobs, :), t(iobs, :), q(iobs, :), layer_base, layer_top) ! thickness is a function
        END DO
        ! PRINT *, 'thickness_200_50 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
        ! CASE ('thickness_50_5')
        !   layer_base = 50.0
        !   layer_top = 5.0
        !   pred(:, ipred) = Calc_Pred_thickness(nlevs, nprofs, pres, t, layer_base, layer_top)
        !   CALL Calc_Pred_thickness(nlevs, nprofs, pres, t, q, layer_base, layer_top, pred(:, ipred)) ! thickness is a function
        !   PRINT *, 'thickness_50_5 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
        !     CASE ('wind_speed_10m')
        !       pred(:, ipred) = sqrt( u10m**2 + v10m **2 )
      CASE ('scan_position_order1')
        ! For polar satellites, the scan predictor is scan position
        ! scan_predictor = Y%scanpos
        ! scanpos = ifov
        ! For geostationary satellites, zenith angle is used
        ! IF (SIZE(Y%ObsFields,1) > 0) THEN
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_ScanAngle(satzenith(iobs), 1) ! ScanAgnle is a function
        END DO
        ! PRINT *, 'scan_position_order1 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
      CASE ('scan_position_order2')
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_ScanAngle(satzenith(iobs), 2) ! ScanAgnle is a function
        END DO
        ! PRINT *, 'scan_position_order2 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
      CASE ('scan_position_order3')
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_ScanAngle(satzenith(iobs), 3) ! ScanAgnle is a function
        END DO
        ! PRINT *, 'scan_position_order3 values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
        !     CASE ('lapse_rate')
        !       pred(:, ipred) = Calc_Pred_T_lapse() ! T_lapse is a function
      CASE ('skin_temperature')
        pred(:, ipred) = tskin
        PRINT *, 'skin_temperature values: ', MAXVAL(pred(:, ipred)), MINVAL(pred(:, ipred))
      CASE ('TPW_guess')
        DO iobs = 1, nobs
          IF (MAXVAL(pres(iobs, :)) - MINVAL(pres(iobs, :)) < 1.0D0) THEN
            ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
            CYCLE
          END IF
          pred(:, ipred) = Calc_Pred_qv2tpw(nlevs, pres(iobs, :), q(iobs, :)) ! qv2tpw is a function
        END DO
        ! PRINT *, 'TPW_guess values: ', maxval(pred(:, ipred)), minval(pred(:, ipred))
        !     CASE ('TPW_obs')
        !       pred(:, ipred) = Y%TPW ! Calculated from O instead of B
        !     CASE ('LPW_guess')
        !       pred(:, ipred) = Calc_Pred_qv2lpw(X%qvapor) ! qv2lpw is a function
        !     CASE ('LPW_obs')
        !       pred(:, ipred) = Y%LPW ! Calculated from O instead of B

        !     CASE ('clw_guess')
        !       pred(:, ipred) = Calc_Pred_clw_guess(this)
        !     CASE ('emissivity')
        !       pred(:, ipred) = HX%emis
        !     CASE ('sin_of_latitude')
        !       pred(:, ipred) = sin(this%latitude*degree2radian)
        !     CASE ('orbital_angle')
        !       CALL Calc_Pred_OrbitalAngle(pred(:, ipred)) ! OrbitalAngle is a function
        !     CASE ('zenith_angle')
        !       pred(:, ipred) = Y%zenangle
        !     CASE ('latitude')
        !       pred(:, ipred) = Y%latitude
        !     CASE ('observed_BT')
        !       pred(:, ipred) = Y%BT
        !     CASE ('model_BT')
        !       pred(:, ipred) = HX%BT
        !     CASE ('cloud_top_height')
        !       pred(:, ipred) = Y%CTH
      CASE default
        pred(:, ipred) = 1.0
      END SELECT
    END DO

    !   ! Normalization
    !   DO ipred = 2, npred
    !     pred_mean(ipred) = sum(pred(:, ipred))/real(npred)
    !     pred_std(ipred) = sqrt(pred(:, ipred) - pred_mean(ipred))/real(npred)
    !     pred(ipred) = (pred(ipred) - pred_mean(ipred))/pred_std(ipred)
    !   END DO

    !   ! Read in coefficients for each predictor
    !   ! Currently, we only have constant, thickness_troposphere, thickness_tropopause, scan_position_order1,
    !   ! skin_temperature, total precipitable water.
    !   ! OPEN ()
    !   ! CLOSE ()

  END SUBROUTINE BC_predictors

  ! FUNCTION Calc_Pred_T_lapse() RESULT(lapse)
  !   REAL(r_kind), INTENT(in) ::
  !   REAL(r_kind) :: lapse

  ! END RUNCTION Calc_Pred_T_lapse

  ! FUNCTION Calc_Pred_clw_guess(this) RESULT(clw_guess)
  !   ! clw or cloud content is for MW radiance only
  !   REAL(r_kind), INTENT(in) :: cloud
  !   REAL(r_kind) :: cloud_cont ! cloud content
  !   REAL(r_kind) :: clw_guess

  !   clw_guess = ZERO
  !   cloud_cont = qc*kgkg_kgm2
  !   clw_guess = clw_guess + cloud_cont

  ! END RUNCTION Calc_Pred_clw_guess

  FUNCTION Calc_Pred_qv2tpw(nlevs, pres, qv) RESULT(tpw)
    REAL(r_kind), INTENT(in) :: pres(:), qv(:)
    REAL(r_kind) :: tpw
    REAL(r_kind), ALLOCATABLE :: qm(:), dp(:), pres1(:)
    INTEGER(i_kind) :: nlev, k, iprof

    tpw = ZERO

    IF (.NOT. ALL(qv .LT. 0.0001D0)) THEN
      pres1 = pres(:) * Pa2hPa
      dp = pres1; dp = ZERO; qm = dp
      dp(1:nlevs - 1) = pres1(1:nlevs - 1) - pres1(2:nlevs)

      ! Full-level to half-level
      DO k = 1, nlevs - 1
        qm(k) = 0.5D0 * (qv(k) + qv(k + 1))
      END DO

      tpw = 100.0 / g * SUM(qm(1:nlevs - 1) * dp(1:nlevs - 1))
      ! PRINT *, 'tpw guess = ', tpw

    END IF

    IF (ALLOCATED(pres1)) DEALLOCATE (pres1, qm, dp)

  END FUNCTION Calc_Pred_qv2tpw

  ! FUNCTION Calc_Pred_qv2lpw(qv) RESULT(lpw)
  !   REAL(r_kind), INTENT(in) :: qv
  !   REAL(r_kind) :: lpw

  ! END RUNCTION Calc_Pred_qv2lpw

  FUNCTION Calc_Pred_thickness(nlevs, pres, t, q, layer_base, layer_top) RESULT(layer_thickness)
    ! Algorithm: hydrostatic equation
    ! For simplication, equation of state of wet air is calculated with the same form as that of dry air,
    ! but with T repaced by T_v
    ! Equations:
    ! (1) \frac{\partial p}{\partial z}=-rho*g
    ! (2) P = rho * R * T_v
    ! (3) T_v = ( 1 + 0.608*q ) * T
    ! With the above three equations, we have
    ! delta_z = -R * T_v * ( lnP2 - ln P1 ) / g
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: nlevs
    REAL(r_kind), INTENT(IN) :: layer_base, layer_top
    REAL(r_kind), INTENT(IN) :: pres(:), t(:), q(:)
    REAL(r_kind) :: layer_thickness
    REAL(r_kind) :: model_top, model_bot
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: qm, tv, dlp, pres1
    INTEGER(i_kind) :: k, k_tmp(1)
    INTEGER(i_kind) :: k_top(1), k_base(1)
    INTEGER(i_kind) :: iprof
    REAL(r_kind) :: add_thk

    ! PRINT *, 'Calc_Pred_thickness bewteen ', layer_base, ' and ', layer_top

    ALLOCATE (pres1(nlevs), qm(nlevs), tv(nlevs), dlp(nlevs))
    layer_thickness = ZERO

    pres1 = pres(:) * Pa2hPa  ! pres(1) is the bottom, pres(nlevs) is the top
    model_top = MINVAL(pres1)
    model_bot = MAXVAL(pres1)
    ! Note: rttov uses hPa, MOTOR-DA uses Pa

    IF (.NOT. ALL(pres1 .LT. 0.1D0)) THEN

      k_tmp(1:1) = MINLOC(ABS(pres1 - layer_base)) ! Example: 1000 hPa
      k_base = k_tmp
      k_tmp(1:1) = MINLOC(ABS(pres1 - layer_top)) ! Example: 300 hPa
      k_top = k_tmp

      DO K = 1, nlevs - 1
        qm(k) = 0.5D0 * (q(k) + q(k + 1))
        tv(k) = 0.5D0 * (t(k) + t(k + 1)) * (1.0D0 + 0.608 * qm(k))
      END DO
      dlp(1:nlevs - 1) = LOG(pres1(1:nlevs - 1)) - LOG(pres1(2:nlevs))

      IF ((1000.0 > pres1(1)) .AND. layer_base > 900.0) THEN
        add_thk = dry_air_gas_const * (tv(1) * (LOG(pres1(1)) - LOG(1000.0))) / g  ! approximation for levels below 1000.0 hPa
      ELSE
        add_thk = 0.0D0
      END IF

      ! PRINT *, 'pres1 = ', pres1
      layer_thickness = dry_air_gas_const * km2m * SUM(tv(k_base(1):k_top(1) + 1) * dlp(k_base(1):k_top(1) + 1)) / g
      ! PRINT *, 'layer_thickness = ', layer_thickness
    END IF

    DEALLOCATE (pres1, qm, tv, dlp)

  END FUNCTION Calc_Pred_thickness

  FUNCTION Calc_Pred_ScanAngle(scan_angle, order) RESULT(sa)
    ! the scan bias b^scan is a function of latitude as well as scan position
    REAL(r_kind), INTENT(in) :: scan_angle
    INTEGER(i_kind), INTENT(in) :: order
    REAL(r_kind) :: sa

    sa = (scan_angle * radian2degree)**order
    ! PRINT *, 'scan_angle = ', sa

  END FUNCTION Calc_Pred_ScanAngle

  SUBROUTINE Calc_Pred_OrbitalAngle(orbital_angle, order, cos_oa)
    REAL(r_kind), INTENT(in) :: orbital_angle(:)
    INTEGER(i_kind), INTENT(in) :: order
    REAL(r_kind), INTENT(INOUT) :: cos_oa(:)

    ! NOTE: Fortran uses radian as the unit for trigonometric functions
    cos_oa = COS(orbital_angle * order)

  END SUBROUTINE Calc_Pred_OrbitalAngle

  SUBROUTINE GETSCORR(SCORR, LAT, vmnrlb, JSCAN)

    IMPLICIT NONE

    ! subroutine arguments
    !-------------------------
    !REAL(r_kind),            INTENT(OUT) :: SCORR(JPCHAN)
    REAL(r_kind), INTENT(OUT) :: SCORR
    REAL(r_kind), INTENT(IN)  :: LAT
    !REAL(r_kind),            INTENT(INOUT)  :: vmnrlb(JPCHAN,JPSCAN,JBAND)
    REAL(r_kind), INTENT(IN)  :: vmnrlb(:, :)
    INTEGER(i_kind), INTENT(IN)  :: JSCAN

    ! local variables
    !-------------------
    INTEGER :: sband, i, JPSCAN, JBAND
    REAL    :: BLAT

    INTEGER :: BOFF
    REAL    :: BDIV

    BOFF = JBAND / 2 + 1
    BDIV = 180.0 / JBAND + 0.0001
    sband = FLOOR(LAT / BDIV) + BOFF
    BLAT = FLOOR(LAT / BDIV) * BDIV

    IF (LAT >= BLAT + BDIV / 2) THEN

      IF (sband < JBAND) THEN
        !SCORR(1:JPCHAN) = (LAT -   (BLAT+BDIV/2)) * vmnrlb(1:JPCHAN,JSCAN,sband+1) / BDIV &
        !               + ((BLAT+3*BDIV/2) - LAT) * vmnrlb(1:JPCHAN,JSCAN,sband  ) / BDIV
        SCORR = (LAT - (BLAT + BDIV / 2)) * vmnrlb(JSCAN, sband + 1) / BDIV &
                + ((BLAT + 3 * BDIV / 2) - LAT) * vmnrlb(JSCAN, sband) / BDIV
      ELSE
        !SCORR(1:JPCHAN) = vmnrlb(1:JPCHAN,JSCAN,sband)
        SCORR = vmnrlb(JSCAN, sband)
      END IF

    ELSEIF (LAT < BLAT + BDIV / 2) THEN

      IF (sband > 1) THEN
        !SCORR(1:JPCHAN) = (LAT - (BLAT-BDIV/2)) * vmnrlb(1:JPCHAN,JSCAN,sband  ) / BDIV &
        !               + ((BLAT+BDIV/2) - LAT) * vmnrlb(1:JPCHAN,JSCAN,sband-1) / BDIV
        !ELSE
        !SCORR(1:JPCHAN) = vmnrlb(1:JPCHAN,JSCAN,sband)
        SCORR = (LAT - (BLAT - BDIV / 2)) * vmnrlb(JSCAN, sband) / BDIV &
                + ((BLAT + BDIV / 2) - LAT) * vmnrlb(JSCAN, sband - 1) / BDIV
      ELSE
        SCORR = vmnrlb(JSCAN, sband)
      END IF

    END IF

  END SUBROUTINE GETSCORR

  SUBROUTINE Biasprep(X, PassDomainCheck, pres, t, q, tskin, pres_BC, t_BC, q_BC, tskin_BC, angles_BC, &
    ca_mean, tb_obs, ca_mean_BC, tb_obs_BC, num_rad, TotalObs, npred, &
    njplev, olatlon, obsHght, obsTime, Angles, norm_pred, zpred)
    ! USE ObsGIIRS_m, ONLY: ObsGIIRS_t
    IMPLICIT NONE
    ! ACCUMULATE STATS FROM FEEDBACK FILES FOR BIAS CORRECTION PROGS
    !--------------------------------------------------------------------------
    TYPE(State_t), INTENT(IN)  :: X
    LOGICAL, INTENT(IN)      :: PassDomainCheck(:)
    REAL(r_kind), INTENT(IN)      :: pres(:, :), t(:, :), q(:, :), tskin(:)
    REAL(r_kind), INTENT(IN)      :: pres_BC(:, :), t_BC(:, :), q_BC(:, :), tskin_BC(:)
    REAL(r_kind), INTENT(IN), OPTIONAL   :: ca_mean(:), tb_obs(:), ca_mean_BC(:), tb_obs_BC(:)
    INTEGER(i_kind), INTENT(IN)  :: num_rad, TotalObs, npred, njplev
    REAL(r_kind), INTENT(IN)  :: olatlon(:, :), obsHght(:)
    INTEGER(i_kind), INTENT(IN)  :: obsTime(:)
    REAL(r_kind), INTENT(IN)  :: angles_BC(:, :), Angles(:,:)
    LOGICAL, INTENT(IN)  :: norm_pred
    REAL(r_kind), INTENT(OUT) :: zpred(num_rad, npred)
    !---------------------------------------------------------------------------
    ! INITIALIZE CONSTANTS AND VARIABLES.
    INTEGER            :: i, j, jx, jch, ivar, nvalid, num_obs_total
    REAL(r_kind), ALLOCATABLE :: zpred_BC(:, :), pred_BC(:)
    REAL(r_kind), ALLOCATABLE   :: pred(:), mean_pred(:), std_pred(:)
    INTEGER, ALLOCATABLE :: valid_num(:)

    INCLUDE "mpif.h"

    ALLOCATE (pred(1:npred))
    ALLOCATE (mean_pred(npred), std_pred(npred))
    ALLOCATE (zpred_BC(TotalObs, npred))
    ALLOCATE (pred_BC(1:npred))
    ALLOCATE (valid_num(1:npred))
    mean_pred = ZERO
    std_pred = ZERO
    zpred = ZERO
    zpred_BC = ZERO
    pred = ZERO
    pred_BC = ZERO

    WRITE (*, 1) num_rad
1   FORMAT('Biasprep - num_rad: ', I8)

    nvalid = 0
    DO j = 1, num_rad
      ! IF (MAXVAL(pres(j, :)) - MINVAL(pres(j, :)) < 1.0D0) THEN
      !   ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
      !   CYCLE
      ! END IF
      IF (PassDomainCheck(j)) THEN
      ! PRINT *, 'Check pres before calling ', maxval(pres), minval(pres), maxval(pres(j,:)), minval(pres(j,:))
        IF (npred .EQ. 3) THEN
          CALL PRED_CONV_CMAGFS(pred(1:npred), npred, njplev, t(j, 1:njplev), q(j, 1:njplev), tskin(j), Angles(1, j), pres(j, 1:njplev)) ! For zhanghua's predictors and coeffs
        ELSE
          IF (PRESENT(ca_mean) .AND. PRESENT(tb_obs)) THEN
            CALL calc_predictors(pred(1:npred), npred, njplev, &
                                t(j, 1:njplev), q(j, 1:njplev), tskin(j), Angles(1, j), pres(j, 1:njplev), &
                                ca_mean=ca_mean(j), tb_obs=tb_obs(j))
          ELSE
            CALL calc_predictors(pred(1:npred), npred, njplev, &
                                t(j, 1:njplev), q(j, 1:njplev), tskin(j), Angles(1, j), pres(j, 1:njplev))
          END IF
        END IF
        zpred(j, 1:npred) = pred(1:npred)
        nvalid = nvalid + 1
      END IF
    END DO

    valid_num = ZERO
    DO j = 1, TotalObs
      IF (npred .EQ. 3) THEN
        CALL PRED_CONV_CMAGFS(pred_BC(1:npred), npred, njplev, t_BC(j, 1:njplev), q_BC(j, 1:njplev), tskin_BC(j), angles_BC(1, j), pres_BC(j, 1:njplev))
      ELSE
        IF (PRESENT(ca_mean) .AND. PRESENT(tb_obs)) THEN
          CALL calc_predictors(pred_BC(1:npred), npred, njplev, &
                              t_BC(j, 1:njplev), q_BC(j, 1:njplev), tskin_BC(j), angles_BC(1, j), pres_BC(j, 1:njplev), &
                              ca_mean=ca_mean_BC(j), tb_obs=tb_obs_BC(j))
        ELSE
          CALL calc_predictors(pred_BC(1:npred), npred, njplev, &
                              t_BC(j, 1:njplev), q_BC(j, 1:njplev), tskin_BC(j), angles_BC(1, j), pres_BC(j, 1:njplev))
        END IF
      END IF
      zpred_BC(j, 1:npred) = pred_BC(1:npred)
      DO i = 1, npred
        IF (ABS(zpred_BC(j,i) - missing) > 1.0) THEN
          valid_num(i) = valid_num(i) + 1
          mean_pred(i) = mean_pred(i) + zpred_BC(j, i)
        END IF
      END DO
    END DO

    ! calculate mean and std of predictors
    mean_pred = mean_pred / REAL(valid_num)
    DO j = 1, TotalObs
      DO i = 1, npred
        IF (ABS(zpred_BC(j,i) - missing) > 1.0) THEN
          std_pred(i) = std_pred(i) + (zpred_BC(j, i) - mean_pred(i))**2
        END IF
      END DO
    END DO
    std_pred = SQRT(std_pred / REAL(valid_num - 1))
    ! PRINT *, 'valid_num = ', valid_num

    mean_pred = (anint(mean_pred*1000))/1000
    std_pred = (anint(std_pred*1000))/1000
    ! PRINT *, 'mean_pred = ', mean_pred
    ! PRINT *, 'std_pred = ', std_pred

    ! For parallel testing only
    ! PRINT *, 'check mean_pred: ', mean_pred
    ! PRINT *, 'check std_pred: ', std_pred
    ! PRINT *, 'valided number: ', nvalid, ' out of ', num_rad
    ! PRINT *, 'max pres: ', MAXVAL(pres(:,1)),MAXVAL(t(:,2)),MAXVAL(q(:,3))
    ! PRINT *, 'max zpred: ', MAXVAL(zpred(:,1)),MAXVAL(zpred(:,2)),MAXVAL(zpred(:,3))
    ! PRINT *, 'min zpred: ', MINVAL(zpred(:,1)),MINVAL(zpred(:,2)),MINVAL(zpred(:,3))

    ! FOR AGRI, we use the coefficients from normalized predictors
    ! FOR GIIRS, predictors use their own units
    IF (norm_pred) THEN
      DO j = 1, num_rad

        ! IF (MAXVAL(pres(j, :)) - MINVAL(pres(j, :)) < 1.0D0) THEN
        !   ! PRINT *, 'The pressure profile is invalid; SKIP NOW'
        !   CYCLE
        ! END IF
        IF (PassDomainCheck(j)) THEN
          DO i = 2, npred
            IF (ABS(zpred(j,i) - missing) > 1.0) THEN
              zpred(j,i) = ( zpred(j,i) - mean_pred(i) ) / std_pred(i)
            ELSE
              zpred(j,i) = 0.0D0
            END IF 
          END DO
        END IF

      END DO
    END IF
    ! PRINT *, 'Check normalization of mean_pred: ', mean_pred
    ! PRINT *, 'Check normalization of std_pred: ', std_pred

    DEALLOCATE (pred, mean_pred, std_pred, zpred_BC, pred_BC, valid_num)

  END SUBROUTINE Biasprep

  SUBROUTINE calc_predictors(pred,npred,nlevs,temp,hum,tskin,satzen,pres,ca_mean,tb_obs)
  IMPLICIT NONE   
  ! pred(1) - gas constant (1)
  ! pred(2) - 1000-300 hPa thickness
  ! pred(3) - 200-50 hPa thickness
  ! pred(4) - tskin
  ! pred(5) - total column precipitable water
  ! pred(6) - satzen
  ! pred(7) - satzen ** 2
  ! pred(8) - satzen ** 3
  
    INTEGER, INTENT(IN)  :: npred, nlevs
    REAL(r_kind), INTENT(IN)     :: temp(nlevs), hum(nlevs), pres(nlevs)
    REAL(r_kind), INTENT(IN), OPTIONAL     :: satzen, tskin, ca_mean, tb_obs
    REAL(r_kind), INTENT(OUT)    :: pred(npred)
    REAL(r_kind) :: layer_base, layer_top, layer_thickness

    pred(1) = 1.0D0

    ! pred(2) - 1000-300 hPa thickness
    layer_base = 1000.0
    layer_top = 300.0
    layer_thickness = Calc_Pred_thickness(nlevs, pres, temp, hum, layer_base, layer_top)
    pred(2) = layer_thickness

    ! pred(3) - 200-50 hPa thickness
    layer_base = 200.0
    layer_top = 50.0
    layer_thickness = Calc_Pred_thickness(nlevs, pres, temp, hum, layer_base, layer_top)
    pred(3) = layer_thickness

    ! ! pred(3) - cos(satzen) (wuyl comment: satzen has already been in the radian unit)
    ! pred(3) = cos(satzen)

    ! pred(4) - tskin
    pred(4) = tskin

    ! pred(5) - total column precipitable water
    pred(5) = Calc_Pred_qv2tpw(nlevs, pres, hum)

    ! pred(6) - satzen
    pred(6) = Calc_Pred_ScanAngle(satzen, 1)

    ! pred(7) - satzen ** 2
    pred(7) =  Calc_Pred_ScanAngle(satzen, 2)

    ! pred(8) - satzen ** 3
    pred(8) =  Calc_Pred_ScanAngle(satzen, 3)

    ! pred(9) - ca_mean
    IF (PRESENT(ca_mean)) pred(9) =  ca_mean

    ! pred(10) - tb_obs
    IF (PRESENT(tb_obs)) pred(10) =  tb_obs
  
  END SUBROUTINE calc_predictors

  SUBROUTINE PRED_CONV_CMAGFS(pred, npred, JPRTLEV, temp, hum, tskin, satzen, pres)
    IMPLICIT NONE

    ! temp - model level temperatures
    ! hum  - model level moistures
    ! T-skin - model skin temperature

    ! pred(1) - 1000-300 hPa thickness
    ! pred(2) - 200-50 hPa thickness
    ! pred(3) - cos(satzen)

    REAL, PARAMETER :: Kth = 287.0597 * 0.5 / 9.806650
    REAL, PARAMETER :: Kpc = 100.0 * 0.5 / 9.806650

    INTEGER, INTENT(IN)  :: npred, JPRTLEV
    REAL(r_kind), INTENT(IN)     :: temp(JPRTLEV), hum(JPRTLEV), pres(JPRTLEV)
    REAL(r_kind), INTENT(IN)     :: satzen, tskin
    REAL(r_kind), INTENT(OUT)    :: pred(npred)

    LOGICAL, SAVE :: FIRST = .TRUE.

    REAL, ALLOCATABLE   :: tv(:)
    REAL, ALLOCATABLE   :: DLP(:)
    REAL, ALLOCATABLE   :: DP(:)

    INTEGER :: IST, ILAY, I

    REAL, PARAMETER :: p300 = 30000.0
    REAL, PARAMETER :: p1000 = 100000.0
    REAL, PARAMETER :: p50 = 5000.0
    REAL, PARAMETER :: p200 = 20000.0
    REAL, PARAMETER :: p10 = 1000.0
    REAL, PARAMETER :: p2 = 200.0
    REAL, PARAMETER :: p1 = 100.0
    REAL, PARAMETER :: p5 = 500.0D0
    REAL, PARAMETER :: pi = 3.1415926D0

    INTEGER :: ip300, ip1000, ip50, ip200, ip10, ip2, ip1, ip5
    INTEGER :: lloc(1)

    ! PRINT *, 'Look for pres levels from ', maxval(pres), ' to ', minval(pres)
    ! Locating the number of level for specifed pressure level
    lloc = MINLOC(ABS(pres - p300))
    ip300 = lloc(1)
    lloc = MINLOC(ABS(pres - p1000))
    ip1000 = lloc(1)
    lloc = MINLOC(ABS(pres - p50))
    ip50 = lloc(1)
    lloc = MINLOC(ABS(pres - p200))
    ip200 = lloc(1)
    lloc = MINLOC(ABS(pres - p10))
    ip10 = lloc(1)
    lloc = MINLOC(ABS(pres - p2))
    ip2 = lloc(1)
    lloc = MINLOC(ABS(pres - P1))
    ip1 = lloc(1)
    lloc = MINLOC(ABS(pres - P5))
    ip5 = lloc(1)

    ! If first time compute log P(i)/P(i-1) and P(i)-P(i-1)
    ALLOCATE (DLP(JPRTLEV - 1))
    ALLOCATE (DP(JPRTLEV - 1))
    ALLOCATE (tv(JPRTLEV))
    !  IF (FIRST) THEN
    ! DLP(1:JPRTLEV-1) = LOG(pres(2:JPRTLEV)) - LOG(pres(1:JPRTLEV-1))
    ! DP(1:JPRTLEV-1) = pres(2:JPRTLEV) - pres(1:JPRTLEV-1)
    DLP(1:JPRTLEV - 1) = LOG(pres(1:JPRTLEV - 1)) - LOG(pres(2:JPRTLEV))
    DP(1:JPRTLEV - 1) = pres(1:JPRTLEV - 1) - pres(2:JPRTLEV)

    ! 1.0 Convert all temperatures to virtual temperatures
    ! ----------------------------------------------------
    tv = temp / (1.0 - 0.6 * hum)

    ! 2.0 Construct averages for NESDIS thick layers
    ! ----------------------------------------------
    ! print *,'PRED_CONV:',SIZE(tv),ip300,ip1000
    IF (ip300 .LE. 1) THEN
      WRITE (*, 1) ip300, ip1000, MINVAL(pres), MAXVAL(pres)
1     FORMAT('PRED_CONV - invalid ip300/ip1000 value: ', 2I5, &
             ' Pressure max-min: ', 2D12.4, ' check pressure and rerun!', /,/)
      STOP
    END IF
    ! The original pred(1:2) assume the pressure is arranged from top to bottom: e.g., ip300 < ip1000
    ! the integration formula is wrong by off 2 times, as the middle point formula
    ! pred(1) = Kth*SUM((tv(ip300:ip1000)+tv(ip300-1:ip1000-1))*DLP(ip300-1:ip1000-1))
    ! pred(2) = Kth*SUM((tv(ip50:ip200)+tv(ip50-1:ip200-1))*DLP(ip50-1:ip200-1))
    ! Yuanfu Xie changed these pred(1:2) calculation assuming pressure arranged from bottom to top
    pred(1) = 0.5D0 * Kth * SUM((tv(ip1000:ip300 - 1) + tv(ip1000 + 1:ip300)) * DLP(ip1000:ip300 - 1))
    pred(2) = 0.5D0 * Kth * SUM((tv(ip200:ip50 - 1) + tv(ip200 + 1:ip50)) * DLP(ip200:ip50 - 1))

    !   sa = (scan_angle*degree2radian)**order ?
    !pred(3) = cos(satzen*pi/180.)
    pred(3) = COS(satzen)

    DEALLOCATE (tv)
    DEALLOCATE (DP)
    DEALLOCATE (DLP)

    RETURN

  END SUBROUTINE PRED_CONV_CMAGFS

END MODULE BcPredictors_m
