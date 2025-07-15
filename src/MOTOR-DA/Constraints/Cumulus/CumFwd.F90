!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

MODULE CumFwd_m

  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE Interp1D_m, ONLY: Interp1D
  USE Tendency_m, ONLY: Tendency
  USE divergence_m, ONLY: divergence
  USE parameters_m, ONLY: g, cp, Lv, k_d, r_k_d, epsw_r, gamma_d, r_d, TK

CONTAINS

  SUBROUTINE CumFwd(uwnd, vwnd, wwnd, pres, &
                    theta, qvapor, X, ztop, &
                    par_theta, par_q, precip)

    IMPLICIT NONE
    TYPE(State_t), INTENT(IN) :: X

    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: uwnd, vwnd, wwnd, theta, &
                                                 qvapor, pres
    REAL(r_kind), DIMENSION(:, :), INTENT(INOUT) :: par_theta, par_q
    REAL(r_kind), DIMENSION(:), INTENT(INOUT) :: precip
    REAL(r_kind), INTENT(IN) :: ztop

    INTEGER(i_kind) :: i, k, k_ref, k_bot, &
                       i_e, i_w, i_n, i_s, keps, kgph, size_3m, vLevel1d

    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: rhod, qu, qv, qw, temp

    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: uwnd1d, vwnd1d, wwnd1d, theta1d, temp1d, &
                                               qvapor1d, pres1d, gph1d, C_theta, C_q, &
                                               Tlcl, theta_e, temp_c, theta_c, q_c, dtheta, &
                                               theta_down, beta_theta, &
                                               theta_con, Q_cont, Q_con, &
                                               divq, psat, qsat, act_Ht, act_Ht_m, &
                                               dtheta_dz, dtheta_m, gph1dm, gph1dm2, dtheta_m2, dtheta_dz2, &
                                               Pm, Pm2t, Pm2, act_dthz_ct, act_dth_ct, act_dthz_lfc, act_dth_lfc, &
                                               act_p500_ct, act_p500_lfc, Act_ct, Act_lfc, dgphm, &
                                               Act_lfc2, Act_Pt, act_cin, act_dth_cin, z_u, act_zu, act_dth_cape, &
                                               z_u2, act_zu2, act_dth_D, act_dth_m, qu_e, qu_w, qv_n, qv_s, qw_c, &
                                               theta1d_e, theta1d_w, theta1d_n, theta1d_s, q1d_n, q1d_s, q1d_e, q1d_w, &
                                               theta1d_c, uwnd1d_c, vwnd1d_c, wwnd1d_c, q1d_c, d_theta_con, &
                                               d_theta_con_m, Q_con_m, qu_c, qv_c, rhod1d

    REAL(r_kind) :: A, qvapor1dmean, pres1dmean, gph1dmean, &
                    temp1dmean, Tb, Tb1, Hb, tempt1, tempt2, &
                    est1, est2, gamma_m1, gamma_m2, cape, cin, avg_dtheta, &
                    std_dtheta, max_dtheta, D, &
                    act_dz, act_dthetae, act_P500, lfct, Htt, lfc, Ht, Pt, Act_sum_ct, &
                    act_cape, Tlcl1, fTs, dfTs, &
                    RH, act_RH, b_param, Ir, act_Ir, &
                    sum_act_zu2, sum_act_zu2_p2, dth_m, dth_std, q1d_sum, qsat_sum, &
                    C_star, d_C_star, d_3m

    ASSOCIATE (cell_type => X%sg%cell_type, &
               cell_stcl => X%sg%cell_stcl, &
               cell_dist => X%sg%cell_dist, &
               vLevel => X%sg%vLevel, &
               num_cell => X%sg%num_cell, &
               num_icell => X%sg%num_icell, &
               gph => X%sg%zHght, &
               sigma => X%sg%sigma, &
               topo => X%sg%topo &
               )

      vLevel1d = vLevel - 1

      size_3m = 10000

      ALLOCATE (gph1dm2(size_3m), dtheta_m2(size_3m), dtheta_dz2(size_3m), &
                Pm2t(size_3m), Pm2(size_3m), act_dthz_ct(size_3m), &
                act_dth_ct(size_3m), act_dthz_lfc(size_3m), &
                act_dth_lfc(size_3m), &
                act_p500_ct(size_3m), &
                act_p500_lfc(size_3m), Act_ct(size_3m), Act_lfc(size_3m), &
                dgphm(size_3m - 1), Act_Pt(size_3m), Act_lfc2(size_3m))

      ALLOCATE (rhod(vLevel, num_cell), qu(vLevel, num_cell), &
                qv(vLevel, num_cell), qw(vLevel, num_cell), temp(vLevel, num_cell), &
                uwnd1d(vLevel1d), vwnd1d(vLevel1d), wwnd1d(vLevel1d), &
                theta1d(vLevel1d), temp1d(vLevel1d), qvapor1d(vLevel1d), &
                pres1d(vLevel1d), gph1d(vLevel1d), C_theta(vLevel), C_q(vLevel), &
                Tlcl(vLevel1d), beta_theta(vLevel1d), &
                theta_e(vLevel1d), temp_c(vLevel1d), theta_c(vLevel1d), q_c(vLevel1d), &
                dtheta(vLevel1d), theta_down(vLevel1d), &
                theta_con(vLevel1d), Q_con(vLevel1d), Q_cont(vLevel1d), &
                gph1dm(vLevel1d - 1), Pm(vLevel1d - 1), &
                psat(vLevel1d), qsat(vLevel1d), act_Ht(vLevel1d), &
                dtheta_dz(vLevel1d - 1), dtheta_m(vLevel1d - 1), act_cin(vLevel1d), act_dth_cin(vLevel1d), &
                z_u(vLevel1d), act_zu(vLevel1d), act_dth_cape(vLevel1d), z_u2(vLevel1d), act_zu2(vLevel1d), &
                act_dth_D(vLevel1d), act_dth_m(vLevel1d), divq(vLevel), &
                qu_e(vLevel), qu_w(vLevel), qv_s(vLevel), qv_n(vLevel), qw_c(vLevel), &
                theta1d_e(vLevel), theta1d_w(vLevel), theta1d_n(vLevel), theta1d_s(vLevel), theta1d_c(vLevel), &
                q1d_n(vLevel), q1d_s(vLevel), q1d_e(vLevel), q1d_w(vLevel), q1d_c(vLevel), &
                uwnd1d_c(vLevel), vwnd1d_c(vLevel), wwnd1d_c(vLevel), d_theta_con(vLevel1d), &
                act_Ht_m(vLevel1d - 1), d_theta_con_m(vLevel1d - 1), Q_con_m(vlevel1d - 1), &
                qu_c(vLevel), qv_c(vLevel), rhod1d(vLevel1d))

      A = 2.4304D02 - TK
      temp = theta * (pres / 1.0D05)**k_d
      rhod = pres / temp / r_d / (1 + qvapor)
      qu = rhod * qvapor * uwnd
      qv = rhod * qvapor * vwnd
      qw = rhod * qvapor * wwnd

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          uwnd1d = uwnd(2:vLevel, i)
          vwnd1d = vwnd(2:vLevel, i)
          wwnd1d = wwnd(2:vLevel, i)
          pres1d = pres(2:vLevel, i) / 1.0D02
          theta1d = theta(2:vLevel, i)
          temp1d = temp(2:vLevel, i)
          qvapor1d = qvapor(2:vLevel, i)
          gph1d = gph(2:vLevel, i)
          rhod1d = rhod(2:vLevel, i)

          gph1dm = (gph1d(2:vLevel1d) + gph1d(1:vLevel1d - 1)) / 2.0D0

          d_3m = (gph1dm(vLevel1d - 1) - gph1dm(1)) / size_3m

          DO kgph = 1, size_3m
            gph1dm2(kgph) = gph1dm(1) + d_3m * (kgph - 1)
          END DO

          C_theta = 0.0D0
          C_q = 0.0D0

          Tlcl = temp1d
          DO k = 1, vLevel1d
            DO keps = 1, 5
              fTs = (Tlcl(k) + A) / k_d * dlog(Tlcl(k)) &
                    + Tlcl(k) * (dlog(qvapor1d(k)) + dlog(pres1d(k)) - 1.0D0 / k_d * dlog(temp1d(k)) &
                                 - dlog((epsw_r + qvapor1d(k)) * 6.1094D0) - 17.625D0) &
                    + A * (dlog(qvapor1d(k)) + dlog(pres1d(k)) - 1.0D0 / k_d * dlog(temp1d(k)) &
                           - dlog((epsw_r + qvapor1d(k)) * 6.1094D0)) + 17.625D0 * TK
              dfTs = dlog(Tlcl(k)) / k_d + (Tlcl(k) + A) / k_d / Tlcl(k) &
                     + (dlog(qvapor1d(k)) + dlog(pres1d(k)) &
                        - 1.0D0 / k_d * dlog(temp1d(k)) - dlog((epsw_r + qvapor1d(k)) * 6.1094D0) - 17.625D0)
              Tlcl1 = Tlcl(k) - fTs / dfTs
              Tlcl(k) = Tlcl1
            END DO
          END DO

          theta_e = theta1d * dexp(Lv * qvapor1d / cp / Tlcl)
          k = 1
          k_ref = vLevel / 2
          act_dthetae = (TANH((theta_e(1) - theta_e(k_ref) - 1.0D0) * 2.0D0) + 1.0D0) / 2.0D0 ! The first activation function
          qvapor1dmean = SUM(qvapor1d(1:3)) / 3.0D0
          pres1dmean = SUM(pres1d(1:3)) / 3.0D0
          temp1dmean = SUM(temp1d(1:3)) / 3.0D0
          gph1dmean = SUM(gph1d(1:3)) / 3.0D0

          Tb = temp1dmean
          DO keps = 1, 5
            fTs = (Tb + A) * r_k_d * dlog(Tb) &
                  + Tb * (dlog(qvapor1dmean) + dlog(pres1dmean) - r_k_d * dlog(temp1dmean) - dlog((epsw_r + qvapor1dmean) * 6.1094D0) - 1.7625D01) &
                  + A * (dlog(qvapor1dmean) + dlog(pres1dmean) - r_k_d * dlog(temp1dmean) - dlog((epsw_r + qvapor1dmean) * 6.1094D0)) + 1.7625D01 * TK
            dfTs = dlog(Tb) * r_k_d + (Tb + A) * r_k_d / Tb &
                   + (dlog(qvapor1dmean) + dlog(pres1dmean) - r_k_d * dlog(temp1dmean) - dlog((epsw_r + qvapor1dmean) * 6.1094D0) - 1.7625D01)
            Tb1 = Tb - fTs / dfTs
            Tb = Tb1
          END DO

          Hb = (temp1dmean - Tb) / g * cp + gph1dmean
          k_bot = MINLOC(dabs(gph1d - Hb), 1) ! Cloud base found

          q_c = qvapor1d
          theta_c(1:k_bot) = theta1d(1)
          temp_c(k_bot) = Tb
          DO k = k_bot, vLevel1d - 1
            tempt1 = temp_c(k)
            IF (tempt1 .LE. 3.011D01) THEN
              gamma_m1 = -gamma_d
              gamma_m2 = -gamma_d
              q_c(k) = 0.0D0
            ELSE
              est1 = 6.1094D0 * EXP(1.7625D01 * (tempt1 - TK) / (tempt1 - 3.011D01))
              q_c(k) = est1 / (pres1d(k) - est1) * epsw_r
              gamma_m1 = -gamma_d * ((1.0D0 + Lv * est1 * epsw_r / pres1d(k) / r_d / tempt1) &
                                     / (1.0D0 + Lv**2.0D0 * epsw_r**2.0D0 * est1 / cp / pres1d(k) / r_d / tempt1**2.0D0))
              tempt2 = tempt1 + gamma_m1 * (gph1d(k + 1) - gph1d(k))
              IF (tempt2 .GT. 3.011D01) THEN
                est2 = 6.1094D0 * EXP(1.7625D01 * (tempt2 - TK) / (tempt2 - 3.011D01))
                gamma_m2 = -gamma_d * ((1.0D0 + Lv * est2 * epsw_r / pres1d(k + 1) / r_d / tempt2) &
                                       / (1.0D0 + Lv**2.0D0 * epsw_r**2.0D0 * est2 / cp / pres1d(k + 1) / r_d / tempt2**2.0D0))
              ELSE
                gamma_m2 = -gamma_d
              END IF
              temp_c(k + 1) = tempt1 + (gamma_m1 + gamma_m2) * (gph1d(k + 1) - gph1d(k)) / 2.0D0
              theta_c(k + 1) = temp_c(k + 1) * (1.0D03 / pres1d(k + 1))**k_d
            END IF
          END DO

          dtheta = theta_c - theta1d
          dtheta_dz = (dtheta(2:vLevel1d) - dtheta(1:vLevel1d - 1)) &
                      / (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))
          dtheta_m = (dtheta(2:vLevel1d) + dtheta(1:vLevel1d - 1)) / 2.0D0
          Pm = (dlog(pres1d(2:vLevel1d)) + dlog(pres1d(1:vLevel1d - 1))) / 2.0D0

          CALL Interp1D(gph1dm, dtheta_m, gph1dm2, dtheta_m2)
          CALL Interp1D(gph1dm, dtheta_dz, gph1dm2, dtheta_dz2)
          CALL Interp1D(gph1dm, Pm, gph1dm2, Pm2t)

          Pm2 = EXP(pm2t)
          act_p500_ct = (TANH((5.0D02 - Pm2) / 8.0D01) + 1.0D0) / 2.0D0
          act_p500_lfc = (TANH((-5.0D02 + Pm2) / 8.0D01) + 1.0D0) / 2.0D0

          act_dthz_ct = (-TANH(dtheta_dz2 * 3.0D04) + 1.0D0) / 2.0D0
          act_dth_ct = 1.0D0 / (1.0D0 + 5.0D03 * dtheta_m2**4.0D0)
          act_dthz_lfc = (TANH(dtheta_dz2 * 1.0D04) + 1.0D0) / 2.0D0
          act_dth_lfc = 1.0D0 / (1.0D0 + 5.0D03 * dtheta_m2**4.0D0)
          Act_ct = act_p500_ct * act_dthz_ct * act_dth_ct

          Act_sum_ct = (1.0D0 / (1.0D0 + EXP(-1.0D01 * SUM(Act_ct)))) * 2.0D0 - 1.0D0
          Act_lfc = act_p500_lfc * act_dthz_lfc * act_dth_lfc
          dgphm = (gph1dm2(2:size_3m) - gph1dm2(1:size_3m - 1))
          Htt = Act_sum_ct * SUM((Act_ct(1:size_3m - 1) * gph1dm2(1:size_3m - 1) &
                                  + Act_ct(2:size_3m) * gph1dm2(2:size_3m)) * dgphm) / &
                SUM((Act_ct(1:size_3m - 1) + Act_ct(2:size_3m)) * dgphm)

          lfct = Act_sum_ct * SUM((Act_lfc(1:size_3m - 1) * gph1dm2(1:size_3m - 1) + &
                                   Act_lfc(2:size_3m) * gph1dm2(2:size_3m)) * dgphm) / &
                 SUM((Act_lfc(1:size_3m - 1) + Act_lfc(2:size_3m)) * dgphm)
          Act_Pt = 1.0D0 / (1.0D0 + (gph1dm2 - Htt)**4.0D0 / 1.0D08)
          Act_lfc2 = 1.0D0 / (1.0D0 + (gph1dm2 - lfct)**4.0D0 / 1.0D08)
          Pt = SUM((Act_Pt(1:size_3m - 1) * Pm2(1:size_3m - 1) + &
                    Act_Pt(2:size_3m) * Pm2(2:size_3m)) * dgphm) / &
               SUM((Act_Pt(1:size_3m - 1) + Act_Pt(2:size_3m)) * dgphm)
          Ht = SUM((Act_Pt(1:size_3m - 1) * gph1dm2(1:size_3m - 1) + &
                    Act_Pt(2:size_3m) * gph1dm2(2:size_3m)) * dgphm) / &
               SUM((Act_Pt(1:size_3m - 1) + Act_Pt(2:size_3m)) * dgphm)
          lfc = SUM((Act_lfc2(1:size_3m - 1) * gph1dm2(1:size_3m - 1) + &
                     Act_lfc2(2:size_3m) * gph1dm2(2:size_3m)) * dgphm) / &
                SUM((Act_lfc2(1:size_3m - 1) + Act_lfc2(2:size_3m)) * dgphm)

          act_dz = (TANH((Ht - lfc - 3.0D03) / 5.0D02) + 1.0D0) / 2.0D0
          act_P500 = (TANH((5.0D02 - Pt) / 5.0D01) + 1.0D0) / 2.0D0 ! Third activation function
          act_cin = (TANH((lfc - gph1d) / 5.0D02) + 1.0D0) / 2.0D0
          act_dth_cin = act_cin * dtheta

          cin = SUM((act_dth_cin(2:vLevel1d) + act_dth_cin(1:vLevel1d - 1)) &
                    * (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))) / 2.0D0
          z_u = (2.0D0 * gph1d - Ht - lfc) / (Ht - lfc)

          act_zu = 1.0D0 / (1.0D0 + z_u**8.0D0)
          act_dth_cape = act_zu * dtheta
          cape = SUM((act_dth_cape(2:vLevel1d) + act_dth_cape(1:vLevel1d - 1)) &
                     * (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))) / 2.0D0

          act_cape = (TANH((cape + 1.0D0 * cin) / 1.0D03 - 2.5D0) + 1.0D0) / 2.0D0 ! Fourth activation function

          z_u2 = (2.0D0 * gph1d - Ht - Hb) / (Ht - Hb)

          act_zu2 = 1.0D0 / (1.0D0 + z_u2**8.0D0)

          sum_act_zu2 = SUM(act_zu2)
          sum_act_zu2_p2 = SUM(act_zu2**2.0D0)
          act_dth_D = act_zu2 * dtheta
          dth_m = SUM(act_dth_D) / sum_act_zu2
          act_dth_m = dth_m * act_zu2
          dth_std = SQRT(SUM((act_dth_D - act_dth_m)**2.0D0) / sum_act_zu2_p2)
          max_dtheta = dth_m + 1.5D0 * dth_std

          D = 1.0D0 / (1.0D0 + EXP(-max_dtheta + 2.5D0)) * 5.0D0

          theta_down = -EXP(-4.0D0 * SQRT(gph1d - gph1d(1)) / SQRT(Ht)) * D + theta1d

          beta_theta = EXP(-5.0D0 * SQRT(gph1d - gph1d(1)) / SQRT(Ht))
          theta_con = beta_theta * theta_down + (1.0D0 - beta_theta) * theta_c

          Q_cont = q_c - qvapor1d
          Q_con = Q_cont * act_zu2

          i_e = cell_stcl(6, i)
          i_w = cell_stcl(4, i)
          i_n = cell_stcl(8, i)
          i_s = cell_stcl(2, i)
          qu_c = qu(:, i)
          qu_e = qu(:, i_e)
          qu_w = qu(:, i_w)
          qv_c = qv(:, i)
          qv_n = qv(:, i_n)
          qv_s = qv(:, i_s)
          qw_c = qw(:, i)
          CALL divergence(qu_c, qu_w, qu_e, &
                          qv_c, qv_s, qv_n, &
                          qw_c, sigma, ztop, &
                          topo(i), topo(i_w), &
                          topo(i_e), topo(i_s), &
                          topo(i_n), cell_dist(4, i), &
                          cell_dist(2, i), divq)

          act_Ht = (TANH((-gph1d + Ht) / 1.0D02) + 1.0D0) / 2.0D0
          Ir = -SUM((divq(2:vLevel - 1) * act_Ht(1:vLevel1d - 1) + &
                     divq(3:vLevel) * act_Ht(2:vLevel1d)) * &
                    (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))) / 2.0D0

          psat = 6.1094D0 * EXP(1.7625D01 * (temp1d - TK) / (temp1d - 3.011D01))
          qsat = psat / (pres1d - psat) * epsw_r

          q1d_sum = SUM(qvapor1d * act_Ht)
          qsat_sum = SUM(qsat * act_Ht)
          RH = q1d_sum / qsat_sum
          act_RH = (TANH((RH - 6.8D-1) * 4.0D01) + 1.0D0) / 2.0D0
          b_param = act_RH * (1.0D0 - RH) * 5.0D-1 + (1.0D0 - act_RH)

          act_Ht_m = (act_Ht(2:vLevel1d) + act_Ht(1:vLevel1d - 1)) / 2.0D0
          d_theta_con = theta_con - theta1d
          d_theta_con_m = (d_theta_con(2:vLevel1d) + d_theta_con(1:vLevel1d - 1)) / 2.0D0
          Q_con_m = (Q_con(2:vLevel1d) + Q_con(1:vLevel1d - 1)) / 2.0D0

          act_Ir = (TANH((Ir * 1.0D06 - 1.5D0) * 2.0D0) + 1.0D0) / 2.0D0
          C_star = act_Ir * act_dthetae * act_dz * act_P500 * act_cape * Ir

          d_C_star = b_param * C_star
          precip(i) = C_star - d_C_star
          C_theta(2:vLevel) = (C_star - d_C_star) / &
                              SUM(act_Ht_m * d_theta_con_m * (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))) &
                              * act_Ht * d_theta_con * &
                              Lv / cp * (1.0D03 / pres1d)**k_d
          C_q(2:vLevel) = d_C_star &
                          / SUM(act_Ht_m * Q_con_m * (gph1d(2:vLevel1d) - gph1d(1:vLevel1d - 1))) &
                          * act_Ht * Q_con / rhod1d

          theta1d_e = theta(:, i_e)
          theta1d_w = theta(:, i_w)
          theta1d_n = theta(:, i_n)
          theta1d_s = theta(:, i_s)
          theta1d_c = theta(:, i)
          uwnd1d_c = uwnd(:, i)
          vwnd1d_c = vwnd(:, i)
          wwnd1d_c = wwnd(:, i)
          CALL Tendency(uwnd1d_c, vwnd1d_c, wwnd1d_c, &
                        theta1d_c, theta1d_w, theta1d_e, theta1d_s, theta1d_n, &
                        sigma, ztop, topo(i), topo(i_w), topo(i_e), &
                        topo(i_s), topo(i_n), cell_dist(4, i), cell_dist(2, i), &
                        C_theta, par_theta(:, i))

          q1d_e = qvapor(:, i_e)
          q1d_w = qvapor(:, i_w)
          q1d_n = qvapor(:, i_n)
          q1d_s = qvapor(:, i_s)
          q1d_c = qvapor(:, i)

          CALL Tendency(uwnd1d_c, vwnd1d_c, wwnd1d_c, &
                        q1d_c, q1d_w, q1d_e, q1d_s, q1d_n, &
                        sigma, ztop, topo(i), topo(i_w), topo(i_e), &
                        topo(i_s), topo(i_n), cell_dist(4, i), &
                        cell_dist(2, i), C_q, par_q(:, i))
        END IF

      END DO

      DEALLOCATE (gph1dm2, dtheta_m2, dtheta_dz2, &
                  Pm2t, Pm2, act_dthz_ct, &
                  act_dth_ct, act_dthz_lfc, &
                  act_dth_lfc, &
                  act_p500_ct, &
                  act_p500_lfc, Act_ct, Act_lfc, &
                  dgphm, Act_Pt, Act_lfc2)

      DEALLOCATE (rhod, qu, &
                  qv, qw, temp, &
                  uwnd1d, vwnd1d, wwnd1d, &
                  theta1d, temp1d, qvapor1d, &
                  pres1d, gph1d, C_theta, C_q, &
                  Tlcl, beta_theta, &
                  theta_e, temp_c, theta_c, q_c, &
                  dtheta, theta_down, &
                  theta_con, Q_con, Q_cont, &
                  gph1dm, Pm, &
                  psat, qsat, act_Ht, &
                  dtheta_dz, dtheta_m, act_cin, act_dth_cin, &
                  z_u, act_zu, act_dth_cape, z_u2, act_zu2, &
                  act_dth_D, act_dth_m, divq, &
                  qu_e, qu_w, qv_s, qv_n, qw_c, &
                  theta1d_e, theta1d_w, theta1d_n, theta1d_s, theta1d_c, &
                  q1d_n, q1d_s, q1d_e, q1d_w, q1d_c, &
                  uwnd1d_c, vwnd1d_c, wwnd1d_c, d_theta_con, &
                  act_Ht_m, d_theta_con_m, Q_con_m, &
                  qu_c, qv_c)

    END ASSOCIATE
  END SUBROUTINE CumFwd
END MODULE CumFwd_m
