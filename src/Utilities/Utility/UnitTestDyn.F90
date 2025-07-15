MODULE UnitTestDyn_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: Omega, EarthRadius, machineEps, pi, g
  USE RossbyHaurwitz_m, ONLY: RossbyHaurwitzSphere1_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  ! Rossby-Haurwitz parameters:
  REAL(r_kind), PARAMETER :: &
    rh_kapa = 7.848D-6, &
    rh_omga = 7.848D-6, &
    rh_high = 8.0D3, &
    rh_bigR = 4.0D0, &
    rh_angl = &
    (rh_bigR * (3.0D0 + rh_bigR) * rh_omga / 1.0D6 - &
     2.0D0 * Omega) / (1.0D0 + rh_bigR) / (2.0D0 + rh_bigR), &
    omega_t = 2 * pi / 36000.0D0

  TYPE :: UnitTestDyn_t
    INTEGER(i_kind), POINTER :: num_cell, vLevel
    REAL(r_kind), POINTER :: sigma(:), z_top, f(:, :)
    REAL(r_kind), ALLOCATABLE :: &
      sigma3d(:, :), zHght(:, :), DM(:, :), Hz(:, :), z_s(:), &
      psi1st(:, :, :), zta1st(:, :, :), & ! 1st order derivatives of psi
      psi1st_dclat(:, :, :), psi1st_dclatp2(:, :, :), psi1st_dclatp3(:, :, :), & !derivatives devided by cos(lat)
      zta1st_dclat(:, :, :), & ! 1st order derivatives devided by cos(lat)
      psi2nd(:, :, :), & ! 2nd order derivatives of psi
      psi2nd_dclat(:, :, :), psi2nd_dclatp2(:, :, :), &
      psiThetaLambda_dclat(:, :), & ! 2nd order of psi_theta_lambda devided by cos(lat)
      psi2ndLambda_dclat(:, :), & ! 2nd order of psi_lambda devided by cos(lat)
      psi3rd(:, :, :), & ! 3rd order derivatives of psi
      psi3rd_dclat(:, :, :), &
      chi1st(:, :, :), del1st(:, :, :), & ! 1st order derivatives of chi
      chi1st_dclat(:, :, :), chi1st_dclatp2(:, :, :), chi1st_dclatp3(:, :, :), &
      del1st_dclat(:, :, :), & ! 1st order derivatives devided by cos(lat)
      chi2nd(:, :, :), & ! 2nd order derivatives of chi
      chi2nd_dclat(:, :, :), chi2nd_dclatp2(:, :, :), &
      chiThetaLambda_dclat(:, :), & ! 2nd order of chi_theta_lambda divided by cos(lat)
      chi2ndLambda_dclat(:, :), & ! 2nd order of chi_lambda devided by cos(lat)
      chi3rd(:, :, :), & ! 3rd order derivatives of chi
      chi3rd_dclat(:, :, :), &
      zs1sttheta(:), & ! 1st order of z_s, just function with lat
      zs1sttheta_dclat(:), & ! 1st order z_s divided by cos(lat)
      zs2ndtheta(:), & ! 2end order of z_s, just function with lat
      zs2ndtheta_dclat(:), & ! 2end order z_s divided by cos(lat)
      z1sttheta(:, :), z1sttheta_dclat(:, :), & ! 1st order derivatives of zhight, just function with lat
      d_z1stthetadclat_dtheta(:, :), & ! 1st order derivative for z1sttheta_dclat, used for DM calculation
      z2ndtheta(:, :), z2ndtheta_dclat(:, :), & ! 2end order derivatives of zhight, just function with lat
      Hz1sttheta(:, :), Hz1sttheta_dclat(:, :), & ! 1st order derivatives of Hz
      DM1st(:, :, :), DM1stLambda_dclat(:, :), & ! 1st order derivatives of DM,
      zfunct(:, :, :), & ! saving the derivative coefs in vertical direction
      pres1sttheta(:, :), pres1sttheta_dclat(:, :), &
      pres2ndtheta(:, :), pres2ndtheta_dclat(:, :), &
      pres1stsigma(:, :), pres2ndsigmatheta(:, :), &
      rho1sttheta(:, :), rho1stsigma(:, :), & !1st order derivatives of rho
      Tp1sttheta(:, :), Tp1stsigma(:, :), & ! 1st order derivatives of potential temperature
      q1sttheta(:, :), q1stsigma(:, :) ! 1st order derivatives of qvapor

  CONTAINS
    PROCEDURE :: initial
    PROCEDURE :: destroy
    PROCEDURE :: RH_psi_func
    PROCEDURE :: RH_chi_func
    PROCEDURE :: Hght_func
    PROCEDURE :: DM_func
    PROCEDURE :: VorTen
    PROCEDURE :: pres_func
    PROCEDURE :: rho_func
    PROCEDURE :: Tp_func
    PROCEDURE :: q_func
    PROCEDURE :: DivTen
    PROCEDURE :: KE_func
    PROCEDURE :: RhoTen
    PROCEDURE :: wTen
    PROCEDURE :: TpTen
    PROCEDURE :: qTen
    PROCEDURE :: righthandsHalf
    PROCEDURE :: righthandsFull
  END TYPE

CONTAINS
  SUBROUTINE initial(this, sg)
    CLASS(UnitTestDyn_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    INTEGER(i_kind), POINTER :: num_cell, vLevel
    INTEGER(i_kind) :: i

    this%num_cell => sg%num_cell
    this%vLevel => sg%vLevel
    this%sigma => sg%sigma
    this%z_top => sg%ztop
    num_cell => sg%num_cell
    vLevel => sg%vLevel
    this%f => sg%f

    ! PRINT *, 'ztop in unitestdyn is: ', this%z_top

    ALLOCATE (this%psi1st(3, vLevel, num_cell), this%psi1st_dclat(2, vLevel, num_cell), &
              this%psi1st_dclatp2(2, vLevel, num_cell), this%psi1st_dclatp3(2, vLevel, num_cell), &
              this%zta1st(3, vLevel, num_cell), this%zta1st_dclat(2, vLevel, num_cell), &
              this%psi2nd(3, vLevel, num_cell), &
              this%psi2nd_dclat(3, vLevel, num_cell), this%psi2nd_dclatp2(3, vLevel, num_cell), &
              this%psiThetaLambda_dclat(vLevel, num_cell), &
              this%psi2ndLambda_dclat(vLevel, num_cell), &
              this%psi3rd(4, vLevel, num_cell), &
              this%psi3rd_dclat(4, vLevel, num_cell), &
              this%chi1st(3, vLevel, num_cell), this%chi1st_dclat(2, vLevel, num_cell), &
              this%chi1st_dclatp2(2, vLevel, num_cell), this%chi1st_dclatp3(2, vLevel, num_cell), &
              this%del1st(3, vLevel, num_cell), this%del1st_dclat(2, vLevel, num_cell), &
              this%chi2nd(3, vLevel, num_cell), &
              this%chi2nd_dclat(3, vLevel, num_cell), this%chi2nd_dclatp2(3, vLevel, num_cell), &
              this%chiThetaLambda_dclat(vLevel, num_cell), &
              this%chi2ndLambda_dclat(vLevel, num_cell), &
              this%chi3rd(4, vLevel, num_cell), &
              this%chi3rd_dclat(4, vLevel, num_cell), &
              this%z1sttheta(vLevel, num_cell), this%z1sttheta_dclat(vLevel, num_cell), &
              this%d_z1stthetadclat_dtheta(vLevel, num_cell), &
              this%z2ndtheta(vLevel, num_cell), this%z2ndtheta_dclat(vLevel, num_cell), &
              this%Hz1sttheta(vLevel, num_cell), this%Hz1sttheta_dclat(vLevel, num_cell), &
              this%DM1st(2, vLevel, num_cell), this%DM1stLambda_dclat(vLevel, num_cell), &
              this%zfunct(2, vLevel, num_cell), &
              this%zs1sttheta(num_cell), this%zs1sttheta_dclat(num_cell), &
              this%zs2ndtheta(num_cell), this%zs2ndtheta_dclat(num_cell), &
              this%pres1sttheta(vLevel, num_cell), this%pres1sttheta_dclat(vLevel, num_cell), &
              this%pres2ndtheta(vLevel, num_cell), this%pres2ndtheta_dclat(vLevel, num_cell), &
              this%pres1stsigma(vLevel, num_cell), &
              this%rho1sttheta(vLevel, num_cell), this%rho1stsigma(vLevel, num_cell), &
              this%pres2ndsigmatheta(vLevel, num_cell), &
              this%Tp1sttheta(vLevel, num_cell), this%Tp1stsigma(vLevel, num_cell), &
              this%q1sttheta(vLevel, num_cell), this%q1stsigma(vLevel, num_cell), &
              this%sigma3d(vLevel, num_cell), this%zHght(vLevel, num_cell), &
              this%DM(vLevel, num_cell), this%Hz(vLevel, num_cell), &
              this%z_s(num_cell))
    DO i = 1, sg%num_cell
      this%sigma3d(:, i) = sg%sigma
    END DO

    NULLIFY (num_cell)
    NULLIFY (vLevel)
  END SUBROUTINE initial

  SUBROUTINE destroy(this)
    CLASS(UnitTestDyn_t) :: this

    IF (ASSOCIATED(this%vLevel)) THEN
      NULLIFY (this%vLevel)
      NULLIFY (this%num_cell)
      NULLIFY (this%sigma)
      NULLIFY (this%z_top)
      NULLIFY (this%f)
    END IF

    IF (ALLOCATED(this%psi1st)) THEN
      DEALLOCATE (this%psi1st, this%zta1st, &
                  this%psi1st_dclat, &
                  this%psi1st_dclatp2, this%psi1st_dclatp3, &
                  this%zta1st_dclat, &
                  this%psi2nd, &
                  this%psi2nd_dclat, this%psi2nd_dclatp2, &
                  this%psiThetaLambda_dclat, &
                  this%psi2ndLambda_dclat, &
                  this%psi3rd, &
                  this%psi3rd_dclat, &
                  this%chi1st, this%del1st, &
                  this%chi1st_dclat, this%del1st_dclat, &
                  this%chi1st_dclatp2, this%chi1st_dclatp3, &
                  this%chi2nd, &
                  this%chi2nd_dclat, this%chi2nd_dclatp2, &
                  this%chiThetaLambda_dclat, &
                  this%chi2ndLambda_dclat, &
                  this%chi3rd, &
                  this%chi3rd_dclat, &
                  this%z1sttheta, this%z1sttheta_dclat, &
                  this%d_z1stthetadclat_dtheta, &
                  this%z2ndtheta, this%z2ndtheta_dclat, &
                  this%Hz1sttheta, this%Hz1sttheta_dclat, &
                  this%DM1st, this%DM1stLambda_dclat, &
                  this%zfunct, &
                  this%zs1sttheta, this%zs1sttheta_dclat, &
                  this%zs2ndtheta, this%zs2ndtheta_dclat, &
                  this%pres1sttheta, this%pres1sttheta_dclat, &
                  this%pres2ndtheta, this%pres2ndtheta_dclat, &
                  this%pres1stsigma, &
                  this%rho1sttheta, this%rho1stsigma, &
                  this%pres2ndsigmatheta, &
                  this%Tp1sttheta, this%Tp1stsigma, &
                  this%q1sttheta, this%q1stsigma, &
                  this%sigma3d, this%zHght, &
                  this%DM, this%Hz, &
                  this%z_s)

    END IF

  END SUBROUTINE destroy

  SUBROUTINE RH_psi_func(this, attime, latlon, zhght, stream, vortct, par_vor_t)
    IMPLICIT NONE

    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  attime, zhght(:, :), latlon(:, :)

    REAL(r_kind), INTENT(OUT) :: stream(:, :), &
                                 vortct(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_vor_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: R0, clat, slat, tlat, clon, slon, clatR(0:5), &
                       szhght, czhght, z_ratio, zfun, zfun2, parzfun_sigma, pi_ratio

    z_ratio = 0.1D0
    pi_ratio = 0.1D0
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        R0 = EarthRadius + zhght(k, i)
        szhght = DSIN(zhght(k, i) / this%z_top * pi_ratio * pi)
        czhght = DCOS(zhght(k, i) / this%z_top * pi_ratio * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = (2.0D0 / R0 + &
                 z_ratio * czhght * pi_ratio * pi / this%z_top / zfun)
        parzfun_sigma = z_ratio * czhght * pi_ratio * pi / this%z_top
        this%zfunct(1, k, i) = zfun
        this%zfunct(2, k, i) = zfun2

        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        tlat = DTAN(latlon(1, i))
        clon = DCOS(rh_bigR * (latlon(2, i) + rh_angl * attime))
        slon = DSIN(rh_bigR * (latlon(2, i) + rh_angl * attime))

        clatR(0) = clat**(rh_bigR - 3.0D0)  ! c**(R-3)
        clatR(1) = clatR(0) * clat          ! c**(R-2)
        clatR(2) = clatR(1) * clat          ! c**(R-1)
        clatR(3) = clatR(2) * clat          ! c** R
        clatR(4) = clatR(3) * clat          ! c**(R+1)
        clatR(5) = clatR(4) * clat          ! c**(R+2)

        ! ! Stream function:
        ! stream(k, i) = &
        !    -R0*R0*zfun*slat*(rh_omga - rh_kapa*clatR(3)*clon)

        ! ! Relative vorticity:
        ! vortct(k, i) = zfun*slat*(2.0D0*rh_omga - rh_kapa*clatR(3)* &
        !                           (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon)

        ! Stream function:
        stream(k, i) = &
          R0 * R0 * zfun * slat * (rh_kapa * clatR(3) * clon)

        ! Relative vorticity:
        vortct(k, i) = -zfun * slat * (rh_kapa * clatR(3) * &
                                       (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon)
        IF (PRESENT(par_vor_t)) THEN
          par_vor_t(k, i) = rh_bigR * rh_angl * zfun * slat * (rh_kapa * clatR(3) * &
                                                               (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon)
        END IF

    !!! Note: these derivatives are derived in the document of
    !!!       Rossby_HaurwitzTest under
    !!!         /Users/xieyuanfu/developments/models/square/doc
    !!!       they need to double check the doc and coding here.
    !!! For testing the Poisson solver on a sphere, I currently
    !!! use the above stream function and vorticity only

        ! The first order derivatives: for the last dimension, psi: (1) sigma, (2) lat, (3) lon; zta: (1) lat, (2) lon, (3)sigma
        ! *_dclat are parameters divived by cos(lat)
        this%psi1st(1, k, i) = zfun2 &
                               * stream(k, i)
        this%psi1st(2, k, i) = R0 * R0 * zfun * &
                               rh_kapa * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * clon
        this%psi1st(3, k, i) = -R0 * R0 * zfun &
                               * rh_kapa * rh_bigR * clatR(3) * slat * slon
        this%psi1st_dclat(1, k, i) = R0 * R0 * zfun * &
                                     (rh_kapa * &
                                      ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * clon)
        this%psi1st_dclat(2, k, i) = -R0 * R0 * zfun &
                                     * rh_kapa * rh_bigR * clatR(2) * slat * slon

        this%psi1st_dclatp2(1, k, i) = R0 * R0 * zfun * &
                                       (rh_kapa * &
                                        ((rh_bigR + 1.0D0) * clatR(2) - rh_bigR * clatR(0)) * clon)
        this%psi1st_dclatp2(2, k, i) = -R0 * R0 * zfun &
                                       * rh_kapa * rh_bigR * clatR(1) * slat * slon

        this%psi1st_dclatp3(1, k, i) = R0 * R0 * zfun * &
                                       (rh_kapa * &
                                        ((rh_bigR + 1.0D0) * clatR(1) - rh_bigR) * clon)
        this%psi1st_dclatp3(2, k, i) = -R0 * R0 * zfun &
                                       * rh_kapa * rh_bigR * clatR(0) * slat * slon

        this%zta1st(1, k, i) = -(rh_kapa * &
                                 ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                                 (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                               * zfun
        this%zta1st(2, k, i) = rh_kapa * rh_bigR * slat * clatR(3) * &
                               (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                               * zfun
        this%zta1st(3, k, i) = parzfun_sigma / zfun * vortct(k, i)

        this%zta1st_dclat(1, k, i) = -(rh_kapa * &
                                       ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * &
                                       (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                                     * zfun
        this%zta1st_dclat(2, k, i) = rh_kapa * rh_bigR * slat * clatR(2) * &
                                     (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                                     * zfun

        ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
        this%psi2nd(1, k, i) = R0 * R0 * slat &
                               * zfun &
                               * rh_kapa * &
                               ((rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                 (rh_bigR + 1.0D0)**2.0D0 * clatR(3))) * clon
        this%psi2nd(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                               zfun * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * slon
        this%psi2nd(3, k, i) = -R0 * R0 * rh_kapa * rh_bigR**2.0D0 * clatR(3) * slat * clon &
                               * zfun

        this%psi2nd_dclat(1, k, i) = R0 * R0 * slat &
                                     * zfun &
                                     * (rh_kapa * &
                                        ((rh_bigR * (rh_bigR - 1.0D0) * clatR(0) - &
                                          (rh_bigR + 1.0D0)**2.0D0 * clatR(2))) * clon)
        this%psi2nd_dclat(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                                     zfun * &
                                     ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon
        this%psi2nd_dclat(3, k, i) = -R0 * R0 * rh_kapa * rh_bigR**2.0D0 * clatR(2) * slat * clon &
                                     * zfun

        this%psi2nd_dclatp2(1, k, i) = R0 * R0 * slat &
                                       * zfun &
                                       * (rh_kapa * &
                                          ((rh_bigR * (rh_bigR - 1.0D0) - &
                                            (rh_bigR + 1.0D0)**2.0D0 * clatR(1))) * clon)
        this%psi2nd_dclatp2(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                                       zfun * &
                                       ((rh_bigR + 1.0D0) * clatR(2) - rh_bigR * clatR(0)) * slon
        this%psi2nd_dclatp2(3, k, i) = -R0 * R0 * rh_kapa * rh_bigR**2.0D0 * clatR(1) * slat * clon &
                                       * zfun

        this%psiThetaLambda_dclat(k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                                          ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon &
                                          * zfun
        this%psi2ndLambda_dclat(k, i) = -R0 * R0 * rh_kapa * rh_bigR**2.0D0 * clatR(2) * slat * clon &
                                        * zfun

        ! the third order derivatives of psi: (1) lat, (2)latp2_lon, (3)lonp2_lat, (4)lon
        this%psi3rd(1, k, i) = R0 * R0 * rh_kapa * (rh_bigR * (rh_bigR + 1.0D0)**2.0D0 * clatR(2) &
                                                    - rh_bigR * (rh_bigR - 1.0D0) * (rh_bigR - 2.0D0) * clatR(0)) &
                               * slat * slat * clon * zfun &
                               + R0 * R0 * rh_kapa * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) &
                                                      - (rh_bigR + 1.0D0)**2.0D0 * clatR(3)) * clat * clon * zfun

        this%psi3rd(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR &
                               * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) &
                                  - (rh_bigR + 1.0D0)**2.0D0 * clatR(3)) * slat * slon * zfun
        this%psi3rd(3, k, i) = R0 * R0 * rh_kapa * rh_bigR**2.0D0 * (rh_bigR * clatR(2) - (rh_bigR + 1.0D0) * clatR(4)) * clon * zfun
        this%psi3rd(4, k, i) = R0 * R0 * rh_kapa * rh_bigR**3.0D0 * clatR(3) * slat * slon * zfun

        this%psi3rd_dclat(1, k, i) = R0 * R0 * rh_kapa * (rh_bigR * (rh_bigR + 1.0D0)**2.0D0 * clatR(1) &
                                                          - rh_bigR * (rh_bigR - 1.0D0) * (rh_bigR - 2.0D0)) &
                                     * slat * slat * clon * zfun &
                                     + R0 * R0 * rh_kapa * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) &
                                                            - (rh_bigR + 1.0D0)**2.0D0 * clatR(2)) * clat * clon * zfun
        this%psi3rd_dclat(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR &
                                     * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) &
                                        - (rh_bigR + 1.0D0)**2.0D0 * clatR(2)) * slat * slon * zfun
        this%psi3rd_dclat(3, k, i) = R0 * R0 * rh_kapa * rh_bigR**2.0D0 * (rh_bigR * clatR(1) - (rh_bigR + 1.0D0) * clatR(3)) * clon * zfun
        this%psi3rd_dclat(4, k, i) = R0 * R0 * rh_kapa * rh_bigR**3.0D0 * clatR(2) * slat * slon * zfun

        ! ! The first order derivatives: for the last dimension, psi: (1) sigma, (2) lat, (3) lon; zta: (1) lat, (2) lon, (3)sigma
        ! ! *_dclat are parameters divived by cos(lat)
        ! this%psi1st(1, k, i) = zfun2 &
        !                        *stream(k, i)
        ! this%psi1st(2, k, i) = -R0*R0*zfun* &
        !                        (rh_omga*clat - rh_kapa* &
        !                         ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))*clon)
        ! this%psi1st(3, k, i) = -R0*R0*zfun &
        !                        *rh_kapa*rh_bigR*clatR(3)*slat*slon
        ! this%psi1st_dclat(1, k, i) = -R0*R0*zfun* &
        !                              (rh_omga - rh_kapa* &
        !                               ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*clon)
        ! this%psi1st_dclat(2, k, i) = -R0*R0*zfun &
        !                              *rh_kapa*rh_bigR*clatR(2)*slat*slon

        ! this%psi1st_dclatp2(1, k, i) = -R0*R0*zfun* &
        !                              (rh_omga - rh_kapa* &
        !                               ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*clon)
        ! this%psi1st_dclatp2(2, k, i) = -R0*R0*zfun &
        !                              *rh_kapa*rh_bigR*clatR(2)*slat*slon

        ! this%zta1st(1, k, i) = (2.0D0*rh_omga*clat - rh_kapa* &
        !                         ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))* &
        !                         (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon) &
        !                        *zfun
        ! this%zta1st(2, k, i) = rh_kapa*rh_bigR*slat*clatR(3)* &
        !                        (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*slon &
        !                        *zfun
        ! this%zta1st(3, k, i) = parzfun_sigma/zfun*vortct(k, i)

        ! this%zta1st_dclat(1, k, i) = (2.0D0*rh_omga - rh_kapa* &
        !                               ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))* &
        !                               (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon) &
        !                              *zfun
        ! this%zta1st_dclat(2, k, i) = rh_kapa*rh_bigR*slat*clatR(2)* &
        !                              (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*slon &
        !                              *zfun

        ! ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
        ! this%psi2nd(1, k, i) = R0*R0*slat &
        !                        *zfun &
        !                        *(rh_omga + rh_kapa* &
        !                          ((rh_bigR*(rh_bigR - 1.0D0)*clatR(1) - &
        !                            (rh_bigR + 1.0D0)**2*clatR(3)))*clon)
        ! this%psi2nd(2, k, i) = -R0*R0*rh_kapa*rh_bigR* &
        !                        zfun* &
        !                        ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))*slon
        ! this%psi2nd(3, k, i) = -R0*R0*rh_kapa*rh_bigR**2*clatR(3)*slat*clon &
        !                        *zfun

        ! this%psiThetaLambda_dclat(k, i) = -R0*R0*rh_kapa*rh_bigR* &
        !                                   ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*slon &
        !                                   *zfun
        ! this%psi2ndLambda_dclat(k, i) = -R0*R0*rh_kapa*rh_bigR**2*clatR(2)*slat*clon &
        !                                 *zfun

        !  zta2nd(1, i) = slat*(-2.0D0*rh_omga + rh_kapa* &
        !                            ((rh_bigR + 1.0D0)**2*clatR(4) - rh_bigR*(rh_bigR - 1.0D0)*clatR(1))* &
        !                            (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon)
        !  zta2nd(2, i) = rh_kapa*rh_bigR* &
        !                      ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))* &
        !                      (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*slon
        !  zta2nd(3, i) = rh_kapa*rh_bigR**2*clatR(3)* &
        !                      (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon

        ! The third order derivatives: 1 lat3; 2 lat2lon; 3 latlon2; 4 lon3
        !  psi3rd(1, i) = R0*R0*(rh_omga*clat + rh_kapa* &
        !                             (rh_bigR*(rh_bigR - 1.0D0)*clatR(1) - &
        !                              (rh_bigR + 1.0D0)*rh_bigR*clatR(3))*clat*clon - rh_kapa* &
        !                             (rh_bigR*(rh_bigR - 1.0D0)*(rh_bigR - 2.0D0)*clatR(0) - &
        !                              rh_bigR**2*(rh_bigR + 1.0D0)*clatR(2))*slat*slat*clon)
        !  psi3rd(2, i) = -R0*R0*slat*rh_kapa*rh_bigR**2* &
        !                      ((rh_bigR - 1.0D0)*clatR(1) - (rh_bigR + 1.0D0)*clatR(3))*slon
        !  psi3rd(3, i) = -R0*R0*rh_kapa*rh_bigR**2* &
        !                      (-rh_bigR*clatR(2) + (rh_bigR + 1.0D0)*clatR(4))*clon
        !  psi3rd(4, i) = R0*R0*rh_kapa*rh_bigR**3*clatR(3)*slat*slon

        ! Jacobian and flux divergence:
        !  IF (ABS(clat) .GT. machineEps) THEN
        ! Use equation (39) for divergence, (11) or (43) for Jacobian
        ! in the document formula
        ! jacobi(k, i) = (this%zta1st_dclat(2, k, i)*this%psi1st(2, k, i) - &
        !                 (this%zta1st(1, k, i) + 2.0D0*Omega*clat)*this%psi1st_dclat(2, k, i)) &
        !                /R0/R0
        ! flxdiv(k, i) = ((this%zta1st(1, k, i) + 2.0D0*Omega*clat)*this%psi1st(2, k, i) + &
        !                 this%zta1st_dclat(2, k, i)*this%psi1st_dclat(2, k, i))/R0/R0 + &
        !                (vortct(k, i) + 2.0D0*Omega*slat)*vortct(k, i)
        ! !  END IF
      END DO
    END DO

  END SUBROUTINE RH_psi_func

  SUBROUTINE RH_chi_func(this, attime, latlon, zhght, potential, divergence, par_div_t)
    IMPLICIT NONE

    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  attime, zhght(:, :), latlon(:, :)
    REAL(r_kind), INTENT(OUT) :: potential(:, :), &
                                 divergence(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_div_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: R0, clat, slat, tlat, clon, slon, clatR(0:5), &
                       rh_kapa_t, rh_omga_t, rh_angl_t, &
                       szhght, czhght, z_ratio, zfun, zfun2, parzfun_sigma, pi_ratio

    rh_kapa_t = rh_kapa / 1000.0D0
    rh_omga_t = rh_omga / 1000.0D0
    rh_angl_t = &
      (rh_bigR * (3.0D0 + rh_bigR) * rh_omga_t - &
       2.0D0 * Omega) / (1.0D0 + rh_bigR) / (2.0D0 + rh_bigR)

    z_ratio = 0.1D0
    pi_ratio = 0.1D0
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        R0 = EarthRadius + zhght(k, i)
        szhght = DSIN(zhght(k, i) / this%z_top * pi_ratio * pi)
        czhght = DCOS(zhght(k, i) / this%z_top * pi_ratio * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = (2.0D0 / R0 + &
                 z_ratio * czhght * pi_ratio * pi / this%z_top / zfun)
        parzfun_sigma = z_ratio * czhght * pi_ratio * pi / this%z_top
        this%zfunct(1, k, i) = zfun
        this%zfunct(2, k, i) = zfun2

        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        tlat = DTAN(latlon(1, i))
        clon = DCOS(rh_bigR * (latlon(2, i) + rh_angl_t * attime))
        slon = DSIN(rh_bigR * (latlon(2, i) + rh_angl_t * attime))

        clatR(0) = clat**(rh_bigR - 3.0D0)  ! c**(R-3)
        clatR(1) = clatR(0) * clat          ! c**(R-2)
        clatR(2) = clatR(1) * clat          ! c**(R-1)
        clatR(3) = clatR(2) * clat          ! c** R
        clatR(4) = clatR(3) * clat          ! c**(R+1)
        clatR(5) = clatR(4) * clat          ! c**(R+2)

        ! potential function:
        potential(k, i) = &
          R0 * R0 * zfun * slat * (rh_kapa_t * clatR(3) * clon)

        ! divergence:
        divergence(k, i) = -zfun * slat * (rh_kapa_t * clatR(3) * &
                                           (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon)
        IF (PRESENT(par_div_t)) THEN
          par_div_t(k, i) = rh_bigR * rh_angl * zfun * slat * (rh_kapa_t * clatR(3) * &
                                                               (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon)
        END IF

    !!! Note: these derivatives are derived in the document of
    !!!       Rossby_HaurwitzTest under
    !!!         /Users/xieyuanfu/developments/models/square/doc
    !!!       they need to double check the doc and coding here.
    !!! For testing the Poisson solver on a sphere, I currently
    !!! use the above stream function and vorticity only

        ! The first order derivatives: for the last dimension, chi: (1) sigma, (2) lat, (3) lon; del: (1) lat, (2) lon, (3)sigma
        ! *_dclat are parameters divived by cos(lat)
        this%chi1st(1, k, i) = zfun2 &
                               * potential(k, i)
        this%chi1st(2, k, i) = R0 * R0 * zfun * &
                               rh_kapa_t * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * clon
        this%chi1st(3, k, i) = -R0 * R0 * zfun &
                               * rh_kapa_t * rh_bigR * clatR(3) * slat * slon
        this%chi1st_dclat(1, k, i) = R0 * R0 * zfun * &
                                     (rh_kapa_t * &
                                      ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * clon)
        this%chi1st_dclat(2, k, i) = -R0 * R0 * zfun &
                                     * rh_kapa_t * rh_bigR * clatR(2) * slat * slon

        this%chi1st_dclatp2(1, k, i) = R0 * R0 * zfun * &
                                       (rh_kapa_t * &
                                        ((rh_bigR + 1.0D0) * clatR(2) - rh_bigR * clatR(0)) * clon)
        this%chi1st_dclatp2(2, k, i) = -R0 * R0 * zfun &
                                       * rh_kapa_t * rh_bigR * clatR(1) * slat * slon

        this%chi1st_dclatp3(1, k, i) = R0 * R0 * zfun * &
                                       (rh_kapa_t * &
                                        ((rh_bigR + 1.0D0) * clatR(1) - rh_bigR) * clon)
        this%chi1st_dclatp3(2, k, i) = -R0 * R0 * zfun &
                                       * rh_kapa_t * rh_bigR * clatR(0) * slat * slon

        this%del1st(1, k, i) = -(rh_kapa_t * &
                                 ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                                 (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                               * zfun
        this%del1st(2, k, i) = rh_kapa_t * rh_bigR * slat * clatR(3) * &
                               (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                               * zfun
        this%del1st(3, k, i) = parzfun_sigma / zfun * divergence(k, i)

        this%del1st_dclat(1, k, i) = -(rh_kapa_t * &
                                       ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * &
                                       (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                                     * zfun
        this%del1st_dclat(2, k, i) = rh_kapa_t * rh_bigR * slat * clatR(2) * &
                                     (rh_bigR**2.0D0 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                                     * zfun

        ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
        this%chi2nd(1, k, i) = R0 * R0 * slat &
                               * zfun &
                               * rh_kapa_t * &
                               ((rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                 (rh_bigR + 1.0D0)**2.0D0 * clatR(3))) * clon
        this%chi2nd(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * &
                               zfun * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * slon
        this%chi2nd(3, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * clatR(3) * slat * clon &
                               * zfun

        this%chi2nd_dclat(1, k, i) = R0 * R0 * slat &
                                     * zfun &
                                     * (rh_kapa_t * &
                                        ((rh_bigR * (rh_bigR - 1.0D0) * clatR(0) - &
                                          (rh_bigR + 1.0D0)**2.0D0 * clatR(2))) * clon)
        this%chi2nd_dclat(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * &
                                     zfun * &
                                     ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon
        this%chi2nd_dclat(3, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * clatR(2) * slat * clon &
                                     * zfun

        this%chi2nd_dclatp2(1, k, i) = R0 * R0 * slat &
                                       * zfun &
                                       * (rh_kapa_t * &
                                          ((rh_bigR * (rh_bigR - 1.0D0) - &
                                            (rh_bigR + 1.0D0)**2.0D0 * clatR(1))) * clon)
        this%chi2nd_dclatp2(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * &
                                       zfun * &
                                       ((rh_bigR + 1.0D0) * clatR(2) - rh_bigR * clatR(0)) * slon
        this%chi2nd_dclatp2(3, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * clatR(1) * slat * clon &
                                       * zfun

        this%chiThetaLambda_dclat(k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * &
                                          ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon &
                                          * zfun
        this%chi2ndLambda_dclat(k, i) = -R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * clatR(2) * slat * clon &
                                        * zfun

        ! the third order derivatives of chi: (1) lat, (2)latp2_lon, (3)lonp2_lat, (4)lon
        this%chi3rd(1, k, i) = R0 * R0 * rh_kapa_t * (rh_bigR * (rh_bigR + 1.0D0)**2.0D0 * clatR(2) &
                                                      - rh_bigR * (rh_bigR - 1.0D0) * (rh_bigR - 2.0D0) * clatR(0)) &
                               * slat * slat * clon * zfun &
                               + R0 * R0 * rh_kapa_t * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) &
                                                        - (rh_bigR + 1.0D0)**2.0D0 * clatR(3)) * clat * clon * zfun

        this%chi3rd(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR &
                               * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) &
                                  - (rh_bigR + 1.0D0)**2.0D0 * clatR(3)) * slat * slon * zfun
        this%chi3rd(3, k, i) = R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * (rh_bigR * clatR(2) - (rh_bigR + 1.0D0) * clatR(4)) * clon * zfun
        this%chi3rd(4, k, i) = R0 * R0 * rh_kapa_t * rh_bigR**3.0D0 * clatR(3) * slat * slon * zfun

        this%chi3rd_dclat(1, k, i) = R0 * R0 * rh_kapa_t * (rh_bigR * (rh_bigR + 1.0D0)**2.0D0 * clatR(1) &
                                                            - rh_bigR * (rh_bigR - 1.0D0) * (rh_bigR - 2.0D0)) &
                                     * slat * slat * clon * zfun &
                                     + R0 * R0 * rh_kapa_t * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) &
                                                              - (rh_bigR + 1.0D0)**2.0D0 * clatR(2)) * clat * clon * zfun
        this%chi3rd_dclat(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR &
                                     * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) &
                                        - (rh_bigR + 1.0D0)**2.0D0 * clatR(2)) * slat * slon * zfun
        this%chi3rd_dclat(3, k, i) = R0 * R0 * rh_kapa_t * rh_bigR**2.0D0 * (rh_bigR * clatR(1) - (rh_bigR + 1.0D0) * clatR(3)) * clon * zfun
        this%chi3rd_dclat(4, k, i) = R0 * R0 * rh_kapa_t * rh_bigR**3.0D0 * clatR(2) * slat * slon * zfun

        ! this%chi1st(2, k, i) = -R0*R0*zfun*(rh_omga_t*clat - rh_kapa_t* &
        !                                     ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))*clon)
        ! this%chi1st(3, k, i) = -R0*R0*zfun*rh_kapa_t*rh_bigR*clatR(3)*slat*slon
        ! this%chi1st_dclat(1, k, i) = -R0*R0 &
        !                              *zfun*(rh_omga_t - rh_kapa_t* &
        !                                     ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*clon)
        ! this%chi1st_dclat(2, k, i) = -R0*R0*zfun*rh_kapa_t*rh_bigR*clatR(2)*slat*slon

        ! this%del1st(1, k, i) = (2.0D0*rh_omga_t*clat - rh_kapa_t* &
        !                         ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))* &
        !                         (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon)*zfun
        ! this%del1st(2, k, i) = zfun*rh_kapa_t*rh_bigR*slat*clatR(3)* &
        !                        (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*slon
        ! this%del1st(3, k, i) = parzfun_sigma/zfun*divergence(k, i)

        ! this%del1st_dclat(1, k, i) = (2.0D0*rh_omga_t - rh_kapa_t* &
        !                               ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))* &
        !                               (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*clon)*zfun
        ! this%del1st_dclat(2, k, i) = zfun*rh_kapa_t*rh_bigR*slat*clatR(2)* &
        !                              (rh_bigR**2 + 3.0D0*rh_bigR + 2.0D0)*slon

        ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:

        ! this%chi2nd(1, k, i) = R0*R0*slat &
        !                        *zfun &
        !                        *(rh_kapa_t* &
        !                          ((rh_bigR*(rh_bigR - 1.0D0)*clatR(1) - &
        !                            (rh_bigR + 1.0D0)**2.0D0*clatR(3)))*clon)
        ! this%chi2nd(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR* &
        !                        zfun* &
        !                        ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))*slon
        ! this%chi2nd(3, k, i) = -R0*R0*rh_kapa_t*rh_bigR**2*clatR(3)*slat*clon &
        !                        *zfun

        ! this%chi2nd_dclat(1, k, i) = R0*R0*slat &
        !                              *zfun &
        !                              *(rh_kapa_t* &
        !                                ((rh_bigR*(rh_bigR - 1.0D0)*clatR(0) - &
        !                                  (rh_bigR + 1.0D0)**2*clatR(2)))*clon)
        ! this%chi2nd_dclat(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR* &
        !                              zfun* &
        !                              ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*slon
        ! this%chi2nd_dclat(3, k, i) = -R0*R0*rh_kapa_t*rh_bigR**2*clatR(2)*slat*clon &
        !                              *zfun

        ! this%chi2nd_dclatp2(1, k, i) = R0*R0*slat &
        !                                *zfun &
        !                                *(rh_kapa_t* &
        !                                  ((rh_bigR*(rh_bigR - 1.0D0) - &
        !                                    (rh_bigR + 1.0D0)**2*clatR(1)))*clon)
        ! this%chi2nd_dclatp2(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR* &
        !                                zfun* &
        !                                ((rh_bigR + 1.0D0)*clatR(2) - rh_bigR*clatR(0))*slon
        ! this%chi2nd_dclatp2(3, k, i) = -R0*R0*rh_kapa_t*rh_bigR**2*clatR(1)*slat*clon &
        !                                *zfun

        ! ! this%chi2nd(1, k, i) = R0*R0*zfun*slat*(rh_omga_t + rh_kapa_t* &
        ! !                                         ((rh_bigR*(rh_bigR - 1.0D0)*clatR(1) - &
        ! !                                           (rh_bigR + 1.0D0)**2*clatR(3)))*clon)
        ! ! this%chi2nd(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR*zfun* &
        ! !                        ((rh_bigR + 1.0D0)*clatR(4) - rh_bigR*clatR(2))*slon
        ! ! this%chi2nd(3, k, i) = -R0*R0*zfun*rh_kapa_t*rh_bigR**2*clatR(3)*slat*clon

        ! this%chiThetaLambda_dclat(k, i) = -R0*R0*rh_kapa_t*rh_bigR*zfun* &
        !                                   ((rh_bigR + 1.0D0)*clatR(3) - rh_bigR*clatR(1))*slon
        ! this%chi2ndLambda_dclat(k, i) = -R0*R0*zfun*rh_kapa_t*rh_bigR**2*clatR(2)*slat*clon

        ! ! the third order derivatives of psi: (1) lat, (2)latp2_lon, (3)lonp2_lat, (4)lon
        ! this%chi3rd(1, k, i) = R0*R0*rh_kapa_t*(-rh_bigR*(rh_bigR - 1.0D0)*(rh_bigR - 2.0D0)*clatR(0) &
        !                                         + rh_bigR*(rh_bigR + 1.0D0)**2.0D0*clatR(2))*slat*slat*clon*zfun &
        !                        + R0*R0*rh_kapa_t*(rh_bigR*(rh_bigR - 1.0D0)*clatR(1) &
        !                                           - (rh_bigR + 1.0D0)**2.0D0*clatR(3))*clat*clon*zfun
        ! this%chi3rd(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR &
        !                        *(rh_bigR*(rh_bigR - 1.0D0)*clatR(1) &
        !                          - (rh_bigR + 1.0D0)**2.0D0*clatR(3))*slat*slon*zfun
        ! this%chi3rd(3, k, i) = R0*R0*rh_kapa_t*rh_bigR**2.0D0*(rh_bigR*clatR(2) - (rh_bigR + 1.0D0)*clatR(4))*clon*zfun
        ! this%chi3rd(4, k, i) = R0*R0*rh_kapa_t*rh_bigR**3.0D0*clatR(3)*slat*slon*zfun

        ! this%chi3rd_dclat(1, k, i) = R0*R0*rh_kapa_t*(-rh_bigR*(rh_bigR - 1.0D0)*(rh_bigR - 2.0D0) &
        !                                               + rh_bigR*(rh_bigR + 1.0D0)**2.0D0*clatR(1))*slat*slat*clon*zfun &
        !                              + R0*R0*rh_kapa_t*(rh_bigR*(rh_bigR - 1.0D0)*clatR(0) &
        !                                                 - (rh_bigR + 1.0D0)**2.0D0*clatR(2))*clat*clon*zfun
        ! this%chi3rd_dclat(2, k, i) = -R0*R0*rh_kapa_t*rh_bigR &
        !                              *(rh_bigR*(rh_bigR - 1.0D0)*clatR(0) &
        !                                - (rh_bigR + 1.0D0)**2.0D0*clatR(2))*slat*slon*zfun
        ! this%chi3rd_dclat(3, k, i) = R0*R0*rh_kapa_t*rh_bigR**2.0D0*(rh_bigR*clatR(1) - (rh_bigR + 1.0D0)*clatR(3))*clon*zfun
        ! this%chi3rd_dclat(4, k, i) = R0*R0*rh_kapa_t*rh_bigR**3.0D0*clatR(2)*slat*slon*zfun

        ! Jacobian and flux divergence:
        !  IF (ABS(clat) .GT. machineEps) THEN
        ! Use equation (39) for divergence, (11) or (43) for Jacobian
        ! in the document formula
        ! jacobi(k, i) = (this%del1st_dclat(2, k, i)*this%chi1st(2, k, i) - &
        !                 this%del1st(1, k, i)*this%chi1st_dclat(2, k, i)) &
        !                /R0/R0
        ! flxdiv(k, i) = (this%del1st(1, k, i)*this%chi1st(2, k, i) + &
        !                 this%del1st_dclat(2, k, i)*this%chi1st_dclat(2, k, i))/R0/R0 + &
        !                divergence(k, i)*divergence(k, i)
        !  END IF
      END DO
    END DO

  END SUBROUTINE RH_chi_func

  SUBROUTINE Hght_func(this, latlon, z, z_s, Hz)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) :: latlon(:, :)
    REAL(r_kind), INTENT(OUT) :: z(:, :), z_s(:), Hz(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, R0, clat, slat, slat_p2, clon, &
                       slon, clatR(0:5), Rclat_t

    omg = 1500.0D0

    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        slat_p2 = slat**2.0D0
        clon = DCOS(rh_bigR * latlon(2, i))
        slon = DSIN(rh_bigR * latlon(2, i))

        clatR(0) = clat**(rh_bigR - 3.0D0)  ! c**(R-3)
        clatR(1) = clatR(0) * clat          ! c**(R-2)
        clatR(2) = clatR(1) * clat          ! c**(R-1)
        clatR(3) = clatR(2) * clat          ! c** R
        clatR(4) = clatR(3) * clat          ! c**(R+1)
        clatR(5) = clatR(4) * clat          ! c**(R+2)

        Rclat_t = (rh_bigR + 2.0D0) * clat**2 - rh_bigR

        ! IF (k .EQ. 1) THEN
        !    z_s(i) = omg*clatR(3)*slat_p2 + 1.0D0
        !    this%zs1sttheta(i) = omg*slat*clatR(2)*Rclat_t
        !    this%zs2ndtheta(i) = -(rh_bigR - 1.0D0)*omg*clatR(1)*slat_p2*Rclat_t &
        !                         + omg*clatR(3)*Rclat_t &
        !                         - 2*omg*(rh_bigR + 2.0D0)*clatR(3)*slat_p2
        !    this%zs1sttheta_dclat(i) = omg*slat*clatR(1)*Rclat_t
        !    this%zs2ndtheta_dclat(i) = -(rh_bigR - 1.0D0)*omg*clatR(0)*slat_p2*Rclat_t &
        !                               + omg*clatR(2)*Rclat_t &
        !                               - 2*omg*(rh_bigR + 2.0D0)*clatR(2)*slat_p2
        ! END IF

        IF (k .EQ. 1) THEN
          z_s(i) = omg * clatR(3)
          this%zs1sttheta(i) = -omg * rh_bigR * slat * clatR(2)
          this%zs2ndtheta(i) = rh_bigR * (rh_bigR - 1.0D0) * omg * clatR(1) * slat_p2 - &
                               rh_bigR * omg * clatR(3)
          this%zs1sttheta_dclat(i) = -omg * rh_bigR * slat * clatR(1)
          this%zs2ndtheta_dclat(i) = rh_bigR * (rh_bigR - 1.0D0) * omg * clatR(0) * slat_p2 - &
                                     rh_bigR * omg * clatR(2)
        END IF

        z(k, i) = this%sigma(k) * (this%z_top - z_s(i)) / this%z_top + z_s(i)! + (1.0D0 - this%sigma(k)/this%z_top)*z_s(i)
        Hz(k, i) = this%z_top / (this%z_top - z_s(i))

        this%z1sttheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * this%zs1sttheta(i)

        this%z1sttheta_dclat(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * this%zs1sttheta_dclat(i)

        this%d_z1stthetadclat_dtheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * &
                                             (omg * rh_bigR * (rh_bigR - 2.0D0) * slat_p2 * clatR(0) - &
                                              omg * rh_bigR * clatR(2))

        this%z2ndtheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * this%zs2ndtheta(i)
        this%z2ndtheta_dclat(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * this%zs2ndtheta_dclat(i)

        this%Hz1sttheta(k, i) = this%z_top / (this%z_top - z_s(i))**2.0D0 * this%zs1sttheta(i)
        this%Hz1sttheta_dclat(k, i) = this%z_top / (this%z_top - z_s(i))**2.0D0 * this%zs1sttheta_dclat(i)
      END DO
    END DO

  END SUBROUTINE Hght_func

  SUBROUTINE pres_func(this, attime, latlon, z, z_s, pres, par_pres_t)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), z(:, :), z_s(:), attime
    REAL(r_kind), INTENT(OUT) :: pres(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_pres_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, P0, alpha, R0, clat, slat, slat_p2, clon, &
                       slon, clatR(0:5), Rclat_t, expz, zratio, ct, st

    omg = 1.5D-2
    P0 = 1.0D5
    alpha = 1.073D-4 / 5.0D0

    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        ! R0 = EarthRadius + z(k, i)
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        slat_p2 = slat**2.0D0
        clon = DCOS(rh_bigR * latlon(2, i))
        slon = DSIN(rh_bigR * latlon(2, i))
        ct = DCOS(rh_angl * attime)
        st = DSIN(rh_angl * attime)

        clatR(0) = clat**(rh_bigR - 3.0D0)  ! c**(R-3)
        clatR(1) = clatR(0) * clat          ! c**(R-2)
        clatR(2) = clatR(1) * clat          ! c**(R-1)
        clatR(3) = clatR(2) * clat          ! c** R
        clatR(4) = clatR(3) * clat          ! c**(R+1)
        clatR(5) = clatR(4) * clat          ! c**(R+2)

        expz = DEXP(-alpha * z(k, i))
        zratio = (1.0D0 - this%sigma(k) / this%z_top)

        pres(k, i) = P0 * (1.0D0 + omg * clatR(3) * slat) * expz * ct
        this%pres1stsigma(k, i) = -alpha * pres(k, i) * (this%z_top - z_s(i)) / this%z_top
        this%pres1sttheta(k, i) = P0 * omg * (-rh_bigR * clatR(2) + (rh_bigR + 1.0D0) * clatR(4)) * expz * ct &
                                  - alpha * pres(k, i) * zratio * this%zs1sttheta(i)
        this%pres1sttheta_dclat(k, i) = P0 * omg * (-rh_bigR * clatR(1) + (rh_bigR + 1.0D0) * clatR(3)) * expz * ct &
                                        - alpha * pres(k, i) * zratio * this%zs1sttheta_dclat(i)
        this%pres2ndtheta(k, i) = P0 * omg * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) * slat - (rh_bigR + 1.0D0)**2 * clatR(3) * slat) * expz * ct &
                                  - alpha * P0 * omg * (-rh_bigR * clatR(2) + (rh_bigR + 1.0D0) * clatR(4)) * expz * zratio * this%zs1sttheta(i) * ct &
                                  - alpha * zratio * this%pres1sttheta(k, i) * this%zs1sttheta(i) &
                                  - alpha * pres(k, i) * zratio * this%zs2ndtheta(i)
        this%pres2ndtheta_dclat(k, i) = P0 * omg * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) * slat - (rh_bigR + 1.0D0)**2.0D0 * clatR(2) * slat) * expz * ct &
                                        - alpha * P0 * omg * (-rh_bigR * clatR(1) + (rh_bigR + 1.0D0) * clatR(3)) * expz * zratio * this%zs1sttheta(i) * ct &
                                        - alpha * zratio * this%pres1sttheta(k, i) * this%zs1sttheta_dclat(i) &
                                        - alpha * pres(k, i) * zratio * this%zs2ndtheta_dclat(i)
        this%pres2ndsigmatheta(k, i) = -alpha * (this%z_top - z_s(i)) / this%z_top * this%pres1sttheta(k, i) &
                                       + alpha * pres(k, i) / this%z_top * this%zs1sttheta(i)
        IF (PRESENT(par_pres_t)) THEN
          par_pres_t(k, i) = -P0 * (1.0D0 + omg * clatR(3) * slat) * expz * st * rh_angl
        END IF

      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE rho_func(this, attime, latlon, zhight, z_s, rho, par_rho_t)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), zhight(:, :), z_s(:), attime
    REAL(r_kind), INTENT(OUT) :: rho(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_rho_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, rho0, alpha, R0, clat, slat, slat_p2, clon, &
                       slon, clatR(0:5), Rclat_t, expz, z_ratio, ct, st, &
                       zfun, zfun2, zfun2_h, szhght, czhght, pi_ratio

    omg = 1.0D-1
    rho0 = 0.125D0
    alpha = omg

    z_ratio = 0.01D0
    pi_ratio = 0.1D0

    ! DO k = 1, this%vLevel
    !    DO i = 1, this%num_cell
    !       ! Trigonometry of latlon:
    !       ! R0 = EarthRadius + z(k, i)

    !       szhght = DSIN(this%sigma3d(k, i)/this%z_top*2.0D0*pi)
    !       czhght = DCOS(this%sigma3d(k, i)/this%z_top*2.0D0*pi)
    !       zfun = (1.0D0 + z_ratio*szhght)
    !       zfun2 = z_ratio*czhght*2.0D0*pi &
    !               /this%z_top

    !       clat = DCOS(latlon(1, i))
    !       slat = DSIN(latlon(1, i))
    !       slat_p2 = slat**2.0D0
    !       clon = DCOS(rh_bigR*latlon(2, i))
    !       slon = DSIN(rh_bigR*latlon(2, i))
    !       ct = DCOS(rh_angl*attime)
    !       st = DSIN(rh_angl*attime)

    !       clatR(0) = clat**(rh_bigR - 3.0D0)  ! c**(R-3)
    !       clatR(1) = clatR(0)*clat          ! c**(R-2)
    !       clatR(2) = clatR(1)*clat          ! c**(R-1)
    !       clatR(3) = clatR(2)*clat          ! c** R
    !       clatR(4) = clatR(3)*clat          ! c**(R+1)
    !       clatR(5) = clatR(4)*clat          ! c**(R+2)

    !       rho(k, i) = rho0*(1.0D0 + omg*slat)*zfun*ct
    !       this%rho1stsigma(k, i) = rho0*(1.0D0 + omg*slat)*zfun2*ct
    !       this%rho1sttheta(k, i) = rho0*omg*clat*zfun*ct

    !       IF (PRESENT(par_rho_t)) THEN
    !          par_rho_t(k, i) = -rho0*(1.0D0 + omg*slat)*zfun*st*rh_angl
    !       END IF

    !    END DO

    ! END DO
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        szhght = DSIN(zhight(k, i) / this%z_top * pi_ratio * pi)
        czhght = DCOS(zhight(k, i) / this%z_top * pi_ratio * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = z_ratio * czhght * pi_ratio * pi &
                / this%z_top**2.0D0 * (this%z_top - this%z_s(i))
        zfun2_h = z_ratio * czhght * pi_ratio * pi &
                  / this%z_top * this%z1sttheta(k, i)

        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        ct = DCOS(rh_angl * attime)
        st = DSIN(rh_angl * attime)

        rho(k, i) = rho0 * (1.0D0 + omg * slat) * zfun * ct
        this%rho1sttheta(k, i) = rho0 * omg * clat * zfun * ct + rho0 * (1.0D0 + omg * slat) * zfun2_h * ct
        this%rho1stsigma(k, i) = rho0 * (1.0D0 + omg * slat) * zfun2 * ct

        IF (PRESENT(par_rho_t)) THEN
          par_rho_t(k, i) = -rho0 * (1.0D0 + omg * slat) * zfun * st * rh_angl
        END IF

      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE Tp_func(this, attime, latlon, zhight, Tp, par_Tp_t)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), zhight(:, :), attime
    REAL(r_kind), INTENT(OUT) :: Tp(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_tp_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, Tp0, zratio, clat, slat, ct, st, &
                       szhght, czhght, zfun, zfun2, zfun2_h, hfun, &
                       pi_ratio

    omg = 1.5D-2
    Tp0 = 1.0D0
    zratio = 0.01D0
    pi_ratio = 0.1D0

    ! DO k = 1, this%vLevel
    !    DO i = 1, this%num_cell
    !       szhght = DSIN(this%sigma3d(k, i)/this%z_top*2.0D0*pi)
    !       czhght = DCOS(this%sigma3d(k, i)/this%z_top*2.0D0*pi)
    !       clat = DCOS(latlon(1, i))
    !       slat = DSIN(latlon(1, i))
    !       ct = DCOS(rh_angl*attime)
    !       st = DSIN(rh_angl*attime)
    !       hfun = 1.0D0 + omg*clat
    !       zfun = 1.0D0 + zratio*szhght

    !       Tp(k, i) = Tp0*zfun*hfun*ct
    !       this%Tp1sttheta(k, i) = -Tp0*zfun*omg*slat*ct
    !       this%Tp1stsigma(k, i) = zratio*Tp0*hfun*2.0D0*pi/this%z_top*czhght*ct

    !       IF (PRESENT(par_tp_t)) THEN
    !          par_tp_t(k, i) = -Tp0*zfun*hfun*st*rh_angl
    !       END IF

    !    END DO
    ! END DO
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        szhght = DSIN(zhight(k, i) / this%z_top * pi_ratio * pi)
        czhght = DCOS(zhight(k, i) / this%z_top * pi_ratio * pi)
        zfun = (1.0D0 + zratio * szhght)
        zfun2 = zratio * czhght * pi_ratio * pi &
                / this%z_top**2.0D0 * (this%z_top - this%z_s(i))
        zfun2_h = zratio * czhght * pi_ratio * pi &
                  / this%z_top * this%z1sttheta(k, i)

        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        ct = DCOS(rh_angl * attime)
        st = DSIN(rh_angl * attime)

        Tp(k, i) = Tp0 * (1.0D0 + omg * slat) * zfun * ct
        this%Tp1sttheta(k, i) = Tp0 * omg * clat * zfun * ct + Tp0 * (1.0D0 + omg * slat) * zfun2_h * ct
        this%Tp1stsigma(k, i) = Tp0 * (1.0D0 + omg * slat) * zfun2 * ct

        IF (PRESENT(par_Tp_t)) THEN
          par_Tp_t(k, i) = -Tp0 * (1.0D0 + omg * slat) * zfun * st * rh_angl
        END IF

      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE q_func(this, attime, latlon, zhight, q, par_q_t)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), zhight(:, :), attime
    REAL(r_kind), INTENT(OUT) :: q(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: par_q_t(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, q0, alpha, clat, slat, z_ratio, ct, st, &
                       zfun, zfun2, zfun2_h, szhght, czhght, pi_ratio

    omg = 1.0D-1
    q0 = 2.0D-2
    alpha = omg

    z_ratio = 0.01D0
    pi_ratio = 0.1D0

    ! DO k = 1, this%vLevel
    !    DO i = 1, this%num_cell
    !       szhght = DSIN(this%sigma(k)/this%z_top*2.0D0*pi)
    !       czhght = DCOS(this%sigma(k)/this%z_top*2.0D0*pi)
    !       zfun = (1.0D0 + z_ratio*szhght)
    !       zfun2 = z_ratio*czhght*2.0D0*pi &
    !               /this%z_top

    !       clat = DCOS(latlon(1, i))
    !       slat = DSIN(latlon(1, i))
    !       ct = DCOS(rh_angl*attime)
    !       st = DSIN(rh_angl*attime)

    !       q(k, i) = q0*(1.0D0 + omg*slat)*zfun*ct
    !       this%q1sttheta(k, i) = q0*omg*clat*zfun*ct
    !       this%q1stsigma(k, i) = q0*(1.0D0 + omg*slat)*zfun2*ct

    !       IF (PRESENT(par_q_t)) THEN
    !          par_q_t(k, i) = -q0*(1.0D0 + omg*slat)*zfun*st*rh_angl
    !       END IF

    !    END DO
    ! END DO

    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        szhght = DSIN(zhight(k, i) / this%z_top * pi_ratio * pi)
        czhght = DCOS(zhight(k, i) / this%z_top * pi_ratio * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = z_ratio * czhght * pi_ratio * pi &
                / this%z_top**2.0D0 * (this%z_top - this%z_s(i))
        zfun2_h = z_ratio * czhght * pi_ratio * pi &
                  / this%z_top * this%z1sttheta(k, i)

        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))
        ct = DCOS(rh_angl * attime)
        st = DSIN(rh_angl * attime)

        q(k, i) = q0 * (1.0D0 + omg * slat) * zfun * ct
        this%q1sttheta(k, i) = q0 * omg * clat * zfun * ct + q0 * (1.0D0 + omg * slat) * zfun2_h * ct
        this%q1stsigma(k, i) = q0 * (1.0D0 + omg * slat) * zfun2 * ct

        IF (PRESENT(par_q_t)) THEN
          par_q_t(k, i) = -q0 * (1.0D0 + omg * slat) * zfun * st * rh_angl
        END IF

      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE DM_func(this, Hz, DM)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  Hz(:, :)
    REAL(r_kind), INTENT(OUT) :: DM(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind) :: R0

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        DM(k, i) = Hz(k, i) / R0**2 * (this%z1sttheta(k, i) * this%chi1st(2, k, i) &
                                       + this%z1sttheta_dclat(k, i) * this%psi1st(3, k, i))

        ! 1st order of Dm, (1) theta; (2) Lambda
        this%DM1st(1, k, i) = 1.0D0 / R0**2 * (this%Hz1sttheta(k, i) * this%z1sttheta(k, i) * this%chi1st(2, k, i) &
                                               + Hz(k, i) * this%z2ndtheta(k, i) * this%chi1st(2, k, i) &
                                               + Hz(k, i) * this%z1sttheta(k, i) * this%chi2nd(1, k, i)) &
                              + 1.0D0 / R0**2 * (this%Hz1sttheta(k, i) * this%z1sttheta_dclat(k, i) * this%psi1st(3, k, i) &
                                                 + Hz(k, i) * this%z1sttheta_dclat(k, i) * this%psi2nd(2, k, i) &
                                                 + Hz(k, i) * this%d_z1stthetadclat_dtheta(k, i) * this%psi1st(3, k, i))
        this%DM1st(2, k, i) = Hz(k, i) / R0**2 * (this%z1sttheta(k, i) * this%chi2nd(2, k, i) &
                                                  + this%z1sttheta_dclat(k, i) * this%psi2nd(3, k, i))

        this%DM1stLambda_dclat(k, i) = Hz(k, i) / R0**2 * (this%z1sttheta_dclat(k, i) * this%chi2nd(2, k, i) &
                                                           + this%z1sttheta_dclat(k, i) * this%psi2ndLambda_dclat(k, i))
      END DO
    END DO

  END SUBROUTINE DM_func

  SUBROUTINE KE_func(this, latlon, psi, chi, LKE, KE)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  psi(:, :), chi(:, :), latlon(:, :)
    REAL(r_kind), INTENT(OUT) :: LKE(:, :)
    REAL(r_kind), INTENT(INOUT), OPTIONAL :: KE(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind) :: KEt(this%vLevel, this%num_cell)
    REAL(r_kind):: R0, slat, slat_p2, clat, &
                   F_psi_psi, F_chi_chi, J_psi_chi, &
                   KE2ndTheta, &
                   KE1stTheta_dclat, KE2ndLambda_dclatp2

    ! ALLOCATE (KE1stTheta(this%vLevel, this%num_cell), &
    !           KE2ndTheta(this%vLevel, this%num_cell), &
    !           KE1stTheta_dclat(this%vLevel, this%num_cell), &
    !           KE2ndLambda_dclatp2(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        slat = DSIN(latlon(1, i))
        slat_p2 = slat * slat
        clat = DCOS(latlon(1, i))
        F_psi_psi = (this%psi1st(2, k, i) * this%psi1st(2, k, i) + &
                     this%psi1st_dclat(2, k, i) * this%psi1st_dclat(2, k, i)) / R0 / R0

        F_chi_chi = (this%chi1st(2, k, i) * this%chi1st(2, k, i) + &
                     this%chi1st_dclat(2, k, i) * this%chi1st_dclat(2, k, i)) / R0 / R0

        J_psi_chi = (this%psi1st_dclat(2, k, i) * this%chi1st(2, k, i) - &
                     this%psi1st(2, k, i) * this%chi1st_dclat(2, k, i)) &
                    / R0 / R0

        KEt(k, i) = (F_psi_psi + F_chi_chi) / 2.0D0 + J_psi_chi

        ! KE1stTheta = (this%psi1st(2, k, i)*this%psi2nd(1, k, i) &
        !               + this%psi1st_dclatp2(2, k, i)*this%psi2nd(2, k, i) &
        !               + slat*clat*this%psi1st_dclatp2(2, k, i)**2.0D0 &
        !               + this%chi1st(2, k, i)*this%chi2nd(1, k, i) &
        !               + this%chi1st_dclatp2(2, k, i)*this%chi2nd(2, k, i) &
        !               + slat*clat*this%chi1st_dclatp2(2, k, i)**2.0D0 &
        !               + this%psi2nd(2, k, i)*this%chi1st_dclat(1, k, i) &
        !               + this%psi1st_dclat(2, k, i)*this%chi2nd(1, k, i) &
        !               - this%psi2nd(1, k, i)*this%chi1st_dclat(2, k, i) &
        !               - this%psi1st_dclat(1, k, i)*this%chi2nd(2, k, i) &
        !               + slat*(this%psi1st_dclatp2(2, k, i)*this%chi1st(2, k, i) &
        !                       - this%psi1st_dclatp2(1, k, i)*this%chi1st(3, k, i)))/R0/R0

        KE1stTheta_dclat = (this%psi1st_dclat(1, k, i) * this%psi2nd(1, k, i) &
                            + this%psi1st_dclatp3(2, k, i) * this%psi2nd(2, k, i) &
                            + slat * this%psi1st_dclatp2(2, k, i)**2.0D0 &
                            + this%chi1st_dclat(1, k, i) * this%chi2nd(1, k, i) &
                            + this%chi1st_dclatp3(2, k, i) * this%chi2nd(2, k, i) &
                            + slat * this%chi1st_dclatp2(2, k, i)**2.0D0 &
                            + this%psi2nd(2, k, i) * this%chi1st_dclatp2(1, k, i) &
                            + this%psi1st_dclatp2(2, k, i) * this%chi2nd(1, k, i) &
                            - this%psi2nd(1, k, i) * this%chi1st_dclatp2(2, k, i) &
                            - this%psi1st_dclatp2(1, k, i) * this%chi2nd(2, k, i) &
                            + slat * (this%psi1st_dclatp3(2, k, i) * this%chi1st(2, k, i) &
                                      - this%psi1st_dclatp3(1, k, i) * this%chi1st(3, k, i))) / R0 / R0
        KE2ndTheta = (this%psi2nd(1, k, i)**2.0D0 &
                      + this%psi1st(2, k, i) * this%psi3rd(1, k, i) &
                      + this%psi2nd_dclat(2, k, i)**2.0D0 &
                      + this%psi1st_dclatp2(2, k, i) * this%psi3rd(2, k, i) &
                      + 2.0D0 * slat * this%psi1st_dclatp3(2, k, i) * this%psi2nd(2, k, i) &
                      + 3.0D0 * slat_p2 * this%psi1st_dclatp2(2, k, i)**2.0D0 &
                      + this%psi1st_dclat(2, k, i)**2.0D0 &
                      + 2.0D0 * slat * this%psi1st_dclatp3(2, k, i) * this%psi2nd(2, k, i) &
                      + this%chi2nd(1, k, i)**2.0D0 &
                      + this%chi1st(2, k, i) * this%chi3rd(1, k, i) &
                      + this%chi2nd_dclat(2, k, i)**2.0D0 &
                      + this%chi1st_dclatp2(2, k, i) * this%chi3rd(2, k, i) &
                      + 2.0D0 * slat * this%chi1st_dclatp3(2, k, i) * this%chi2nd(2, k, i) &
                      + 3.0D0 * slat_p2 * this%chi1st_dclatp2(2, k, i)**2.0D0 &
                      + this%chi1st_dclat(2, k, i)**2.0D0 &
                      + 2.0D0 * slat * this%chi1st_dclatp3(2, k, i) * this%chi2nd(2, k, i) &
                      + this%psi3rd(2, k, i) * this%chi1st_dclat(1, k, i) &
                      + 2.0D0 * this%psi2nd(2, k, i) * this%chi2nd_dclat(1, k, i) &
                      + this%psi1st_dclat(2, k, i) * this%chi3rd(1, k, i) &
                      - this%psi3rd(1, k, i) * this%chi1st_dclat(2, k, i) &
                      - 2 * this%psi2nd_dclat(1, k, i) * this%chi2nd(2, k, i) &
                      - this%psi1st_dclat(1, k, i) * this%chi3rd(2, k, i) &
                      + slat * (this%psi2nd(2, k, i) * this%chi1st_dclatp2(1, k, i) &
                                + this%psi1st_dclatp2(2, k, i) * this%chi2nd(1, k, i) &
                                - this%psi2nd(1, k, i) * this%chi1st_dclatp2(2, k, i) &
                                - this%psi1st_dclatp2(1, k, i) * this%chi2nd(2, k, i)) &
                      + 2.0D0 * slat_p2 * (this%psi1st_dclatp3(2, k, i) * this%chi1st(2, k, i) &
                                           - this%psi1st_dclatp3(1, k, i) * this%chi1st(3, k, i)) &
                      + this%psi1st_dclat(2, k, i) * this%chi1st(2, k, i) &
                      - this%psi1st_dclat(1, k, i) * this%chi1st(3, k, i) &
                      + slat * (this%psi2nd(2, k, i) * this%chi1st_dclatp2(1, k, i) &
                                + this%psi1st_dclatp2(2, k, i) * this%chi2nd(1, k, i) &
                                - this%psi2nd(1, k, i) * this%chi1st_dclatp2(2, k, i) &
                                - this%psi1st_dclatp2(1, k, i) * this%chi2nd(2, k, i)) &
                      ) / R0 / R0

        KE2ndLambda_dclatp2 = &
          (this%psi2nd_dclat(2, k, i)**2.0D0 &
           + this%psi1st_dclatp2(1, k, i) * this%psi3rd(3, k, i) &
           + this%psi2nd_dclatp2(3, k, i)**2.0D0 &
           + this%psi1st_dclatp3(2, k, i) * this%psi3rd_dclat(4, k, i) &
           + this%chi2nd_dclat(2, k, i)**2.0D0 &
           + this%chi1st_dclatp2(1, k, i) * this%chi3rd(3, k, i) &
           + this%chi2nd_dclatp2(3, k, i)**2.0D0 &
           + this%chi1st_dclatp3(2, k, i) * this%chi3rd_dclat(4, k, i) &
           + this%psi3rd(4, k, i) * this%chi1st_dclatp3(1, k, i) &
           + 2.0D0 * this%psi2nd_dclatp2(3, k, i) * this%chi2nd_dclat(2, k, i) &
           + this%psi1st_dclatp3(2, k, i) * this%chi3rd(3, k, i) &
           - this%psi3rd(3, k, i) * this%chi1st_dclatp3(2, k, i) &
           - 2.0D0 * this%psi2nd_dclat(2, k, i) * this%chi2nd_dclatp2(3, k, i) &
           - this%psi1st_dclatp3(1, k, i) * this%chi3rd(4, k, i) &
           ) / R0 / R0

        LKE(k, i) = (-slat * KE1stTheta_dclat + KE2ndTheta + KE2ndLambda_dclatp2) / R0 / R0
      END DO
    END DO

    IF (PRESENT(KE)) KE = KEt

  END SUBROUTINE KE_func

  SUBROUTINE VorTen(this, latlon, vortct, divergence, psi, chi, DM, tenvor, X1)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                vortct(:, :), divergence(:, :), DM(:, :)
    REAL(r_kind), INTENT(OUT) :: tenvor(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat, J_DM_chi_R0, F_DM_psi_R0

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: J_eta_psi, F_eta_chi, F_DM_psisigma, &
                                                  J_DM_chisigma
    ALLOCATE (J_eta_psi(this%vLevel, this%num_cell), &
              F_eta_chi(this%vLevel, this%num_cell), &
              F_DM_psisigma(this%vLevel, this%num_cell), &
              J_DM_chisigma(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_eta_psi(k, i) = (this%zta1st_dclat(2, k, i) * this%psi1st(2, k, i) - &
                           (this%zta1st(1, k, i) + 2.0D0 * Omega * clat) * this%psi1st_dclat(2, k, i)) &
                          / R0 / R0

        F_eta_chi(k, i) = ((this%zta1st(1, k, i) + 2.0D0 * Omega * clat) * this%chi1st(2, k, i) + &
                           this%zta1st_dclat(2, k, i) * this%chi1st_dclat(2, k, i)) / R0 / R0 + &
                          (vortct(k, i) + 2.0D0 * Omega * slat) * divergence(k, i)

        F_DM_psisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(1, k, i) * this%psi1st(2, k, i) &
                              + this%zfunct(2, k, i) &
                              / R0**2 * this%DM1stLambda_dclat(k, i) * this%psi1st_dclat(2, k, i) &
                              + this%zfunct(2, k, i) &
                              * DM(k, i) * vortct(k, i)

        J_DM_chisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(2, k, i) * this%chi1st_dclat(1, k, i) &
                              - this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(1, k, i) * this%chi1st_dclat(2, k, i)
        J_DM_chi_R0 = (this%DM1st(2, k, i) * this%chi1st_dclat(1, k, i) - &
                       this%DM1st(1, k, i) * this%chi1st_dclat(2, k, i)) &
                      / R0 / R0 / R0
        F_DM_psi_R0 = (this%DM1st(1, k, i) * this%psi1st(2, k, i) + &
                       this%DM1stLambda_dclat(k, i) * this%psi1st_dclat(2, k, i)) / R0 / R0 / R0
        tenvor(k, i) = J_eta_psi(k, i) - F_eta_chi(k, i) + F_DM_psisigma(k, i) + J_DM_chisigma(k, i) - F_DM_psi_R0 - J_DM_chi_R0
      END DO
    END DO

    IF (PRESENT(X1)) THEN
      X1%fields(X1%getVarIdx('J_eta_psit'))%DATA(:, :, 1) = J_eta_psi
      X1%fields(X1%getVarIdx('F_eta_chit'))%DATA(:, :, 1) = F_eta_chi
      X1%fields(X1%getVarIdx('F_DM_psisigmat'))%DATA(:, :, 1) = F_DM_psisigma
      X1%fields(X1%getVarIdx('J_DM_chisigmat'))%DATA(:, :, 1) = J_DM_chisigma
    END IF

    DEALLOCATE (J_eta_psi, F_eta_chi, F_DM_psisigma, &
                J_DM_chisigma)

  END SUBROUTINE

  SUBROUTINE DivTen(this, latlon, divergence, vortict, psi, chi, DM, pres, rho, z, Hz, tendiv, X1)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                divergence(:, :), DM(:, :), pres(:, :), rho(:, :), &
                                Hz(:, :), z(:, :), vortict(:, :)
    REAL(r_kind), INTENT(OUT) :: tendiv(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat, f_theta

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: F_vor_psi, J_vor_chi, KE, LKE, J_DM_psisigma, &
                                                  F_DM_chisigma, J_DM_psi_R0, F_DM_chi_R0, F_invrho_P, &
                                                  F_HzinvRhoPsigma_z, F_f_psi, J_f_chi, &
                                                  J_DM_psi, F_DM_chi
    ALLOCATE (F_vor_psi(this%vLevel, this%num_cell), &
              J_vor_chi(this%vLevel, this%num_cell), &
              KE(this%vLevel, this%num_cell), &
              LKE(this%vLevel, this%num_cell), &
              J_DM_psisigma(this%vLevel, this%num_cell), &
              F_DM_chisigma(this%vLevel, this%num_cell), &
              J_DM_psi_R0(this%vLevel, this%num_cell), &
              F_DM_chi_R0(this%vLevel, this%num_cell), &
              F_invrho_P(this%vLevel, this%num_cell), &
              F_HzinvRhoPsigma_z(this%vLevel, this%num_cell), &
              F_f_psi(this%vLevel, this%num_cell), &
              J_f_chi(this%vLevel, this%num_cell), &
              J_DM_psi(this%vLevel, this%num_cell), &
              F_DM_chi(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_vor_chi(k, i) = (this%zta1st(2, k, i) * this%chi1st_dclat(1, k, i) - &
                           this%zta1st(1, k, i) * this%chi1st_dclat(2, k, i)) &
                          / R0 / R0

        F_vor_psi(k, i) = (this%zta1st(1, k, i) * this%psi1st(2, k, i) + &
                           this%zta1st_dclat(2, k, i) * this%psi1st_dclat(2, k, i)) / R0 / R0 + &
                          vortict(k, i)**2.0D0

        F_DM_chisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2.0D0 * this%DM1st(1, k, i) * this%chi1st(2, k, i) &
                              + this%zfunct(2, k, i) &
                              / R0**2.0D0 * this%DM1stLambda_dclat(k, i) * this%chi1st_dclat(2, k, i) &
                              + this%zfunct(2, k, i) &
                              * DM(k, i) * divergence(k, i)

        J_DM_psisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2.0D0 * this%DM1st(2, k, i) * this%psi1st_dclat(1, k, i) &
                              - this%zfunct(2, k, i) &
                              / R0**2.0D0 * this%DM1st(1, k, i) * this%psi1st_dclat(2, k, i)
        J_DM_psi_R0(k, i) = (this%DM1st(2, k, i) * this%psi1st_dclat(1, k, i) - &
                             this%DM1st(1, k, i) * this%psi1st_dclat(2, k, i)) &
                            / R0 / R0 / R0
        F_DM_chi_R0(k, i) = (this%DM1st(1, k, i) * this%chi1st(2, k, i) + &
                             this%DM1stLambda_dclat(k, i) * this%chi1st_dclat(2, k, i)) / R0 / R0 / R0! &
        !+ DM(k, i)*divergence(k, i)/R0

        J_DM_psi(k, i) = (this%DM1st(2, k, i) * this%psi1st_dclat(1, k, i) - &
                          this%DM1st(1, k, i) * this%psi1st_dclat(2, k, i)) &
                         / R0 / R0
        F_DM_chi(k, i) = (this%DM1st(1, k, i) * this%chi1st(2, k, i) + &
                          this%DM1stLambda_dclat(k, i) * this%chi1st_dclat(2, k, i)) / R0 / R0 &
                         + DM(k, i) * divergence(k, i)

        F_invrho_P(k, i) = 1.0D0 / rho(k, i) * (-slat * this%pres1sttheta_dclat(k, i) &
                                                + this%pres2ndtheta(k, i)) / R0 / R0 &
                           - 1.0D0 / rho(k, i)**2.0D0 * this%rho1sttheta(k, i) * &
                           this%pres1sttheta(k, i) / R0 / R0

        F_HzinvRhoPsigma_z(k, i) = Hz(k, i) / R0 / R0 * (1.0D0 / rho(k, i) * this%pres2ndsigmatheta(k, i) &
                                                         - 1.0D0 / rho(k, i)**2.0D0 * this%pres1stsigma(k, i) &
                                                         * this%rho1sttheta(k, i)) * this%z1sttheta(k, i) &
                                   + 1.0D0 / rho(k, i) / R0 / R0 * this%pres1stsigma(k, i) &
                                   * this%Hz1sttheta(k, i) * this%z1sttheta(k, i) &
                                   + Hz(k, i) * (-slat * this%z1sttheta_dclat(k, i) + this%z2ndtheta(k, i)) &
                                   / rho(k, i) * this%pres1stsigma(k, i) / R0 / R0
        f_theta = 2.0D0 * Omega * clat

        ! PRINT *, 'size of f is ', SIZE(this%f)!, '----', this%f(1:10, 1)
        F_f_psi(k, i) = f_theta / R0 / R0 * this%psi1st(2, k, i) &
                        + this%f(k, i) * vortict(k, i)

        J_f_chi(k, i) = -f_theta / R0 / R0 * this%chi1st_dclat(2, k, i)

      END DO
    END DO

    CALL this%KE_func(latlon, psi, chi, LKE, KE)

    tendiv = F_vor_psi + J_vor_chi - LKE - J_DM_psisigma + J_DM_psi_R0 &
             + F_DM_chisigma - F_DM_chi_R0 &
             - F_invrho_P + F_HzinvRhoPsigma_z &
             + F_f_psi + J_f_chi

    IF (PRESENT(X1)) THEN
      X1%fields(X1%getVarIdx('F_vor_psit'))%DATA(:, :, 1) = F_vor_psi
      X1%fields(X1%getVarIdx('J_vor_chit'))%DATA(:, :, 1) = J_vor_chi
      X1%fields(X1%getVarIdx('LKEt'))%DATA(:, :, 1) = LKE
      X1%fields(X1%getVarIdx('J_DM_psisigmat'))%DATA(:, :, 1) = J_DM_psisigma
      X1%fields(X1%getVarIdx('F_DM_chisigmat'))%DATA(:, :, 1) = F_DM_chisigma
      X1%fields(X1%getVarIdx('J_DM_psit'))%DATA(:, :, 1) = J_DM_psi
      X1%fields(X1%getVarIdx('F_DM_chit'))%DATA(:, :, 1) = F_DM_chi
      X1%fields(X1%getVarIdx('F_invRho_Pt'))%DATA(:, :, 1) = F_invrho_P
      X1%fields(X1%getVarIdx('F_HzinvRhoPsigma_zt'))%DATA(:, :, 1) = F_HzinvRhoPsigma_z
      X1%fields(X1%getVarIdx('F_f_psit'))%DATA(:, :, 1) = F_f_psi
      X1%fields(X1%getVarIdx('J_f_chit'))%DATA(:, :, 1) = J_f_chi
    END IF

    DEALLOCATE (F_vor_psi, J_vor_chi, LKE, J_DM_psisigma, &
                F_DM_chisigma, J_DM_psi_R0, F_DM_chi_R0, F_invrho_P, &
                F_HzinvRhoPsigma_z, F_f_psi, J_f_chi, &
                J_DM_psi, F_DM_chi)

  END SUBROUTINE

  SUBROUTINE RhoTen(this, latlon, rho, divergence, psi, chi, DM, z, Hz, tenrho, X1)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                rho(:, :), divergence(:, :), z(:, :), &
                                DM(:, :), Hz(:, :)
    REAL(r_kind), INTENT(OUT) :: tenrho(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: J_rho_psi, F_rho_chi, F_z_chisigma, &
                                                  J_z_psisigma, J_z_psi_r0, F_z_chi_r0
    ALLOCATE (J_rho_psi(this%vLevel, this%num_cell), &
              F_rho_chi(this%vLevel, this%num_cell), &
              F_z_chisigma(this%vLevel, this%num_cell), &
              J_z_psisigma(this%vLevel, this%num_cell), &
              J_z_psi_r0(this%vLevel, this%num_cell), &
              F_z_chi_r0(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_rho_psi(k, i) = -this%rho1sttheta(k, i) * this%psi1st_dclat(2, k, i) / R0 / R0

        F_rho_chi(k, i) = this%rho1sttheta(k, i) * this%chi1st(2, k, i) / R0 / R0 + rho(k, i) * divergence(k, i)

        F_z_chisigma(k, i) = this%zfunct(2, k, i) &
                             / R0**2.0D0 * this%z1sttheta(k, i) * this%chi1st(2, k, i)
        J_z_psisigma(k, i) = -this%zfunct(2, k, i) &
                             / R0**2.0D0 * this%z1sttheta(k, i) * this%psi1st_dclat(2, k, i)
        F_z_chi_r0(k, i) = 1.0D0 / R0**3.0D0 * this%z1sttheta(k, i) * this%chi1st(2, k, i)
        J_z_psi_r0(k, i) = -1.0D0 / R0**3.0D0 * this%z1sttheta(k, i) * this%psi1st_dclat(2, k, i)
        tenrho(k, i) = J_rho_psi(k, i) &
                       - F_rho_chi(k, i) &
                       + Hz(k, i) * rho(k, i) &
                       * (F_z_chisigma(k, i) - J_z_psisigma(k, i) &
                          - F_z_chi_r0(k, i) + J_z_psi_r0(k, i)) &
                       + DM(k, i) * this%rho1stsigma(k, i)
      END DO
    END DO
    IF (PRESENT(X1)) THEN
      X1%fields(X1%getVarIdx('prho_sgt'))%DATA(:, :, 1) = this%rho1stsigma
      ! X1%fields(X1%getVarIdx('par_rho_thetat'))%data(:, :, 1) = this%rho1sttheta
      X1%fields(X1%getVarIdx('J_rho_psit'))%DATA(:, :, 1) = J_rho_psi
      X1%fields(X1%getVarIdx('F_rho_chit'))%DATA(:, :, 1) = F_rho_chi
      X1%fields(X1%getVarIdx('F_z_chisgt'))%DATA(:, :, 1) = F_z_chisigma
      X1%fields(X1%getVarIdx('J_z_psisgt'))%DATA(:, :, 1) = J_z_psisigma
      X1%fields(X1%getVarIdx('F_z_chirt'))%DATA(:, :, 1) = F_z_chi_r0
      X1%fields(X1%getVarIdx('J_z_psirt'))%DATA(:, :, 1) = J_z_psi_r0
    END IF

    DEALLOCATE (J_rho_psi, F_rho_chi, F_z_chisigma, &
                J_z_psisigma, J_z_psi_r0, F_z_chi_r0)

  END SUBROUTINE

  SUBROUTINE wTen(this, latlon, rho, Hz, tenw)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    ! TYPE(State_t), INTENT(INOUT) :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), Hz(:, :), &
                                rho(:, :)
    REAL(r_kind), INTENT(OUT) :: tenw(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    tenw = -g - Hz * this%pres1stsigma / rho

    PRINT *, "max Hz and pre1stsigma are: ", MAXVAL(Hz), MAXVAL(this%pres1stsigma), &
      MINVAL(Hz), MINVAL(this%pres1stsigma), MAXVAL(rho), MINVAL(rho)

    ! X1%fields(X1%getVarIdx('J_rho_psit'))%data(:, :, 1) = J_rho_psi
    ! X1%fields(X1%getVarIdx('F_rho_chit'))%data(:, :, 1) = F_rho_chi
    ! X1%fields(X1%getVarIdx('F_z_chisigmat'))%data(:, :, 1) = F_DM_psisigma
    ! X1%fields(X1%getVarIdx('J_z_psisigmat'))%data(:, :, 1) = J_DM_chisigma

  END SUBROUTINE

  SUBROUTINE TpTen(this, latlon, psi, chi, z, Hz, div, tenTp, X)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                z(:, :), Hz(:, :), div(:, :)
    REAL(r_kind), INTENT(OUT) :: tenTp(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: J_Tp_psi, F_Tp_chi, J_z_psi, F_z_chi, F_z_chi0
    ALLOCATE (J_Tp_psi(this%vLevel, this%num_cell), &
              F_Tp_chi(this%vLevel, this%num_cell), &
              F_z_chi(this%vLevel, this%num_cell), &
              F_z_chi0(this%vLevel, this%num_cell), &
              J_z_psi(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_Tp_psi(k, i) = -this%Tp1sttheta(k, i) * this%psi1st_dclat(2, k, i) / R0 / R0

        F_Tp_chi(k, i) = this%Tp1sttheta(k, i) * this%chi1st(2, k, i) / R0 / R0
        F_z_chi(k, i) = 1.0D0 / R0**2.0D0 * this%z1sttheta(k, i) * this%chi1st(2, k, i)
        F_z_chi0(k, i) = F_z_chi(k, i) + z(k, i) * div(k, i)
        J_z_psi(k, i) = -1.0D0 / R0**2.0D0 * this%z1sttheta(k, i) * this%psi1st_dclat(2, k, i)
        tenTp(k, i) = J_Tp_psi(k, i) &
                      - F_Tp_chi(k, i) &
                      + Hz(k, i) * this%Tp1stsigma(k, i) &
                      * (F_z_chi(k, i) - J_z_psi(k, i))
      END DO
    END DO

    IF (PRESENT(X)) THEN
      X%fields(X%getVarIdx('J_z_psiFt'))%DATA(:, :, 1) = J_z_psi
      X%fields(X%getVarIdx('F_z_chiFt'))%DATA(:, :, 1) = F_z_chi0
      X%fields(X%getVarIdx('J_th_psit'))%DATA(:, :, 1) = J_Tp_psi
      X%fields(X%getVarIdx('F_th_chit'))%DATA(:, :, 1) = F_Tp_chi
      X%fields(X%getVarIdx('chit'))%DATA(:, :, 1) = chi
      X%fields(X%getVarIdx('psit'))%DATA(:, :, 1) = psi
    END IF

    ! X1%fields(X1%getVarIdx('J_rho_psit'))%data(:, :, 1) = J_rho_psi
    ! X1%fields(X1%getVarIdx('F_rho_chit'))%data(:, :, 1) = F_rho_chi
    ! X1%fields(X1%getVarIdx('F_z_chisigmat'))%data(:, :, 1) = F_DM_psisigma
    ! X1%fields(X1%getVarIdx('J_z_psisigmat'))%data(:, :, 1) = J_DM_chisigma

    DEALLOCATE (J_Tp_psi, F_Tp_chi, J_z_psi, F_z_chi, F_z_chi0)

  END SUBROUTINE

  SUBROUTINE qTen(this, latlon, psi, chi, z, Hz, div, tenq, X)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                z(:, :), Hz(:, :), div(:, :)
    REAL(r_kind), INTENT(OUT) :: tenq(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: J_q_psi, F_q_chi, J_z_psi, F_z_chi, F_z_chi0
    ALLOCATE (J_q_psi(this%vLevel, this%num_cell), &
              F_q_chi(this%vLevel, this%num_cell), &
              F_z_chi(this%vLevel, this%num_cell), &
              J_z_psi(this%vLevel, this%num_cell), &
              F_z_chi0(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_q_psi(k, i) = -this%q1sttheta(k, i) * this%psi1st_dclat(2, k, i) / R0 / R0

        F_q_chi(k, i) = this%q1sttheta(k, i) * this%chi1st(2, k, i) / R0 / R0
        F_z_chi(k, i) = 1.0D0 / R0**2.0D0 * this%z1sttheta(k, i) * this%chi1st(2, k, i)
        J_z_psi(k, i) = -1.0D0 / R0**2.0D0 * this%z1sttheta(k, i) * this%psi1st_dclat(2, k, i)
        F_z_chi0(k, i) = F_z_chi(k, i) + z(k, i) * div(k, i)
        tenq(k, i) = J_q_psi(k, i) &
                     - F_q_chi(k, i) &
                     + Hz(k, i) * this%q1stsigma(k, i) &
                     * (F_z_chi(k, i) - J_z_psi(k, i))
      END DO
    END DO

    PRINT *, "max of zhght and Hz in qTen are: ", MAXVAL(this%zhght(1, :)), MAXVAL(this%Hz(1, :))

    IF (PRESENT(X)) THEN
      X%fields(X%getVarIdx('J_q_psit'))%DATA(:, :, 1) = J_q_psi
      X%fields(X%getVarIdx('F_q_chit'))%DATA(:, :, 1) = F_q_chi
      X%fields(X%getVarIdx('J_z_psit'))%DATA(:, :, 1) = J_z_psi
      X%fields(X%getVarIdx('F_z_chit'))%DATA(:, :, 1) = F_z_chi0
      X%fields(X%getVarIdx('pq_sgt'))%DATA(:, :, 1) = this%q1stsigma

    END IF

    DEALLOCATE (J_q_psi, F_q_chi, J_z_psi, F_z_chi, F_z_chi0)

  END SUBROUTINE

  SUBROUTINE righthandsHalf(this, t0, XHalf)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XHalf

    REAL(r_kind), INTENT(IN) :: t0

    INTEGER(i_kind) :: div_idxt, vor_idxt, psi_idxt, chi_idxt, q_idxt, rho_idxt, &
                       pres_idxt, w_idxt, theta_idxt, &
                       tdiv_idxt, tvor_idxt, tq_idxt, trho_idxt, &
                       tw_idxt, ttheta_idxt, &
                       pdiv_idxt, pvor_idxt, pq_idxt, prho_idxt, &
                       ptheta_idxt, ppres_idxt

    vor_idxt = XHalf%getVarIdx('vort')
    div_idxt = XHalf%getVarIdx('divt')
    rho_idxt = XHalf%getVarIdx('rhot')
    q_idxt = XHalf%getVarIdx('qvaport')
    pres_idxt = XHalf%getVarIdx('prest')
    psi_idxt = XHalf%getVarIdx('psit')
    chi_idxt = XHalf%getVarIdx('chit')
    ! w_idxt = XFull%getVarIdx('wt')
    theta_idxt = XHalf%getVarIdx('thetat')

    tvor_idxt = XHalf%getVarIdx('tenvort')
    tdiv_idxt = XHalf%getVarIdx('tendivt')
    trho_idxt = XHalf%getVarIdx('tenrhot')
    tq_idxt = XHalf%getVarIdx('tenqvaport')
    ! tw_idxt = XFull%getVarIdx('tenwt')
    ! ttheta_idxt = XFull%getVarIdx('tenthetat')

    pvor_idxt = XHalf%getVarIdx('parvort')
    pdiv_idxt = XHalf%getVarIdx('pardivt')
    prho_idxt = XHalf%getVarIdx('parrhot')
    pq_idxt = XHalf%getVarIdx('parqvaport')
    ptheta_idxt = XHalf%getVarIdx('parthetat')
    ppres_idxt = XHalf%getVarIdx('parprest')

    ASSOCIATE (sg => XHalf%sg, &
               fields => XHalf%fields)

      CALL this%RH_psi_func(t0, sg%cell_cntr, this%sigma3d, &
                            fields(psi_idxt)%DATA(:, :, 1), fields(vor_idxt)%DATA(:, :, 1))
      CALL this%RH_chi_func(t0, sg%cell_cntr, this%sigma3d, &
                            fields(chi_idxt)%DATA(:, :, 1), fields(div_idxt)%DATA(:, :, 1))
      CALL this%pres_func(t0, sg%cell_cntr, this%zHght, this%z_s, &
                          fields(pres_idxt)%DATA(:, :, 1))
      CALL this%rho_func(t0, sg%cell_cntr, this%zHght, this%z_s, &
                         fields(rho_idxt)%DATA(:, :, 1))
      CALL this%q_func(t0, sg%cell_cntr, this%zHght, fields(q_idxt)%DATA(:, :, 1))
      CALL this%Tp_func(t0, sg%cell_cntr, this%zHght, fields(theta_idxt)%DATA(:, :, 1))

      CALL this%DM_func(this%Hz, this%DM)
      CALL this%DivTen(sg%cell_cntr, fields(div_idxt)%DATA(:, :, 1), fields(vor_idxt)%DATA(:, :, 1), &
                       fields(psi_idxt)%DATA(:, :, 1), fields(chi_idxt)%DATA(:, :, 1), &
                       this%DM, fields(pres_idxt)%DATA(:, :, 1), fields(rho_idxt)%DATA(:, :, 1), &
                       this%zHght, this%Hz, fields(tdiv_idxt)%DATA(:, :, 1))
      CALL this%VorTen(sg%cell_cntr, fields(vor_idxt)%DATA(:, :, 1), fields(div_idxt)%DATA(:, :, 1), &
                       fields(psi_idxt)%DATA(:, :, 1), fields(chi_idxt)%DATA(:, :, 1), &
                       this%DM, fields(tvor_idxt)%DATA(:, :, 1))
      CALL this%RhoTen(sg%cell_cntr, fields(rho_idxt)%DATA(:, :, 1), fields(div_idxt)%DATA(:, :, 1), &
                       fields(psi_idxt)%DATA(:, :, 1), fields(chi_idxt)%DATA(:, :, 1), &
                       this%DM, this%zHght, this%Hz, fields(trho_idxt)%DATA(:, :, 1), XHalf)
      CALL this%qTen(sg%cell_cntr, fields(psi_idxt)%DATA(:, :, 1), fields(chi_idxt)%DATA(:, :, 1), &
                     this%zHght, this%Hz, fields(div_idxt)%DATA(:, :, 1), fields(tq_idxt)%DATA(:, :, 1), XHalf)

    END ASSOCIATE

  END SUBROUTINE

  SUBROUTINE righthandsFull(this, t0, XFull)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XFull

    REAL(r_kind), INTENT(IN) :: t0

    INTEGER(i_kind) :: div_idxt, vor_idxt, psi_idxt, chi_idxt, q_idxt, rho_idxt, &
                       pres_idxt, w_idxt, theta_idxt, &
                       tdiv_idxt, tvor_idxt, tq_idxt, trho_idxt, &
                       tw_idxt, ttheta_idxt, &
                       pdiv_idxt, pvor_idxt, pq_idxt, prho_idxt, &
                       ptheta_idxt

    vor_idxt = XFull%getVarIdx('vort')
    div_idxt = XFull%getVarIdx('divt')
    rho_idxt = XFull%getVarIdx('rhot')
    q_idxt = XFull%getVarIdx('qvaport')
    pres_idxt = XFull%getVarIdx('prest')
    psi_idxt = XFull%getVarIdx('psit')
    chi_idxt = XFull%getVarIdx('chit')
    w_idxt = XFull%getVarIdx('wt')
    theta_idxt = XFull%getVarIdx('thetat')

    tw_idxt = XFull%getVarIdx('tenwt')
    ttheta_idxt = XFull%getVarIdx('tenthetat')

    ptheta_idxt = XFull%getVarIdx('parthetat')

    ASSOCIATE (sg => XFull%sg, &
               fields => XFull%fields)

      CALL this%RH_psi_func(t0, sg%cell_cntr, this%sigma3d, &
                            fields(psi_idxt)%DATA(:, :, 1), fields(vor_idxt)%DATA(:, :, 1))
      CALL this%RH_chi_func(t0, sg%cell_cntr, this%sigma3d, &
                            fields(chi_idxt)%DATA(:, :, 1), fields(div_idxt)%DATA(:, :, 1))

      CALL this%rho_func(t0, sg%cell_cntr, this%zHght, this%z_s, &
                         fields(rho_idxt)%DATA(:, :, 1))
      CALL this%pres_func(t0, sg%cell_cntr, this%zHght, this%z_s, &
                          fields(pres_idxt)%DATA(:, :, 1))

      CALL this%Tp_func(t0, sg%cell_cntr, this%zHght, fields(theta_idxt)%DATA(:, :, 1))

      CALL this%wTen(sg%cell_cntr, fields(rho_idxt)%DATA(:, :, 1), this%Hz, fields(tw_idxt)%DATA(:, :, 1))
      CALL this%TpTen(sg%cell_cntr, fields(psi_idxt)%DATA(:, :, 1), fields(chi_idxt)%DATA(:, :, 1), &
                      this%zHght, this%Hz, fields(div_idxt)%DATA(:, :, 1), fields(ttheta_idxt)%DATA(:, :, 1), XFull)
    END ASSOCIATE

  END SUBROUTINE

END MODULE
