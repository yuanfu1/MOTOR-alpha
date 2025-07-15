MODULE UnitTestDyn_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: Omega, EarthRadius, machineEps, pi
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
     2.0D0 * Omega) / (1.0D0 + rh_bigR) / (2.0D0 + rh_bigR)

  TYPE :: UnitTestDyn_t
    INTEGER(i_kind), POINTER :: num_cell, vLevel
    REAL(r_kind), POINTER :: sigma(:), z_top, f(:, :)
    REAL(r_kind), ALLOCATABLE :: &
      psi1st(:, :, :), zta1st(:, :, :), & ! 1st order derivatives of psi
      psi1st_dclat(:, :, :), zta1st_dclat(:, :, :), & ! 1st order derivatives devided by cos(lat)
      psi2nd(:, :, :), & ! 2nd order derivatives of psi
      psiThetaLambda_dclat(:, :), & ! 2nd order of psi_theta_lambda devided by cos(lat)
      psi2ndLambda_dclat(:, :), & ! 2nd order of psi_lambda devided by cos(lat)
      chi1st(:, :, :), del1st(:, :, :), & ! 1st order derivatives of chi
      chi1st_dclat(:, :, :), del1st_dclat(:, :, :), & ! 1st order derivatives devided by cos(lat)
      chi2nd(:, :, :), & ! 2nd order derivatives of chi
      chiThetaLambda_dclat(:, :), & ! 2nd order of chi_theta_lambda divided by cos(lat)
      chi2ndLambda_dclat(:, :), & ! 2nd order of chi_lambda devided by cos(lat)
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
      rho1sttheta(:, :), rho1stsigma(:, :) !1st order derivatives of rho

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
    PROCEDURE :: DivTen
  END TYPE

CONTAINS
  SUBROUTINE initial(this, sg)
    CLASS(UnitTestDyn_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    INTEGER(i_kind), POINTER :: num_cell, vLevel

    this%num_cell => sg%num_cell
    this%vLevel => sg%vLevel
    this%sigma => sg%sigma
    this%z_top => sg%ztop
    num_cell => sg%num_cell
    vLevel => sg%vLevel
    IF (ALLOCATED(sg%f)) this%f => sg%f

    PRINT *, 'ztop in unitestdyn is: ', this%z_top

    ALLOCATE (this%psi1st(3, vLevel, num_cell), this%zta1st(3, vLevel, num_cell), &
              this%psi1st_dclat(2, vLevel, num_cell), this%zta1st_dclat(2, vLevel, num_cell), &
              this%psi2nd(3, vLevel, num_cell), &
              this%psiThetaLambda_dclat(vLevel, num_cell), &
              this%psi2ndLambda_dclat(vLevel, num_cell), &
              this%chi1st(3, vLevel, num_cell), this%del1st(3, vLevel, num_cell), &
              this%chi1st_dclat(2, vLevel, num_cell), this%del1st_dclat(2, vLevel, num_cell), &
              this%chi2nd(3, vLevel, num_cell), &
              this%chiThetaLambda_dclat(vLevel, num_cell), &
              this%chi2ndLambda_dclat(vLevel, num_cell), &
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
              this%pres2ndsigmatheta(vLevel, num_cell))

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
                  this%psi1st_dclat, this%zta1st_dclat, &
                  this%psi2nd, &
                  this%psiThetaLambda_dclat, &
                  this%psi2ndLambda_dclat, &
                  this%chi1st, this%del1st, &
                  this%chi1st_dclat, this%del1st_dclat, &
                  this%chi2nd, &
                  this%chiThetaLambda_dclat, &
                  this%chi2ndLambda_dclat, &
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
                  this%pres2ndsigmatheta)

    END IF

  END SUBROUTINE destroy

  SUBROUTINE RH_psi_func(this, attime, latlon, zhght, stream, vortct)
    IMPLICIT NONE

    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  attime, zhght(:, :), latlon(:, :)
    REAL(r_kind), INTENT(OUT) :: stream(:, :), &
                                 vortct(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: R0, clat, slat, tlat, clon, slon, clatR(0:5), &
                       szhght, czhght, z_ratio, zfun, zfun2, parzfun_sigma

    z_ratio = 0.1D0
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        R0 = EarthRadius + zhght(k, i)
        szhght = DSIN(zhght(k, i) / this%z_top * 2.0D0 * pi)
        czhght = DCOS(zhght(k, i) / this%z_top * 2.0D0 * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = (2.0D0 / R0 + &
                 z_ratio * czhght * 2.0D0 * pi / this%z_top / zfun)
        parzfun_sigma = z_ratio * czhght * 2.0D0 * pi / this%z_top
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

        ! Stream function:
        stream(k, i) = &
          -R0 * R0 * zfun * slat * (rh_omga - rh_kapa * clatR(3) * clon)

        ! Relative vorticity:
        vortct(k, i) = zfun * slat * (2.0D0 * rh_omga - rh_kapa * clatR(3) * &
                                      (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon)

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
        this%psi1st(2, k, i) = -R0 * R0 * zfun * &
                               (rh_omga * clat - rh_kapa * &
                                ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * clon)
        this%psi1st(3, k, i) = -R0 * R0 * zfun &
                               * rh_kapa * rh_bigR * clatR(3) * slat * slon
        this%psi1st_dclat(1, k, i) = -R0 * R0 * zfun * &
                                     (rh_omga - rh_kapa * &
                                      ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * clon)
        this%psi1st_dclat(2, k, i) = -R0 * R0 * zfun &
                                     * rh_kapa * rh_bigR * clatR(2) * slat * slon

        this%zta1st(1, k, i) = (2.0D0 * rh_omga * clat - rh_kapa * &
                                ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                                (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                               * zfun
        this%zta1st(2, k, i) = rh_kapa * rh_bigR * slat * clatR(3) * &
                               (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                               * zfun
        this%zta1st(3, k, i) = parzfun_sigma / zfun * vortct(k, i)

        this%zta1st_dclat(1, k, i) = (2.0D0 * rh_omga - rh_kapa * &
                                      ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * &
                                      (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon) &
                                     * zfun
        this%zta1st_dclat(2, k, i) = rh_kapa * rh_bigR * slat * clatR(2) * &
                                     (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon &
                                     * zfun

        ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
        this%psi2nd(1, k, i) = R0 * R0 * slat &
                               * zfun &
                               * (rh_omga + rh_kapa * &
                                  ((rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                    (rh_bigR + 1.0D0)**2 * clatR(3))) * clon)
        this%psi2nd(2, k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                               zfun * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * slon
        this%psi2nd(3, k, i) = -R0 * R0 * rh_kapa * rh_bigR**2 * clatR(3) * slat * clon &
                               * zfun

        this%psiThetaLambda_dclat(k, i) = -R0 * R0 * rh_kapa * rh_bigR * &
                                          ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon &
                                          * zfun
        this%psi2ndLambda_dclat(k, i) = -R0 * R0 * rh_kapa * rh_bigR**2 * clatR(2) * slat * clon &
                                        * zfun

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

  SUBROUTINE RH_chi_func(this, attime, latlon, zhght, potential, divergence)
    IMPLICIT NONE

    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  attime, zhght(:, :), latlon(:, :)
    REAL(r_kind), INTENT(OUT) :: potential(:, :), &
                                 divergence(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: R0, clat, slat, tlat, clon, slon, clatR(0:5), &
                       rh_kapa_t, rh_omga_t, rh_angl_t, &
                       szhght, czhght, z_ratio, zfun, zfun2, parzfun_sigma

    rh_kapa_t = rh_kapa / 500.0D0
    rh_omga_t = rh_omga / 500.0D0
    rh_angl_t = &
      (rh_bigR * (3.0D0 + rh_bigR) * rh_omga_t - &
       2.0D0 * Omega) / (1.0D0 + rh_bigR) / (2.0D0 + rh_bigR)

    z_ratio = 0.1D0
    ! Williamson et al. eqn (pg. 141)
    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        R0 = EarthRadius + zhght(k, i)
        szhght = DSIN(zhght(k, i) / this%z_top * 2.0D0 * pi)
        czhght = DCOS(zhght(k, i) / this%z_top * 2.0D0 * pi)
        zfun = (1.0D0 + z_ratio * szhght)
        zfun2 = (2.0D0 / R0 + &
                 z_ratio * czhght * 2.0D0 * pi / this%z_top / zfun)
        parzfun_sigma = z_ratio * czhght * 2.0D0 * pi / this%z_top
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
          -R0 * R0 * zfun * slat * (rh_omga_t - rh_kapa_t * clatR(3) * clon)

        ! divergence:
        divergence(k, i) = zfun * slat * (2.0D0 * rh_omga_t - rh_kapa_t * clatR(3) * &
                                          (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon)

        ! The first order derivatives: for the last dimension, chi: (1) sigma, (2) lat, (3) lon; del: (1) lat, (2) lon, (3) sigma
        ! *_dclat are parameters divived by cos(lat)
        this%chi1st(1, k, i) = zfun2 * potential(k, i)
        this%chi1st(2, k, i) = -R0 * R0 * zfun * (rh_omga_t * clat - rh_kapa_t * &
                                                  ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * clon)
        this%chi1st(3, k, i) = -R0 * R0 * zfun * rh_kapa_t * rh_bigR * clatR(3) * slat * slon
        this%chi1st_dclat(1, k, i) = -R0 * R0 &
                                     * zfun * (rh_omga_t - rh_kapa_t * &
                                               ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * clon)
        this%chi1st_dclat(2, k, i) = -R0 * R0 * zfun * rh_kapa_t * rh_bigR * clatR(2) * slat * slon

        this%del1st(1, k, i) = (2.0D0 * rh_omga_t * clat - rh_kapa_t * &
                                ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                                (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon) * zfun
        this%del1st(2, k, i) = zfun * rh_kapa_t * rh_bigR * slat * clatR(3) * &
                               (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon
        this%del1st(3, k, i) = parzfun_sigma / zfun * divergence(k, i)

        this%del1st_dclat(1, k, i) = (2.0D0 * rh_omga_t - rh_kapa_t * &
                                      ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * &
                                      (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon) * zfun
        this%del1st_dclat(2, k, i) = zfun * rh_kapa_t * rh_bigR * slat * clatR(2) * &
                                     (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon

        ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
        this%chi2nd(1, k, i) = R0 * R0 * zfun * slat * (rh_omga_t + rh_kapa_t * &
                                                        ((rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                                          (rh_bigR + 1.0D0)**2 * clatR(3))) * clon)
        this%chi2nd(2, k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * zfun * &
                               ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * slon
        this%chi2nd(3, k, i) = -R0 * R0 * zfun * rh_kapa_t * rh_bigR**2 * clatR(3) * slat * clon

        this%chiThetaLambda_dclat(k, i) = -R0 * R0 * rh_kapa_t * rh_bigR * zfun * &
                                          ((rh_bigR + 1.0D0) * clatR(3) - rh_bigR * clatR(1)) * slon
        this%chi2ndLambda_dclat(k, i) = -R0 * R0 * zfun * rh_kapa_t * rh_bigR**2 * clatR(2) * slat * clon

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

    omg = 1000.0D0

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

        IF (k .EQ. 1) THEN
          z_s(i) = omg * clatR(3) * slat_p2 + 1.0D0
          this%zs1sttheta(i) = omg * slat * clatR(2) * Rclat_t
          this%zs2ndtheta(i) = -(rh_bigR - 1.0D0) * omg * clatR(1) * slat_p2 * Rclat_t &
                               + omg * clatR(3) * Rclat_t &
                               - 2 * omg * (rh_bigR + 2.0D0) * clatR(3) * slat_p2
          this%zs1sttheta_dclat(i) = omg * slat * clatR(1) * Rclat_t
          this%zs2ndtheta_dclat(i) = -(rh_bigR - 1.0D0) * omg * clatR(0) * slat_p2 * Rclat_t &
                                     + omg * clatR(2) * Rclat_t &
                                     - 2 * omg * (rh_bigR + 2.0D0) * clatR(2) * slat_p2
        END IF

        z(k, i) = this%sigma(k) * (this%z_top - z_s(i)) / this%z_top + z_s(i)! + (1.0D0 - this%sigma(k)/this%z_top)*z_s(i)
        Hz(k, i) = this%z_top / (this%z_top - z_s(i))

        this%z1sttheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * omg * slat * clatR(2) * Rclat_t

        this%z1sttheta_dclat(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * omg * slat * clatR(1) * Rclat_t

        this%d_z1stthetadclat_dtheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * omg * &
                                             ((2.0D0 - rh_bigR) * clatR(0) * slat_p2 * Rclat_t + &
                                              clatR(2) * Rclat_t - &
                                              2.0D0 * (rh_bigR + 2.0D0) * clatR(2) * slat_p2)

        this%z2ndtheta(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * omg * &
                               ((1.0D0 - rh_bigR) * clatR(1) * slat_p2 * Rclat_t + &
                                clatR(3) * Rclat_t - &
                                2.0D0 * (rh_bigR + 2.0D0) * clatR(3) * slat_p2)
        this%z2ndtheta_dclat(k, i) = (1.0D0 - this%sigma(k) / this%z_top) * omg * &
                                     ((1.0D0 - rh_bigR) * clatR(0) * slat_p2 * Rclat_t + &
                                      clatR(2) * Rclat_t - &
                                      2.0D0 * (rh_bigR + 2.0D0) * clatR(2) * slat_p2)

        this%Hz1sttheta(k, i) = (this%z_top / (this%z_top - z_s(i))**2 * slat * clatR(2) * &
                                 ((rh_bigR + 2.0) * clat**2 - rh_bigR)) * omg
        this%Hz1sttheta_dclat(k, i) = (this%z_top / (this%z_top - z_s(i))**2 * slat * clatR(2) * &
                                       ((rh_bigR + 2.0) * clat**2 - rh_bigR)) * omg
      END DO
    END DO

  END SUBROUTINE Hght_func

  SUBROUTINE pres_func(this, latlon, z, z_s, pres)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), z(:, :), z_s(:)
    REAL(r_kind), INTENT(OUT) :: pres(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, P0, alpha, R0, clat, slat, slat_p2, clon, &
                       slon, clatR(0:5), Rclat_t, expz, zratio

    omg = 1.5D-2
    P0 = 1.0D5
    alpha = 1.073D-4

    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        ! R0 = EarthRadius + z(k, i)
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

        expz = EXP(-alpha * z(k, i))
        zratio = (1.0D0 - this%sigma(k) / this%z_top)

        pres(k, i) = P0 * (1.0D0 + omg * clatR(3) * slat) * expz
        this%pres1stsigma(k, i) = -alpha * pres(k, i) * (this%z_top - z_s(i)) / this%z_top
        this%pres1sttheta(k, i) = P0 * omg * (-rh_bigR * clatR(2) + (rh_bigR + 1.0D0) * clatR(4)) * expz &
                                  - alpha * pres(k, i) * zratio * this%zs1sttheta(i)
        this%pres1sttheta_dclat(k, i) = P0 * omg * (-rh_bigR * clatR(1) + (rh_bigR + 1.0D0) * clatR(3)) * expz &
                                        - alpha * pres(k, i) * zratio * this%zs1sttheta_dclat(i)
        this%pres2ndtheta(k, i) = P0 * omg * (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) * slat - (rh_bigR + 1.0D0)**2 * clatR(3) * slat) * expz &
                                  - alpha * P0 * omg * (-rh_bigR * clatR(2) + (rh_bigR + 1.0D0) * clatR(4)) * expz * zratio * this%zs1sttheta(i) &
                                  - alpha * zratio * this%pres1sttheta(k, i) * this%zs1sttheta(i) &
                                  - alpha * pres(k, i) * zratio * this%zs2ndtheta(i)
        this%pres2ndtheta_dclat(k, i) = P0 * omg * (rh_bigR * (rh_bigR - 1.0D0) * clatR(0) * slat - (rh_bigR + 1.0D0)**2.0D0 * clatR(2) * slat) * expz &
                                        - alpha * P0 * omg * (-rh_bigR * clatR(1) + (rh_bigR + 1.0D0) * clatR(3)) * expz * zratio * this%zs1sttheta(i) &
                                        - alpha * zratio * this%pres1sttheta(k, i) * this%zs1sttheta_dclat(i) &
                                        - alpha * pres(k, i) * zratio * this%zs2ndtheta_dclat(i)
        this%pres2ndsigmatheta(k, i) = -alpha * (this%z_top - z_s(i)) / this%z_top * this%pres1sttheta(k, i) &
                                       + alpha * pres(k, i) / this%z_top * this%zs1sttheta(i)

      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE rho_func(this, latlon, z, z_s, rho)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    REAL(r_kind), INTENT(IN) ::  latlon(:, :), z(:, :), z_s(:)
    REAL(r_kind), INTENT(OUT) :: rho(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k
    REAL(r_kind)    :: omg, rho0, alpha, R0, clat, slat, slat_p2, clon, &
                       slon, clatR(0:5), Rclat_t, expz, zratio

    omg = 1.5D-2
    rho0 = 1.25D0
    alpha = 1.073D-4

    DO k = 1, this%vLevel
      DO i = 1, this%num_cell
        ! Trigonometry of latlon:
        ! R0 = EarthRadius + z(k, i)
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

        expz = EXP(-alpha * z(k, i))
        zratio = (1.0D0 - this%sigma(k) / this%z_top)
        rho(k, i) = rho0 * (1.0D0 + omg * slat) * expz
        this%rho1stsigma(k, i) = -alpha * rho(k, i) * (this%z_top - z_s(i)) / this%z_top
        this%rho1sttheta(k, i) = rho0 * omg * clat * expz - alpha * rho(k, i) * zratio * this%zs1sttheta(i)

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

  SUBROUTINE VorTen(this, latlon, vortct, divergence, psi, chi, DM, tenvor, X1)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                vortct(:, :), divergence(:, :), DM(:, :)
    REAL(r_kind), INTENT(OUT) :: tenvor(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat

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
                              + this%zta1st(3, k, i) &
                              * DM(k, i)

        J_DM_chisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(2, k, i) * this%chi1st_dclat(1, k, i) &
                              - this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(1, k, i) * this%chi1st_dclat(2, k, i)
        tenvor(k, i) = J_eta_psi(k, i) - F_eta_chi(k, i) + F_DM_psisigma(k, i) + J_DM_chisigma(k, i)
      END DO
    END DO

    X1%fields(X1%getVarIdx('J_eta_psit'))%DATA(:, :, 1) = J_eta_psi
    X1%fields(X1%getVarIdx('F_eta_chit'))%DATA(:, :, 1) = F_eta_chi
    X1%fields(X1%getVarIdx('F_DM_psisigmat'))%DATA(:, :, 1) = F_DM_psisigma
    X1%fields(X1%getVarIdx('J_DM_chisigmat'))%DATA(:, :, 1) = J_DM_chisigma

    DEALLOCATE (J_eta_psi, F_eta_chi, F_DM_psisigma, &
                J_DM_chisigma)

  END SUBROUTINE

  SUBROUTINE DivTen(this, latlon, divergence, vortict, psi, chi, DM, pres, rho, z, Hz, tendiv, X1)
    IMPLICIT NONE
    CLASS(UnitTestDyn_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X1
    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
                                divergence(:, :), DM(:, :), pres(:, :), rho(:, :), &
                                Hz(:, :), z(:, :), vortict(:, :)
    REAL(r_kind), INTENT(OUT) :: tendiv(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i, k

    REAL(r_kind) :: R0, clat, slat, f_theta

    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: J_del_psi, F_del_chi, J_DM_psisigma, &
                                                  F_DM_chisigma, F_invrho_P, &
                                                  F_HzinvRhoPsigma_z, F_f_psi, J_f_chi
    ALLOCATE (J_del_psi(this%vLevel, this%num_cell), &
              F_del_chi(this%vLevel, this%num_cell), &
              J_DM_psisigma(this%vLevel, this%num_cell), &
              F_DM_chisigma(this%vLevel, this%num_cell), &
              F_invrho_P(this%vLevel, this%num_cell), &
              F_HzinvRhoPsigma_z(this%vLevel, this%num_cell), &
              F_f_psi(this%vLevel, this%num_cell), &
              J_f_chi(this%vLevel, this%num_cell))

    DO k = 1, this%vLevel
      R0 = EarthRadius + this%sigma(k)
      DO i = 1, this%num_cell
        clat = DCOS(latlon(1, i))
        slat = DSIN(latlon(1, i))

        J_del_psi(k, i) = (this%del1st_dclat(2, k, i) * this%psi1st(2, k, i) - &
                           this%del1st(1, k, i) * this%psi1st_dclat(2, k, i)) &
                          / R0 / R0

        F_del_chi(k, i) = (this%del1st(1, k, i) * this%chi1st(2, k, i) + &
                           this%del1st_dclat(2, k, i) * this%chi1st_dclat(2, k, i)) / R0 / R0 + &
                          divergence(k, i)**2.0D0
        F_DM_chisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(1, k, i) * this%chi1st(2, k, i) &
                              + this%zfunct(2, k, i) &
                              / R0**2 * this%DM1stLambda_dclat(k, i) * this%chi1st_dclat(2, k, i) &
                              + this%zfunct(2, k, i) &
                              * DM(k, i) * divergence(k, i)

        J_DM_psisigma(k, i) = this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(2, k, i) * this%psi1st_dclat(1, k, i) &
                              - this%zfunct(2, k, i) &
                              / R0**2 * this%DM1st(1, k, i) * this%psi1st_dclat(2, k, i)

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
        F_f_psi(k, i) = f_theta / R0 / R0 * this%psi1st(2, k, i) &
                        + this%f(k, i) * vortict(k, i)

        J_f_chi(k, i) = -f_theta / R0 / R0 * this%chi1st_dclat(2, k, i)

      END DO
    END DO

    tendiv = J_del_psi - F_del_chi - J_DM_psisigma + F_DM_chisigma &
             - F_invrho_P + F_HzinvRhoPsigma_z &
             + F_f_psi + J_f_chi

    X1%fields(X1%getVarIdx('J_del_psit'))%DATA(:, :, 1) = J_del_psi
    X1%fields(X1%getVarIdx('F_del_chit'))%DATA(:, :, 1) = F_del_chi
    X1%fields(X1%getVarIdx('J_DM_psisigmat'))%DATA(:, :, 1) = J_DM_psisigma
    X1%fields(X1%getVarIdx('F_DM_chisigmat'))%DATA(:, :, 1) = F_DM_chisigma
    X1%fields(X1%getVarIdx('F_invRho_Pt'))%DATA(:, :, 1) = F_invrho_P
    X1%fields(X1%getVarIdx('F_HzinvRhoPsigma_zt'))%DATA(:, :, 1) = F_HzinvRhoPsigma_z
    X1%fields(X1%getVarIdx('F_f_psit'))%DATA(:, :, 1) = F_f_psi
    X1%fields(X1%getVarIdx('J_f_chit'))%DATA(:, :, 1) = J_f_chi

    DEALLOCATE (J_del_psi, &
                F_del_chi, &
                J_DM_psisigma, &
                F_DM_chisigma, &
                F_invrho_P, &
                F_HzinvRhoPsigma_z, &
                F_f_psi, &
                J_f_chi)

  END SUBROUTINE

  ! SUBROUTINE RhoTen(this, latlon, rho, divergence, psi, chi, DM, tenrho, X1)
  !    IMPLICIT NONE
  !    CLASS(UnitTestDyn_t) :: this
  !    TYPE(State_t), INTENT(INOUT) :: X1
  !    REAL(r_kind), INTENT(IN) :: latlon(:, :), psi(:, :), chi(:, :), &
  !                                rho(:, :), divergence(:, :), DM(:, :)
  !    REAL(r_kind), INTENT(OUT) :: tenrho(:, :)

  !    ! Local variables:
  !    INTEGER(i_kind) :: i, k

  !    REAL(r_kind) :: R0, clat, slat

  !    REAL(r_kind), ALLOCATABLE, dimension(:, :) :: J_rho_psi, F_rho_chi, F_z_chisigma, &
  !                                                  J_z_psisigma
  !    ALLOCATE (J_rho_psi(this%vLevel, this%num_cell), &
  !              F_rho_chi(this%vLevel, this%num_cell), &
  !              F_z_chisigma(this%vLevel, this%num_cell), &
  !              J_z_psisigma(this%vLevel, this%num_cell))

  !    DO k = 1, this%vLevel
  !       R0 = EarthRadius + this%sigma(k)
  !       DO i = 1, this%num_cell
  !          clat = DCOS(latlon(1, i))
  !          slat = DSIN(latlon(1, i))

  !          J_rho_psi(k, i) = (this%rho1sttheta(k, i)*this%psi1st_dclat(2, k, i) - &
  !                             (this%zta1st(1, k, i) + 2.0D0*Omega*clat)*this%psi1st_dclat(2, k, i)) &
  !                            /R0/R0

  !          F_eta_chi(k, i) = ((this%zta1st(1, k, i) + 2.0D0*Omega*clat)*this%chi1st(2, k, i) + &
  !                             this%zta1st_dclat(2, k, i)*this%chi1st_dclat(2, k, i))/R0/R0 + &
  !                            (vortct(k, i) + 2.0D0*Omega*slat)*divergence(k, i)
  !          F_DM_psisigma(k, i) = this%zfunct(2, k, i) &
  !                                /R0**2*this%DM1st(1, k, i)*this%psi1st(2, k, i) &
  !                                + this%zfunct(2, k, i) &
  !                                /R0**2*this%DM1stLambda_dclat(k, i)*this%psi1st_dclat(2, k, i) &
  !                                + this%zfunct(2, k, i) &
  !                                *DM(k, i)*vortct(k, i)

  !          J_DM_chisigma(k, i) = this%zfunct(2, k, i) &
  !                                /R0**2*this%DM1st(2, k, i)*this%chi1st_dclat(1, k, i) &
  !                                - this%zfunct(2, k, i) &
  !                                /R0**2*this%DM1st(1, k, i)*this%chi1st_dclat(2, k, i)
  !          tenvor(k, i) = J_eta_psi(k, i) - F_eta_chi(k, i) + F_DM_psisigma(k, i) + J_DM_chisigma(k, i)
  !       END DO
  !    END DO

  !    X1%fields(X1%getVarIdx('J_eta_psit'))%data(:, :, 1) = J_eta_psi
  !    X1%fields(X1%getVarIdx('F_eta_chit'))%data(:, :, 1) = F_eta_chi
  !    X1%fields(X1%getVarIdx('F_DM_psisigmat'))%data(:, :, 1) = F_DM_psisigma
  !    X1%fields(X1%getVarIdx('J_DM_chisigmat'))%data(:, :, 1) = J_DM_chisigma

  !    DEALLOCATE (J_eta_psi, F_eta_chi, F_DM_psisigma, &
  !                J_DM_chisigma)

  ! END SUBROUTINE

END MODULE
