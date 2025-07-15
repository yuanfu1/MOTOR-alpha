SUBROUTINE RHV_sphereFunc(this, attime, latlon, stream, vortct, &
                          thicks, derivS, derivV, jacobi, flxdiv)

  !>
  !!=================================================================
  !!  This routine is a member of RossbyHaurwitzSphere1_t data type.
  !!
  !!  Inputs:
  !!    attime: time in seconds
  !!    latlon: latitude and longitude of all grid cells
  !!
  !!  Output:
  !!    values: in num_cell x num_var
  !!            num_var
  !!            1. stream:    streamfunction
  !!            2. vortct:    vorticity
  !!            3. thicks:    thickness
  !!            4. derivS:    partial derivatives of stream function
  !!                          in terms of latitude and longitude
  !!            5. derivV:    partial derivatives of vorticity
  !!            6. jacobi:    jacobian of eta and psi, J(eta,psi)
  !!            7. flxdiv:    flux divergence of eta and psi,
  !!               \nabla dot (eta \nabla psi)
  !!
  !!  \author Yuanfu Xie
  !!  \b History
  !!  \b 2019-12 created by Yuanfu Xie
  !!  \b 2021-02 modified by Yuanfu Xie for the header context with
  !!             more detailed output variable names
  !!=================================================================
  !
  IMPLICIT NONE

  CLASS(RossbyHaurwitzSphere1_t) :: this
  REAL(r_kind), INTENT(IN) :: attime, latlon(2, this%num_cell)
  REAL(r_kind), INTENT(OUT) :: stream(this%num_cell), &
                               vortct(this%num_cell), &
                               thicks(this%num_cell), &
                               derivS(2, this%num_cell), &
                               derivV(2, this%num_cell), &
                               jacobi(this%num_cell), &
                               flxdiv(this%num_cell)

  ! Local variables:
  INTEGER(i_kind) :: i
  REAL(r_kind)    :: clat, slat, tlat, clon, slon, clatR(0:5)

  ! Williamson et al. eqn (pg. 141)
  DO i = 1, this%num_cell

    ! Trigonometry of latlon:
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
    stream(i) = &
      -EarthRadius * EarthRadius * slat * (rh_omga - rh_kapa * clatR(3) * clon)

    ! Relative vorticity:
    vortct(i) = slat * (2.0D0 * rh_omga - rh_kapa * clatR(3) * &
                        (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon)

    thicks(i) = 0.0D4

    !!! Note: these derivatives are derived in the document of
    !!!       Rossby_HaurwitzTest under
    !!!         /Users/xieyuanfu/developments/models/square/doc
    !!!       they need to double check the doc and coding here.
    !!! For testing the Poisson solver on a sphere, I currently
    !!! use the above stream function and vorticity only

    ! The first order derivatives: 1: lat; 2: lon;
    this%psi1st(1, i) = -EarthRadius * EarthRadius * (rh_omga * clat - rh_kapa * &
                                                      ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * clon)
    this%psi1st(2, i) = -EarthRadius * EarthRadius * rh_kapa * rh_bigR * clatR(3) * slat * slon

    this%zta1st(1, i) = 2.0D0 * rh_omga * clat - rh_kapa * &
                        ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                        (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon
    this%zta1st(2, i) = rh_kapa * rh_bigR * slat * clatR(3) * &
                        (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon

    ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
    this%psi2nd(1, i) = EarthRadius * EarthRadius * slat * (rh_omga + rh_kapa * &
                                                            ((rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                                              (rh_bigR + 1.0D0)**2 * clatR(3))) * clon)
    this%psi2nd(2, i) = -EarthRadius * EarthRadius * rh_kapa * &
                        ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * slon
    this%psi2nd(3, i) = -EarthRadius * EarthRadius * rh_kapa * rh_bigR**2 * clatR(3) * slat * clon

    this%zta2nd(1, i) = slat * (-2.0D0 * rh_omga + rh_kapa * &
                                ((rh_bigR + 1.0D0)**2 * clatR(4) - rh_bigR * (rh_bigR - 1.0D0) * clatR(1)) * &
                                (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon)
    this%zta2nd(2, i) = rh_kapa * rh_bigR * &
                        ((rh_bigR + 1.0D0) * clatR(4) - rh_bigR * clatR(2)) * &
                        (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * slon
    this%zta2nd(3, i) = rh_kapa * rh_bigR**2 * clatR(3) * &
                        (rh_bigR**2 + 3.0D0 * rh_bigR + 2.0D0) * clon

    ! The third order derivatives: 1 lat3; 2 lat2lon; 3 latlon2; 4 lon3
    this%psi3rd(1, i) = EarthRadius * EarthRadius * (rh_omga * clat + rh_kapa * &
                                                     (rh_bigR * (rh_bigR - 1.0D0) * clatR(1) - &
                                                      (rh_bigR + 1.0D0) * rh_bigR * clatR(3)) * clat * clon - rh_kapa * &
                                                     (rh_bigR * (rh_bigR - 1.0D0) * (rh_bigR - 2.0D0) * clatR(0) - &
                                                      rh_bigR**2 * (rh_bigR + 1.0D0) * clatR(2)) * slat * slat * clon)
    this%psi3rd(2, i) = -EarthRadius * EarthRadius * slat * rh_kapa * rh_bigR**2 * &
                        ((rh_bigR - 1.0D0) * clatR(1) - (rh_bigR + 1.0D0) * clatR(3)) * slon
    this%psi3rd(3, i) = -EarthRadius * EarthRadius * rh_kapa * rh_bigR**2 * &
                        (-rh_bigR * clatR(2) + (rh_bigR + 1.0D0) * clatR(4)) * clon
    this%psi3rd(4, i) = EarthRadius * EarthRadius * rh_kapa * rh_bigR**3 * clatR(3) * slat * slon

    ! Jacobian and flux divergence:
    IF (ABS(clat) .GT. machineEps) THEN
      ! Use equation (39) for divergence, (11) or (43) for Jacobian
      ! in the document formula
      jacobi(i) = (this%zta1st(2, i) * this%psi1st(1, i) - &
                   (this%zta1st(1, i) + 2.0D0 * Omega * clat) * this%psi1st(2, i)) &
                  / clat / EarthRadius / EarthRadius
      flxdiv(i) = ((this%zta1st(1, i) + 2.0D0 * Omega * clat) * this%psi1st(1, i) + &
                   (vortct(i) + 2.0D0 * Omega * slat) * this%psi2nd(1, i) - &
                   (vortct(i) + 2.0D0 * Omega * slat) * this%psi1st(1, i) * tlat + &
                   (this%zta1st(2, i) * this%psi1st(2, i) + &
                    (vortct(i) + 2.0D0 * Omega * slat) * this%psi2nd(3, i)) / clat / clat) &
                  / EarthRadius / EarthRadius
    END IF

  END DO
  ! Output the derivatives of stream function and vorticity:
  derivS = this%psi1st
  derivV = this%zta1st

END SUBROUTINE RHV_sphereFunc
