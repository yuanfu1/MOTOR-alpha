MODULE RossbyHaurwitzSphere_m

  !>
  !!=================================================================
  !!  This module defines a Rossby Haurwitz wave
  !!
  !!  Note: Use Rossby Haurwitz function for stream function but
  !!  calculate the vorticity and thickness different from Williamson
  !!  et al. For details, see Yuanfu's note, Rossby_HaurwitzTest
  !!  under: /Users/xieyuanfu/developments/models/square/doc on
  !!  xieyuanfudeiMac, my desktop.
  !!
  !!  uthor Yuanfu Xie
  !!   History 2019-6 created.
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: Omega,g

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RossbyHaurwitzSphere_t

  TYPE   :: RossbyHaurwitzSphere_t
    REAL(r_kind) :: rh_K, rh_o, rh_R, rh_0, angular

    REAL(r_kind) :: ae ! debugging for Ning's ml_itmesh code. Deleted and used the one in parameters_m

    REAL(r_kind), ALLOCATABLE :: psi0th(:, :), zta0th(:, :), & ! 1st order
                                 psi1st(:, :), zta1st(:, :), & ! 1st order
                                 psi2nd(:, :), zta2nd(:, :), & ! 2nd order
                                 psi3rd(:, :)                ! 3rd order
  CONTAINS
    PROCEDURE, PUBLIC  :: initializ
    PROCEDURE, PUBLIC  :: deallocat
    PROCEDURE, PUBLIC  :: updateVar
    PROCEDURE, PUBLIC  :: Thickness
    PROCEDURE, PUBLIC  :: SWForcing
    PROCEDURE, PUBLIC  :: DivAdvect

    PROCEDURE, PRIVATE :: StreamFunDerivatives
    PROCEDURE, PRIVATE :: VorticityDerivatives
  END TYPE RossbyHaurwitzSphere_t

CONTAINS

  SUBROUTINE Initializ(this, numCel)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH stream function values
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this
    INTEGER(i_kind) :: numCel

    ! Debugging for Ning's ml_itmesh:
    this%ae = 6371.220D3

    this%rh_K = 7.848D-6  ! Rossby Haurwitz: K
    this%rh_o = 7.848D-6  ! Rossby Haurwitz: omega
    this%rh_R = 4.0D0     ! Rossby Haurwitz: R
    this%rh_0 = 5.0D3     ! Rossby Haurwitz: h0
    ! Angular velocity: (R(3+R)omega-2 Omega)/(1+R)/(2+R)
    this%angular = this%rh_R * (3.0D0 + this%rh_R) * this%rh_o - 2.0D0 * Omega
    this%angular = this%angular / (1.0D0 + this%rh_R) / (2.0D0 + this%rh_R)

    ! Internal variables saved:
    ALLOCATE (this%psi0th(1, numCel), this%zta0th(1, numCel), &
              this%psi1st(2, numCel), this%zta1st(2, numCel), &
              this%psi2nd(3, numCel), this%zta2nd(3, numCel), &
              this%psi3rd(4, numCel))

  END SUBROUTINE Initializ

  SUBROUTINE Deallocat(this)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH stream function values
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this

    IF (ALLOCATED(this%psi0th)) DEALLOCATE (this%psi0th)
    IF (ALLOCATED(this%zta0th)) DEALLOCATE (this%zta0th)

    IF (ALLOCATED(this%psi1st)) DEALLOCATE (this%psi1st)
    IF (ALLOCATED(this%zta1st)) DEALLOCATE (this%zta1st)
    IF (ALLOCATED(this%psi2nd)) DEALLOCATE (this%psi2nd)
    IF (ALLOCATED(this%zta2nd)) DEALLOCATE (this%zta2nd)
    IF (ALLOCATED(this%psi3rd)) DEALLOCATE (this%psi3rd)

  END SUBROUTINE Deallocat

  SUBROUTINE UpdateVar(this, numCel, latlon, atimes, stream, vortic, height,rhs)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH stream function and
  !!  relative vorticity values and the needed derivatives of them
  !!  for later thickness and its derivatives because of the need of
  !!  Poisson solutions.
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numCel
    REAL(r_kind), INTENT(IN) :: latlon(2, numCel), atimes
    REAL(r_kind), INTENT(OUT) :: stream(numCel), vortic(numCel), height(numCel)
    REAL(r_kind), INTENT(OUT) :: rhs(numCel,3)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind)    :: clat, slat, clon, slon, clatR(0:5)

    ! Debugging Ning's ml_itmesh code:
    REAL(r_kind) :: ae

    ae = this%ae

    ! Williamson et al. eqn (141):
    DO i = 1, numCel

      ! Trigonometry of latlon:
      clat = DCOS(latlon(1, i))
      slat = DSIN(latlon(1, i))
      clon = DCOS(this%rh_R * latlon(2, i) + this%angular * atimes)
      slon = DSIN(this%rh_R * latlon(2, i) + this%angular * atimes)

      clatR(0) = clat**(this%rh_R - 3.0D0)  ! c**(R-3)
      clatR(1) = clatR(0) * clat            ! c**(R-2)
      clatR(2) = clatR(1) * clat            ! c**(R-1)
      clatR(3) = clatR(2) * clat            ! c** R
      clatR(4) = clatR(3) * clat            ! c**(R+1)
      clatR(5) = clatR(4) * clat            ! c**(R+2)

      ! Stream function:
      this%psi0th(1, i) = &
        -ae * ae * slat * (this%rh_o - this%rh_K * clatR(3) * clon)

      ! Relative vorticity:
      this%zta0th(1, i) = slat * (2.0D0 * this%rh_o - this%rh_K * clatR(3) * &
                                  (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * clon)

        !!! Note: these derivatives are derived in the document of
        !!!       Rossby_HaurwitzTest under
        !!!         /Users/xieyuanfu/developments/models/square/doc
        !!!       they need to double check the doc and coding here.
        !!! For testing the Poisson solver on a sphere, I currently
        !!!
        !!! use the above stream function and vorticity only
        !!!

      ! The first order derivatives: 1: lat; 2: lon;
      this%psi1st(1, i) = -ae * ae * (this%rh_o * clat - this%rh_K * &
                                      ((this%rh_R + 1.0D0) * clatR(4) - this%rh_R * clatR(2)) * clon)
      this%psi1st(2, i) = -ae * ae * this%rh_K * this%rh_R * clatR(3) * slat * slon

      this%zta1st(1, i) = -2.0D0 * this%rh_o * clat - this%rh_K * &
                          ((this%rh_R + 1.0D0) * clatR(4) - this%rh_R * clatR(2)) * &
                          (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * clon
      this%zta1st(2, i) = this%rh_K * this%rh_R * slat * clatR(3) * &
                          (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * slon

      ! The second order derivatives: 1 lat_lat; 2 lat_lon; 3 lon_lon:
      this%psi2nd(1, i) = ae * ae * slat * (this%rh_o + this%rh_K * &
                                            ((this%rh_R * (this%rh_R - 1.0D0) * clatR(1) - &
                                              this%rh_R * (this%rh_R + 1.0D0) * clatR(3))) * clon)
      this%psi2nd(2, i) = -ae * ae * this%rh_K * &
                          ((this%rh_R + 1.0D0) * clat - this%rh_R * clatR(2)) * slon
      this%psi2nd(3, i) = -ae * ae * this%rh_K * this%rh_R**2 * clatR(3) * slat * clon

      this%zta2nd(1, i) = slat * (2.0D0 * this%rh_o + this%rh_K * &
                                  ((this%rh_R + 1.0D0)**2 * clatR(4) - this%rh_R * (this%rh_R - 1.0D0) * clatR(1)) * &
                                  (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * clon)
      this%zta2nd(2, i) = this%rh_K * this%rh_R * &
                          ((this%rh_R + 1.0D0) * clatR(4) - this%rh_R * clatR(2)) * &
                          (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * slon
      this%zta2nd(3, i) = this%rh_K * this%rh_R**2 * clatR(3) * &
                          (this%rh_R**2 + 3.0D0 * this%rh_R + 2.0D0) * clon

      ! The third order derivatives: 1 lat3; 2 lat2lon; 3 latlon2; 4 lon3
      this%psi3rd(1, i) = ae * ae * (this%rh_o * clat + this%rh_K * &
                                     (this%rh_R * (this%rh_R - 1.0D0) * clatR(1) - &
                                      (this%rh_R + 1.0D0) * this%rh_R * clatR(3)) * clat * clon - this%rh_K * &
                                     (this%rh_R * (this%rh_R - 1.0D0) * (this%rh_R - 2.0D0) * clatR(0) - &
                                      this%rh_R**2 * (this%rh_R + 1.0D0) * clatR(2)) * slat * slat * clon)
      this%psi3rd(2, i) = -ae * ae * slat * this%rh_K * this%rh_R**2 * &
                          ((this%rh_R - 1.0D0) * clatR(1) - (this%rh_R + 1.0D0) * clatR(3)) * slon
      this%psi3rd(3, i) = -ae * ae * this%rh_K * this%rh_R**2 * &
                          (-this%rh_R * clatR(2) + (this%rh_R + 1.0D0) * clatR(4)) * clon
      this%psi3rd(4, i) = ae * ae * this%rh_K * this%rh_R**3 * clatR(3) * slat * slon

      ! Right hand sides of a Z-grid shallow water equation:
      rhs(i,1) = ((this%zta1st(1,i)+2.0D0*DCOS(clat))*this%psi1st(2,i) - &
        this%zta1st(2,i)*this%psi1st(1,i))/ae/ae/DCOS(clat)
      rhs(i,2) = ((this%zta1st(1,i)+2.0D0*DCOS(clat))*this%psi1st(1,i) + &
        this%zta0th(1,i)*this%psi2nd(1,i))/ae/ae
      ! Temporary use of stream and vortic to save -h_theta and -h_lambda:
      stream(i) = (this%psi1st(2,i)*this%psi2nd(3,i)+this%psi1st(1,i)*this%psi2nd(2,i))/g/ae/ae
      vortic(i) = (this%psi1st(2,i)*this%psi2nd(2,i)+this%psi1st(1,i)*this%psi2nd(1,i))/g/ae/ae
      ! Height is not derived from zero divergence equation but canceling kinetic energy term:
      ! See Rossby_HaurwitzTest.docx
      height(i) = -0.5D0*(this%psi1st(2,i)**2+this%psi1st(1,i)**2)/g/ae/ae

      rhs(i,3) = (vortic(i)*this%psi1st(2,i)-stream(i)*this%psi1st(1,i))/ae/ae/clat
    END DO

    ! Save stream function and relative vorticity:
    stream = this%psi0th(1, 1:numCel)
    vortic = this%zta0th(1, 1:numCel)

  END SUBROUTINE UpdateVar

  SUBROUTINE Thickness(this, numCel, latlon, divAdv, thicks)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH thickness values for a
  !!  given gridded Laplacian inverse values of the divergence
  !!  advection term
  !!
  !!  Inputs:
  !!    divAdv: real array holding the
  !!
  !! abla^ (-2)
  !!  abla\dot(eta*nabla psi)
  !!  Output:
  !!    thicks: real array holding: h = -h_s - (K - divAdv)/g
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numCel
    REAL(r_kind), INTENT(IN) :: latlon(2, numCel), divAdv(numCel)
    REAL(r_kind), INTENT(OUT) :: thicks(numCel)
  END SUBROUTINE Thickness

  SUBROUTINE SWForcing(this, numCel, latlon, forces)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH forcing functions
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numCel
    REAL(r_kind), INTENT(IN) :: latlon(2, numCel)
    REAL(r_kind), INTENT(OUT) :: forces(2, numCel)
  END SUBROUTINE SWForcing

  SUBROUTINE DivAdvect(this, numCel, latlon, divAdv, thiLat, thiLon)
    !>
  !!=================================================================
  !!  This routine calculates the gridded RH divergence advection
  !!  values of divAdv, divAdv_lat and divAdv_lon as right hand sides
  !!  of Laplacian operator
  !!=================================================================
    !
    CLASS(RossbyHaurwitzSphere_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numCel
    REAL(r_kind), INTENT(IN) :: latlon(2, numCel)
    REAL(r_kind), INTENT(OUT) :: divAdv(numCel), &
                                 thiLat(numCel), thiLon(numCel)
  END SUBROUTINE DivAdvect

  ! Private procedures:
  SUBROUTINE StreamFunDerivatives(this)
    CLASS(RossbyHaurwitzSphere_t) :: this
  END SUBROUTINE StreamFunDerivatives

  SUBROUTINE VorticityDerivatives(this)
    CLASS(RossbyHaurwitzSphere_t) :: this
  END SUBROUTINE VorticityDerivatives

END MODULE RossbyHaurwitzSphere_m
