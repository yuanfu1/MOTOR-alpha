MODULE RossbyHaurwitz_m

  !>
  !!=================================================================
  !!  This module defines a Rossby Haurwitz wave
  !!
  !!  Note: Use Rossby Haurwitz function for stream function but
  !!  calculate the vorticity and thickness different from Williamson
  !!  et al. For details, see Yuanfu's note, Rossby_HaurwitzTest
  !!  under: /Users/xieyuanfu/developments/models/square/doc on
  !!  xieyuanfudeiMac, my desktop
  !!
  !!  \author Yuanfu Xie
  !!  \b History 2019-6 created.
  !!  \b Yuanfu Xie 2019-12 Modified to make it an abstract for both
  !!                                 rectangle and sphere
  !!  \b Yuanfu Xie 2020-09 Modified to add outputs of derivatives of
  !!                                 streamfunction and vorticity for
  !!                                 general debugging purposes.
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: Omega, EarthRadius, machineEps

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RossbyHaurwitz_t, RossbyHaurwitzSphere1_t

  ! Rossby-Haurwitz parameters:
  REAL(r_kind), PARAMETER :: &
    rh_kapa = 7.848D-6, &
    rh_omga = 7.848D-6, &
    rh_high = 8.0D3, &
    rh_bigR = 4.0D0, &
    rh_angl = &
    (rh_bigR * (3.0D0 + rh_bigR) * rh_omga - &
     2.0D0 * Omega) / (1.0D0 + rh_bigR) / (2.0D0 + rh_bigR)

  TYPE, ABSTRACT :: RossbyHaurwitz_t
    INTEGER(i_kind) :: num_cell
    REAL(r_kind), ALLOCATABLE :: &
      psi1st(:, :), zta1st(:, :), & ! 1st order derivatives
      psi2nd(:, :), zta2nd(:, :), & ! 2nd order derivatives
      psi3rd(:, :)                ! 3rd order derivatives

  CONTAINS
    PROCEDURE(get_values), PUBLIC, DEFERRED :: RHV_func

    PROCEDURE, PUBLIC :: RH_initial
    PROCEDURE, PUBLIC :: RH_destroy
  END TYPE RossbyHaurwitz_t

  ABSTRACT INTERFACE
    SUBROUTINE get_values(this, attime, latlon, stream, vortct, thicks, &
                          derivS, derivV, jacobi, flxdiv)
      IMPORT RossbyHaurwitz_t, r_kind
      CLASS(RossbyHaurwitz_t) :: this
      REAL(r_kind), INTENT(IN) :: attime, latlon(2, this%num_cell)
      REAL(r_kind), INTENT(OUT) :: stream(this%num_cell), &
                                   vortct(this%num_cell), &
                                   thicks(this%num_cell), &
                                   derivS(2, this%num_cell), &
                                   derivV(2, this%num_cell), &
                                   jacobi(this%num_cell), &
                                   flxdiv(this%num_cell)
    END SUBROUTINE get_values
  END INTERFACE

  ! This type is essentially identical to RossbyHaurwitzSphere_t
  ! but as a derived type for abstraction for including rectangles
  TYPE, EXTENDS(RossbyHaurwitz_t) :: RossbyHaurwitzSphere1_t
  CONTAINS
    PROCEDURE :: RHV_func => RHV_sphereFunc
  END TYPE RossbyHaurwitzSphere1_t

CONTAINS
  !>
    !!===============================================================
    !!  These routines are public for any general grids
    !!===============================================================
  !
  SUBROUTINE RH_initial(this, ncells)
    CLASS(RossbyHaurwitz_t) :: this
    INTEGER(i_kind), INTENT(IN) :: ncells

    this%num_cell = ncells

    ! Memory allocation:
    ALLOCATE (this%psi1st(2, ncells), this%zta1st(2, ncells), &
              this%psi2nd(3, ncells), this%zta2nd(3, ncells), &
              this%psi3rd(4, ncells))
  END SUBROUTINE RH_initial

  SUBROUTINE RH_destroy(this)
    CLASS(RossbyHaurwitz_t) :: this

    ! Memory deallocation:
    DEALLOCATE (this%psi1st, this%zta1st, &
                this%psi2nd, this%zta2nd, &
                this%psi3rd)
  END SUBROUTINE RH_destroy

  INCLUDE 'RossbyHaurwitzSphere1.F90'

END MODULE RossbyHaurwitz_m
