!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.possionSolver_Kp.psSolutn_adjoint
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTHOR(S)         : [Your Name]
! VERSION           : V 1.0
! HISTORY           :
!   Created by [Your Name], [Date], @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the adjoint data type for model states.
!! It extends the original psSolutn_t type by adding adjoint variables.
!! @note
!! Since the original module does not contain computational subroutines that require adjoint counterparts,
!! the adjoint module focuses on managing adjoint variables for `rights` and `solutn`.
MODULE psSolutn_adjoint_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE psSolutn_m, ONLY: psSolutn_t

  IMPLICIT NONE

  TYPE, EXTENDS(psSolutn_t) :: psSolutn_adj_t
    ! Adjoint variables corresponding to `rights` and `solutn`
    REAL(r_kind), ALLOCATABLE :: adj_rights(:, :)
    REAL(r_kind), ALLOCATABLE :: adj_solutn(:, :)
  CONTAINS
    PROCEDURE :: constructor_adj
    FINAL     :: destructor_adj
  END TYPE psSolutn_adj_t

CONTAINS

  !> @brief
  !! Adjoint constructor for psSolutn_adj_t.
  !! Allocates and initializes adjoint and forward variables.
  SUBROUTINE constructor_adj(this, vLevel, num_cell)
    CLASS(psSolutn_adj_t), INTENT(INOUT) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, num_cell

    ! Allocate forward variables (assuming `rights` and `solutn` are in the base type)
    ALLOCATE (this%rights(vLevel, num_cell))
    ALLOCATE (this%solutn(vLevel, num_cell))

    ! Initialize forward variables to zero
    this%rights = 0.0_R_KIND
    this%solutn = 0.0_R_KIND

    ! Allocate adjoint variables
    ALLOCATE (this%adj_rights(vLevel, num_cell))
    ALLOCATE (this%adj_solutn(vLevel, num_cell))

    ! Initialize adjoint variables to zero
    this%adj_rights = 0.0_R_KIND
    this%adj_solutn = 0.0_R_KIND

  END SUBROUTINE constructor_adj

  !> @brief
  !! Adjoint destructor for psSolutn_adj_t.
  !! Deallocates adjoint and forward variables.
  IMPURE ELEMENTAL SUBROUTINE destructor_adj(this)
    TYPE(psSolutn_adj_t), INTENT(INOUT) :: this

    ! Deallocate adjoint variables
    IF (ALLOCATED(this%adj_rights)) DEALLOCATE (this%adj_rights)
    IF (ALLOCATED(this%adj_solutn)) DEALLOCATE (this%adj_solutn)

    ! Deallocate forward variables
    IF (ALLOCATED(this%rights)) DEALLOCATE (this%rights)
    IF (ALLOCATED(this%solutn)) DEALLOCATE (this%solutn)

  END SUBROUTINE destructor_adj

END MODULE psSolutn_adjoint_m
