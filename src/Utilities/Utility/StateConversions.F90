!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-Utilities-Utility
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/23, @GBA-MWF, Shenzhen
!     A module for model states conversion.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains tranformation of model state variables.
!
MODULE stateConversions_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  TYPE :: stateConversions_t
    INTEGER(i_kind) :: nlvl
    REAL(r_kind), ALLOCATABLE :: coefDz(:), coef_sigma(:, :)

    ! QR factorization matrices:
    DOUBLE PRECISION, ALLOCATABLE :: hydroA(:, :)

  CONTAINS
    PROCEDURE, PUBLIC :: initial
    PROCEDURE, PUBLIC :: hydroPoisson

  END TYPE stateConversions_t

CONTAINS

  ! Initialization:
  SUBROUTINE initial(this, n, h)
    CLASS(stateConversions_t) :: this
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: h(n)

    ! Local variables:
    INTEGER(i_kind) :: i

    this%nlvl = n

    ! Allocate memory:
    ALLOCATE (this%hydroA(n - 2, n - 2), this%coefDz(n - 2), this%coef_sigma(3, n - 2))

    ! First order interpolation coefficients for an uneven grid:
    DO i = 2, n
      this%coef_sigma(1, i) = (h(i + 1) - h(i)) / (h(i - 1) - h(i)) / (h(i + 1) - h(i - 1))
      this%coef_sigma(2, i) = (h(i + 1) + h(i - 1) - 2.0D0 * h(i)) / (h(i + 1) - h(i)) / (h(i) - h(i - 1))
      this%coef_sigma(3, i) = (h(i) - h(i - 1)) / (h(i + 1) - h(i)) / (h(i + 1) - h(i - 1))
    END DO

  END SUBROUTINE initial

  ! A hydrostatic conversion through a Poisson solver with pressures are known at top and bottom
  SUBROUTINE hydroPoisson(this)
    CLASS(stateConversions_t) :: this
  END SUBROUTINE hydroPoisson

  ! The destructor:
  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(stateConversions_t), INTENT(INOUT) :: this

    ! Deallocate memory:
    IF (ALLOCATED(this%hydroA)) DEALLOCATE (this%hydroA, this%coefDz, this%coef_sigma)

    WRITE (*, 1)
1   FORMAT('End of destructor of stateConversions_t.')
  END SUBROUTINE destructor

END MODULE stateConversions_m
