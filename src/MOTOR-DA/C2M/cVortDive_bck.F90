!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/4/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines a control variable case using vorticity and divergence as control variables.
!! Control variables: $\zeta$, $\delta$, $\ln(P)$, $temp$ and $rh$, where rh is scaled mixing ratio
!! by an exponential function profile.
!!
MODULE cVortDive_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE C2MBase_m, ONLY: C2MBase_t
  USE State_m, ONLY: State_t
  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t):: cVortDive_t
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => cVortDiveForward
    PROCEDURE, PUBLIC :: bckwardOpr_s => cVortDiveBckward
    PROCEDURE, PUBLIC :: tangentOpr_s => cVortDiveTangent
    PROCEDURE, PUBLIC :: adjointOpr_s => cVortDiveAdjoint
    PROCEDURE, PUBLIC :: bckAdjoint_s => cVortDiveBckAdjoint
  END TYPE cVortDive_t

CONTAINS

  SUBROUTINE initialize_s(this, mgStart, mgEnd, X)

    IMPLICIT NONE

    CLASS(cVortDive_t) :: this
    INTEGER(i_kind) :: mgStart, mgEnd
    TYPE(State_t), INTENT(IN) :: X(mgStart:mgEnd)

    ! Local variables:
    INTEGER(i_kind) :: i, img
    TYPE(SingleGrid_t) :: sg  ! Assuming sg refers to a SingleGrid type.

    ! Set analysis to the multigrid model states as default:
    ALLOCATE (this%analysis(mgStart:mgEnd))
    this%analysis = X

    ! Loop over the grid points
    DO img = mgStart, mgEnd
      ! Get the single grid information for the current analysis point
      sg = this%analysis(img)%sg

      ! Loop over the fields and check their names
      DO i = 1, UBOUND(this%analysis(img)%fields, 1)
        PRINT *, 'Module state: ', i, this%analysis(img)%fields(i)%Get_Name()

        SELECT CASE (TRIM(this%analysis(img)%fields(i)%Get_Name()))
        CASE ('vorticity', 'divergence', 'ln(P)', 'temp', 'rh')
          ! Set as control variables
          this%analysis(img)%fields(i)%numTotalMask = sg%num_icell_global
          this%analysis(img)%fields(i)%maskHorizontal = 1
          this%analysis(img)%fields(i)%maskTemporal = 0
          this%analysis(img)%fields(i)%maskTemporal(1) = 1
          this%analysis(img)%fields(i)%maskVertical = 0
        CASE DEFAULT
          ! Non-control variables
          this%analysis(img)%fields(i)%numTotalMask = 0
          this%analysis(img)%fields(i)%maskHorizontal = 0
          this%analysis(img)%fields(i)%maskTemporal = 0
          this%analysis(img)%fields(i)%maskVertical = 0
        END SELECT

      END DO
    END DO

  END SUBROUTINE initialize_s

  FUNCTION cVortDiveForward(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION cVortDiveForward

  FUNCTION cVortDiveBckward(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION cVortDiveBckward

  FUNCTION cVortDiveTangent(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION cVortDiveTangent

  FUNCTION cVortDiveAdjoint(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION cVortDiveAdjoint

  FUNCTION cVortDiveBckAdjoint(this, X) RESULT(Y)
    CLASS(cVortDive_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y
  END FUNCTION cVortDiveBckAdjoint
END MODULE cVortDive_m
