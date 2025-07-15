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
!! This module defines a 3DVar control variable case using U and V as control variables.
!! Control variables: $u$, $v$, $\ln(P)$, $temp$ and $rh$, where rh is scaled mixing ratio
!! by an exponential function profile.
!!
MODULE cDefault_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE C2MBase_m, ONLY: C2MBase_t
  USE State_m, ONLY: State_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE

  TYPE, EXTENDS(C2MBase_t):: cDefault_t
    REAL(r_kind), ALLOCATABLE :: verticalProfiles(:, :, :)  ! for scaling pressure and water vapor
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s
    PROCEDURE, PUBLIC :: forwardOpr_s => cDefaultForward
    PROCEDURE, PUBLIC :: tangentOpr_s => cDefaultTangent
    PROCEDURE, PUBLIC :: adjointOpr_s => cDefaultAdjoint
    PROCEDURE, PUBLIC :: bckwardOpr_s => cDefaultBckwardOpr
  END TYPE cDefault_t

CONTAINS
  SUBROUTINE initialize_s(this, sg, X, dyc)

    IMPLICIT NONE

    CLASS(cDefault_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), INTENT(IN) :: X
    CLASS(dyCoresBase_t), TARGET, OPTIONAL :: dyc

    ! Local variables:
    CHARACTER(LEN=1024) :: yamlFile
    INTEGER(i_kind) :: i

    ! Set analysis to the multigrid model states as default:
    this%analysis = X
    this%sg => sg

    ! Get yaml file:
    ! CALL GETARG(1, yamlFile)

    ! Check the state variables:
    ! ALLOCATE(this%verticalProfiles(this%analysis(mgEnd)%sg%vLevel, &
    !   UBOUND(this%analysis(mgEnd)%fields,1),mgStart:mgEnd))

    ! Use of Yali's water vapor scaling in MGOpts module:
    DO i = LBOUND(this%analysis%fields, 1), UBOUND(this%analysis%fields, 1)

      this%analysis%fields(i)%numTotalMask = &
        sg%num_icell_global * sg%vLevel * sg%tSlots   ! Default value

      IF (this%analysis%fields(i)%Get_Name() .EQ. 'qvapor' .OR. &
          this%analysis%fields(i)%Get_Name() .EQ. 'qvapor_ctl') THEN ! Temporarily leave the ctl name
        this%analysis%fields(i)%lowBound = 1
        ! The current analysis does not analyze the vapor fields until the last two levels
        ! Yuanfu Xie 2025-06-25 temporarily set the mask for level coarser than 7 to zero
        ! replace the original except the last two levels for the new design of MOTOR-variableDesign on 2025-06-25
        IF (this%sg%glevel .LT. 7) THEN
          this%analysis%fields(i)%maskVertical = 0
        END IF

        this%analysis%fields(i)%maskVertical(1) = 1  ! Surface qvapor is analyzed in the current
        this%analysis%fields(i)%numTotalMask = sg%num_icell * sg%tSlots
      END IF
      
      ! Pressure: Mask setting -- only surface pressures are control variables;
      ! Pressure at other level is computed using hydrostatic relation
      IF (this%analysis%fields(i)%Get_Name() .EQ. 'pres' .OR. &
          this%analysis%fields(i)%Get_Name() .EQ. 'psl_ctl') THEN
        this%analysis%fields(i)%maskVertical = 0
        this%analysis%fields(i)%maskVertical(1) = 1
        this%analysis%fields(i)%numTotalMask = sg%num_icell_global * sg%tSlots
      END IF
      IF (this%analysis%fields(i)%Get_Name() .EQ. 'temp') &
        this%analysis%fields(i)%lowBound = 1
    END DO

  END SUBROUTINE initialize_s

  FUNCTION cDefaultForward(this, X) RESULT(Y)
    CLASS(cDefault_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = X
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA * 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) * X%sg%SVapor
        END DO
      END SELECT
    END DO
  END FUNCTION cDefaultForward

  FUNCTION cDefaultBckwardOpr(this, X) RESULT(Y)
    CLASS(cDefault_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = this%analysis
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) / X%sg%SVapor
        END DO
      CASE DEFAULT
        Y%fields(j)%DATA = X%fields(j)%DATA
      END SELECT
    END DO
  END FUNCTION cDefaultBckwardOpr

  FUNCTION cDefaultTangent(this, X) RESULT(Y)
    CLASS(cDefault_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    Y = X ! Identical tangent
  END FUNCTION cDefaultTangent

  FUNCTION cDefaultAdjoint(this, X, dt) RESULT(Y)
    CLASS(cDefault_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(IN), OPTIONAL :: dt
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = this%analysis
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA * 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) * X%sg%SVapor
        END DO
      CASE DEFAULT
        Y%fields(j)%DATA = X%fields(j)%DATA
      END SELECT
    END DO
  END FUNCTION cDefaultAdjoint

  FUNCTION cDefaultBckAdjoint(this, X) RESULT(Y)
    CLASS(cDefault_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j

    Y = X ! Identity matrix transpose
    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      SELECT CASE (TRIM(X%fields(j)%Get_Name()))
      CASE ('pres', 'psl')
        Y%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
      CASE ('qvapor')
        DO i = 1, X%sg%tSlots
          Y%fields(j)%DATA(:, :, i) = X%fields(j)%DATA(:, :, i) / X%sg%SVapor
        END DO
      END SELECT
    END DO
  END FUNCTION cDefaultBckAdjoint
END MODULE cDefault_m
