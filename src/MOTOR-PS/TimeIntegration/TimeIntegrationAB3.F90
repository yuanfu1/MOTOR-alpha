!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.TimeIntegration.TimeIntegrationAB3
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-12   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides a Adam-Bashford 3rd order time integration scheme.
MODULE TimeIntegrationAB3_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  ! USE C2MBase_m, ONLY: C2MBase_t
  USE obsMG_m, ONLY: obsMG_t
  ! USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  USE ObsSetArray_m, ONLY: ObsSetArray_t

  TYPE, EXTENDS(TimeIntegrationBase_t):: TimeIntegrationAB3_t
    REAL(r_kind) :: coeff(3)
  CONTAINS
    PROCEDURE, PUBLIC :: initialization_s => initializationAB3_s
    PROCEDURE, PUBLIC :: TimeIntegrationBase_fwd => TimeIntegrationAB3_fwd
    PROCEDURE, PUBLIC :: TimeIntegrationBase_adj => TimeIntegrationAB3_adj
  END TYPE TimeIntegrationAB3_t

CONTAINS
  SUBROUTINE initializationAB3_s(this, rhs, X, yamlFile)
    CLASS(TimeIntegrationAB3_t) :: this
    CLASS(rhsBase_t), TARGET :: rhs
    TYPE(State_t), INTENT(IN) :: X
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile

    this%rhs => rhs
    this%yamlFile = yamlFile
    this%coeff = (/23.0D0 / 12.0D0, 16.0D0 / 12.0D0, 5.0D0 / 12.0D0/)

    CALL this%ab3_adj%initialize(this%yamlFile, X%sg, 3)
    this%ab3_adj = this%ab3_adj%zeroCopy()  !< Adjoint initial conditions
  END SUBROUTINE initializationAB3_s

  SUBROUTINE TimeIntegrationAB3_fwd(this, dt, int, it, iC, X, yk)
    IMPLICIT NONE
    CLASS(TimeIntegrationAB3_t) :: this
    REAL(r_kind), INTENT(IN) :: dt
    INTEGER(i_kind), INTENT(IN) :: int, it
    TYPE(State_t), INTENT(INOUT) :: iC
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(State_t), OPTIONAL, INTENT(INOUT) :: yk

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, l, num_vars

    num_vars = UBOUND(X%fields, 1)
    ! Check if there are three time steps of initial condition:
    IF (UBOUND(iC%fields(1)%DATA, 3) .LT. 3) THEN
      WRITE (*, 2) UBOUND(iC%fields(1)%DATA, 3)
2     FORMAT('TimeIntegrationAB3_fwd - Requires 3 IC steps but you provided: ', I2)
      STOP
    END IF

    ! Check if the right hand sides at previous time steps are available:
    IF (.NOT. PRESENT(yk)) THEN
      WRITE (*, 3)
3     FORMAT('For this AB3 time stepping, the RHS of the previous 2 time steps are required, rerun!')
      STOP
    END IF

    ! Get the right hand side values for the step:
    WRITE(*,1) int, it, dt
1   FORMAT('TimeIntegrationAB3_fwd - at intermedium step: ',I3,' and model frame: ',I4,' with dt: ',E14.6)
    CALL this%rhs%rightHandSide(3, iC, yk)

    ! AB3 update:
    DO l = 1, num_vars
      X%fields(l)%DATA(:, :, it) = iC%fields(l)%DATA(:, :, 3) + &
                                   dt * (this%coeff(1) * yk%fields(l)%DATA(:, :, 3) + &
                                         this%coeff(2) * yk%fields(l)%DATA(:, :, 2) + &
                                         this%coeff(3) * yk%fields(l)%DATA(:, :, 1))
    END DO
    PRINT*,'Updating forecast states...',it, (yk .DOT. yk), (iC .DOT. iC), (X .DOT. X)

    ! Shift the right hand sides and iC:
    DO i = 1, 2
      DO l = 1, num_vars
        iC%fields(l)%DATA(:, :, i) = iC%fields(l)%DATA(:, :, i + 1)
        yk%fields(l)%DATA(:, :, i) = yk%fields(l)%DATA(:, :, i + 1)
      END DO
    END DO
    ! Update the last step:
    DO l = 1, num_vars
      iC%fields(l)%DATA(:, :, 3) = X%fields(l)%DATA(:, :, it)
    END DO
  END SUBROUTINE TimeIntegrationAB3_fwd

  !<
  !< We have assume the forward model state stored in X:
  !< There may be several time steps between X time frames but the states are not saved
  SUBROUTINE TimeIntegrationAB3_adj(this, dt, it, X, Xadj)
    IMPLICIT NONE
    CLASS(TimeIntegrationAB3_t) :: this
    REAL(r_kind), INTENT(IN) :: dt
    INTEGER(i_kind), INTENT(IN) :: it
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Xadj

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, n, num_vars, intval
    REAL(r_kind) :: dTime

    num_vars = UBOUND(X%fields, 1)

    ! Ab3 adjoint needs the sum of adjoint variables at previous 3 time steps:
    DO j = 1, num_vars
      Xadj%fields(j)%DATA(:, :, it) = &
        this%ab3_adj%fields(j)%DATA(:, :, 1) + &
        this%ab3_adj%fields(j)%DATA(:, :, 2) + &
        this%ab3_adj%fields(j)%DATA(:, :, 3)
    END DO

    ! Calculate the right hand side:
    PRINT*,'AB3_adj: rhs is calculating...'
    CALL this%rhs%rightHandSide(it, X, Xadj)
  END SUBROUTINE TimeIntegrationAB3_adj
END MODULE TimeIntegrationAB3_m
