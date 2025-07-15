!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.TimeIntegration.TimeIntegrationRK4
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-11   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides a Runge-Kutta 4th order time integration scheme.
MODULE TimeIntegrationRK4_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  ! USE C2MBase_m, ONLY: C2MBase_t
  USE obsMG_m, ONLY: obsMG_t
  ! USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  IMPLICIT NONE

  TYPE, EXTENDS(TimeIntegrationBase_t) :: TimeIntegrationRK4_t
    REAL(r_kind) :: steps(4)  !< RK4 time steps
    REAL(r_kind) :: coeff(4)  !< RK4 coefficients
  CONTAINS
    PROCEDURE, PUBLIC :: initialization_s => initializationRK4_s
    PROCEDURE, PUBLIC :: TimeIntegrationBase_fwd => TimeIntegrationRK4_fwd
    PROCEDURE, PUBLIC :: TimeIntegrationBase_adj => TimeIntegrationRK4_adj
  END TYPE TimeIntegrationRK4_t

CONTAINS
  SUBROUTINE initializationRK4_s(this, rhs, X, yamlFile)
    CLASS(TimeIntegrationRK4_t) :: this
    CLASS(rhsBase_t), TARGET :: rhs
    TYPE(State_t), INTENT(IN) :: X
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile

    this%rhs => rhs
    this%yamlFile = yamlFile
    this%steps = (/0.0D0, 0.5D0, 0.5D0, 1.0D0/)
    this%coeff = (/1.0D0 / 6.0D0, 2.0D0 / 6.0D0, 2.0D0 / 6.0D0, 1.0D0 / 6.0D0/)
  END SUBROUTINE initializationRK4_s

  ! This step forward routine integrates the right hand side one RK4 step forward
  ! and save the result in iC and X in the it frame even though it may not be the
  ! full time step of X. Noted by Yuanfu Xie 2024-10-10
  SUBROUTINE TimeIntegrationRK4_fwd(this, dt, int, it, IC, X, yk)

    IMPLICIT NONE

    CLASS(TimeIntegrationRK4_t) :: this
    REAL(r_kind), INTENT(IN) :: dt
    INTEGER(i_kind), INTENT(IN) :: int, it
    TYPE(State_t), INTENT(INOUT) :: iC
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(State_t), OPTIONAL, INTENT(INOUT) :: yk

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, l, num_vars
    TYPE(State_t) :: yk4, rk4

    ! For RK4, pass the initial condition to X' first time frame:

    num_vars = UBOUND(X%fields, 1)
    ! yk4 and rk4 need 4 time frames only:
    CALL yk4%initialize(this%yamlFile, X%sg, 4)
    CALL rk4%initialize(this%yamlFile, X%sg, 4)

    ! An integration step:
    rk4 = rk4%zeroCopy()

    ! Runge-Kutta steps:
    DO i = 1, 4
      ! Assign the yk value:
      IF (i .EQ. 20) STOP ! Temporary debugging option
      DO l = 1, num_vars
        yk4%fields(l)%DATA(:, :, i) = iC%fields(l)%DATA(:, :, 1) + &
                                      this%steps(i) * rk4%fields(l)%DATA(:, :, i)
      END DO

      ! Calculate the right hand side:
      CALL this%rhs%rightHandSide(i, yk4, rk4)
    END DO

    ! RK update:
    DO l = 1, num_vars
      X%fields(l)%DATA(:, :, it) = iC%fields(l)%DATA(:, :, 1) + dt * ( &
                                   this%coeff(1) * rk4%fields(l)%DATA(:, :, 1) + &
                                   this%coeff(2) * rk4%fields(l)%DATA(:, :, 2) + &
                                   this%coeff(3) * rk4%fields(l)%DATA(:, :, 3) + &
                                   this%coeff(4) * rk4%fields(l)%DATA(:, :, 4))
    END DO

    ! Save the model integrated at state time frame:
    DO l = 1, num_vars
      iC%fields(l)%DATA(:, :, k) = X%fields(l)%DATA(:, :, it)
    END DO

  END SUBROUTINE TimeIntegrationRK4_fwd

  SUBROUTINE TimeIntegrationRK4_adj(this, dt, it, X, Xadj)
    CLASS(TimeIntegrationRK4_t) :: this
    REAL(r_kind), INTENT(IN) :: dt
    INTEGER(i_kind), INTENT(IN) :: it
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(INOUT) :: Xadj

    PRINT *, '+************************ ABORT ***************************+'
    WRITE (*, 1)
1   FORMAT(' | RK adjoint requires calling the RHS adjoints 10 times    |', /, &
           ' | It is too expensive to run so it has not implemented yet |')
    PRINT *, '+**********************************************************+'
    STOP
  END SUBROUTINE TimeIntegrationRK4_adj
END MODULE TimeIntegrationRK4_m
